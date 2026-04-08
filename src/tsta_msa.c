#define _GNU_SOURCE

#include "tsta_common.h"
#include "tsta_graph.h"
#include "tsta_threadpool.h"

#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MIN -120

static int tsta_sequence_view_make(tsta_sequence_view_t *view,
                                   const char *sequence, int length) {
  if (!view || !sequence || length < 0) {
    return -1;
  }

  view->sequence = sequence;
  view->length = length;
  return 0;
}

static int tsta_msa_final_score(tsta_graph_t *graph);
static int tsta_msa_run_alignment_step(tsta_msa_aligner *aligner,
                                       tsta_graph_t *graph,
                                       const tsta_sequence_view_t *sequence,
                                       int num, int sum);
static void tsta_msa_block_alignment(void *pa);

struct tsta_msa_aligner {
  tsta_config config;
  tsta_msa_state state;
  tsta_graph_t *graph;
  tsta_threadpool_t *pool;
  size_t sequence_count;
  packed_sequence_array_t packed_sequences;
};

static void tsta_msa_cache_reset(tsta_msa_aligner *aligner) {
  if (!aligner) {
    return;
  }

  tsta_packed_sequence_array_reset_mm(&aligner->packed_sequences);
}

static int tsta_msa_cache_reserve_slot(tsta_msa_aligner *aligner, size_t slot) {
  if (!aligner) {
    return -1;
  }

  return tsta_packed_sequence_array_reserve_slot(&aligner->packed_sequences,
                                                 slot, 4);
}

typedef struct msa_block_task {
  int diagonal_index;
  int diagonal_block_count;
  int diagonal_block_offset;
  int sequence_index;
  char *packed_sequence;
  tsta_graph_t *graph;
  tsta_msa_state *state;
} msa_block_task;

#define TSTA_MSA_PASSING_SEQ_SET(node, index)                                  \
  do {                                                                         \
    (node)->passing_seq[(size_t)(index) >> 3] |=                               \
        (unsigned char)(1u << ((size_t)(index) & 7u));                         \
  } while (0)

#define TSTA_MSA_PASSING_SEQ_GET(node, index)                                  \
  (((node)->passing_seq[(size_t)(index) >> 3] >> ((size_t)(index) & 7u)) & 1u)

static void tsta_msa_apply_config(tsta_msa_aligner *aligner,
                                  const tsta_config *config) {
  if (config) {
    aligner->config = *config;
  } else {
    tsta_config_default(&aligner->config);
  }

  tsta_init_msa_state(&aligner->state, &aligner->config, block);
  aligner->state.z = 0;
}

static int tsta_msa_build_alignment_rows(tsta_graph_t *graph,
                                         size_t sequence_count, char ***rows,
                                         size_t *aln_length) {
  char **result_rows;
  size_t len = 0;

  if (!graph || !rows || !aln_length || sequence_count == 0) {
    return -1;
  }

  result_rows = (char **)calloc(sequence_count, sizeof(char *));
  if (!result_rows) {
    return -1;
  }

  for (size_t i = 0; i < sequence_count; i++) {
    result_rows[i] = (char *)malloc((size_t)graph->len + 1);
    if (!result_rows[i]) {
      for (size_t j = 0; j < i; j++) {
        free(result_rows[j]);
      }
      free(result_rows);
      return -1;
    }
    memset(result_rows[i], '-', (size_t)graph->len);
  }

  for (int i = 0; i < graph->len; i++) {
    if (tsta_graph_sort_node(graph, i)->node_status == 0) {
      for (size_t j = 0; j < sequence_count; j++) {
        if (TSTA_MSA_PASSING_SEQ_GET(tsta_graph_sort_node(graph, i), j)) {
          result_rows[j][len] = tsta_graph_sort_node(graph, i)->base;
        }
      }
      for (int l = 0; l < tsta_graph_sort_node(graph, i)->mismatch_num; l++) {
        for (size_t j = 0; j < sequence_count; j++) {
          if (TSTA_MSA_PASSING_SEQ_GET(
                  tsta_graph_node(
                      graph, tsta_graph_sort_node(graph, i)->mismatch_node[l]),
                  j)) {
            result_rows[j][len] =
                tsta_graph_node(
                    graph, tsta_graph_sort_node(graph, i)->mismatch_node[l])
                    ->base;
            tsta_graph_node(graph,
                            tsta_graph_sort_node(graph, i)->mismatch_node[l])
                ->node_status = 1;
          }
        }
      }
      len++;
    }
  }

  for (size_t i = 0; i < sequence_count; i++) {
    result_rows[i][len] = '\0';
  }

  *rows = result_rows;
  *aln_length = len;
  return 0;
}

static int tsta_msa_capture_result(tsta_graph_t *graph, size_t sequence_count,
                                   tsta_msa_result_t *result) {
  char **rows = NULL;
  size_t aln_length = 0;

  if (!result) {
    return -1;
  }

  if (tsta_msa_build_alignment_rows(graph, sequence_count, &rows,
                                    &aln_length) != 0) {
    return -1;
  }

  tsta_msa_result_free(result);
  result->score = tsta_msa_final_score(graph);
  result->aln = rows;
  result->sequence_count = sequence_count;
  result->aln_length = aln_length;
  return 0;
}

static void tsta_msa_release_node_trace_storage(tsta_node_t *node) {
  if (!node) {
    return;
  }

  if (node->source_store) {
    tsta_trace_block_store_destroy(node->source_store);
    node->source_store = NULL;
  }
  if (node->esource_store) {
    tsta_trace_block_store_destroy(node->esource_store);
    node->esource_store = NULL;
  }
  if (node->fsource_store) {
    tsta_trace_block_store_destroy(node->fsource_store);
    node->fsource_store = NULL;
  }
}

static void tsta_msa_release_trace_storage(tsta_graph_t *graph) {
  if (!graph || !graph->sort.data) {
    return;
  }

  for (int i = 0; i < graph->len; i++) {
    tsta_msa_release_node_trace_storage(tsta_graph_sort_node(graph, i));
  }
}

static void tsta_msa_release_alignment_work_buffers(tsta_graph_t *graph) {
  if (!graph || !graph->nodes.data) {
    return;
  }

  for (size_t i = 0; i < graph->nodes.length; i++) {
    tsta_node_t *node = &graph->nodes.data[i];

    if (!node) {
      continue;
    }

    free(node->simple_sorce);
    node->simple_sorce = NULL;
  }
}

static size_t tsta_msa_estimate_trace_bytes(size_t sequence_count,
                                            size_t trace_length) {
  size_t per_node;

  if (sequence_count == 0 || trace_length == 0) {
    return 0;
  }

  if (trace_length > SIZE_MAX / 3) {
    return SIZE_MAX;
  }
  per_node = trace_length * 3;
  if (sequence_count > SIZE_MAX / per_node) {
    return SIZE_MAX;
  }
  return sequence_count * per_node;
}

static int tsta_msa_should_spill_trace(const tsta_msa_aligner *aligner,
                                       size_t sequence_count,
                                       size_t trace_length) {
  size_t estimated_bytes;

  if (!aligner) {
    return 0;
  }

  if (aligner->config.msa_trace_dump_threshold == 0) {
    return 0;
  }

  estimated_bytes = tsta_msa_estimate_trace_bytes(sequence_count, trace_length);
  return estimated_bytes >= aligner->config.msa_trace_dump_threshold;
}

static int tsta_msa_prepare_trace_storage(tsta_graph_t *graph,
                                          size_t trace_length,
                                          size_t chunk_size, int compress,
                                          int spill_immediately) {
  if (!graph || !graph->sort.data || trace_length == 0) {
    return -1;
  }

  for (int i = 0; i < graph->len; i++) {
    tsta_node_t *node = tsta_graph_sort_node(graph, i);

    if (!node) {
      continue;
    }

    tsta_msa_release_node_trace_storage(node);
    node->source_store =
        tsta_trace_block_store_create(trace_length, chunk_size, compress);
    if (spill_immediately && node->source_store &&
        tsta_trace_block_store_spill(node->source_store) != 0) {
      tsta_msa_release_trace_storage(graph);
      return -1;
    }
    node->esource_store =
        tsta_trace_block_store_create(trace_length, chunk_size, compress);
    if (spill_immediately && node->esource_store &&
        tsta_trace_block_store_spill(node->esource_store) != 0) {
      tsta_msa_release_trace_storage(graph);
      return -1;
    }
    node->fsource_store =
        tsta_trace_block_store_create(trace_length, chunk_size, compress);
    if (spill_immediately && node->fsource_store &&
        tsta_trace_block_store_spill(node->fsource_store) != 0) {
      tsta_msa_release_trace_storage(graph);
      return -1;
    }
    if (!node->source_store || !node->esource_store || !node->fsource_store) {
      tsta_msa_release_trace_storage(graph);
      return -1;
    }
  }

  return 0;
}

static void tsta_msa_spill_trace_storage(tsta_graph_t *graph) {
  if (!graph || !graph->sort.data) {
    return;
  }

  for (int i = 0; i < graph->len; i++) {
    tsta_node_t *node = tsta_graph_sort_node(graph, i);

    if (!node) {
      continue;
    }
    if (node->source_store) {
      tsta_trace_block_store_spill(node->source_store);
    }
    if (node->esource_store) {
      tsta_trace_block_store_spill(node->esource_store);
    }
    if (node->fsource_store) {
      tsta_trace_block_store_spill(node->fsource_store);
    }
  }
}

static int tsta_msa_final_score(tsta_graph_t *graph) {
  int score = INT_MIN;

  if (!graph || !graph->sort.data) {
    return 0;
  }

  for (int i = 0; i < graph->len; i++) {
    if (tsta_graph_sort_node(graph, i) &&
        tsta_graph_sort_node(graph, i)->out == 0 &&
        score < tsta_graph_sort_node(graph, i)->lastsorce) {
      score = tsta_graph_sort_node(graph, i)->lastsorce;
    }
  }

  return score == INT_MIN ? 0 : score;
}

static void tsta_msa_session_reset(tsta_msa_aligner *aligner) {
  if (!aligner) {
    return;
  }

  if (aligner->pool) {
    tsta_threadpool_destroy(aligner->pool);
    aligner->pool = NULL;
  }
  if (aligner->graph) {
    tsta_msa_free_graph(aligner->graph);
    aligner->graph = NULL;
  }
  aligner->sequence_count = 0;
}

tsta_msa_aligner *tsta_msa_aligner_create(const tsta_config *config) {
  tsta_msa_aligner *aligner = (tsta_msa_aligner *)calloc(1, sizeof(*aligner));

  if (!aligner) {
    return NULL;
  }

  tsta_msa_apply_config(aligner, config);
  return aligner;
}

void tsta_msa_aligner_destroy(tsta_msa_aligner *aligner) {
  if (!aligner) {
    return;
  }
  tsta_msa_session_reset(aligner);
  tsta_msa_cache_reset(aligner);
  free(aligner);
}

int tsta_msa_aligner_begin(tsta_msa_aligner *aligner,
                           const char *initial_sequence,
                           int initial_sequence_length) {
  tsta_msa_state *state;
  tsta_sequence_view_t initial_view;

  if (!aligner || tsta_sequence_view_make(&initial_view, initial_sequence,
                                          initial_sequence_length) != 0) {
    return -1;
  }

  tsta_msa_session_reset(aligner);
  state = &aligner->state;
  tsta_init_msa_state(state, &aligner->config, block);
  state->z = 0;

  aligner->pool = tsta_threadpool_create(
      aligner->config.threads > 0 ? aligner->config.threads : 10, 100,
      sizeof(msa_block_task));
  if (!aligner->pool) {
    return -1;
  }

  aligner->graph = (tsta_graph_t *)calloc(1, sizeof(tsta_graph_t));
  if (!aligner->graph) {
    tsta_msa_session_reset(aligner);
    return -1;
  }

  if (!tsta_graph_build_initial(aligner->graph, &initial_view, 1,
                                &aligner->state)) {
    tsta_msa_session_reset(aligner);
    return -1;
  }
  aligner->graph->p = 1;
  aligner->sequence_count = 1;
  return 0;
}

int tsta_msa_aligner_add(tsta_msa_aligner *aligner, const char *sequence,
                         int sequence_length) {
  size_t new_count;
  tsta_sequence_view_t sequence_view;
  tsta_graph_t *graph;

  if (!aligner || !aligner->graph || !aligner->pool ||
      aligner->sequence_count == 0 ||
      tsta_sequence_view_make(&sequence_view, sequence, sequence_length) != 0) {
    return -1;
  }

  new_count = aligner->sequence_count + 1;
  if (tsta_msa_expand_passing_seq(aligner->graph, aligner->sequence_count,
                                  new_count) != 0) {
    tsta_msa_session_reset(aligner);
    return -1;
  }

  graph = aligner->graph;
  if (tsta_msa_run_alignment_step(aligner, graph, &sequence_view,
                                  (int)aligner->sequence_count,
                                  (int)new_count) != 0) {
    tsta_msa_free_graph(graph);
    aligner->graph = NULL;
    tsta_msa_session_reset(aligner);
    return -1;
  }
  aligner->graph = tsta_graph_update_impl(
      graph, &sequence_view, (int)aligner->sequence_count, (int)new_count,
      aligner->state.maxtag + 1, &aligner->state);
  if (!aligner->graph) {
    tsta_msa_free_graph(graph);
    aligner->graph = NULL;
    tsta_msa_session_reset(aligner);
    return -1;
  }
  tsta_msa_release_alignment_work_buffers(aligner->graph);
  aligner->graph = tsta_graph_sort(aligner->graph, 0);
  if (!aligner->graph) {
    tsta_msa_session_reset(aligner);
    return -1;
  }
  aligner->sequence_count = new_count;
  return 0;
}

int tsta_msa_aligner_get_result(tsta_msa_aligner *aligner,
                                tsta_msa_result_t *result) {
  if (!aligner || !aligner->graph || aligner->sequence_count == 0 || !result) {
    return -1;
  }

  aligner->graph = tsta_graph_sort(aligner->graph, 1);
  if (!aligner->graph) {
    return -1;
  }
  return tsta_msa_capture_result(aligner->graph, aligner->sequence_count,
                                 result);
}

int tsta_msa_aligner_align(tsta_msa_aligner *aligner,
                           const char *const *sequences,
                           const int *sequence_lengths, size_t sequence_count,
                           tsta_msa_result_t *result) {
  if (!aligner || !sequences || !sequence_lengths || sequence_count == 0 ||
      !result) {
    return -1;
  }

  if (tsta_msa_aligner_begin(aligner, sequences[0], sequence_lengths[0]) != 0) {
    return -1;
  }

  for (size_t i = 1; i < sequence_count; i++) {
    if (tsta_msa_aligner_add(aligner, sequences[i], sequence_lengths[i]) != 0) {
      tsta_msa_session_reset(aligner);
      return -1;
    }
  }

  int status = tsta_msa_aligner_get_result(aligner, result);
  tsta_msa_session_reset(aligner);
  return status;
}

int tsta_msa_align(const char *const *sequences, const int *sequence_lengths,
                   size_t sequence_count, const tsta_config *config,
                   tsta_msa_result_t *result) {
  tsta_config effective_config;
  size_t longest_sequence_length = 0;
  tsta_msa_aligner *aligner;

  if (!sequences || !sequence_lengths || sequence_count == 0 || !result) {
    return -1;
  }

  for (size_t i = 0; i < sequence_count; i++) {
    if (sequence_lengths[i] > 0 &&
        (size_t)sequence_lengths[i] > longest_sequence_length) {
      longest_sequence_length = (size_t)sequence_lengths[i];
    }
  }

  effective_config = config ? *config : tsta_config_make_default();
  tsta_config_apply_length_fallback(&effective_config, longest_sequence_length,
                                    sequence_count);

  aligner = tsta_msa_aligner_create(&effective_config);
  if (!aligner) {
    return -1;
  }

  int status = tsta_msa_aligner_align(aligner, sequences, sequence_lengths,
                                      sequence_count, result);
  tsta_msa_aligner_destroy(aligner);
  return status;
}

#define nconvert(n)                                                            \
  (state->maxtag > 0 ? (!!(n / state->maxtag)) * state->maxtag +               \
                           (!(n / state->maxtag)) * (n % state->maxtag)        \
                     : 0)
#define NUM2(num2)                                                             \
  ((num2) / state->L) * state->L +                                             \
      ((((num2) % state->L) % state->W) * state->B +                           \
       (((num2) % state->L) / state->W))

#define TRACE_SOURCE(node, index)                                              \
  tsta_trace_block_store_get_byte((node)->source_store, (size_t)(index))
#define TRACE_SOURCE_SET(node, index, value)                                   \
  tsta_trace_block_store_set_byte((node)->source_store, (size_t)(index),       \
                                  (value))
#define TRACE_ESOURCE(node, index)                                             \
  tsta_trace_block_store_get_byte((node)->esource_store, (size_t)(index))
#define TRACE_ESOURCE_SET(node, index, value)                                  \
  tsta_trace_block_store_set_byte((node)->esource_store, (size_t)(index),      \
                                  (value))
#define TRACE_FSOURCE(node, index)                                             \
  tsta_trace_block_store_get_byte((node)->fsource_store, (size_t)(index))
#define TRACE_FSOURCE_SET(node, index, value)                                  \
  tsta_trace_block_store_set_byte((node)->fsource_store, (size_t)(index),      \
                                  (value))

typedef struct msa_thread_workspace {
  int capacity_num;
  int lane_len;
  int block_width;
  char **f_temp;
  char **VC2;
  char **VC1;
  char **r_s;
  char *h_g;
  char *v0;
  char *vc_1;
  char *vc_2;
  int *pd;
  int *te;
} msa_thread_workspace;

static pthread_key_t msa_workspace_key;
static pthread_once_t msa_workspace_once = PTHREAD_ONCE_INIT;
static int msa_workspace_key_ready = 0;

static void msa_workspace_reset(msa_thread_workspace *ws) {
  if (!ws) {
    return;
  }

  if (ws->f_temp) {
    for (int i = 0; i < ws->capacity_num; i++) {
      mm_free(ws->f_temp[i]);
    }
    free(ws->f_temp);
  }
  if (ws->VC2) {
    for (int i = 0; i < ws->capacity_num; i++) {
      mm_free(ws->VC2[i]);
    }
    free(ws->VC2);
  }
  if (ws->VC1) {
    for (int i = 0; i < ws->capacity_num; i++) {
      mm_free(ws->VC1[i]);
    }
    free(ws->VC1);
  }
  if (ws->r_s) {
    for (int i = 0; i < ws->capacity_num; i++) {
      mm_free(ws->r_s[i]);
    }
    free(ws->r_s);
  }

  mm_free(ws->h_g);
  free(ws->v0);
  free(ws->vc_1);
  free(ws->vc_2);
  free(ws->pd);
  free(ws->te);

  ws->f_temp = NULL;
  ws->VC2 = NULL;
  ws->VC1 = NULL;
  ws->r_s = NULL;
  ws->h_g = NULL;
  ws->v0 = NULL;
  ws->vc_1 = NULL;
  ws->vc_2 = NULL;
  ws->pd = NULL;
  ws->te = NULL;
  ws->capacity_num = 0;
  ws->lane_len = 0;
  ws->block_width = 0;
}

static void msa_workspace_destroy(void *ptr) {
  msa_thread_workspace *ws = (msa_thread_workspace *)ptr;

  msa_workspace_reset(ws);
  free(ws);
}

static void msa_workspace_make_key(void) {
  if (pthread_key_create(&msa_workspace_key, msa_workspace_destroy) == 0) {
    msa_workspace_key_ready = 1;
  }
}

static int msa_workspace_resize(msa_thread_workspace *ws, int num, int lane_len,
                                int block_width) {
  int old_cap;

  if (!ws || num <= 0 || lane_len <= 0 || block_width <= 0) {
    return -1;
  }

  old_cap = ws->capacity_num;
  if (ws->capacity_num < num) {
    char **tmp;

    tmp = (char **)realloc(ws->f_temp, (size_t)num * sizeof(char *));
    if (!tmp) {
      return -1;
    }
    ws->f_temp = tmp;

    tmp = (char **)realloc(ws->VC2, (size_t)num * sizeof(char *));
    if (!tmp) {
      return -1;
    }
    ws->VC2 = tmp;

    tmp = (char **)realloc(ws->VC1, (size_t)num * sizeof(char *));
    if (!tmp) {
      return -1;
    }
    ws->VC1 = tmp;

    tmp = (char **)realloc(ws->r_s, (size_t)num * sizeof(char *));
    if (!tmp) {
      return -1;
    }
    ws->r_s = tmp;

    for (int i = old_cap; i < num; i++) {
      ws->f_temp[i] = NULL;
      ws->VC2[i] = NULL;
      ws->VC1[i] = NULL;
      ws->r_s[i] = NULL;
    }
    ws->capacity_num = num;
  }

  if (!ws->v0) {
    ws->v0 = (char *)malloc((size_t)ws->capacity_num * sizeof(char));
  }
  if (!ws->vc_1) {
    ws->vc_1 = (char *)malloc((size_t)ws->capacity_num * sizeof(char));
  }
  if (!ws->vc_2) {
    ws->vc_2 = (char *)malloc((size_t)ws->capacity_num * sizeof(char));
  }
  if (!ws->pd) {
    ws->pd = (int *)malloc((size_t)ws->capacity_num * sizeof(int));
  }
  if (!ws->te) {
    ws->te = (int *)malloc((size_t)ws->capacity_num * sizeof(int));
  }

  if (!ws->v0 || !ws->vc_1 || !ws->vc_2 || !ws->pd || !ws->te) {
    return -1;
  }

  if (ws->lane_len != lane_len || !ws->h_g) {
    char *new_h_g;

    mm_free(ws->h_g);
    ws->h_g = NULL;
    new_h_g = (char *)mm_malloc((size_t)lane_len * sizeof(char));
    if (!new_h_g) {
      return -1;
    }
    ws->h_g = new_h_g;
  }

  if (ws->block_width != block_width) {
    for (int i = 0; i < ws->capacity_num; i++) {
      mm_free(ws->f_temp[i]);
      mm_free(ws->VC2[i]);
      mm_free(ws->VC1[i]);
      mm_free(ws->r_s[i]);
      ws->f_temp[i] = NULL;
      ws->VC2[i] = NULL;
      ws->VC1[i] = NULL;
      ws->r_s[i] = NULL;
    }
  }

  for (int i = 0; i < ws->capacity_num; i++) {
    if (!ws->f_temp[i]) {
      ws->f_temp[i] = (char *)mm_malloc((size_t)block_width * sizeof(char));
    }
    if (!ws->VC2[i]) {
      ws->VC2[i] = (char *)mm_malloc((size_t)block_width * sizeof(char));
    }
    if (!ws->VC1[i]) {
      ws->VC1[i] = (char *)mm_malloc((size_t)block_width * sizeof(char));
    }
    if (!ws->r_s[i]) {
      ws->r_s[i] = (char *)mm_malloc((size_t)block_width * sizeof(char));
    }
    if (!ws->f_temp[i] || !ws->VC2[i] || !ws->VC1[i] || !ws->r_s[i]) {
      return -1;
    }
  }

  ws->lane_len = lane_len;
  ws->block_width = block_width;
  return 0;
}

static msa_thread_workspace *msa_workspace_acquire(int num, int lane_len,
                                                   int block_width) {
  msa_thread_workspace *ws;

  pthread_once(&msa_workspace_once, msa_workspace_make_key);
  if (!msa_workspace_key_ready) {
    return NULL;
  }

  ws = (msa_thread_workspace *)pthread_getspecific(msa_workspace_key);
  if (!ws) {
    ws = (msa_thread_workspace *)calloc(1, sizeof(*ws));
    if (!ws) {
      return NULL;
    }
    if (pthread_setspecific(msa_workspace_key, ws) != 0) {
      free(ws);
      return NULL;
    }
  }

  if (msa_workspace_resize(ws, num, lane_len, block_width) != 0) {
    return NULL;
  }

  return ws;
}

static char *pack_sequence_for_simd(tsta_msa_aligner *aligner,
                                    size_t sequence_slot,
                                    const tsta_sequence_view_t *sequence,
                                    int poal, tsta_msa_state *state) {
  packed_sequence_t *slot_entry;
  char *seq2;
  int a = (int)sequence->length;
  int packed_len;
  int required_capacity;

  if (!aligner || !sequence || !state) {
    return NULL;
  }

  state->length1 = a;
  state->length2 = poal;
  if (a % state->L != 0)
    state->length1 = a + (state->L - a % state->L);
  packed_len = state->length1;

  if (tsta_msa_cache_reserve_slot(aligner, sequence_slot) != 0) {
    return NULL;
  }
  slot_entry = &aligner->packed_sequences.data[sequence_slot];

  required_capacity = state->length1 + 1;
  if (slot_entry->len < required_capacity) {
    mm_free(slot_entry->packed_seq);
    slot_entry->packed_seq =
        (char *)mm_malloc((size_t)required_capacity * sizeof(char));
    if (!slot_entry->packed_seq) {
      slot_entry->len = 0;
      return NULL;
    }
    slot_entry->len = required_capacity;
  }

  seq2 = slot_entry->packed_seq;

  for (int i = 0; i < packed_len; i++) {
    int source_index = i / state->L * state->L +
                       ((i % state->L) % state->B) * state->W +
                       ((i % state->L) / state->B);
    seq2[i] = source_index < a ? sequence->sequence[source_index] : 'N';
  }
  seq2[packed_len] = '\0';
  return seq2;
}

static int tsta_msa_run_alignment_step(tsta_msa_aligner *aligner,
                                       tsta_graph_t *graph,
                                       const tsta_sequence_view_t *sequence,
                                       int num, int sum) {
  unsigned int tsl;
  int j = 0;
  tsta_msa_state *state;
  char *packed_sequence;

  if (!aligner || !graph || !sequence) {
    return -1;
  }

  state = &aligner->state;
  packed_sequence =
      pack_sequence_for_simd(aligner, (size_t)num, sequence, graph->len, state);
  if (!packed_sequence) {
    tsta_msa_release_trace_storage(graph);
    return -1;
  }

  state->length2 = graph->len;
  if (graph->len % state->L != 0) {
    state->length2 = graph->len + (state->L - graph->len % state->L);
  }
  tsl = (unsigned int)((state->length1 + state->length2) / state->L - 1);

  state->trace_dump_compress = aligner->config.msa_trace_dump_compress ? 1 : 0;
  state->trace_dump_enabled = tsta_msa_should_spill_trace(
      aligner, (size_t)graph->len, (size_t)state->length1);

  tsta_msa_release_trace_storage(graph);
  if (tsta_msa_prepare_trace_storage(graph, (size_t)state->length1, 4096,
                                     state->trace_dump_compress,
                                     state->trace_dump_enabled) != 0) {
    tsta_msa_release_trace_storage(graph);
    return -1;
  }

  if (state->length1 >= state->length2) {
    state->fmaxtag = state->length2 / state->L - 1;
    state->lmaxtag = state->length1 / state->L - 1;
  } else {
    state->fmaxtag = state->length1 / state->L - 1;
    state->lmaxtag = state->length2 / state->L - 1;
  }
  state->maxtag = state->length1 / state->L - 1;

  for (int i = 0; i < graph->len; i++) {
    tsta_node_t *node = tsta_graph_sort_node(graph, (size_t)i);

    free(node->simple_sorce);
    node->simple_sorce =
        (int *)malloc((size_t)(state->maxtag + 2) * sizeof(int));
    if (!node->simple_sorce) {
      return -1;
    }
  }

  for (unsigned int i = 0; i < tsl; i++) {
    if (i <= (unsigned int)state->fmaxtag) {
      j++;
    } else if (i <= (unsigned int)state->lmaxtag) {
    } else {
      j--;
    }

    for (int l = 0; l < j; l++) {
      msa_block_task task = {
          .diagonal_index = (int)i,
          .diagonal_block_count = j,
          .diagonal_block_offset = l,
          .sequence_index = num,
          .packed_sequence = packed_sequence,
          .graph = graph,
          .state = state,
      };

      if (tsta_threadpool_submit(aligner->pool, tsta_msa_block_alignment, &task,
                                 sizeof(task)) != 0) {
        tsta_msa_block_alignment(&task);
      }
    }
    tsta_threadpool_wait(aligner->pool);
  }

  if (state->trace_dump_enabled) {
    tsta_msa_spill_trace_storage(graph);
  }

  return 0;
}

void block_line_alignment(tsta_graph_t *graph, int block_i, int block_j,
                          int block_l, tsta_node_t *row, char *seq, int nv,
                          int pc2, tsta_msa_state *state) {
  int m1, m2, m3;
  m1 = m2 = m3 = 0;
  char logo = -6;
  char Logo1 = 60;
  char Logo = 100;
  short reduce = 0;
  int original_pre_num = row->in;
  int pre_num = original_pre_num > 0 ? original_pre_num : 1;
  size_t row_width = (size_t)state->B;
  size_t row_alignment = (size_t)block;
  if (row_width < row_alignment) {
    row_width = row_alignment;
  } else if (row_width % row_alignment != 0) {
    row_width += row_alignment - (row_width % row_alignment);
  }
  size_t trace_chunk_index = (size_t)pc2 / (size_t)state->W;
  char *source_chunk =
      row->source_store
          ? tsta_trace_block_store_chunk(row->source_store, trace_chunk_index)
          : NULL;
  char *esource_chunk =
      row->esource_store
          ? tsta_trace_block_store_chunk(row->esource_store, trace_chunk_index)
          : NULL;
  char *fsource_chunk =
      row->fsource_store
          ? tsta_trace_block_store_chunk(row->fsource_store, trace_chunk_index)
          : NULL;

  char *h_g = NULL;
  char(*f_temp)[row_width] = NULL;
  char(*VC2)[row_width] = NULL;
  char(*VC1)[row_width] = NULL;
  char(*r_s)[row_width] = NULL;
  char *v0 = NULL;
  char *vc_1 = NULL;
  char *vc_2 = NULL;
  int *pd = NULL;
  int *te = NULL;

  if (!source_chunk || !esource_chunk || !fsource_chunk) {
    return;
  }

  h_g = (char *)tsta_aligned_malloc((size_t)state->L * sizeof(char), 64);
  f_temp = (char(*)[row_width])tsta_aligned_malloc(
      (size_t)pre_num * row_width * sizeof(char), 64);
  VC2 = (char(*)[row_width])tsta_aligned_malloc(
      (size_t)pre_num * row_width * sizeof(char), 64);
  VC1 = (char(*)[row_width])tsta_aligned_malloc(
      (size_t)pre_num * row_width * sizeof(char), 64);
  r_s = (char(*)[row_width])tsta_aligned_malloc(
      (size_t)pre_num * row_width * sizeof(char), 64);
  v0 = (char *)malloc((size_t)pre_num * sizeof(char));
  vc_1 = (char *)malloc((size_t)pre_num * sizeof(char));
  vc_2 = (char *)malloc((size_t)pre_num * sizeof(char));
  pd = (int *)malloc((size_t)pre_num * sizeof(int));
  te = (int *)malloc((size_t)pre_num * sizeof(int));
  if (!h_g || !f_temp || !VC2 || !VC1 || !r_s || !v0 || !vc_1 || !vc_2 || !pd ||
      !te) {
    tsta_aligned_free(h_g);
    tsta_aligned_free(f_temp);
    tsta_aligned_free(VC2);
    tsta_aligned_free(VC1);
    tsta_aligned_free(r_s);
    free(v0);
    free(vc_1);
    free(vc_2);
    free(pd);
    free(te);
    return;
  }

  if (original_pre_num == 0) {
    if (block_i == 0) {
      tsta_graph_node(graph, row->prev[0])->sorce[0] = state->O + state->E;
      tsta_graph_node(graph, row->prev[0])->esorce[0] =
          2 * (state->O + state->E);
    } else {
      tsta_graph_node(graph, row->prev[0])->sorce[0] = state->E;
      tsta_graph_node(graph, row->prev[0])->esorce[0] =
          state->E + state->O + state->E;
    }
    row->frist_col_sorce = row->simple_sorce[0] = state->E + state->O;
  }

  for (int i = 0; i < pre_num; i++)
    pd[i] = (tsta_graph_node(graph, row->prev[i])->node_status / 3) * pc2;
  int pc1 = (row->node_status / 3) * pc2;

  if (block_i <= state->lmaxtag && block_l == block_j - 1 && row->in != 0) {
    row->frist_col_sorce =
        tsta_graph_node(graph, row->prev[0])->frist_col_sorce + state->E;
    for (int i = 1; i < pre_num; i++)
      if (row->frist_col_sorce <
          tsta_graph_node(graph, row->prev[i])->frist_col_sorce + state->E)
        row->frist_col_sorce =
            tsta_graph_node(graph, row->prev[i])->frist_col_sorce + state->E;
    row->simple_sorce[0] = row->frist_col_sorce;
    for (int i = 0; i < pre_num; i++) {
      te[i] = row->frist_col_sorce -
              tsta_graph_node(graph, row->prev[i])->frist_col_sorce;
      if (te[i] > Logo) {
        v0[i] = Logo;
        if (te[i] - Logo > 127) {
          vc_2[i] = VC2[i][0] =
              (te[i] - Logo - 127) > 127 ? 127 : (te[i] - Logo - 127);
          vc_1[i] = VC1[i][0] = 127;
        } else {
          vc_2[i] = VC2[i][0] = 0;
          vc_1[i] = VC1[i][0] = te[i] - Logo;
        }
      } else {
        v0[i] = te[i];
        vc_2[i] = VC2[i][0] = 0;
        vc_1[i] = VC1[i][0] = 0;
      }
    }
  } else {
    if (tsta_graph_node(graph, row->prev[0])->sub == -1) {
      v0[0] = row->simple_sorce[nv] -
              (nv * state->L * state->E + (nv > 0 ? state->O : 0));
      vc_2[0] = VC2[0][0] = 0;
      vc_1[0] = VC1[0][0] = 0;
    } else {
      for (int i = 0; i < pre_num; i++) {
        te[i] = row->simple_sorce[nv] -
                tsta_graph_node(graph, row->prev[i])->simple_sorce[nv];
        if (te[i] > Logo) {
          v0[i] = Logo;
          if (te[i] - Logo > 127) {
            vc_2[i] = VC2[i][0] =
                (te[i] - Logo - 127) > 127 ? 127 : (te[i] - Logo - 127);
            vc_1[i] = VC1[i][0] = 127;
          } else {
            vc_2[i] = VC2[i][0] = 0;
            vc_1[i] = VC1[i][0] = te[i] - Logo;
          }
        } else {
          v0[i] = te[i];
          vc_2[i] = VC2[i][0] = 0;
          vc_1[i] = VC1[i][0] = 0;
        }
      }
    }
  }

  if (block_i <= state->lmaxtag && block_l == block_j - 1 &&
      block_i < state->length2 / state->L) {
    if (row->in == 0)
      row->f0[0] = v0[0] + state->E + state->O;
    else
      for (int i = 0; i < row->in; i++)
        row->f0[i] = v0[i] + state->E + state->O;
  }
  __mxxxi zero, top, Smin, h, max, emax, eumax, source, source_num, esource,
      esource_num, fsource, s0, s1, s2, s3, s4, s5, s6, temp, temp1, temp2, mat,
      mis, egap, ogap, base, N, z;
  __mask mask, mask1, mask2, mask3, mask4, mask5, SN, SM, SX;
  __mxxxi y[pre_num], vc2[pre_num], vc1[pre_num], vc0[pre_num], diff[pre_num],
      t[pre_num], e[pre_num], eu[pre_num], ev[pre_num], f[pre_num], fv[pre_num],
      q[pre_num], v[pre_num], sum[pre_num];
  z = mm_set1_epi8(Logo1);
  zero = mm_setzero();
  top = mm_set1_epi8(127);
  Smin = mm_set1_epi8(MIN);
  for (int i = 0; i < pre_num; i++) {
    sum[i] = mm_setzero();
    for (int j = 0; j < state->W; j++) {
      s1 = mm_load((__mxxxi *)tsta_graph_node(graph, row->prev[i])->sorce +
                   pd[i] + j);
      sum[i] = mm_add_epi8(sum[i], s1);
    }
    mm_store((__mxxxi *)r_s[i], sum[i]);
  }

  if (pre_num != 1) {
    for (int i = 0; i < pre_num; i++)
      f_temp[i][0] = v0[i];
    for (int j = 1; j < state->B; j++) {
      for (int i = 0; i < pre_num; i++) {
        te[i] = te[i] - r_s[i][j - 1] + state->W * state->E;
      }
      m1 = te[0];
      for (int s = 1; s < pre_num; s++) {
        if (te[s] < m1)
          m1 = te[s];
      }
      m2 = logo - m1;
      for (int i = 0; i < pre_num; i++) {
        if (te[i] + m2 > Logo) {
          f_temp[i][j] = Logo;
          if (te[i] + m2 - Logo > 127) {
            VC2[i][j] = (te[i] + m2 - Logo - 127) > 127
                            ? 127
                            : (te[i] + m2 - Logo - 127);
            VC1[i][j] = 127;
          } else {
            VC2[i][j] = 0;
            VC1[i][j] = te[i] + m2 - Logo;
          }
        } else {
          f_temp[i][j] = te[i] + m2;
          VC2[i][j] = 0;
          VC1[i][j] = 0;
        }
      }
    }
    for (int i = 0; i < pre_num; i++)
      v[i] = mm_load((__mxxxi *)f_temp[i]);
  } else {
    for (int j = 0; j < state->B; j++) {
      VC2[0][j] = 0;
      VC1[0][j] = 0;
    }
    vc_1[0] = vc_2[0] = 0;
    v[0] = mm_set1_epi8(state->E);
    v[0] = mm_insert_epi8(v[0], v0[0], 0);
  }
  mat = mm_set1_epi8(state->M);
  mis = mm_set1_epi8(state->X);
  egap = mm_set1_epi8(state->E);
  ogap = mm_set1_epi8(state->O + state->E);
  base = mm_set1_epi8(row->base);
  N = mm_set1_epi8('N');
  for (int j = 0; j < pre_num; j++) {
    vc2[j] = mm_load((__mxxxi *)VC2[j]);
    vc1[j] = mm_load((__mxxxi *)VC1[j]); ///
    f[j] = Smin;
    f[j] = mm_insert_epi8(f[j], row->f0[j], 0);
  }
  for (int i = 0; i < state->W; i++) {
    h = mm_load((__mxxxi *)seq + pc2 + i);
    mask = mm_cmpeq_epi8(h, base);
    h = mm_blendv_epi8(mis, mat, mask);
    mm_store((__mxxxi *)h_g + i, h);
    s1 = Smin;
    for (int j = 0; j < pre_num; j++) {
      t[j] = mm_load((__mxxxi *)tsta_graph_node(graph, row->prev[j])->sorce +
                     pd[j] + i);
      e[j] = mm_load((__mxxxi *)tsta_graph_node(graph, row->prev[j])->esorce +
                     pd[j] + i);
      temp = mm_max_epi8(f[j], h);
      temp = mm_max_epi8(e[j], temp);
      temp = mm_subs_epi8(temp, v[j]);
      mask4 = mm_cmpgt_epi8(v[j], z);
      temp = mm_blendv_epi8(temp, ogap, mask4);
      s1 = mm_max_epi8(s1, temp);
    }
    for (int j = 0; j < pre_num; j++) {
      temp = mm_sub_epi8(t[j], egap);
      temp = mm_subs_epi8(f[j], temp);
      temp1 = mm_adds_epi8(s1, ogap);
      temp1 = mm_subs_epi8(temp1, t[j]);
      temp1 = mm_adds_epi8(v[j], temp1);
      f[j] = mm_max_epi8(temp, temp1);

      temp1 = mm_subs_epi8(s1, t[j]);
      vc0[j] = mm_adds_epi8(v[j], temp1);

      mask4 = mm_cmpgt_epi8(temp1, zero);
      temp1 = mm_blendv_epi8(zero, temp1, mask4);
      temp2 = mm_subs_epi8(top, v[j]);
      y[j] = mm_subs_epu8(temp1, temp2);

      v[j] = mm_adds_epi8(vc0[j], vc1[j]);

      mask5 = mm_cmpeq_epi8(vc1[j], zero);
      temp2 = mm_subs_epu8(top, vc0[j]);
      diff[j] = mm_blendv_epi8(temp2, zero, mask5);

      temp2 = vc1[j];
      vc1[j] = mm_subs_epu8(vc1[j], diff[j]);
      vc1[j] = mm_adds_epi8(vc1[j], vc2[j]);
      temp2 = mm_subs_epu8(vc1[j], temp2);
      vc2[j] = mm_subs_epu8(vc2[j], diff[j]);
      vc2[j] = mm_adds_epi8(vc2[j], y[j]);
      vc2[j] = mm_subs_epu8(vc2[j], temp2);
    }
  }

  for (int j = 0; j < pre_num; j++) {
    mm_store((__mxxxi *)f_temp[j], f[j]);
    te[j] = f_temp[j][0];
    for (int x = 1; x < state->B - 1; x++) {
      te[j] = te[j] - r_s[j][x] + state->W * state->E;
      if (te[j] > f_temp[j][x] && te[j] > 125) {
        f_temp[j][x] = 125;
      } else if (te[j] > f_temp[j][x] && te[j] <= 125) {
        f_temp[j][x] = te[j];
      } else if (te[j] <= f_temp[j][x] && f_temp[j][x] > 125) {
        te[j] = f_temp[j][x];
        f_temp[j][x] = 125;
      } else {
        te[j] = f_temp[j][x];
      }
    }
    f[j] = mm_load((__mxxxi *)f_temp[j]);
    temp1 = mm_subs_epi8(f[j], egap);
    f[j] = mm_slli(f[j]);
    f[j] = mm_insert_epi8(f[j], row->f0[j], 0);

    vc0[j] = mm_max_epi8(temp1, v[j]);
    vc0[j] = mm_slli(vc0[j]);
    vc0[j] = mm_insert_epi8(vc0[j], v0[j], 0);

    vc1[j] = mm_slli(vc1[j]);
    vc1[j] = mm_insert_epi8(vc1[j], vc_1[j], 0);
    v[j] = mm_adds_epi8(vc0[j], vc1[j]);

    vc2[j] = mm_slli(vc2[j]);
    vc2[j] = mm_insert_epi8(vc2[j], vc_2[j], 0);
  }

  sum[0] = zero;
  s0 = mm_set1_epi8(42);
  s2 = mm_add_epi8(s0, s0);
  s3 = mm_add_epi8(s2, s0);
  s4 = mm_set1_epi8(1);
  s5 = mm_add_epi8(s0, s4);
  s6 = mm_add_epi8(s4, s4);

  for (int i = 0; i < state->W; i++) {
    h = mm_load((__mxxxi *)seq + pc2 + i);
    SN = mm_cmpeq_epi8(h, N);
    h = mm_load((__mxxxi *)h_g + i);
    SM = mm_cmpeq_epi8(mat, h);
    SX = mm_cmpeq_epi8(mis, h);
    max = eumax = Smin;
    for (int j = 0; j < pre_num; j++) {
      t[j] = mm_load((__mxxxi *)tsta_graph_node(graph, row->prev[j])->sorce +
                     pd[j] + i);
      e[j] = mm_load((__mxxxi *)tsta_graph_node(graph, row->prev[j])->esorce +
                     pd[j] + i);
      fv[j] = mm_subs_epi8(f[j], v[j]);
      eu[j] = mm_subs_epi8(e[j], v[j]);
      q[j] = mm_subs_epi8(h, v[j]);
      temp = mm_max_epi8(fv[j], eu[j]);
      temp = mm_max_epi8(temp, q[j]);
      mask4 = mm_cmpgt_epi8(v[j], z);
      temp = mm_blendv_epi8(temp, ogap, mask4);
      max = mm_max_epi8(max, temp);
      ev[j] = mm_subs_epi8(e[j], t[j]);
      eumax = mm_max_epi8(eumax, eu[j]);
    }
    max = mm_blendv_epi8(max, zero, SN);
    sum[0] = mm_add_epi8(sum[0], max);

    // source
    source = s3;
    source_num = zero;
    for (int j = pre_num - 1; j >= 0; j--) {
      mask = mm_cmpeq_epi8(max, eu[j]);
      source = mm_blendv_epi8(source, zero, mask);
      source_num = mm_blendv_epi8(source_num, mm_set1_epi8(j), mask);
    }
    for (int j = pre_num - 1; j >= 0; j--) {
      mask = mm_and_epi8(mm_cmpeq_epi8(max, q[j]), SX);
      source = mm_blendv_epi8(source, s2, mask);
      source_num = mm_blendv_epi8(source_num, mm_set1_epi8(j), mask);
    }
    for (int j = pre_num - 1; j >= 0; j--) {
      mask = mm_and_epi8(mm_cmpeq_epi8(max, q[j]), SM);
      source = mm_blendv_epi8(source, s0, mask);
      source_num = mm_blendv_epi8(source_num, mm_set1_epi8(j), mask);
    }
    source = mm_add_epi8(source, source_num); // now:0/42/84/126
    mm_store((__mxxxi *)source_chunk + i, source);
    mm_store((__mxxxi *)row->sorce + pc1 + i, max);

    // esource+fsource
    esource = fsource = s4;
    esource_num = zero;
    temp = mm_adds_epi8(max, ogap);
    emax = Smin;
    for (int j = pre_num - 1; j >= 0; j--) {
      f[j] = mm_adds_epi8(f[j], egap);
      s1 = mm_adds_epi8(temp, v[j]);
      mask1 = mm_cmpeq_epi8(f[j], s1);
      f[j] = mm_max_epi8(f[j], s1);
      f[j] = mm_subs_epi8(f[j], t[j]);
      mask = mm_cmpeq_epi8(fv[j], ogap);
      fsource = mm_blendv_epi8(fsource, s6, mask);

      e[j] = mm_adds_epi8(e[j], egap);
      e[j] = mm_subs_epi8(e[j], v[j]);
      mask2 = mm_cmpeq_epi8(temp, e[j]);
      temp1 = mm_max_epi8(temp, e[j]);
      emax = mm_max_epi8(emax, temp1);

      mask3 = mm_cmpeq_epi8(eu[j], eumax);
      esource_num = mm_blendv_epi8(esource_num, mm_set1_epi8(j), mask3);
      mask = mm_cmpeq_epi8(ev[j], ogap);
      temp1 = mm_blendv_epi8(s4, s5, mask);
      temp1 = mm_add_epi8(temp1, esource_num);
      esource = mm_blendv_epi8(esource, temp1, mask3);
      temp1 = mm_sub_epi8(zero, esource);
      mask = mm_and_epi8(mask3, mask2);
      esource = mm_blendv_epi8(esource, temp1, mask);

      temp1 = mm_subs_epi8(max, t[j]);
      vc0[j] = mm_adds_epi8(v[j], temp1);

      mask4 = mm_cmpgt_epi8(temp1, zero);
      temp1 = mm_blendv_epi8(zero, temp1, mask4);
      temp2 = mm_subs_epi8(top, v[j]);
      y[j] = mm_subs_epu8(temp1, temp2);

      v[j] = mm_adds_epi8(vc0[j], vc1[j]);

      mask5 = mm_cmpeq_epi8(vc1[j], zero);
      temp2 = mm_subs_epu8(top, vc0[j]);
      diff[j] = mm_blendv_epi8(temp2, zero, mask5);

      temp2 = vc1[j];
      vc1[j] = mm_subs_epu8(vc1[j], diff[j]);
      vc1[j] = mm_adds_epi8(vc1[j], vc2[j]);
      temp2 = mm_subs_epu8(vc1[j], temp2);
      vc2[j] = mm_subs_epu8(vc2[j], diff[j]);
      vc2[j] = mm_adds_epi8(vc2[j], y[j]);
      vc2[j] = mm_subs_epu8(vc2[j], temp2);
      /*temp1 = mm_subs_epi8(max, t[j]);
      v[j] = mm_adds_epi8(v[j], temp1);*/
    }
    temp1 = mm_sub_epi8(zero, fsource);
    fsource = mm_blendv_epi8(fsource, temp1, mask1);
    mm_store((__mxxxi *)fsource_chunk + i, fsource);  //-
    mm_store((__mxxxi *)esource_chunk + i, esource);  //-
    mm_store((__mxxxi *)row->esorce + pc1 + i, emax); //-
  }
  tsta_trace_block_store_mark_dirty(row->source_store, trace_chunk_index);
  tsta_trace_block_store_mark_dirty(row->esource_store, trace_chunk_index);
  tsta_trace_block_store_mark_dirty(row->fsource_store, trace_chunk_index);
  for (int j = 0; j < pre_num; j++)
    row->f0[j] = mm_extract_epi8(f[j]);
  s1 = mm_reduceto16_epi8(sum[0]);
  s1 = mm_hadd_epi8(s1);
  reduce = mm_reduce_epi8(s1);
  row->simple_sorce[nv + 1] = row->simple_sorce[nv] + reduce;

  if (row->out == 0 && block_i >= state->maxtag && block_l == 0)
    row->lastsorce = row->simple_sorce[nv + 1];

  // cross-block
  int kk = (row->sub / state->L + 1) * state->L;
  for (int i = 0; i < row->out; i++) {
    if (tsta_graph_node(graph, row->next[i])->sub >= kk &&
        row->node_status != 3) {
      char *t_sorce = (char *)mm_malloc(state->length1 * sizeof(char)); //
      memcpy(t_sorce, row->sorce, state->L * sizeof(char));
      mm_free(row->sorce);
      row->sorce = t_sorce;

      char *t_esorce = (char *)mm_malloc(state->length1 * sizeof(char)); //
      memcpy(t_esorce, row->esorce, state->L * sizeof(char));
      mm_free(row->esorce);
      row->esorce = t_esorce;

      row->node_status = 3;
    }
  }
  tsta_aligned_free(h_g);
  tsta_aligned_free(f_temp);
  tsta_aligned_free(VC2);
  tsta_aligned_free(VC1);
  tsta_aligned_free(r_s);
  free(v0);
  free(vc_1);
  free(vc_2);
  free(pd);
  free(te);
}

void tsta_msa_block_alignment(void *pa) {
  msa_block_task *task = (msa_block_task *)pa;
  tsta_msa_state *state = task->state;
  int block_i = task->diagonal_index;
  int block_j = task->diagonal_block_count;
  int block_l = task->diagonal_block_offset;
  int num = task->sequence_index;
  char *seq = task->packed_sequence;
  tsta_graph_t *p = task->graph;

  int a1, a2;
  int nv = nconvert(block_i) - block_l;
  int pc2 = nv * state->L / state->B;

  if (num <= 0) {
    return;
  }

  a1 = (((block_i - state->maxtag) > 0) * (block_i - state->maxtag) + block_l) *
       state->L;
  for (int i = 0; i < state->L; i++) {
    a2 = a1 + i;
    if (a2 >= p->len)
      break;
    block_line_alignment(p, block_i, block_j, block_l,
                         tsta_graph_sort_node(p, a2), seq, nv, pc2, state);
  }
}
