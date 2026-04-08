#define _GNU_SOURCE

#include "tsta_common.h"

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct tsta_psa_aligner {
  tsta_config config;
  tsta_psa_state state;
  packed_sequence_array_t packed_sequences;
  char *packed_sequence_slots[2];
  char *cached_sorce;
  char *cached_esorce;
  int *cached_real;
  char *cached_V;
  char *cached_F;
  size_t cached_main_len;
  size_t cached_aux_len;
};

#include "tsta_threadpool.h"

#define MIN -100
#define I_MIN -2000000000
static inline int tsta_psa_trace_index(const tsta_psa_state *state, int j) {
  int lane_len = state->length[0];

  return (j / lane_len) * lane_len + ((((j % lane_len) % state->W) * state->B) +
                                      ((j % lane_len) / state->W));
}

typedef struct psa_block_task {
  int diagonal_index;
  int diagonal_block_offset;
  tsta_psa_state *state;
} psa_block_task;

typedef struct psa_thread_workspace {
  int lane_len;
  int block_width;
  int *maxsorce;
  char *h_s;
  char *rf;
  char *t_temp;
  char *e_temp;
  char *q_temp;
  int *r_temp;
  char *source[3];
} psa_thread_workspace;

static pthread_key_t psa_workspace_key;
static pthread_once_t psa_workspace_once = PTHREAD_ONCE_INIT;
static int psa_workspace_key_ready = 0;

static void psa_workspace_reset(psa_thread_workspace *ws) {
  if (!ws) {
    return;
  }

  free(ws->maxsorce);
  ws->maxsorce = NULL;
  mm_free(ws->h_s);
  ws->h_s = NULL;
  mm_free(ws->rf);
  ws->rf = NULL;
  mm_free(ws->t_temp);
  ws->t_temp = NULL;
  mm_free(ws->e_temp);
  ws->e_temp = NULL;
  mm_free(ws->q_temp);
  ws->q_temp = NULL;
  mm_free(ws->r_temp);
  ws->r_temp = NULL;
  for (int i = 0; i < 3; i++) {
    mm_free(ws->source[i]);
    ws->source[i] = NULL;
  }
  ws->lane_len = 0;
  ws->block_width = 0;
}

static void psa_workspace_destroy(void *ptr) {
  psa_thread_workspace *ws = (psa_thread_workspace *)ptr;

  psa_workspace_reset(ws);
  free(ws);
}

static void psa_workspace_make_key(void) {
  if (pthread_key_create(&psa_workspace_key, psa_workspace_destroy) == 0) {
    psa_workspace_key_ready = 1;
  }
}

static int psa_workspace_resize(psa_thread_workspace *ws, int lane_len,
                                int block_width) {
  if (!ws || lane_len <= 0 || block_width <= 0) {
    return -1;
  }

  if (ws->lane_len == lane_len && ws->block_width == block_width &&
      ws->maxsorce && ws->h_s && ws->rf && ws->t_temp && ws->e_temp &&
      ws->q_temp && ws->r_temp && ws->source[0] && ws->source[1] &&
      ws->source[2]) {
    return 0;
  }

  psa_workspace_reset(ws);
  ws->maxsorce = (int *)malloc((size_t)lane_len * sizeof(int));
  ws->h_s = (char *)mm_malloc((size_t)lane_len * sizeof(char));
  ws->rf = (char *)mm_malloc((size_t)block_width * sizeof(char));
  ws->t_temp = (char *)mm_malloc((size_t)lane_len * sizeof(char));
  ws->e_temp = (char *)mm_malloc((size_t)lane_len * sizeof(char));
  ws->q_temp = (char *)mm_malloc((size_t)lane_len * sizeof(char));
  ws->r_temp = (int *)mm_malloc((size_t)lane_len * sizeof(int));
  for (int i = 0; i < 3; i++) {
    ws->source[i] = (char *)mm_malloc((size_t)lane_len * sizeof(char));
  }

  if (!ws->maxsorce || !ws->h_s || !ws->rf || !ws->t_temp || !ws->e_temp ||
      !ws->q_temp || !ws->r_temp || !ws->source[0] || !ws->source[1] ||
      !ws->source[2]) {
    psa_workspace_reset(ws);
    return -1;
  }

  ws->lane_len = lane_len;
  ws->block_width = block_width;
  return 0;
}

static psa_thread_workspace *psa_workspace_acquire(int lane_len,
                                                   int block_width) {
  psa_thread_workspace *ws;

  pthread_once(&psa_workspace_once, psa_workspace_make_key);
  if (!psa_workspace_key_ready) {
    return NULL;
  }

  ws = (psa_thread_workspace *)pthread_getspecific(psa_workspace_key);
  if (!ws) {
    ws = (psa_thread_workspace *)calloc(1, sizeof(*ws));
    if (!ws) {
      return NULL;
    }
    if (pthread_setspecific(psa_workspace_key, ws) != 0) {
      free(ws);
      return NULL;
    }
  }

  if (psa_workspace_resize(ws, lane_len, block_width) != 0) {
    return NULL;
  }

  return ws;
}

static inline int mm128_max_reduce(__mxxxi a) {
  __mxxxi b, c;
  b = mm_shuffle_epi32(a, _MM_SHUFFLE(3, 3, 1, 1));
  b = mm_max_epi32(a, b);
  c = mm_shuffle_epi32(b, _MM_SHUFFLE(2, 2, 2, 2));
  c = mm_max_epi32(b, c);
  return mm_cvtsixxx_si32(c);
}

static inline int mm256_max_reduce(__mxxxi e) {
  __mxxxi f, g, h;
  f = mm_shuffle_epi32(e, _MM_SHUFFLE(3, 3, 1, 1));
  f = mm_max_epi32(e, f);
  g = mm_shuffle_epi32(f, _MM_SHUFFLE(2, 2, 2, 2));
  g = mm_max_epi32(f, g);
  h = mm_permute4x64_epi64(g, _MM_SHUFFLE(2, 2, 2, 2));
  h = mm_max_epi32(g, h);
  return mm_cvtsixxx_si32(h);
}

static inline void blockmatrix_init(tsta_psa_state *state) {
  memset(state->sorce, state->E, state->length[0]);
  memset(state->esorce, state->E + state->E + state->O, state->length[0]);
  state->sorce[0] = state->E + state->O;
  state->esorce[0] = 2 * (state->E + state->O);
  for (int i = 0; i < state->length[0]; i++)
    state->real[i] = state->O + ((i / state->L * state->L +
                                  ((i % state->L) % state->B) * state->W +
                                  ((i % state->L) / state->B)) +
                                 1) *
                                    state->E;

  memset(state->V, state->E, state->length[3]);
  memset(state->F, state->E + state->E + state->O, state->length[3]);
  state->V[0] = state->E + state->O;
  state->F[0] = 2 * (state->E + state->O);
}

static inline void row(int *maxsorce, int y, int block_i, int block_l, int pc1,
                       int pc2, int pc4, char *h_s, char *t_temp, char *e_temp,
                       char *q_temp, char *rf, int *r_temp, char **source,
                       tsta_psa_state *state) {
  __mxxxi zero, s1, s2, Smin, h, b2, e, ev, f, fv, t, s, v, mat, mis, egap,
      ogap, v1, h1, trace, etrace, ftrace, temp1;
  __mask mask, mask1;
  int j = 0;
  zero = mm_setzero();
  egap = mm_set1_epi8(state->E);
  ogap = mm_set1_epi8(state->O + state->E);
  mis = mm_set1_epi8(state->X);
  mat = mm_set1_epi8(state->M);
  Smin = mm_set1_epi8(MIN);
  f = Smin;
  f = mm_insert_epi8(f, state->F[pc2], 0);
  b2 = mm_set1_epi8(state->seq[1][pc2]);

  for (int x = 0; x < state->W; x++) {
    h = mm_load(((__mxxxi *)q_temp) + x);
    mask = mm_cmpeq_epi8(h, b2);
    h = mm_blendv_epi8(mis, mat, mask);
    mm_store(((__mxxxi *)h_s) + x, h);
    t = mm_load(((__mxxxi *)t_temp) + x);
    e = mm_load(((__mxxxi *)e_temp) + x);
    s = mm_max_epi8(h, e);
    s = mm_max_epi8(s, f);
    f = mm_add_epi8(f, egap);
    h1 = mm_add_epi8(s, ogap);
    f = mm_max_epi8(f, h1);
    f = mm_sub_epi8(f, t);
  }

  mm_store((__mxxxi *)rf, f);
  for (int x = 1; x < state->B; x++)
    if (rf[x - 1] + state->W * state->E -
            (r_temp[state->L - state->B + x] -
             r_temp[state->L - state->B + x - 1]) >
        rf[x])
      rf[x] = rf[x - 1] + state->W * state->E -
              (r_temp[state->L - state->B + x] -
               r_temp[state->L - state->B + x - 1]);

  f = mm_load((__mxxxi *)rf);
  temp1 = mm_sub_epi8(f, egap);
  f = mm_slli(f);
  f = mm_insert_epi8(f, state->F[pc2], 0);

  v = mm_sub_epi8(s, t);
  v = mm_max_epi8(temp1, v);
  v = mm_slli(v);
  v = mm_insert_epi8(v, state->V[pc2], 0);

  s1 = mm_set1_epi8(1);
  s2 = mm_add_epi8(s1, s1);
  b2 = mm_set1_epi32(I_MIN);
  for (int x = 0; x < state->W; x++) {
    h1 = mm_load(((__mxxxi *)h_s) + x);
    t = mm_load(((__mxxxi *)t_temp) + x);
    e = mm_load(((__mxxxi *)e_temp) + x);
    s = mm_max_epi8(e, f);
    s = mm_max_epi8(s, h1);
    h = mm_sub_epi8(s, v);
    mm_store(((__mxxxi *)t_temp) + x, h);

    trace = s2;
    mask = mm_cmpeq_epi8(s, f);
    trace = mm_blendv_epi8(trace, zero, mask);
    mask = mm_cmpeq_epi8(s, h1);
    trace = mm_blendv_epi8(trace, s1, mask);
    mm_store(((__mxxxi *)source[0]) + x, trace);

    h1 = mm_add_epi8(s, ogap);
    fv = mm_sub_epi8(f, v);
    f = mm_add_epi8(f, egap);
    mask1 = mm_cmpeq_epi8(f, h1);
    f = mm_max_epi8(f, h1);
    f = mm_sub_epi8(f, t);

    mask = mm_cmpeq_epi8(fv, ogap);
    ftrace = mm_blendv_epi8(s1, s2, mask);
    temp1 = mm_sub_epi8(zero, ftrace);
    mask = mm_and_epi8(mask, mask1);
    ftrace = mm_blendv_epi8(ftrace, temp1, mask);
    mm_store(((__mxxxi *)source[1]) + x, ftrace);

    ev = mm_sub_epi8(e, t);
    e = mm_add_epi8(e, egap);
    mask1 = mm_cmpeq_epi8(e, h1);
    e = mm_max_epi8(e, h1);
    e = mm_sub_epi8(e, v);
    mm_store(((__mxxxi *)e_temp) + x, e);

    mask = mm_cmpeq_epi8(ev, ogap);
    etrace = mm_blendv_epi8(s1, s2, mask);
    temp1 = mm_sub_epi8(zero, etrace);
    mask = mm_and_epi8(mask, mask1);
    etrace = mm_blendv_epi8(etrace, temp1, mask);
    mm_store(((__mxxxi *)source[2]) + x, etrace);

    v = mm_sub_epi8(s, t);
    v1 = mm0_epi8cvt32(v);
    h1 = mm_load(((__mxxxi *)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi *)r_temp) + j, h1);
    j++;
    v1 = mm_epi8cvt32(v, 1);
    h1 = mm_load(((__mxxxi *)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi *)r_temp) + j, h1);
    j++;
    v1 = mm_epi8cvt32(v, 2);
    h1 = mm_load(((__mxxxi *)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi *)r_temp) + j, h1);
    j++;
    v1 = mm_epi8cvt32(v, 3);
    h1 = mm_load(((__mxxxi *)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi *)r_temp) + j, h1);
    j++;
  }
  state->F[pc2] = mm_extract_epi8(f);
  state->V[pc2] = mm_extract_epi8(v);
  maxsorce[y] = mm_reduce_max_epi32(b2);
  if (state->trace_enabled && state->trace.matrices.back &&
      state->trace.matrices.fback && state->trace.matrices.eback) {
    memcpy(state->trace.matrices.back[pc2] + pc4, source[0], state->L);
    memcpy(state->trace.matrices.fback[pc2] + pc4, source[1], state->L);
    memcpy(state->trace.matrices.eback[pc2] + pc4, source[2], state->L);
  }
}

static inline void tsta_psa_block_alignment(void *p) {
  psa_block_task *task = (psa_block_task *)p;
  tsta_psa_state *state = task->state;
  int block_i = task->diagonal_index;
  int block_l = task->diagonal_block_offset;
  psa_thread_workspace *ws;
  int use_tls_workspace = 0;

  int *maxsorce = NULL;
  int pc0, pc1, pc2, pc4;
  char *h_s = NULL;
  char *rf = NULL;
  char *t_temp = NULL;
  char *e_temp = NULL;
  char *q_temp = NULL;
  int *r_temp = NULL;
  char *source[3] = {NULL, NULL, NULL};

  ws = psa_workspace_acquire(state->L, block);
  if (ws) {
    use_tls_workspace = 1;
    maxsorce = ws->maxsorce;
    h_s = ws->h_s;
    rf = ws->rf;
    t_temp = ws->t_temp;
    e_temp = ws->e_temp;
    q_temp = ws->q_temp;
    r_temp = ws->r_temp;
    source[0] = ws->source[0];
    source[1] = ws->source[1];
    source[2] = ws->source[2];
  } else {
    maxsorce = (int *)malloc(state->L * sizeof(int));
    h_s = (char *)mm_malloc(state->L * sizeof(char));
    rf = (char *)mm_malloc(block * sizeof(char));
    t_temp = (char *)mm_malloc(state->L * sizeof(char));
    e_temp = (char *)mm_malloc(state->L * sizeof(char));
    q_temp = (char *)mm_malloc(state->L * sizeof(char));
    r_temp = (int *)mm_malloc(state->L * sizeof(int));
    for (int i = 0; i < 3; i++) {
      source[i] = (char *)mm_malloc(state->L * sizeof(char));
    }
    if (!maxsorce || !h_s || !rf || !t_temp || !e_temp || !q_temp || !r_temp ||
        !source[0] || !source[1] || !source[2]) {
      goto cleanup;
    }
  }

  if (block_i <= state->lmaxtag)
    pc0 = block_i - block_l;
  else
    pc0 = state->lmaxtag - block_l;
  pc1 = pc0 * state->W;
  pc4 = pc0 * state->L;

  memcpy(t_temp, state->sorce + pc4, state->L);
  memcpy(e_temp, state->esorce + pc4, state->L);
  memcpy(r_temp, state->real + pc4, state->L * sizeof(int));
  for (int i = 0; i < state->L; i++)
    q_temp[i] = state->seq[0][pc4 + (i % state->B) * state->W + (i / state->B)];

  for (int i = 0; i < state->L; i++) {
    if (block_i <= state->lmaxtag)
      pc2 = block_l * state->L + i;
    else
      pc2 = (block_l + block_i - state->lmaxtag) * state->L + i;
    if (pc2 >= state->length[3]) {
      for (int s = i; s < state->L; s++)
        maxsorce[s] = I_MIN;
      break;
    }
    row(maxsorce, i, block_i, block_l, pc1, pc2, pc4, h_s, t_temp, e_temp,
        q_temp, rf, r_temp, source, state);
  }
  memcpy(state->sorce + pc4, t_temp, state->L);
  memcpy(state->esorce + pc4, e_temp, state->L);
  memcpy(state->real + pc4, r_temp, state->L * sizeof(int));

  for (int x = 1; x < state->L; x++)
    if (maxsorce[0] <= maxsorce[x])
      maxsorce[0] = maxsorce[x];

  pthread_mutex_lock(&state->mutex);
  if (state->ms <= maxsorce[0])
    state->ms = maxsorce[0];
  pthread_mutex_unlock(&state->mutex);

cleanup:
  if (!use_tls_workspace) {
    mm_free(t_temp);
    mm_free(e_temp);
    mm_free(r_temp);
    mm_free(q_temp);
    mm_free(h_s);
    mm_free(rf);
    free(maxsorce);
    for (int i = 0; i < 3; i++)
      mm_free(source[i]);
  }
}
static void tsta_psa_trace_result_reset(tsta_psa_trace_result *result) {
  if (!result) {
    return;
  }

  free(result->aligned_seq1);
  free(result->aligned_seq2);
  free(result->formatted_text);
  result->aligned_seq1 = NULL;
  result->aligned_seq2 = NULL;
  result->formatted_text = NULL;
  result->aligned_length = 0;
  result->formatted_length = 0;
}

static void tsta_psa_trace_context_reset(tsta_psa_trace_context *trace) {
  if (!trace) {
    return;
  }

  if (trace->matrices.back) {
    for (int i = 0; i < trace->matrices.rows; i++) {
      mm_free(trace->matrices.back[i]);
    }
  }
  if (trace->matrices.eback) {
    for (int i = 0; i < trace->matrices.rows; i++) {
      mm_free(trace->matrices.eback[i]);
    }
  }
  if (trace->matrices.fback) {
    for (int i = 0; i < trace->matrices.rows; i++) {
      mm_free(trace->matrices.fback[i]);
    }
  }

  free(trace->matrices.back);
  free(trace->matrices.eback);
  free(trace->matrices.fback);
  trace->matrices.back = NULL;
  trace->matrices.eback = NULL;
  trace->matrices.fback = NULL;
  trace->matrices.rows = 0;
  trace->matrices.cols = 0;

  tsta_psa_trace_result_reset(&trace->result);
}

static int tsta_psa_trace_context_alloc(tsta_psa_state *state) {
  tsta_psa_trace_matrices *matrices = &state->trace.matrices;

  matrices->rows = state->length[3];
  matrices->cols = state->length[0];
  matrices->back = (char **)calloc((size_t)matrices->rows, sizeof(char *));
  matrices->eback = (char **)calloc((size_t)matrices->rows, sizeof(char *));
  matrices->fback = (char **)calloc((size_t)matrices->rows, sizeof(char *));
  if (!matrices->back || !matrices->eback || !matrices->fback) {
    tsta_psa_trace_context_reset(&state->trace);
    return -1;
  }

  for (int i = 0; i < matrices->rows; i++) {
    matrices->back[i] =
        (char *)mm_malloc((size_t)matrices->cols * sizeof(char));
    matrices->eback[i] =
        (char *)mm_malloc((size_t)matrices->cols * sizeof(char));
    matrices->fback[i] =
        (char *)mm_malloc((size_t)matrices->cols * sizeof(char));
    if (!matrices->back[i] || !matrices->eback[i] || !matrices->fback[i]) {
      tsta_psa_trace_context_reset(&state->trace);
      return -1;
    }
  }

  return 0;
}

static int tsta_psa_trace_build_alignment(tsta_psa_state *state) {
  tsta_psa_trace_result *result = &state->trace.result;
  int i = state->length[3] - 1;
  int j = state->length[2] - 1;
  int insertion_count = 0;
  size_t write_pos;

  tsta_psa_trace_result_reset(result);

  while (i >= 0 && j >= 0) {
    int current = tsta_psa_trace_index(state, j);

    if (state->trace.matrices.back[i][current] == 1) {
      i--;
      j--;
    } else if (state->trace.matrices.back[i][current] == 0) {
      if (j - 1 >= 0 &&
          ((state->trace.matrices.fback[i][current] == 1 ||
            state->trace.matrices.fback[i][current] == -1) ||
           ((state->trace.matrices.fback[i][current] == 2 ||
             state->trace.matrices.fback[i][current] == -2) &&
            state->trace.matrices.fback[i][tsta_psa_trace_index(state, j - 1)] <
                0))) {
        state->trace.matrices.back[i][tsta_psa_trace_index(state, j - 1)] = 0;
      }
      j--;
    } else {
      if (i - 1 >= 0 && ((state->trace.matrices.eback[i][current] == 1 ||
                          state->trace.matrices.eback[i][current] == -1) ||
                         ((state->trace.matrices.eback[i][current] == 2 ||
                           state->trace.matrices.eback[i][current] == -2) &&
                          state->trace.matrices.eback[i - 1][current] < 0))) {
        state->trace.matrices.back[i - 1][current] = 2;
      }
      i--;
      insertion_count++;
    }
  }

  if (i >= 0) {
    insertion_count += i + 1;
  }

  result->aligned_length = (size_t)state->length[2] + (size_t)insertion_count;
  result->aligned_seq1 = (char *)malloc(result->aligned_length + 1);
  result->aligned_seq2 = (char *)malloc(result->aligned_length + 1);
  if (!result->aligned_seq1 || !result->aligned_seq2) {
    tsta_psa_trace_result_reset(result);
    return -1;
  }

  i = state->length[3] - 1;
  j = state->length[2] - 1;
  write_pos = result->aligned_length;
  while (i >= 0 && j >= 0) {
    write_pos--;
    int current = tsta_psa_trace_index(state, j);

    if (state->trace.matrices.back[i][current] == 1) {
      result->aligned_seq1[write_pos] = state->seq[0][j];
      result->aligned_seq2[write_pos] = state->seq[1][i];
      i--;
      j--;
    } else if (state->trace.matrices.back[i][current] == 0) {
      result->aligned_seq1[write_pos] = state->seq[0][j];
      result->aligned_seq2[write_pos] = '-';
      j--;
    } else {
      result->aligned_seq1[write_pos] = '-';
      result->aligned_seq2[write_pos] = state->seq[1][i];
      i--;
    }
  }

  while (j >= 0) {
    write_pos--;
    result->aligned_seq1[write_pos] = state->seq[0][j];
    result->aligned_seq2[write_pos] = '-';
    j--;
  }

  while (i >= 0) {
    write_pos--;
    result->aligned_seq2[write_pos] = state->seq[1][i];
    result->aligned_seq1[write_pos] = '-';
    i--;
  }

  result->aligned_seq1[result->aligned_length] = '\0';
  result->aligned_seq2[result->aligned_length] = '\0';
  return 0;
}

static int tsta_psa_trace_build_text(tsta_psa_state *state) {
  tsta_psa_trace_result *result = &state->trace.result;
  char *cursor;
  size_t total_len;

  if (!result->aligned_seq1 || !result->aligned_seq2) {
    return -1;
  }

  total_len = 3 + result->aligned_length + 4 + result->aligned_length;
  result->formatted_text = (char *)malloc(total_len + 1);
  if (!result->formatted_text) {
    return -1;
  }

  cursor = result->formatted_text;
  memcpy(cursor, ">1\n", 3);
  cursor += 3;
  memcpy(cursor, result->aligned_seq1, result->aligned_length);
  cursor += result->aligned_length;
  memcpy(cursor, "\n>2\n", 4);
  cursor += 4;
  memcpy(cursor, result->aligned_seq2, result->aligned_length);
  cursor += result->aligned_length;
  *cursor = '\0';

  result->formatted_length = total_len;
  return 0;
}

static void tsta_psa_apply_config(tsta_psa_aligner *aligner,
                                  const tsta_config *config) {
  if (config) {
    aligner->config = *config;
  } else {
    tsta_config_default(&aligner->config);
  }

  tsta_init_psa_state(&aligner->state, &aligner->config, block);
}

static int tsta_sequence_view_make(tsta_sequence_view_t *view,
                                   const char *sequence, int length) {
  if (!view || !sequence || length < 0) {
    return -1;
  }

  view->sequence = sequence;
  view->length = length;
  return 0;
}

static char *tsta_psa_make_padded_copy(const char *sequence_data,
                                       size_t sequence_length,
                                       int padded_length) {
  char *copy;
  size_t length = sequence_length;

  copy = (char *)malloc((size_t)padded_length + 1);
  if (!copy) {
    return NULL;
  }

  memset(copy, 'N', (size_t)padded_length);
  if (length > 0) {
    memcpy(copy, sequence_data, length);
  }
  copy[padded_length] = '\0';
  return copy;
}

static int tsta_psa_ensure_cached_sequence(tsta_psa_aligner *aligner, int slot,
                                           int padded_length,
                                           const char *sequence_data,
                                           size_t sequence_length) {
  packed_sequence_t *slot_entry;
  char *buffer;
  int required_capacity;

  if (!aligner || slot < 0 || slot >= 2 || !sequence_data ||
      padded_length < 0) {
    return -1;
  }

  if (tsta_packed_sequence_array_reserve_slot(&aligner->packed_sequences,
                                              (size_t)slot, 2) != 0) {
    return -1;
  }
  slot_entry = &aligner->packed_sequences.data[slot];

  required_capacity = padded_length + 1;
  if (slot_entry->len < required_capacity) {
    buffer = tsta_psa_make_padded_copy(sequence_data, sequence_length,
                                       padded_length);
    if (!buffer) {
      return -1;
    }
    free(slot_entry->packed_seq);
    slot_entry->packed_seq = buffer;
    slot_entry->len = required_capacity;
    return 0;
  }

  buffer = slot_entry->packed_seq;
  memset(buffer, 'N', (size_t)padded_length);
  if (sequence_length > 0) {
    memcpy(buffer, sequence_data, sequence_length);
  }
  buffer[padded_length] = '\0';
  return 0;
}

static int
tsta_psa_prepare_sequences(const char *sequence1_data, size_t sequence1_length,
                           const char *sequence2_data, size_t sequence2_length,
                           tsta_psa_aligner *aligner, tsta_psa_state *state) {
  const char *longer_data = sequence1_data;
  const char *shorter_data = sequence2_data;
  size_t longer_length = sequence1_length;
  size_t shorter_length = sequence2_length;

  if (sequence2_length > sequence1_length) {
    longer_data = sequence2_data;
    shorter_data = sequence1_data;
    longer_length = sequence2_length;
    shorter_length = sequence1_length;
  }

  state->length[2] = (int)longer_length;
  state->length[3] = (int)shorter_length;
  state->length[0] = (int)longer_length;
  state->length[1] = (int)shorter_length;
  if (state->length[0] % state->L != 0) {
    state->length[0] += state->L - state->length[0] % state->L;
  }
  if (state->length[1] % state->L != 0) {
    state->length[1] += state->L - state->length[1] % state->L;
  }

  if (tsta_psa_ensure_cached_sequence(aligner, 0, state->length[0], longer_data,
                                      longer_length) != 0) {
    return -1;
  }
  if (tsta_psa_ensure_cached_sequence(aligner, 1, state->length[1],
                                      shorter_data, shorter_length) != 0) {
    return -1;
  }

  aligner->packed_sequence_slots[0] =
      aligner->packed_sequences.data[0].packed_seq;
  aligner->packed_sequence_slots[1] =
      aligner->packed_sequences.data[1].packed_seq;
  state->seq = aligner->packed_sequence_slots;

  return 0;
}

static void tsta_psa_detach_runtime_buffers(tsta_psa_aligner *aligner) {
  if (!aligner) {
    return;
  }

  aligner->state.seq = NULL;
  aligner->state.sorce = NULL;
  aligner->state.esorce = NULL;
  aligner->state.real = NULL;
  aligner->state.F = NULL;
  aligner->state.V = NULL;
}

static void tsta_psa_free_runtime_buffers(tsta_psa_aligner *aligner,
                                          int release_cache) {
  if (!aligner) {
    return;
  }

  tsta_psa_detach_runtime_buffers(aligner);
  tsta_psa_trace_context_reset(&aligner->state.trace);

  if (!release_cache) {
    return;
  }

  free(aligner->cached_sorce);
  free(aligner->cached_esorce);
  free(aligner->cached_real);
  free(aligner->cached_F);
  free(aligner->cached_V);
  aligner->cached_sorce = NULL;
  aligner->cached_esorce = NULL;
  aligner->cached_real = NULL;
  aligner->cached_F = NULL;
  aligner->cached_V = NULL;
  aligner->cached_main_len = 0;
  aligner->cached_aux_len = 0;
}

static int tsta_psa_ensure_runtime_buffers(tsta_psa_aligner *aligner,
                                           size_t main_len, size_t aux_len) {
  char *new_sorce;
  char *new_esorce;
  int *new_real;
  char *new_F;
  char *new_V;
  char *old_sorce;
  char *old_esorce;
  int *old_real;
  char *old_F;
  char *old_V;

  if (!aligner || main_len == 0 || aux_len == 0) {
    return -1;
  }

  if (aligner->cached_main_len >= main_len &&
      aligner->cached_aux_len >= aux_len && aligner->cached_sorce &&
      aligner->cached_esorce && aligner->cached_real && aligner->cached_F &&
      aligner->cached_V) {
    aligner->state.sorce = aligner->cached_sorce;
    aligner->state.esorce = aligner->cached_esorce;
    aligner->state.real = aligner->cached_real;
    aligner->state.F = aligner->cached_F;
    aligner->state.V = aligner->cached_V;
    return 0;
  }

  new_sorce = (char *)malloc(main_len * sizeof(char));
  new_esorce = (char *)malloc(main_len * sizeof(char));
  new_real = (int *)malloc(main_len * sizeof(int));
  new_F = (char *)malloc(aux_len * sizeof(char));
  new_V = (char *)malloc(aux_len * sizeof(char));
  if (!new_sorce || !new_esorce || !new_real || !new_F || !new_V) {
    free(new_sorce);
    free(new_esorce);
    free(new_real);
    free(new_F);
    free(new_V);
    return -1;
  }

  old_sorce = aligner->cached_sorce;
  old_esorce = aligner->cached_esorce;
  old_real = aligner->cached_real;
  old_F = aligner->cached_F;
  old_V = aligner->cached_V;

  aligner->cached_sorce = new_sorce;
  aligner->cached_esorce = new_esorce;
  aligner->cached_real = new_real;
  aligner->cached_F = new_F;
  aligner->cached_V = new_V;
  aligner->cached_main_len = main_len;
  aligner->cached_aux_len = aux_len;

  free(old_sorce);
  free(old_esorce);
  free(old_real);
  free(old_F);
  free(old_V);

  aligner->state.sorce = aligner->cached_sorce;
  aligner->state.esorce = aligner->cached_esorce;
  aligner->state.real = aligner->cached_real;
  aligner->state.F = aligner->cached_F;
  aligner->state.V = aligner->cached_V;
  return 0;
}

static void tsta_psa_free_cached_sequences(tsta_psa_aligner *aligner) {
  if (!aligner) {
    return;
  }

  tsta_packed_sequence_array_reset_malloc(&aligner->packed_sequences);

  for (int i = 0; i < 2; i++) {
    aligner->packed_sequence_slots[i] = NULL;
  }
}

static int tsta_psa_build_cigar(const tsta_psa_trace_result *trace,
                                tsta_psa_result_t *result) {
  size_t i;
  size_t write_pos = 0;
  char prev_op = '\0';
  unsigned int run_len = 0;
  size_t capacity;
  char *cigar;

  if (!trace || !result || !trace->aligned_seq1 || !trace->aligned_seq2) {
    return -1;
  }

  capacity = trace->aligned_length > 0 ? trace->aligned_length * 12 + 1 : 2;
  cigar = (char *)malloc(capacity);
  if (!cigar) {
    return -1;
  }

  for (i = 0; i < trace->aligned_length; i++) {
    char a = trace->aligned_seq1[i];
    char b = trace->aligned_seq2[i];
    char op;

    if (a == '-' && b != '-') {
      result->ins_count++;
      op = 'I';
    } else if (a != '-' && b == '-') {
      result->del_count++;
      op = 'D';
    } else if (a == b) {
      result->match_count++;
      op = 'M';
    } else {
      result->mismatch_count++;
      op = 'X';
    }

    if (prev_op == '\0') {
      prev_op = op;
      run_len = 1;
      continue;
    }

    if (op == prev_op) {
      run_len++;
      continue;
    }

    write_pos += (size_t)snprintf(cigar + write_pos, capacity - write_pos,
                                  "%u%c", run_len, prev_op);
    prev_op = op;
    run_len = 1;
  }

  if (prev_op != '\0') {
    write_pos += (size_t)snprintf(cigar + write_pos, capacity - write_pos,
                                  "%u%c", run_len, prev_op);
  } else {
    cigar[0] = '0';
    cigar[1] = 'M';
    write_pos = 2;
  }

  cigar[write_pos] = '\0';
  result->cigar = cigar;
  return 0;
}

static int tsta_psa_capture_result(tsta_psa_result_t *result,
                                   tsta_psa_state *state) {
  tsta_psa_trace_result *trace;

  if (!result || !state) {
    return -1;
  }

  if (tsta_psa_trace_build_alignment(state) != 0) {
    return -1;
  }

  trace = &state->trace.result;
  tsta_psa_result_free(result);
  result->score = state->ms;
  result->aln_length = trace->aligned_length;
  result->aln[0] = (char *)malloc(trace->aligned_length + 1);
  result->aln[1] = (char *)malloc(trace->aligned_length + 1);
  if (!result->aln[0] || !result->aln[1]) {
    tsta_psa_result_free(result);
    return -1;
  }

  memcpy(result->aln[0], trace->aligned_seq1, trace->aligned_length + 1);
  memcpy(result->aln[1], trace->aligned_seq2, trace->aligned_length + 1);

  if (tsta_psa_build_cigar(trace, result) != 0) {
    tsta_psa_result_free(result);
    return -1;
  }

  return 0;
}

static int tsta_psa_aligner_align_internal(
    tsta_psa_aligner *aligner, const char *sequence1, int sequence1_length,
    const char *sequence2, int sequence2_length, tsta_psa_result_t *result) {
  tsta_threadpool_t *pool = NULL;
  unsigned int tsl;
  int j = 0;
  int status = -1;
  tsta_config effective;
  tsta_psa_state *state;
  int mutex_initialized = 0;
  tsta_sequence_view_t left;
  tsta_sequence_view_t right;

  if (!aligner) {
    return -1;
  }
  if (tsta_sequence_view_make(&left, sequence1, sequence1_length) != 0 ||
      tsta_sequence_view_make(&right, sequence2, sequence2_length) != 0) {
    return -1;
  }

  effective = aligner->config;
  state = &aligner->state;
  tsta_init_psa_state(state, &effective, block);
  if (result) {
    state->trace_enabled = 1;
  }

  state->ms = MIN;
  if (tsta_psa_prepare_sequences(left.sequence, (size_t)left.length,
                                 right.sequence, (size_t)right.length, aligner,
                                 state) != 0) {
    tsta_psa_free_runtime_buffers(aligner, 0);
    return -1;
  }

  if (pthread_mutex_init(&state->mutex, NULL) != 0) {
    tsta_psa_free_runtime_buffers(aligner, 0);
    return -1;
  }
  mutex_initialized = 1;

  tsl = (unsigned int)((state->length[0] + state->length[1]) / state->L - 1);
  state->fmaxtag = state->length[1] / state->L - 1;
  state->lmaxtag = state->length[0] / state->L - 1;

  if (tsta_psa_ensure_runtime_buffers(aligner, (size_t)state->length[0],
                                      (size_t)state->length[3]) != 0) {
    goto cleanup;
  }

  if (state->trace_enabled) {
    if (tsta_psa_trace_context_alloc(state) != 0) {
      goto cleanup;
    }
  }

  pool = tsta_threadpool_create(effective.threads > 0 ? effective.threads : 10,
                                100, sizeof(psa_block_task));
  if (!pool) {
    goto cleanup;
  }

  blockmatrix_init(state);
  j = 0;
  for (unsigned int i = 0; i < tsl; i++) {
    if ((int)i <= state->fmaxtag) {
      j++;
    } else if ((int)i > state->lmaxtag) {
      j--;
    }
    for (int l = 0; l < j; l++) {
      psa_block_task task = {
          .diagonal_index = (int)i,
          .diagonal_block_offset = l,
          .state = state,
      };
      if (tsta_threadpool_submit(pool, tsta_psa_block_alignment, &task,
                                 sizeof(task)) != 0) {
        tsta_psa_block_alignment(&task);
      }
    }
    tsta_threadpool_wait(pool);
  }

  status = 0;
  if (state->trace_enabled && tsta_psa_capture_result(result, state) != 0) {
    status = -1;
  } else if (!state->trace_enabled) {
    status = -1;
  }

cleanup:
  if (pool) {
    tsta_threadpool_destroy(pool);
  }
  tsta_psa_free_runtime_buffers(aligner, 0);
  if (mutex_initialized) {
    pthread_mutex_destroy(&state->mutex);
  }
  return status;
}

tsta_psa_aligner *tsta_psa_aligner_create(const tsta_config *config) {
  tsta_psa_aligner *aligner = (tsta_psa_aligner *)calloc(1, sizeof(*aligner));

  if (!aligner) {
    return NULL;
  }

  tsta_psa_apply_config(aligner, config);
  return aligner;
}

void tsta_psa_aligner_destroy(tsta_psa_aligner *aligner) {
  if (!aligner) {
    return;
  }
  tsta_psa_free_runtime_buffers(aligner, 1);
  tsta_psa_free_cached_sequences(aligner);
  free(aligner);
}

int tsta_psa_aligner_align(tsta_psa_aligner *aligner, const char *sequence1,
                           int sequence1_length, const char *sequence2,
                           int sequence2_length, tsta_psa_result_t *result) {
  if (!result) {
    return -1;
  }

  return tsta_psa_aligner_align_internal(aligner, sequence1, sequence1_length,
                                         sequence2, sequence2_length, result);
}

int tsta_psa_align(const char *sequence1, int sequence1_length,
                   const char *sequence2, int sequence2_length,
                   const tsta_config *config, tsta_psa_result_t *result) {
  tsta_config effective_config;
  size_t longest_sequence_length;
  tsta_psa_aligner *aligner;

  effective_config = config ? *config : tsta_config_make_default();
  longest_sequence_length =
      (size_t)(sequence1_length > sequence2_length ? sequence1_length
                                                   : sequence2_length);
  tsta_config_apply_length_fallback(&effective_config, longest_sequence_length,
                                    2);

  aligner = tsta_psa_aligner_create(&effective_config);
  if (!aligner) {
    return -1;
  }

  int status = tsta_psa_aligner_align(aligner, sequence1, sequence1_length,
                                      sequence2, sequence2_length, result);
  tsta_psa_aligner_destroy(aligner);
  return status;
}
