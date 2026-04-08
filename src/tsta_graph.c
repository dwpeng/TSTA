#define _GNU_SOURCE

#include "tsta_graph.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TSTA_GRAPH_NODE(graph, index) (&(graph)->nodes.data[(size_t)(index)])
#define TSTA_GRAPH_SORT_INDEX(graph, index)                                    \
  ((graph)->sort.data[(size_t)(index)])
#define TSTA_GRAPH_UNSORT_INDEX(graph, index)                                  \
  ((graph)->unsort.data[(size_t)(index)])
#define TSTA_GRAPH_SORT_NODE(graph, index)                                     \
  TSTA_GRAPH_NODE((graph), TSTA_GRAPH_SORT_INDEX((graph), (index)))
#define TSTA_GRAPH_UNSORT_NODE(graph, index)                                   \
  TSTA_GRAPH_NODE((graph), TSTA_GRAPH_UNSORT_INDEX((graph), (index)))

#define TSTA_GRAPH_INVALID_INDEX UINT32_MAX

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

#define TSTA_TRACE_PARENT_SLOT(trace)                                          \
  ((size_t)((trace) < 0 ? 0 : ((trace) % 42)))

#define TSTA_MSA_PASSING_SEQ_SET(node, index)                                  \
  do {                                                                         \
    (node)->passing_seq[(size_t)(index) >> 3] |=                               \
        (unsigned char)(1u << ((size_t)(index) & 7u));                         \
  } while (0)

size_t tsta_msa_passing_seq_bytes(size_t sequence_count) {
  return (sequence_count + 7) / 8;
}

static int tsta_msa_mismatch_list_reserve(tsta_node_t *node,
                                          size_t mismatch_num) {
  tsta_node_index_t *new_list;

  if (!node) {
    return -1;
  }

  if (mismatch_num == 0) {
    free(node->mismatch_node);
    node->mismatch_node = NULL;
    return 0;
  }

  new_list = (tsta_node_index_t *)realloc(node->mismatch_node,
                                          mismatch_num * sizeof(*new_list));
  if (!new_list) {
    return -1;
  }

  node->mismatch_node = new_list;
  return 0;
}

static int tsta_graph_activate_sequence_node(tsta_node_t *node,
                                             size_t sequence_count,
                                             int sequence_index,
                                             tsta_msa_state *state,
                                             int edge_weight_value) {
  size_t passing_bytes;

  if (!node || !state) {
    return -1;
  }

  passing_bytes = tsta_msa_passing_seq_bytes(sequence_count);
  if (!node->sorce) {
    node->sorce = (char *)mm_malloc(state->L * sizeof(char));
  }
  if (!node->esorce) {
    node->esorce = (char *)mm_malloc(state->L * sizeof(char));
  }
  if (!node->passing_seq) {
    node->passing_seq =
        (unsigned char *)calloc(passing_bytes, sizeof(unsigned char));
  }
  if (!node->edge_weight) {
    node->edge_weight = (int *)malloc(sizeof(int));
  }
  if (!node->f0) {
    node->f0 = (char *)malloc(sizeof(char));
  }
  if (!node->sorce || !node->esorce || !node->passing_seq ||
      !node->edge_weight || !node->f0) {
    return -1;
  }

  node->edge_weight[0] = edge_weight_value;
  TSTA_MSA_PASSING_SEQ_SET(node, (size_t)sequence_index);
  return 0;
}

void free_node(tsta_node_t *n) {
  if (!n) {
    return;
  }

  if (n->prev) {
    free(n->prev);
    n->prev = NULL;
  }
  if (n->next) {
    free(n->next);
    n->next = NULL;
  }
  mm_free(n->sorce);
  n->sorce = NULL;
  mm_free(n->esorce);
  n->esorce = NULL;
  if (n->source_store) {
    tsta_trace_block_store_destroy(n->source_store);
    n->source_store = NULL;
  }
  if (n->esource_store) {
    tsta_trace_block_store_destroy(n->esource_store);
    n->esource_store = NULL;
  }
  if (n->fsource_store) {
    tsta_trace_block_store_destroy(n->fsource_store);
    n->fsource_store = NULL;
  }
  free(n->passing_seq);
  n->passing_seq = NULL;
  free(n->mismatch_node);
  n->mismatch_node = NULL;
  free(n->simple_sorce);
  n->simple_sorce = NULL;
  free(n->f0);
  n->f0 = NULL;
  free(n->edge_weight);
  n->edge_weight = NULL;
}

void tsta_msa_free_graph(tsta_graph_t *graph) {
  if (!graph) {
    return;
  }

  if (graph->nodes.data) {
    for (size_t i = 0; i < graph->nodes.length; i++) {
      free_node(&graph->nodes.data[i]);
    }
  }

  free(graph->sort.data);
  free(graph->unsort.data);
  free(graph->nodes.data);

  graph->p = TSTA_GRAPH_INVALID_INDEX;
  graph->seed = TSTA_GRAPH_INVALID_INDEX;
  graph->nodes.data = NULL;
  graph->nodes.length = 0;
  graph->nodes.capacity = 0;
  graph->unsort.data = NULL;
  graph->unsort.length = 0;
  graph->unsort.capacity = 0;
  graph->sort.data = NULL;
  graph->sort.length = 0;
  graph->sort.capacity = 0;
  graph->len = 0;
  graph->last_node_num = 0;
  graph->sequence_count = 0;
  graph->passing_seq_bytes = 0;
  free(graph);
}

int tsta_msa_expand_passing_seq(tsta_graph_t *graph, size_t old_count,
                                size_t new_count) {
  if (!graph || new_count <= old_count) {
    return -1;
  }

  size_t old_bytes = tsta_msa_passing_seq_bytes(old_count);
  size_t new_bytes = tsta_msa_passing_seq_bytes(new_count);
  size_t node_count = graph->nodes.length;

  for (size_t i = 0; i < node_count; i++) {
    tsta_node_t *node = &graph->nodes.data[i];
    unsigned char *new_buf;

    if (!node || !node->passing_seq) {
      continue;
    }

    new_buf = (unsigned char *)realloc(node->passing_seq, new_bytes);
    if (!new_buf) {
      return -1;
    }
    if (new_bytes > old_bytes) {
      memset(new_buf + old_bytes, 0, new_bytes - old_bytes);
    }
    node->passing_seq = new_buf;
  }

  graph->sequence_count = new_count;
  graph->passing_seq_bytes = new_bytes;
  return 0;
}

tsta_node_t *tsta_graph_build_initial(tsta_graph_t *p,
                                      const tsta_sequence_view_t *sequence,
                                      int sum, tsta_msa_state *state) {
  size_t len_a;
  const char *a;
  tsta_node_index_t head_index;
  tsta_node_index_t seed_index;
  tsta_node_index_t end_index;

  if (!p || !sequence || !sequence->sequence || sequence->length <= 0) {
    return NULL;
  }

  memset(p, 0, sizeof(*p));

  len_a = (size_t)sequence->length;
  a = sequence->sequence;
  seed_index = 0;
  head_index = 1;
  end_index = 1;

  p->len = (int)len_a;
  p->last_node_num = 1;
  p->sequence_count = (size_t)sum;
  p->passing_seq_bytes = tsta_msa_passing_seq_bytes((size_t)sum);
  p->seed = seed_index;
  p->p = head_index;

  if (tsta_node_array_reserve(&p->nodes, len_a + 1) != 0 ||
      tsta_node_index_array_reserve(&p->unsort, len_a) != 0 ||
      tsta_node_index_array_reserve(&p->sort, len_a) != 0) {
    free(p->nodes.data);
    free(p->unsort.data);
    free(p->sort.data);
    p->nodes.data = NULL;
    p->nodes.length = 0;
    p->nodes.capacity = 0;
    p->unsort.data = NULL;
    p->unsort.length = 0;
    p->unsort.capacity = 0;
    p->sort.data = NULL;
    p->sort.length = 0;
    p->sort.capacity = 0;
    return NULL;
  }

  {
    tsta_node_t *seed = tsta_graph_node(p, seed_index);
    memset(seed, 0, sizeof(*seed));
    seed->frist_col_sorce = 0;
    seed->sorce = (char *)mm_malloc(state->L * sizeof(char));
    seed->esorce = (char *)mm_malloc(state->L * sizeof(char));
    if (!seed->sorce || !seed->esorce) {
      goto fail;
    }
    memset(seed->sorce, state->E, state->L);
    memset(seed->esorce, state->E + state->E + state->O, state->L);
    seed->simple_sorce = (int *)malloc(sizeof(int));
    if (!seed->simple_sorce) {
      goto fail;
    }
    seed->simple_sorce[0] = 0;
    seed->sub = -1;
    seed->node_status = -1;
    seed->in = -1;
    seed->out = -1;
    seed->base = 'N';
    seed->node_sorce = 0;
    seed->node_base_len = 0;
    seed->prev = NULL;
    seed->next = NULL;
    seed->edge_weight = NULL;
    seed->f0 = NULL;
    seed->id = seed_index;
  }

  {
    tsta_node_t *head = tsta_graph_node(p, head_index);
    memset(head, 0, sizeof(*head));
    head->prev = (tsta_node_index_t *)malloc(sizeof(tsta_node_index_t));
    head->sorce = (char *)mm_malloc(state->L * sizeof(char));
    head->esorce = (char *)mm_malloc(state->L * sizeof(char));
    if (!head->prev || !head->sorce || !head->esorce) {
      goto fail;
    }
    head->simple_sorce = NULL;
    head->prev[0] = seed_index;
    head->base = a[0];
    head->passing_seq = (unsigned char *)calloc(
        tsta_msa_passing_seq_bytes((size_t)sum), sizeof(unsigned char));
    if (!head->passing_seq) {
      goto fail;
    }
    TSTA_MSA_PASSING_SEQ_SET(head, 0);
    head->in = 0;
    head->out = len_a > 1 ? 1 : 0;
    head->sub = 0;
    head->node_status = 0;
    head->node_sorce = 0;
    head->node_base_len = 1;
    head->edge_weight = (int *)malloc(sizeof(int));
    head->f0 = (char *)malloc(sizeof(char));
    head->mismatch_num = 0;
    head->id = head_index;
    if (!head->edge_weight || !head->f0) {
      goto fail;
    }
    head->edge_weight[0] = 0;
  }

  p->unsort.data[0] = head_index;
  p->sort.data[0] = head_index;

  for (size_t i = 1; i < len_a; i++) {
    tsta_node_t *previous = tsta_graph_node(p, end_index);
    tsta_node_t *node = tsta_graph_node(p, (tsta_node_index_t)(i + 1));

    memset(node, 0, sizeof(*node));
    previous->next = (tsta_node_index_t *)calloc(1, sizeof(tsta_node_index_t));
    node->prev = (tsta_node_index_t *)malloc(sizeof(tsta_node_index_t));
    node->sorce = (char *)mm_malloc(state->L * sizeof(char));
    node->esorce = (char *)mm_malloc(state->L * sizeof(char));
    node->passing_seq = (unsigned char *)calloc(
        tsta_msa_passing_seq_bytes((size_t)sum), sizeof(unsigned char));
    node->edge_weight = (int *)malloc(sizeof(int));
    node->f0 = (char *)malloc(sizeof(char));
    if (!previous->next || !node->prev || !node->sorce || !node->esorce ||
        !node->passing_seq || !node->edge_weight || !node->f0) {
      goto fail;
    }

    previous->next[0] = (tsta_node_index_t)(i + 1);
    node->prev[0] = end_index;
    node->base = a[i];
    node->simple_sorce = NULL;
    node->node_status = 0;
    node->node_sorce = 0;
    node->in = 1;
    node->out = 1;
    node->sub = (int)i;
    node->mismatch_num = 0;
    node->id = (uint32_t)(i + 1);
    node->edge_weight[0] = 1;
    TSTA_MSA_PASSING_SEQ_SET(node, 0);

    p->unsort.data[i] = (tsta_node_index_t)(i + 1);
    p->sort.data[i] = (tsta_node_index_t)(i + 1);
    end_index = (tsta_node_index_t)(i + 1);
  }

  tsta_graph_node(p, end_index)->out = 0;
  p->nodes.length = len_a + 1;
  p->unsort.length = len_a;
  p->sort.length = len_a;
  return tsta_graph_node(p, head_index);

fail:
  if (p->nodes.data) {
    for (size_t i = 0; i < p->nodes.capacity; i++) {
      free_node(&p->nodes.data[i]);
    }
    free(p->nodes.data);
    p->nodes.data = NULL;
    p->nodes.length = 0;
    p->nodes.capacity = 0;
  }
  free(p->unsort.data);
  p->unsort.data = NULL;
  p->unsort.length = 0;
  p->unsort.capacity = 0;
  free(p->sort.data);
  p->sort.data = NULL;
  p->sort.length = 0;
  p->sort.capacity = 0;
  p->p = TSTA_GRAPH_INVALID_INDEX;
  p->seed = TSTA_GRAPH_INVALID_INDEX;
  p->len = 0;
  p->last_node_num = 0;
  return NULL;
}

typedef struct tsta_dfs_frame {
  tsta_node_index_t node;
  int stage;
  int edge_index;
} tsta_dfs_frame;

static inline int tp1_iterative(tsta_graph_t *s, tsta_node_index_t start,
                                int subs, tsta_dfs_frame *stack,
                                int stack_cap) {
  int top = 0;

  if (stack_cap <= 0) {
    return subs;
  }

  stack[top++] = (tsta_dfs_frame){.node = start, .stage = 0, .edge_index = 0};
  while (top > 0) {
    tsta_dfs_frame *frame = &stack[top - 1];
    tsta_node_t *p = tsta_graph_node(s, frame->node);

    if (frame->stage == 0) {
      int max = 0;
      int max_i = 0;

      s->sort.data[subs] = frame->node;
      p->node_status = 0;
      p->sub = subs;

      for (int i = 0; i < p->in; i++) {
        tsta_node_t *pre = tsta_graph_node(s, p->prev[i]);

        if (pre->node_sorce >= 0) {
          if (max < p->edge_weight[i]) {
            max = p->edge_weight[i];
            max_i = i;
          } else if (max == p->edge_weight[i] &&
                     tsta_graph_node(s, p->prev[max_i])->node_sorce <=
                         pre->node_sorce) {
            max = p->edge_weight[i];
            max_i = i;
          }
        }
      }

      p->node_sorce = tsta_graph_node(s, p->prev[max_i])->node_sorce + max;
      p->node_base_len = tsta_graph_node(s, p->prev[max_i])->node_base_len + 1;
      p->node_sorce_source = tsta_graph_node(s, p->prev[max_i])->sub;
      p->in_temp = -1;
      subs++;

      frame->stage = 1;
      frame->edge_index = 0;
      continue;
    }

    if (frame->edge_index >= p->out) {
      top--;
      continue;
    }

    tsta_node_t *next = tsta_graph_node(s, p->next[frame->edge_index++]);
    next->in_temp--;

    if (next->in_temp == 0 && next->mismatch_num == 0 && next->passing != 2) {
      if (top < stack_cap) {
        stack[top++] = (tsta_dfs_frame){.node = tsta_graph_node_index(s, next),
                                        .stage = 0,
                                        .edge_index = 0};
      }
      continue;
    }

    if (next->in_temp == 0 && next->mismatch_num > 0 && next->passing != 2) {
      int ready = 0;
      for (int i = 0; i < next->mismatch_num; i++) {
        if (tsta_graph_node(s, next->mismatch_node[i])->in_temp == 0) {
          ready++;
        }
      }
      if (ready == next->mismatch_num) {
        for (int i = next->mismatch_num - 1; i >= 0; i--) {
          if (tsta_graph_node(s, next->mismatch_node[i])->in_temp == 0 &&
              top < stack_cap) {
            stack[top++] = (tsta_dfs_frame){
                .node = next->mismatch_node[i], .stage = 0, .edge_index = 0};
          }
        }
        if (top < stack_cap) {
          stack[top++] =
              (tsta_dfs_frame){.node = tsta_graph_node_index(s, next),
                               .stage = 0,
                               .edge_index = 0};
        }
      }
    }
  }

  return subs;
}

static inline tsta_graph_t *toposort1(tsta_graph_t *s) {
  int s1 = 0;
  for (size_t i = 0; i < s->unsort.length; i++) {
    tsta_node_t *node = tsta_graph_unsort_node(s, i);

    node->in_temp = node->in;
    node->passing = 0;
    if (node->out == 0 && node->mismatch_num > 0) {
      for (int j = 0; j < node->mismatch_num; j++)
        if (tsta_graph_node(s, node->mismatch_node[j])->out != 0)
          s1 = 1;
      if (s1 != 1)
        node->passing = 2;
    }
    s1 = 0;
  }

  int c = 0;
  int subs = 0;
  int stack_cap = s->len * 2 + 8;
  tsta_dfs_frame *stack =
      (tsta_dfs_frame *)malloc((size_t)stack_cap * sizeof(tsta_dfs_frame));
  if (!stack) {
    return s;
  }

  while (subs < s->len)
    for (size_t i = 0; i < s->unsort.length; i++) {
      tsta_node_t *node = tsta_graph_unsort_node(s, i);

      if (node->in_temp == 0) {
        if (node->mismatch_num == 0) {
          subs = tp1_iterative(s, s->unsort.data[i], subs, stack, stack_cap);
          break;
        } else if (node->in_temp == 0 && node->mismatch_num > 0) {
          for (int j = 0; j < node->mismatch_num; j++)
            if (tsta_graph_node(s, node->mismatch_node[j])->in_temp == 0)
              c++;
          if (c == node->mismatch_num) {
            c = 0;
            subs = tp1_iterative(s, s->unsort.data[i], subs, stack, stack_cap);
            for (int ss = 0; ss < node->mismatch_num; ss++)
              if (tsta_graph_node(s, node->mismatch_node[ss])->in_temp == 0)
                subs = tp1_iterative(s, node->mismatch_node[ss], subs, stack,
                                     stack_cap);
            break;
          }
          c = 0;
        }
      }
    }

  free(stack);

  return s;
}

static inline tsta_graph_t *modify(tsta_graph_t *p) {
  while (1) {
    int max = INT_MIN;
    int max_i = 0;

    for (size_t i = 0; i < p->sort.length; i++) {
      tsta_node_t *node = tsta_graph_sort_node(p, i);

      if (max <= node->node_sorce) {
        max = node->node_sorce;
        max_i = node->sub;
      }
    }

    tsta_node_t *best = tsta_graph_sort_node(p, (size_t)max_i);

    if (best->out == 0) {
      return p;
    }

    for (int i = 0; i < best->out; i++) {
      tsta_node_t *next = tsta_graph_node(p, best->next[i]);

      for (int j = 0; j < next->in; j++) {
        tsta_node_t *prev = tsta_graph_node(p, next->prev[j]);

        if (prev->node_sorce < best->node_sorce && prev->node_sorce > 0) {
          prev->node_sorce = -prev->node_sorce;
        }
      }
      next->node_status = 4;
    }

    for (int i = max_i + 1; i < p->len; i++) {
      int best = 0;
      int best_i = 0;
      tsta_node_t *node = tsta_graph_sort_node(p, (size_t)i);

      if (node->node_sorce >= 0 || node->node_status == 4) {
        for (int j = 0; j < node->in; j++) {
          tsta_node_t *prev = tsta_graph_node(p, node->prev[j]);

          if (prev->node_sorce >= 0) {
            if (best < node->edge_weight[j]) {
              best = node->edge_weight[j];
              best_i = j;
            } else if (best == node->edge_weight[j] &&
                       tsta_graph_node(p, node->prev[best_i])->node_sorce <=
                           prev->node_sorce) {
              best = node->edge_weight[j];
              best_i = j;
            }
          }
        }
        node->node_sorce =
            tsta_graph_node(p, node->prev[best_i])->node_sorce + best;
        node->node_base_len =
            tsta_graph_node(p, node->prev[best_i])->node_base_len + 1;
        node->node_sorce_source = tsta_graph_node(p, node->prev[best_i])->sub;
        node->node_status = 0;
      }
    }
  }
}

static inline int tp_iterative(tsta_graph_t *s, tsta_node_index_t start,
                               int subs, tsta_dfs_frame *stack, int stack_cap) {
  int top = 0;

  if (stack_cap <= 0) {
    return subs;
  }

  stack[top++] = (tsta_dfs_frame){.node = start, .stage = 0, .edge_index = 0};
  while (top > 0) {
    tsta_dfs_frame *frame = &stack[top - 1];
    tsta_node_t *p = tsta_graph_node(s, frame->node);

    if (frame->stage == 0) {
      s->sort.data[subs] = frame->node;
      p->node_status = 0;
      p->sub = subs;
      p->in_temp = -1;
      subs++;

      frame->stage = 1;
      frame->edge_index = 0;
      continue;
    }

    if (frame->stage == 1) {
      if (frame->edge_index >= p->out) {
        frame->stage = 2;
        frame->edge_index = 0;
        continue;
      }

      tsta_node_t *next = tsta_graph_node(s, p->next[frame->edge_index++]);
      if (next->out == 0 && next->passing == 1 && next->in_temp - 1 == 0) {
        next->in_temp--;
        if (next->in_temp == 0 && top < stack_cap) {
          stack[top++] =
              (tsta_dfs_frame){.node = tsta_graph_node_index(s, next),
                               .stage = 0,
                               .edge_index = 0};
        }
      }
      continue;
    }

    if (frame->edge_index >= p->out) {
      top--;
      continue;
    }

    tsta_node_t *next = tsta_graph_node(s, p->next[frame->edge_index++]);
    next->in_temp--;
    if (next->in_temp == 0 && next->passing != 2 && top < stack_cap) {
      stack[top++] = (tsta_dfs_frame){
          .node = tsta_graph_node_index(s, next), .stage = 0, .edge_index = 0};
    }
  }

  return subs;
}

static inline tsta_graph_t *toposort(tsta_graph_t *s) {
  int s1 = 0;
  for (size_t i = 0; i < s->unsort.length; i++) {
    tsta_node_t *node = tsta_graph_unsort_node(s, i);

    node->in_temp = node->in;
    node->passing = 0;
    if (node->out == 0 && node->mismatch_num > 0) {
      for (int j = 0; j < node->mismatch_num; j++)
        if (tsta_graph_node(s, node->mismatch_node[j])->out != 0) {
          node->passing = 1;
          s1 = 1;
        }
      if (s1 != 1)
        node->passing = 2;
    }
    s1 = 0;
  }

  int subs = 0;
  int stack_cap = s->len * 2 + 8;
  tsta_dfs_frame *stack =
      (tsta_dfs_frame *)malloc((size_t)stack_cap * sizeof(tsta_dfs_frame));
  if (!stack) {
    return s;
  }

  while (subs < s->len)
    for (size_t i = 0; i < s->unsort.length; i++)
      if (tsta_graph_unsort_node(s, i)->in_temp == 0) {
        subs = tp_iterative(s, s->unsort.data[i], subs, stack, stack_cap);
        if (subs + s->last_node_num == s->len) {
          for (size_t i1 = 0; i1 < s->unsort.length; i1++)
            if (tsta_graph_unsort_node(s, i1)->in_temp == 0)
              subs =
                  tp_iterative(s, s->unsort.data[i1], subs, stack, stack_cap);
        }
        break;
      }

  free(stack);
  return s;
}

tsta_graph_t *tsta_graph_sort(tsta_graph_t *g, int num) {
  tsta_graph_t *s;
  g->last_node_num = 0;
  for (size_t i = 0; i < g->unsort.length; i++)
    if (tsta_graph_unsort_node(g, i)->out == 0)
      g->last_node_num++;

  if (tsta_node_index_array_reserve(&g->sort, (size_t)g->len) != 0 ||
      tsta_node_index_array_reserve(&g->unsort, (size_t)g->len) != 0) {
    return NULL;
  }
  g->sort.length = (size_t)g->len;

  if (num != 1)
    s = toposort(g);
  else {
    s = toposort1(g);
    s = modify(s);
  }

  for (size_t i = 0; i < s->sort.length; i++) {
    s->unsort.data[i] = s->sort.data[i];
  }
  s->unsort.length = s->sort.length;
  return s;
}

tsta_graph_t *tsta_graph_update_impl(tsta_graph_t *n,
                                     const tsta_sequence_view_t *sequence,
                                     int num, int sum, int last,
                                     tsta_msa_state *state) {
  size_t base_nodes;
  size_t base_unsort;
  size_t len_b;
  tsta_node_t **seq = NULL;
  unsigned char *live_flags = NULL;
  int num1;
  int num2;
  int cont = 0;
  size_t new_node_count = 0;
  int s1 = INT_MIN;
  int s2 = 0;
  int s4 = 0;
  char s5 = 0;

  if (!n || !sequence || !sequence->sequence || sequence->length <= 0) {
    return n;
  }

  len_b = (size_t)sequence->length;
  base_nodes = n->nodes.length;
  base_unsort = n->unsort.length;

  if (tsta_node_array_reserve(&n->nodes, base_nodes + len_b) != 0) {
    return NULL;
  }
  if (tsta_node_index_array_reserve(&n->unsort, base_unsort + len_b) != 0) {
    return NULL;
  }

  seq = (tsta_node_t **)calloc(len_b, sizeof(tsta_node_t *));
  live_flags = (unsigned char *)calloc(len_b, sizeof(unsigned char));
  if (!seq || !live_flags) {
    free(seq);
    free(live_flags);
    return NULL;
  }

  for (size_t i = 0; i < len_b; i++) {
    seq[i] = tsta_graph_node(n, (tsta_node_index_t)(base_nodes + i));
    memset(seq[i], 0, sizeof(*seq[i]));
    seq[i]->prev = (tsta_node_index_t *)malloc(sizeof(tsta_node_index_t));
    if (!seq[i]->prev) {
      goto fail;
    }
    seq[i]->sorce = seq[i]->esorce = NULL;
    seq[i]->simple_sorce = NULL;
    seq[i]->node_status = 0;
    seq[i]->node_sorce = 0;
    seq[i]->node_base_len = 1;
    seq[i]->edge_weight = NULL;
    seq[i]->passing_seq = NULL;
    seq[i]->mismatch_num = 0;
    seq[i]->f0 = NULL;
    seq[i]->id = (uint32_t)(base_nodes + i);
  }

  seq[0]->base = sequence->sequence[0];
  seq[0]->prev[0] = n->seed;
  seq[0]->in = 0;
  seq[0]->out = len_b > 1 ? 1 : 0;
  seq[0]->sub = -1;
  seq[0]->node_status = 0;

  for (size_t i = 1; i < len_b; i++) {
    seq[i]->base = sequence->sequence[i];
    seq[i]->prev[0] = (tsta_node_index_t)(base_nodes + i - 1);
    seq[i - 1]->next =
        (tsta_node_index_t *)calloc(1, sizeof(tsta_node_index_t));
    if (!seq[i - 1]->next) {
      goto fail;
    }
    seq[i - 1]->next[0] = (tsta_node_index_t)(base_nodes + i);
    seq[i]->in = 1;
    seq[i]->out = 1;
    seq[i]->sub = -1;
  }
  seq[len_b - 1]->out = 0;

  num1 = n->len - 1;
  num2 = (int)len_b - 1;

  for (int i = n->len - 1; i > 0; i--) {
    tsta_node_t *node = tsta_graph_sort_node(n, (size_t)i);

    if (node->out == 0) {
      if (s1 <= node->lastsorce) {
        s1 = node->lastsorce;
        num1 = node->sub;
      }
      s2++;
    }
    if (s2 >= n->last_node_num) {
      break;
    }
  }

  while (num1 != -1 && num2 != -1) {
    tsta_node_t *current = tsta_graph_sort_node(n, (size_t)num1);
    int trace = TRACE_SOURCE(current, NUM2(num2));
    size_t trace_parent_index = TSTA_TRACE_PARENT_SLOT(trace);

    if (trace / 42 == 3) {
      tsta_node_t *seq_node = seq[num2];

      cont = 0;
      if (tsta_graph_activate_sequence_node(seq_node, (size_t)sum, num, state,
                                            1) != 0) {
        goto fail;
      }
      seq_node->sub = (int)n->unsort.length;
      n->unsort.data[n->unsort.length++] =
          (tsta_node_index_t)(base_nodes + (size_t)num2);
      new_node_count++;
      live_flags[num2] = 1;
      if (num2 > 0 && current->fsource_store && NUM2(num2 - 1) > 0 &&
          ((TRACE_FSOURCE(current, NUM2(num2)) == 1 ||
            TRACE_FSOURCE(current, NUM2(num2)) == -1) ||
           ((TRACE_FSOURCE(current, NUM2(num2)) == 2 ||
             TRACE_FSOURCE(current, NUM2(num2)) == -2) &&
            TRACE_FSOURCE(current, NUM2(num2 - 1)) < 0))) {
        TRACE_SOURCE_SET(current, NUM2(num2 - 1), 126);
      }
      num2--;
      continue;
    }

    if (trace / 42 == 0) {
      tsta_node_index_t parent_index = current->prev[trace_parent_index];
      tsta_node_t *parent = tsta_graph_node(n, parent_index);

      cont = 3;
      if (parent->sub > 0 &&
          ((TRACE_ESOURCE(current, NUM2(num2)) <= 42 &&
            TRACE_ESOURCE(current, NUM2(num2)) >= -42) ||
           ((TRACE_ESOURCE(current, NUM2(num2)) > 42 ||
             TRACE_ESOURCE(current, NUM2(num2)) < -42) &&
            TRACE_ESOURCE(tsta_graph_node(n, parent->prev[trace_parent_index]),
                          NUM2(num2)) < 0))) {
        s5 = TRACE_ESOURCE(tsta_graph_node(n, parent->prev[trace_parent_index]),
                           NUM2(num2)) %
             42;
        s5 = (s5 >= 0 ? s5 : -s5) - 1;
        TRACE_SOURCE_SET(tsta_graph_node(n, parent->prev[trace_parent_index]),
                         NUM2(num2), s5);
      }
      num1 = parent->sub;
      continue;
    }

    if (trace / 42 == 1) {
      if (num2 == 0) {
        if (cont == 1 || cont == 5) {
          if (live_flags[num2]) {
            live_flags[num2] = 0;
          }
          free_node(seq[num2]);
          seq[num2] = current;
        } else if (len_b > 1) {
          tsta_node_index_t next_index = (tsta_node_index_t)(base_nodes + 1);

          current->out++;
          current->next = (tsta_node_index_t *)realloc(
              current->next, (size_t)current->out * sizeof(tsta_node_index_t));
          if (!current->next) {
            goto fail;
          }
          current->next[current->out - 1] = next_index;
          seq[1]->prev[seq[1]->in - 1] = current->id;
        }
      } else if (num2 == (int)len_b - 1) {
        tsta_node_index_t prev_index = current->prev[trace_parent_index];
        tsta_node_t *prev = tsta_graph_node(n, prev_index);

        if (TRACE_SOURCE(prev, NUM2(num2 - 1)) / 42 == 1) {
          current->edge_weight[trace_parent_index]++;
          if (live_flags[num2]) {
            live_flags[num2] = 0;
          }
          free_node(seq[num2]);
        } else {
          current->in++;
          current->prev = (tsta_node_index_t *)realloc(
              current->prev, (size_t)current->in * sizeof(tsta_node_index_t));
          if (!current->prev) {
            goto fail;
          }
          current->prev[current->in - 1] = seq[num2 - 1]->id;

          current->edge_weight = (int *)realloc(
              current->edge_weight, (size_t)current->in * sizeof(int));
          if (!current->edge_weight) {
            goto fail;
          }
          current->edge_weight[current->in - 1] = 1;
          current->f0 =
              (char *)realloc(current->f0, (size_t)current->in * sizeof(char));
          if (!current->f0) {
            goto fail;
          }
          seq[num2 - 1]->next[seq[num2 - 1]->out - 1] = current->id;
          if (live_flags[num2]) {
            live_flags[num2] = 0;
          }
          free_node(seq[num2]);
          seq[num2] = current;
        }
      } else {
        tsta_node_index_t prev_index = current->prev[trace_parent_index];
        tsta_node_t *prev = tsta_graph_node(n, prev_index);

        if (prev->sub != -1 &&
            TRACE_SOURCE(tsta_graph_node(n, prev->id), NUM2(num2 - 1)) / 42 ==
                1) {
          current->edge_weight[trace_parent_index]++;
        } else {
          current->in++;
          current->prev = (tsta_node_index_t *)realloc(
              current->prev, (size_t)current->in * sizeof(tsta_node_index_t));
          if (!current->prev) {
            goto fail;
          }
          current->prev[current->in - 1] = seq[num2 - 1]->id;

          current->edge_weight = (int *)realloc(
              current->edge_weight, (size_t)current->in * sizeof(int));
          if (!current->edge_weight) {
            goto fail;
          }
          current->edge_weight[current->in - 1] = 1;
          current->f0 =
              (char *)realloc(current->f0, (size_t)current->in * sizeof(char));
          if (!current->f0) {
            goto fail;
          }
          seq[num2 - 1]->next[seq[num2 - 1]->out - 1] = current->id;
        }

        if (cont != 1 && cont != 5) {
          current->out++;
          current->next = (tsta_node_index_t *)realloc(
              current->next, (size_t)current->out * sizeof(tsta_node_index_t));
          if (!current->next) {
            goto fail;
          }
          current->next[current->out - 1] = seq[num2 + 1]->id;
          seq[num2 + 1]->prev[seq[num2 + 1]->in - 1] = current->id;
        }

        if (live_flags[num2]) {
          live_flags[num2] = 0;
        }
        free_node(seq[num2]);
        seq[num2] = current;
      }

      cont = 1;
      TSTA_MSA_PASSING_SEQ_SET(current, (size_t)num);
      num1 = tsta_graph_node(n, current->prev[trace_parent_index])->sub;
      num2--;
      continue;
    }

    s4 = 0;
    for (int s = 0; s < current->mismatch_num; s++) {
      tsta_node_index_t mismatch_index = current->mismatch_node[s];
      tsta_node_t *mismatch = tsta_graph_node(n, mismatch_index);

      if (seq[num2]->base == mismatch->base) {
        if (num2 != 0) {
          tsta_node_index_t parent_index = current->prev[trace_parent_index];
          tsta_node_t *parent = tsta_graph_node(n, parent_index);

          if (parent->sub != -1) {
            if (TRACE_SOURCE(parent, NUM2(num2 - 1)) / 42 == 1) {
              for (int ss = 0; ss < mismatch->in; ss++) {
                if (mismatch->prev[ss] == parent_index) {
                  mismatch->edge_weight[ss]++;
                  s2 = -1;
                }
              }
            }
          }
          if (s2 != -1) {
            mismatch->in++;
            mismatch->prev = (tsta_node_index_t *)realloc(
                mismatch->prev,
                (size_t)mismatch->in * sizeof(tsta_node_index_t));
            if (!mismatch->prev) {
              goto fail;
            }
            mismatch->prev[mismatch->in - 1] = seq[num2 - 1]->id;

            mismatch->edge_weight = (int *)realloc(
                mismatch->edge_weight, (size_t)mismatch->in * sizeof(int));
            if (!mismatch->edge_weight) {
              goto fail;
            }
            mismatch->edge_weight[mismatch->in - 1] = 1;
            seq[num2 - 1]->next[seq[num2 - 1]->out - 1] = mismatch_index;
          }
        }

        s4 = 1;
        if (cont == 1 || cont == 4) {
          for (int ss = 0; ss < seq[num2 + 1]->in; ss++) {
            if (seq[num2 + 1]->prev[ss] == mismatch_index) {
              s4 = 2;
              seq[num2 + 1]->edge_weight[ss]++;
              seq[num2 + 1]->in--;
              seq[num2 + 1]->prev = (tsta_node_index_t *)realloc(
                  seq[num2 + 1]->prev,
                  (size_t)seq[num2 + 1]->in * sizeof(tsta_node_index_t));
              if (!seq[num2 + 1]->prev) {
                goto fail;
              }
              seq[num2 + 1]->edge_weight =
                  (int *)realloc(seq[num2 + 1]->edge_weight,
                                 (size_t)seq[num2 + 1]->in * sizeof(int));
              if (!seq[num2 + 1]->edge_weight) {
                goto fail;
              }
              seq[num2 + 1]->f0 = (char *)realloc(
                  seq[num2 + 1]->f0, (size_t)seq[num2 + 1]->in * sizeof(char));
              if (!seq[num2 + 1]->f0) {
                goto fail;
              }
            }
          }
        }

        if (s4 == 1 && num2 != (int)len_b - 1) {
          seq[num2 + 1]->prev[seq[num2 + 1]->in - 1] = mismatch_index;
          mismatch->out++;
          mismatch->next = (tsta_node_index_t *)realloc(
              mismatch->next,
              (size_t)mismatch->out * sizeof(tsta_node_index_t));
          if (!mismatch->next) {
            goto fail;
          }
          mismatch->next[mismatch->out - 1] = seq[num2 + 1]->id;
        }

        TSTA_MSA_PASSING_SEQ_SET(mismatch, (size_t)num);
        if (s2 == -1) {
          cont = 5;
        } else {
          cont = 4;
        }
        s2 = 0;
        if (live_flags[num2]) {
          live_flags[num2] = 0;
        }
        free_node(seq[num2]);
        seq[num2] = mismatch;
      }
    }

    if (s4 == 0) {
      tsta_node_t *seq_node = seq[num2];

      cont = 2;
      if (tsta_graph_activate_sequence_node(seq_node, (size_t)sum, num, state,
                                            1) != 0) {
        goto fail;
      }

      seq_node->sub = (int)n->unsort.length;
      n->unsort.data[n->unsort.length++] =
          (tsta_node_index_t)(base_nodes + (size_t)num2);
      new_node_count++;
      live_flags[num2] = 1;

      size_t mismatch_count = (size_t)current->mismatch_num;
      size_t new_mismatch_count = mismatch_count + 1;
      int mismatch_failed = 0;

      if (tsta_msa_mismatch_list_reserve(current, new_mismatch_count) != 0 ||
          tsta_msa_mismatch_list_reserve(seq_node, new_mismatch_count) != 0) {
        if (n->unsort.length > 0) {
          n->unsort.data[n->unsort.length - 1] = TSTA_GRAPH_INVALID_INDEX;
          n->unsort.length--;
        }
        if (new_node_count > 0) {
          new_node_count--;
        }
        if (live_flags[num2]) {
          live_flags[num2] = 0;
        }
        free_node(seq_node);
        seq[num2] = current;
        num1 = tsta_graph_node(n, current->prev[trace_parent_index])->sub;
        num2--;
        continue;
      }

      for (size_t s = 0; s < mismatch_count; s++) {
        tsta_node_index_t old_mismatch_index = current->mismatch_node[s];
        tsta_node_t *old_mismatch = tsta_graph_node(n, old_mismatch_index);

        if (tsta_msa_mismatch_list_reserve(
                old_mismatch, (size_t)old_mismatch->mismatch_num + 1) != 0) {
          if (n->unsort.length > 0) {
            n->unsort.data[n->unsort.length - 1] = TSTA_GRAPH_INVALID_INDEX;
            n->unsort.length--;
          }
          if (new_node_count > 0) {
            new_node_count--;
          }
          if (live_flags[num2]) {
            live_flags[num2] = 0;
          }
          free_node(seq_node);
          seq[num2] = current;
          num1 = tsta_graph_node(n, current->prev[trace_parent_index])->sub;
          num2--;
          mismatch_failed = 1;
          goto mismatch_done;
        }
      }

      current->mismatch_num++;
      current->mismatch_node[current->mismatch_num - 1] =
          (tsta_node_index_t)(base_nodes + (size_t)num2);
      seq_node->mismatch_num = current->mismatch_num;
      seq_node->mismatch_node[seq_node->mismatch_num - 1] = current->id;
      for (size_t s = 0; s < mismatch_count; s++) {
        tsta_node_index_t old_mismatch_index = current->mismatch_node[s];
        tsta_node_t *old_mismatch = tsta_graph_node(n, old_mismatch_index);

        old_mismatch->mismatch_num++;
        old_mismatch->mismatch_node[old_mismatch->mismatch_num - 1] =
            (tsta_node_index_t)(base_nodes + (size_t)num2);
        seq_node->mismatch_node[s] = old_mismatch_index;
      }

    mismatch_done:
      if (mismatch_failed) {
        continue;
      }
    }

    num1 = tsta_graph_node(n, current->prev[trace_parent_index])->sub;
    num2--;
  }

  while (num2 > -1) {
    tsta_node_t *seq_node = seq[num2];

    if (tsta_graph_activate_sequence_node(seq_node, (size_t)sum, num, state,
                                          1) != 0) {
      goto fail;
    }
    seq_node->sub = (int)n->unsort.length;
    n->unsort.data[n->unsort.length++] =
        (tsta_node_index_t)(base_nodes + (size_t)num2);
    new_node_count++;
    live_flags[num2] = 1;
    num2--;
  }

  {
    size_t live_count = 0;

    for (size_t i = 0; i < len_b; i++) {
      if (live_flags[i]) {
        live_count++;
      }
    }

    (void)live_count;

    n->len += (int)new_node_count;
    n->sequence_count = (size_t)sum;
    n->passing_seq_bytes = tsta_msa_passing_seq_bytes((size_t)sum);
    n->nodes.length = base_nodes + len_b;
  }

  free(seq);
  free(live_flags);
  return n;

fail:
  n->unsort.length = base_unsort;
  n->nodes.length = base_nodes;
  if (seq) {
    for (size_t i = 0; i < len_b; i++) {
      if (seq[i]) {
        free_node(seq[i]);
      }
    }
  }
  free(seq);
  free(live_flags);
  return NULL;
}
