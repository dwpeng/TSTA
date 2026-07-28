#ifndef TSTA_GRAPH_H
#define TSTA_GRAPH_H

#include "tsta_common.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef uint32_t tsta_node_index_t;

/* ── Graph node (internal) ──────────────────────────────────────────── */

typedef struct tsta_node_t {
    tsta_node_index_t *prev;
    tsta_node_index_t *next;
    char *sorce;
    char *esorce;
    tsta_trace_block_store_t *source_store;
    tsta_trace_block_store_t *esource_store;
    tsta_trace_block_store_t *fsource_store;
    unsigned char *passing_seq;
    tsta_node_index_t *mismatch_node;
    int *simple_sorce;
    char *f0;
    int *edge_weight;
    int mismatch_num;
    int sub;
    int frist_col_sorce;
    int in_temp;
    int in;
    int out;
    int lastsorce;
    int node_sorce;
    int node_sorce_source;
    int node_base_len;
    uint32_t id;
    char base;
    int8_t node_status;
    int8_t passing;
    uint8_t reserved_flags;
} tsta_node_t;

define_array(tsta_node_t, tsta_node_array_t, tsta_node_array);
define_array(tsta_node_index_t, tsta_node_index_array_t, tsta_node_index_array);

/* ── POA graph (internal) ───────────────────────────────────────────── */

typedef struct tsta_graph_t {
    int len;
    int last_node_num;
    tsta_node_index_t p;
    tsta_node_index_t seed;
    tsta_node_array_t nodes;
    tsta_node_index_array_t unsort;
    tsta_node_index_array_t sort;
    size_t sequence_count;
    size_t passing_seq_bytes;
} tsta_graph_t;

/* ── Node access helpers ────────────────────────────────────────────── */

static inline tsta_node_t *tsta_graph_node(tsta_graph_t *graph,
                                           tsta_node_index_t index) {
    return &graph->nodes.data[index];
}

static inline const tsta_node_t *
tsta_graph_cnode(const tsta_graph_t *graph, tsta_node_index_t index) {
    return &graph->nodes.data[index];
}

static inline tsta_node_t *tsta_graph_sort_node(tsta_graph_t *graph,
                                                size_t index) {
    return tsta_graph_node(graph, graph->sort.data[index]);
}

static inline tsta_node_t *tsta_graph_unsort_node(tsta_graph_t *graph,
                                                  size_t index) {
    return tsta_graph_node(graph, graph->unsort.data[index]);
}

static inline tsta_node_index_t
tsta_graph_node_index(const tsta_graph_t *graph, const tsta_node_t *node) {
    return (tsta_node_index_t)(node - graph->nodes.data);
}

/* ── Graph operations ───────────────────────────────────────────────── */

size_t tsta_msa_passing_seq_bytes(size_t sequence_count);
int tsta_msa_expand_passing_seq(tsta_graph_t *graph, size_t old_count,
                                size_t new_count);
void tsta_msa_free_graph(tsta_graph_t *graph);
void free_node(tsta_node_t *n);
tsta_node_t *tsta_graph_build_initial(tsta_graph_t *p,
                                      const tsta_sequence_view_t *sequence,
                                      int sum, tsta_msa_state *state);
tsta_graph_t *tsta_graph_sort(tsta_graph_t *g, int num);
tsta_graph_t *tsta_graph_update_impl(tsta_graph_t *n,
                                     const tsta_sequence_view_t *sequence,
                                     int num, int sum, int last,
                                     tsta_msa_state *state);

#ifdef __cplusplus
}
#endif

#endif
