#ifndef TSTA_INTERNAL_H
#define TSTA_INTERNAL_H

#include "tsta_common.h"
#include "tsta_graph.h"

/* ── Shared alignment macros ────────────────────────────────────────── */

#define nconvert(n)                                                           \
  (state->maxtag > 0 ? (!!(n / state->maxtag)) * state->maxtag                \
                           + (!(n / state->maxtag)) * (n % state->maxtag)     \
                     : 0)

#define NUM2(num2)                                                            \
  ((num2) / state->L) * state->L                                              \
      + ((((num2) % state->L) % state->W) * state->B                          \
         + (((num2) % state->L) / state->W))

/* ── Trace buffer access (direct uint8_t * indexing) ─────────────────── */

#define TRACE_SOURCE(node, index)                                             \
  tsta_trace_block_store_get_byte((node)->source_store, (size_t)(index))

#define TRACE_SOURCE_SET(node, index, value)                                  \
  tsta_trace_block_store_set_byte((node)->source_store, (size_t)(index),      \
                                  (char)(value))

#define TRACE_ESOURCE(node, index)                                            \
  tsta_trace_block_store_get_byte((node)->esource_store, (size_t)(index))

#define TRACE_ESOURCE_SET(node, index, value)                                 \
  tsta_trace_block_store_set_byte((node)->esource_store, (size_t)(index),     \
                                  (char)(value))

#define TRACE_FSOURCE(node, index)                                            \
  tsta_trace_block_store_get_byte((node)->fsource_store, (size_t)(index))

#define TRACE_FSOURCE_SET(node, index, value)                                 \
  tsta_trace_block_store_set_byte((node)->fsource_store, (size_t)(index),     \
                                  (char)(value))

#define TSTA_TRACE_PARENT_SLOT(trace)                                         \
  ((size_t)((trace) < 0 ? 0 : ((trace) % 42)))

/* ── Passing sequence bit array ──────────────────────────────────────── */

#define TSTA_MSA_PASSING_SEQ_SET(node, index)                                 \
  do {                                                                        \
    (node)->passing_seq[(size_t)(index) >> 3] |=                              \
        (unsigned char)(1u << ((size_t)(index) & 7u));                        \
  } while (0)

#define TSTA_MSA_PASSING_SEQ_GET(node, index)                                 \
  (((node)->passing_seq[(size_t)(index) >> 3] >> ((size_t)(index) & 7u)) & 1u)

#endif
