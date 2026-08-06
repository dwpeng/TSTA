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

/* ── Trace buffer access (contiguous per-node block, direct indexing) ──
 * traces holds three planes of trace_width bytes: source, esource, fsource.
 * NULL traces (seed / not-yet-created nodes) read as 0. */

#define TRACE_SOURCE(node, index)                                             \
  ((node)->traces ? (char)((node)->traces[(size_t)(index)]) : 0)

#define TRACE_SOURCE_SET(node, index, value)                                  \
  do {                                                                        \
    if ((node)->traces)                                                       \
      (node)->traces[(size_t)(index)] = (char)(value);                        \
  } while (0)

#define TRACE_ESOURCE(node, index)                                            \
  ((node)->traces                                                             \
       ? (char)((node)->traces[(node)->trace_width + (size_t)(index)])        \
       : 0)

#define TRACE_ESOURCE_SET(node, index, value)                                 \
  do {                                                                        \
    if ((node)->traces)                                                       \
      (node)->traces[(node)->trace_width + (size_t)(index)] = (char)(value);  \
  } while (0)

#define TRACE_FSOURCE(node, index)                                            \
  ((node)->traces                                                             \
       ? (char)((node)->traces[2 * (node)->trace_width + (size_t)(index)])    \
       : 0)

#define TRACE_FSOURCE_SET(node, index, value)                                 \
  do {                                                                        \
    if ((node)->traces)                                                       \
      (node)->traces[2 * (node)->trace_width + (size_t)(index)] =             \
          (char)(value);                                                      \
  } while (0)

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
