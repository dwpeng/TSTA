#ifndef TSTA_COMMON_H
#define TSTA_COMMON_H

#include "tsta.h"
#include "tsta_array.h"
#include "tsta_simd.h"

#include <pthread.h>
#include <stddef.h>
#include <stdint.h>

/* ── Sequence view (internal) ────────────────────────────────────────── */

typedef struct tsta_sequence_view_t {
    const char *sequence;
    int length;
} tsta_sequence_view_t;

int tsta_sequence_view_make(tsta_sequence_view_t *view,
                            const char *sequence, int length);

/* ── Packed sequence cache ──────────────────────────────────────────── */

typedef struct packed_sequence_t {
    char *packed_seq;
    int len;
} packed_sequence_t;

define_array(packed_sequence_t, packed_sequence_array_t,
             tsta_packed_sequence_array);

/* ── Simplified trace block store (memory only) ──────────────────────── */

typedef struct tsta_trace_block_store {
    uint8_t *buffer;
    size_t length;
    size_t chunk_size;
    size_t chunk_count;
} tsta_trace_block_store_t;

tsta_trace_block_store_t *tsta_trace_block_store_create(size_t length,
                                                        size_t chunk_size);
void tsta_trace_block_store_destroy(tsta_trace_block_store_t *store);

static inline uint8_t *
tsta_trace_block_store_chunk(tsta_trace_block_store_t *store,
                             size_t chunk_index) {
    if (!store || chunk_index >= store->chunk_count)
        return NULL;
    return store->buffer + chunk_index * store->chunk_size;
}

static inline char
tsta_trace_block_store_get_byte(tsta_trace_block_store_t *store,
                                size_t index) {
    size_t ci = index / store->chunk_size;
    size_t off = index % store->chunk_size;
    char *chunk = (char *)tsta_trace_block_store_chunk(store, ci);
    return chunk ? chunk[off] : 0;
}

static inline int
tsta_trace_block_store_set_byte(tsta_trace_block_store_t *store,
                                size_t index, char value) {
    size_t ci = index / store->chunk_size;
    size_t off = index % store->chunk_size;
    char *chunk = (char *)tsta_trace_block_store_chunk(store, ci);
    if (!chunk)
        return -1;
    chunk[off] = value;
    return 0;
}

/* ── PSA internal state ─────────────────────────────────────────────── */

typedef struct tsta_psa_trace_matrices {
    char **back;
    char **eback;
    char **fback;
    int rows;
    int cols;
} tsta_psa_trace_matrices;

typedef struct tsta_psa_trace_result {
    char *aligned_seq1;
    char *aligned_seq2;
    size_t aligned_length;
    char *formatted_text;
    size_t formatted_length;
} tsta_psa_trace_result;

typedef struct tsta_psa_trace_context {
    tsta_psa_trace_matrices matrices;
    tsta_psa_trace_result result;
} tsta_psa_trace_context;

typedef struct tsta_psa_state {
    int bS;
    int M;
    int X;
    int E;
    int O;
    int L;
    int W;
    int B;
    int *real;
    tsta_psa_trace_context trace;
    char *sorce;
    char *esorce;
    char *V;
    char *F;
    char **seq;
    int lmaxtag;
    int fmaxtag;
    int length[4];
    int trace_enabled;
    volatile int ms;
    pthread_mutex_t mutex;
    volatile int lock;
} tsta_psa_state;

/* ── MSA internal state ─────────────────────────────────────────────── */

typedef struct tsta_msa_state {
    int bS;
    int M;
    int X;
    int E;
    int O;
    int L;
    int W;
    int B;
    int s_len;
    int maxtag;
    int fmaxtag;
    int lmaxtag;
    int length1;
    int length2;
    char z;
    pthread_mutex_t mutex;
    volatile int lock;
} tsta_msa_state;

/* ── Internal init / helpers ────────────────────────────────────────── */

void tsta_init_psa_state(tsta_psa_state *state, const tsta_config *config,
                         int simd_block);
void tsta_init_msa_state(tsta_msa_state *state, const tsta_config *config,
                         int simd_block);
void tsta_config_apply_length_fallback(tsta_config *config,
                                       size_t longest_sequence_length,
                                       size_t sequence_count);

int tsta_packed_sequence_array_reserve_slot(packed_sequence_array_t *array,
                                            size_t slot,
                                            size_t initial_capacity);
void tsta_packed_sequence_array_reset_malloc(packed_sequence_array_t *array);
void tsta_packed_sequence_array_reset_mm(packed_sequence_array_t *array);

void *tsta_aligned_malloc(size_t size, int base);
void tsta_aligned_free(void *buffer);

#endif
