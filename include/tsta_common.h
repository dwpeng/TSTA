#ifndef TSTA_COMMON_H
#define TSTA_COMMON_H

#include "tsta.h"
#include "tsta_array.h"
#include "tsta_common_simd.h"

#include <pthread.h>
#include <stddef.h>

typedef struct tsta_sequence_view_t {
  const char *sequence;
  int length;
} tsta_sequence_view_t;

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

typedef struct tsta_trace_block_store tsta_trace_block_store_t;

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
  int trace_dump_enabled;
  int trace_dump_compress;
  char z;
  pthread_mutex_t mutex;
  volatile int lock;
} tsta_msa_state;

typedef struct packed_sequence_t {
  char *packed_seq;
  int len;
} packed_sequence_t;

define_array(packed_sequence_t, packed_sequence_array_t,
             tsta_packed_sequence_array);

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

tsta_trace_block_store_t *
tsta_trace_block_store_create(size_t length, size_t chunk_size, int compress);
void tsta_trace_block_store_destroy(tsta_trace_block_store_t *store);
int tsta_trace_block_store_spill(tsta_trace_block_store_t *store);
char *tsta_trace_block_store_chunk(tsta_trace_block_store_t *store,
                                   size_t chunk_index);
char tsta_trace_block_store_get_byte(tsta_trace_block_store_t *store,
                                     size_t index);
int tsta_trace_block_store_set_byte(tsta_trace_block_store_t *store,
                                    size_t index, char value);
int tsta_trace_block_store_mark_dirty(tsta_trace_block_store_t *store,
                                      size_t chunk_index);

void *tsta_aligned_malloc(size_t size, int base);
void tsta_aligned_free(void *buffer);

#endif
