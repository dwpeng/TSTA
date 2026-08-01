#include "tsta_common.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Sequence view ──────────────────────────────────────────────────── */

int
tsta_sequence_view_make(tsta_sequence_view_t* view,
                        const char* sequence,
                        int length)
{
  if (!view || !sequence || length < 0)
    return -1;
  view->sequence = sequence;
  view->length = length;
  return 0;
}

/* ── Trace block store (memory-only) ────────────────────────────────── */

tsta_trace_block_store_t*
tsta_trace_block_store_create(size_t length, size_t chunk_size)
{
  tsta_trace_block_store_t* store;

  if (length == 0 || chunk_size == 0)
    return NULL;

  store = (tsta_trace_block_store_t*)calloc(1, sizeof(*store));
  if (!store)
    return NULL;

  store->length = length;
  store->chunk_size = chunk_size;
  store->chunk_count = (length + chunk_size - 1) / chunk_size;

  store->buffer = (uint8_t*)mm_malloc(store->chunk_count * chunk_size);
  if (!store->buffer) {
    free(store);
    return NULL;
  }
  memset(store->buffer, 0, store->chunk_count * chunk_size);

  return store;
}

void
tsta_trace_block_store_destroy(tsta_trace_block_store_t* store)
{
  if (!store)
    return;
  mm_free(store->buffer);
  free(store);
}

/* Grow the store buffer in place so it can be reused across alignment
 * steps (avoids per-step create/destroy churn). No-op if already large
 * enough; on growth the new tail is zeroed. */
void
tsta_trace_block_store_ensure(tsta_trace_block_store_t* store,
                              size_t length,
                              size_t chunk_size)
{
  uint8_t* new_buffer;
  size_t new_chunk_count;

  if (!store || length == 0 || chunk_size == 0)
    return;
  if (store->buffer && store->length >= length
      && store->chunk_size == chunk_size)
    return;

  new_chunk_count = (length + chunk_size - 1) / chunk_size;
  new_buffer = (uint8_t*)mm_malloc(new_chunk_count * chunk_size);
  if (!new_buffer)
    return;

  mm_free(store->buffer);
  store->buffer = new_buffer;
  store->length = length;
  store->chunk_size = chunk_size;
  store->chunk_count = new_chunk_count;
  memset(store->buffer, 0, new_chunk_count * chunk_size);
}

/* ── Config ─────────────────────────────────────────────────────────── */

void
tsta_config_default(tsta_config* config)
{
  if (!config)
    return;
  config->match = 2;
  config->mismatch = -5;
  config->gap_extend = -2;
  config->gap_open = -4;
  config->block_size = 10;
  config->threads = 10;
}

tsta_config
tsta_config_make_default(void)
{
  tsta_config config;
  tsta_config_default(&config);
  return config;
}

/* ── Result lifecycle ───────────────────────────────────────────────── */

void
tsta_psa_result_init(tsta_psa_result_t* result)
{
  if (!result)
    return;
  result->score = 0;
  result->del_count = 0;
  result->ins_count = 0;
  result->match_count = 0;
  result->mismatch_count = 0;
  result->aln[0] = NULL;
  result->aln[1] = NULL;
  result->aln_length = 0;
  result->cigar = NULL;
}

tsta_psa_result_t
tsta_psa_result_make(void)
{
  tsta_psa_result_t result;
  tsta_psa_result_init(&result);
  return result;
}

void
tsta_psa_result_free(tsta_psa_result_t* result)
{
  if (!result)
    return;
  free(result->aln[0]);
  free(result->aln[1]);
  free(result->cigar);
  tsta_psa_result_init(result);
}

void
tsta_msa_result_init(tsta_msa_result_t* result)
{
  if (!result)
    return;
  result->score = 0;
  result->aln = NULL;
  result->sequence_count = 0;
  result->aln_length = 0;
}

tsta_msa_result_t
tsta_msa_result_make(void)
{
  tsta_msa_result_t result;
  tsta_msa_result_init(&result);
  return result;
}

void
tsta_msa_result_free(tsta_msa_result_t* result)
{
  if (!result)
    return;
  if (result->aln) {
    for (size_t i = 0; i < result->sequence_count; i++)
      free(result->aln[i]);
  }
  free(result->aln);
  tsta_msa_result_init(result);
}

/* ── State init ─────────────────────────────────────────────────────── */

static void
tsta_apply_base_config(const tsta_config* config,
                       int* match,
                       int* mismatch,
                       int* gap_extend,
                       int* gap_open,
                       int* block_size,
                       int* threads)
{
  tsta_config effective;

  if (config)
    effective = *config;
  else
    tsta_config_default(&effective);

  *match = effective.match;
  *mismatch = effective.mismatch;
  *gap_extend = effective.gap_extend;
  *gap_open = effective.gap_open;
  *block_size = effective.block_size > 0 ? effective.block_size : 10;
  *threads = effective.threads > 0 ? effective.threads : 10;
}

void
tsta_init_psa_state(tsta_psa_state* state,
                    const tsta_config* config,
                    int simd_block)
{
  int threads = 0;

  tsta_apply_base_config(config, &state->M, &state->X, &state->E, &state->O,
                         &state->bS, &threads);
  state->L = state->bS * simd_block;
  state->B = simd_block;
  state->W = (state->L + state->B - 1) / state->B;
  (void)threads;
}

void
tsta_init_msa_state(tsta_msa_state* state,
                    const tsta_config* config,
                    int simd_block)
{
  int threads = 0;

  tsta_apply_base_config(config, &state->M, &state->X, &state->E, &state->O,
                         &state->bS, &threads);
  state->L = state->bS * simd_block;
  state->B = simd_block;
  state->W = (state->L + state->B - 1) / state->B;
  (void)threads;
}

void
tsta_config_apply_length_fallback(tsta_config* config,
                                  size_t longest_sequence_length,
                                  size_t sequence_count)
{
  size_t safe_block_size;
  size_t safe_threads;
  size_t lane_width;

  if (!config)
    return;

  if (config->block_size <= 0)
    config->block_size = 10;
  if (config->threads <= 0)
    config->threads = 10;

  lane_width = (size_t)block;
  if (lane_width == 0)
    lane_width = 16;

  if (longest_sequence_length > 0) {
    safe_block_size = (longest_sequence_length + lane_width - 1) / lane_width;
    if (safe_block_size == 0)
      safe_block_size = 1;
    if (config->block_size > (int)safe_block_size)
      config->block_size = (int)safe_block_size;

    safe_threads = safe_block_size;
    if (config->threads > (int)safe_threads)
      config->threads = (int)safe_threads;
  }

  if (sequence_count > 0 && config->threads > (int)sequence_count)
    config->threads = (int)sequence_count;

  if (config->block_size < 1)
    config->block_size = 1;
  if (config->threads < 1)
    config->threads = 1;
}

/* ── Packed sequence helpers ────────────────────────────────────────── */

int
tsta_packed_sequence_array_reserve_slot(packed_sequence_array_t* array,
                                        size_t slot,
                                        size_t initial_capacity)
{
  size_t old_length;
  size_t required_capacity;

  if (!array)
    return -1;
  if (slot < array->length)
    return 0;

  old_length = array->length;
  required_capacity = slot + 1;
  if (initial_capacity > required_capacity)
    required_capacity = initial_capacity;

  if (tsta_packed_sequence_array_reserve(array, required_capacity) != 0)
    return -1;

  for (size_t i = old_length; i <= slot; i++) {
    array->data[i].packed_seq = NULL;
    array->data[i].len = 0;
  }
  array->length = slot + 1;
  return 0;
}

void
tsta_packed_sequence_array_reset_malloc(packed_sequence_array_t* array)
{
  if (!array)
    return;
  for (size_t i = 0; i < array->length; i++) {
    free(array->data[i].packed_seq);
    array->data[i].packed_seq = NULL;
    array->data[i].len = 0;
  }
  free(array->data);
  array->data = NULL;
  array->length = 0;
  array->capacity = 0;
}

void
tsta_packed_sequence_array_reset_mm(packed_sequence_array_t* array)
{
  if (!array)
    return;
  for (size_t i = 0; i < array->length; i++) {
    mm_free(array->data[i].packed_seq);
    array->data[i].packed_seq = NULL;
    array->data[i].len = 0;
  }
  free(array->data);
  array->data = NULL;
  array->length = 0;
  array->capacity = 0;
}

/* ── Aligned alloc ──────────────────────────────────────────────────── */

void*
tsta_aligned_malloc(size_t size, int base)
{
  uint8_t* p;
  uint8_t* q;

  p = (uint8_t*)malloc(size + (size_t)base);
  if (p == NULL)
    return NULL;
  q = (uint8_t*)(((unsigned long long)(p + base))
                 & (~(((unsigned long long)base) - 1)));
  *(q - 1) = (uint8_t)(q - p);
  return q;
}

void
tsta_aligned_free(void* buffer)
{
  uint8_t* p;
  uint8_t* q;

  if (!buffer)
    return;
  q = (uint8_t*)buffer;
  p = q - *(q - 1);
  free(p);
}
