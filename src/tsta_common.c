
#define _GNU_SOURCE

#include "tsta_common.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#if defined(TSTA_HAVE_ZLIB)
#include <zlib.h>
#endif

struct tsta_trace_block_store {
  size_t length;
  size_t chunk_size;
  size_t chunk_count;
  size_t active_chunk;
  char *buffer;
  char **chunks;
  unsigned char *chunk_dirty;
  char *path;
  int spilled;
  int compress;
};

static int tsta_trace_block_store_write_slot(tsta_trace_block_store_t *store,
                                             size_t chunk_index,
                                             const char *chunk) {
  FILE *fp;
  uint32_t stored_length;
  uint32_t flags;
  size_t slot_offset;

  if (!store || !store->path || !chunk || chunk_index >= store->chunk_count) {
    return -1;
  }

  fp = fopen(store->path, "r+b");
  if (!fp) {
    return -1;
  }

  slot_offset = chunk_index * (sizeof(uint32_t) * 2 + store->chunk_size);
  if (fseeko(fp, (off_t)slot_offset, SEEK_SET) != 0) {
    fclose(fp);
    return -1;
  }

  stored_length = (uint32_t)store->chunk_size;
  flags = 0;

#if defined(TSTA_HAVE_ZLIB)
  if (store->compress) {
    uLongf compressed_capacity = compressBound((uLong)store->chunk_size);
    unsigned char *compressed_buffer =
        (unsigned char *)malloc((size_t)compressed_capacity);

    if (compressed_buffer) {
      uLongf compressed_length = compressed_capacity;

      if (compress2(compressed_buffer, &compressed_length, (const Bytef *)chunk,
                    (uLong)store->chunk_size, Z_BEST_SPEED) == Z_OK &&
          compressed_length < store->chunk_size) {
        stored_length = (uint32_t)compressed_length;
        flags = 1;
        if (fwrite(&stored_length, sizeof(stored_length), 1, fp) != 1 ||
            fwrite(&flags, sizeof(flags), 1, fp) != 1 ||
            fwrite(compressed_buffer, 1, (size_t)compressed_length, fp) !=
                (size_t)compressed_length) {
          free(compressed_buffer);
          fclose(fp);
          return -1;
        }
        free(compressed_buffer);
        fclose(fp);
        return 0;
      }
    }
    free(compressed_buffer);
  }
#else
  (void)flags;
#endif

  if (fwrite(&stored_length, sizeof(stored_length), 1, fp) != 1 ||
      fwrite(&flags, sizeof(flags), 1, fp) != 1 ||
      fwrite(chunk, 1, store->chunk_size, fp) != store->chunk_size) {
    fclose(fp);
    return -1;
  }

  fclose(fp);
  return 0;
}

static char *tsta_trace_block_store_load_slot(tsta_trace_block_store_t *store,
                                              size_t chunk_index) {
  FILE *fp;
  uint32_t stored_length = 0;
  uint32_t flags = 0;
  char *chunk;

  if (!store || !store->path || chunk_index >= store->chunk_count) {
    return NULL;
  }

  chunk = (char *)mm_malloc(store->chunk_size * sizeof(char));
  if (!chunk) {
    return NULL;
  }

  fp = fopen(store->path, "rb");
  if (!fp) {
    mm_free(chunk);
    return NULL;
  }

  if (fseeko(fp,
             (off_t)(chunk_index * (sizeof(uint32_t) * 2 + store->chunk_size)),
             SEEK_SET) != 0 ||
      fread(&stored_length, sizeof(stored_length), 1, fp) != 1 ||
      fread(&flags, sizeof(flags), 1, fp) != 1) {
    fclose(fp);
    mm_free(chunk);
    return NULL;
  }

  if (flags & 1U) {
#if defined(TSTA_HAVE_ZLIB)
    unsigned char *compressed_buffer;
    uLongf output_length = (uLongf)store->chunk_size;

    if (stored_length == 0 || stored_length > store->chunk_size) {
      fclose(fp);
      mm_free(chunk);
      return NULL;
    }

    compressed_buffer = (unsigned char *)malloc((size_t)stored_length);
    if (!compressed_buffer) {
      fclose(fp);
      mm_free(chunk);
      return NULL;
    }

    if (fread(compressed_buffer, 1, (size_t)stored_length, fp) !=
            (size_t)stored_length ||
        uncompress((Bytef *)chunk, &output_length, compressed_buffer,
                   (uLong)stored_length) != Z_OK ||
        output_length != store->chunk_size) {
      free(compressed_buffer);
      fclose(fp);
      mm_free(chunk);
      return NULL;
    }
    free(compressed_buffer);
#else
    fclose(fp);
    mm_free(chunk);
    return NULL;
#endif
  } else {
    if (stored_length != store->chunk_size ||
        fread(chunk, 1, store->chunk_size, fp) != store->chunk_size) {
      fclose(fp);
      mm_free(chunk);
      return NULL;
    }
  }

  fclose(fp);
  return chunk;
}

void tsta_trace_block_store_destroy(tsta_trace_block_store_t *store) {
  if (!store) {
    return;
  }

  if (store->buffer) {
    mm_free(store->buffer);
    store->buffer = NULL;
  } else if (store->chunks) {
    for (size_t i = 0; i < store->chunk_count; i++) {
      if (store->chunks[i]) {
        mm_free(store->chunks[i]);
        store->chunks[i] = NULL;
      }
    }
  }

  free(store->chunks);
  free(store->chunk_dirty);
  if (store->path) {
    unlink(store->path);
    free(store->path);
  }
  free(store);
}

tsta_trace_block_store_t *
tsta_trace_block_store_create(size_t length, size_t chunk_size, int compress) {
  tsta_trace_block_store_t *store;
  size_t buffer_size;

  if (length == 0 || chunk_size == 0) {
    return NULL;
  }

  store = (tsta_trace_block_store_t *)calloc(1, sizeof(*store));
  if (!store) {
    return NULL;
  }

  store->length = length;
  store->chunk_size = chunk_size;
  store->chunk_count = (length + chunk_size - 1) / chunk_size;
  store->active_chunk = (size_t)-1;
  store->compress = compress ? 1 : 0;

  store->chunks = (char **)calloc(store->chunk_count, sizeof(char *));
  store->chunk_dirty =
      (unsigned char *)calloc(store->chunk_count, sizeof(unsigned char));
  if (!store->chunks || !store->chunk_dirty) {
    tsta_trace_block_store_destroy(store);
    return NULL;
  }

  buffer_size = store->chunk_count * store->chunk_size;
  store->buffer = (char *)mm_malloc(buffer_size * sizeof(char));
  if (!store->buffer) {
    tsta_trace_block_store_destroy(store);
    return NULL;
  }
  memset(store->buffer, 0, buffer_size);
  for (size_t i = 0; i < store->chunk_count; i++) {
    store->chunks[i] = store->buffer + i * store->chunk_size;
  }

  return store;
}

int tsta_trace_block_store_spill(tsta_trace_block_store_t *store) {
  size_t *loaded_indices = NULL;
  size_t loaded_count = 0;
  int status = -1;

  if (!store) {
    return -1;
  }

  if (!store->path) {
    char template_path[] = "/tmp/tsta-trace-XXXXXX";
    int fd = mkstemp(template_path);
    size_t total_size =
        store->chunk_count * (sizeof(uint32_t) * 2 + store->chunk_size);

    if (fd < 0) {
      return -1;
    }
    if (ftruncate(fd, (off_t)total_size) != 0) {
      close(fd);
      unlink(template_path);
      return -1;
    }
    close(fd);

    store->path = (char *)malloc(strlen(template_path) + 1);
    if (!store->path) {
      unlink(template_path);
      return -1;
    }
    strcpy(store->path, template_path);
  }

  if (store->buffer) {
    for (size_t i = 0; i < store->chunk_count; i++) {
      if (tsta_trace_block_store_write_slot(
              store, i, store->buffer + i * store->chunk_size) != 0) {
        return -1;
      }
    }

    mm_free(store->buffer);
    store->buffer = NULL;
    for (size_t i = 0; i < store->chunk_count; i++) {
      store->chunks[i] = NULL;
      store->chunk_dirty[i] = 0;
    }
    store->active_chunk = (size_t)-1;
    store->spilled = 1;
    return 0;
  }

  loaded_indices = (size_t *)malloc(store->chunk_count * sizeof(size_t));
  if (!loaded_indices) {
    return -1;
  }

  for (size_t i = 0; i < store->chunk_count; i++) {
    if (store->chunks[i]) {
      loaded_indices[loaded_count++] = i;
    }
  }

  for (size_t i = 0; i < loaded_count; i++) {
    size_t chunk_index = loaded_indices[i];

    if (tsta_trace_block_store_write_slot(store, chunk_index,
                                          store->chunks[chunk_index]) != 0) {
      goto cleanup;
    }
  }

  for (size_t i = 0; i < loaded_count; i++) {
    size_t chunk_index = loaded_indices[i];

    mm_free(store->chunks[chunk_index]);
    store->chunks[chunk_index] = NULL;
    store->chunk_dirty[chunk_index] = 0;
  }
  store->active_chunk = (size_t)-1;
  store->spilled = 1;
  status = 0;

cleanup:
  free(loaded_indices);
  return status;
}

char *tsta_trace_block_store_chunk(tsta_trace_block_store_t *store,
                                   size_t chunk_index) {
  char *new_chunk;

  if (!store || chunk_index >= store->chunk_count) {
    return NULL;
  }

  if (store->buffer) {
    return store->buffer + chunk_index * store->chunk_size;
  }

  if (store->chunks[chunk_index]) {
    store->active_chunk = chunk_index;
    return store->chunks[chunk_index];
  }

  new_chunk = tsta_trace_block_store_load_slot(store, chunk_index);
  if (!new_chunk) {
    return NULL;
  }

  if (store->active_chunk != (size_t)-1 && store->active_chunk != chunk_index &&
      store->chunks[store->active_chunk]) {
    size_t old_chunk = store->active_chunk;

    if (store->chunk_dirty[old_chunk]) {
      if (tsta_trace_block_store_write_slot(store, old_chunk,
                                            store->chunks[old_chunk]) != 0) {
        mm_free(new_chunk);
        return NULL;
      }
    }
    mm_free(store->chunks[old_chunk]);
    store->chunks[old_chunk] = NULL;
    store->chunk_dirty[old_chunk] = 0;
  }

  store->chunks[chunk_index] = new_chunk;
  store->active_chunk = chunk_index;
  return new_chunk;
}

char tsta_trace_block_store_get_byte(tsta_trace_block_store_t *store,
                                     size_t index) {
  size_t chunk_index;
  size_t offset;
  char *chunk;

  if (!store || index >= store->length) {
    return 0;
  }

  chunk_index = index / store->chunk_size;
  offset = index % store->chunk_size;
  chunk = tsta_trace_block_store_chunk(store, chunk_index);
  if (!chunk) {
    return 0;
  }

  return chunk[offset];
}

int tsta_trace_block_store_set_byte(tsta_trace_block_store_t *store,
                                    size_t index, char value) {
  size_t chunk_index;
  size_t offset;
  char *chunk;

  if (!store || index >= store->length) {
    return -1;
  }

  chunk_index = index / store->chunk_size;
  offset = index % store->chunk_size;
  chunk = tsta_trace_block_store_chunk(store, chunk_index);
  if (!chunk) {
    return -1;
  }

  chunk[offset] = value;
  if (store->spilled) {
    store->chunk_dirty[chunk_index] = 1;
  }
  return 0;
}

int tsta_trace_block_store_mark_dirty(tsta_trace_block_store_t *store,
                                      size_t chunk_index) {
  if (!store || chunk_index >= store->chunk_count) {
    return -1;
  }

  if (store->spilled) {
    store->chunk_dirty[chunk_index] = 1;
  }
  return 0;
}

void tsta_config_default(tsta_config *config) {
  if (!config) {
    return;
  }

  config->match = 2;
  config->mismatch = -5;
  config->gap_extend = -2;
  config->gap_open = -4;
  config->block_size = 10;
  config->threads = 10;
  config->trace_enabled = 0;
  config->msa_trace_dump_threshold = 64ULL * 1024ULL * 1024ULL;
  config->msa_trace_dump_compress = 1;
}

tsta_config tsta_config_make_default(void) {
  tsta_config config;

  tsta_config_default(&config);
  return config;
}

void tsta_psa_result_init(tsta_psa_result_t *result) {
  if (!result) {
    return;
  }

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

tsta_psa_result_t tsta_psa_result_make(void) {
  tsta_psa_result_t result;

  tsta_psa_result_init(&result);
  return result;
}

void tsta_psa_result_free(tsta_psa_result_t *result) {
  if (!result) {
    return;
  }

  free(result->aln[0]);
  free(result->aln[1]);
  free(result->cigar);
  tsta_psa_result_init(result);
}

void tsta_msa_result_init(tsta_msa_result_t *result) {
  if (!result) {
    return;
  }

  result->score = 0;
  result->aln = NULL;
  result->sequence_count = 0;
  result->aln_length = 0;
}

tsta_msa_result_t tsta_msa_result_make(void) {
  tsta_msa_result_t result;

  tsta_msa_result_init(&result);
  return result;
}

void tsta_msa_result_free(tsta_msa_result_t *result) {
  if (!result) {
    return;
  }

  if (result->aln) {
    for (size_t i = 0; i < result->sequence_count; i++) {
      free(result->aln[i]);
    }
  }
  free(result->aln);
  tsta_msa_result_init(result);
}

static void tsta_apply_base_config(const tsta_config *config, int *match,
                                   int *mismatch, int *gap_extend,
                                   int *gap_open, int *block_size, int *threads,
                                   int *trace_enabled) {
  tsta_config effective;

  if (config) {
    effective = *config;
  } else {
    tsta_config_default(&effective);
  }

  *match = effective.match;
  *mismatch = effective.mismatch;
  *gap_extend = effective.gap_extend;
  *gap_open = effective.gap_open;
  *block_size = effective.block_size > 0 ? effective.block_size : 10;
  *threads = effective.threads > 0 ? effective.threads : 10;
  if (trace_enabled) {
    *trace_enabled = effective.trace_enabled ? 1 : 0;
  }
}

void tsta_init_psa_state(tsta_psa_state *state, const tsta_config *config,
                         int simd_block) {
  int threads = 0;
  int trace_enabled = 0;

  tsta_apply_base_config(config, &state->M, &state->X, &state->E, &state->O,
                         &state->bS, &threads, &trace_enabled);
  state->L = state->bS * simd_block;
  state->B = simd_block;
  state->W = (state->L + state->B - 1) / state->B;
  state->trace_enabled = trace_enabled;
  (void)threads;
}

void tsta_init_msa_state(tsta_msa_state *state, const tsta_config *config,
                         int simd_block) {
  int threads = 0;

  tsta_apply_base_config(config, &state->M, &state->X, &state->E, &state->O,
                         &state->bS, &threads, NULL);
  state->L = state->bS * simd_block;
  state->B = simd_block;
  state->W = (state->L + state->B - 1) / state->B;
  state->trace_dump_enabled = 0;
  state->trace_dump_compress = 0;
  (void)threads;
}

void tsta_config_apply_length_fallback(tsta_config *config,
                                       size_t longest_sequence_length,
                                       size_t sequence_count) {
  size_t safe_block_size;
  size_t safe_threads;
  size_t lane_width;

  if (!config) {
    return;
  }

  if (config->block_size <= 0) {
    config->block_size = 10;
  }
  if (config->threads <= 0) {
    config->threads = 10;
  }

  lane_width = (size_t)block;
  if (lane_width == 0) {
    lane_width = 16;
  }

  if (longest_sequence_length > 0) {
    safe_block_size = (longest_sequence_length + lane_width - 1) / lane_width;
    if (safe_block_size == 0) {
      safe_block_size = 1;
    }
    if (config->block_size > (int)safe_block_size) {
      config->block_size = (int)safe_block_size;
    }

    safe_threads = safe_block_size;
    if (config->threads > (int)safe_threads) {
      config->threads = (int)safe_threads;
    }
  }

  if (sequence_count > 0 && config->threads > (int)sequence_count) {
    config->threads = (int)sequence_count;
  }

  if (config->block_size < 1) {
    config->block_size = 1;
  }
  if (config->threads < 1) {
    config->threads = 1;
  }
}

int tsta_packed_sequence_array_reserve_slot(packed_sequence_array_t *array,
                                            size_t slot,
                                            size_t initial_capacity) {
  size_t old_length;
  size_t required_capacity;

  if (!array) {
    return -1;
  }
  if (slot < array->length) {
    return 0;
  }

  old_length = array->length;
  required_capacity = slot + 1;
  if (initial_capacity > required_capacity) {
    required_capacity = initial_capacity;
  }

  if (tsta_packed_sequence_array_reserve(array, required_capacity) != 0) {
    return -1;
  }

  for (size_t i = old_length; i <= slot; i++) {
    array->data[i].packed_seq = NULL;
    array->data[i].len = 0;
  }

  array->length = slot + 1;
  return 0;
}

void tsta_packed_sequence_array_reset_malloc(packed_sequence_array_t *array) {
  if (!array) {
    return;
  }

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

void tsta_packed_sequence_array_reset_mm(packed_sequence_array_t *array) {
  if (!array) {
    return;
  }

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

void *tsta_aligned_malloc(size_t size, int base) {
  uint8_t *p;
  uint8_t *q;

  p = (uint8_t *)malloc(size + (size_t)base);
  if (p == NULL) {
    return NULL;
  }
  q = (uint8_t *)(((unsigned long long)(p + base)) &
                  (~(((unsigned long long)base) - 1)));
  *(q - 1) = (uint8_t)(q - p);
  return q;
}

void tsta_aligned_free(void *buffer) {
  uint8_t *p;
  uint8_t *q;

  if (!buffer) {
    return;
  }
  q = (uint8_t *)buffer;
  p = q - *(q - 1);
  free(p);
}
