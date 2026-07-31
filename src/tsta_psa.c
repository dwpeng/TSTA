#define _GNU_SOURCE

#include "tsta_common.h"
#include "tsta_psa_core.h"
#include "tsta_threadpool.h"

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct tsta_psa_aligner {
  tsta_config config;
  tsta_psa_state state;
  packed_sequence_array_t packed_sequences;
  char* packed_sequence_slots[2];
  tsta_psa_runtime_buffers buffers;
  tsta_threadpool_t* pool; /* created lazily, reused across aligns */
};

typedef struct psa_block_task {
  int diagonal_index;
  int diagonal_block_offset;
  tsta_psa_state* state;
} psa_block_task;

/* ── Per-thread workspace (threaded driver only) ─────────────────────── */

static pthread_key_t psa_workspace_key;
static pthread_once_t psa_workspace_once = PTHREAD_ONCE_INIT;
static int psa_workspace_key_ready = 0;

static void
psa_workspace_destroy(void* ptr)
{
  tsta_psa_workspace* ws = (tsta_psa_workspace*)ptr;
  if (!ws)
    return;
  tsta_psa_workspace_reset(ws);
  free(ws);
}

static void
psa_workspace_make_key(void)
{
  if (pthread_key_create(&psa_workspace_key, psa_workspace_destroy) == 0)
    psa_workspace_key_ready = 1;
}

static tsta_psa_workspace*
psa_workspace_acquire(int lane_len, int block_width)
{
  tsta_psa_workspace* ws;

  pthread_once(&psa_workspace_once, psa_workspace_make_key);
  if (!psa_workspace_key_ready)
    return NULL;

  ws = (tsta_psa_workspace*)pthread_getspecific(psa_workspace_key);
  if (!ws) {
    ws = (tsta_psa_workspace*)calloc(1, sizeof(*ws));
    if (!ws)
      return NULL;
    if (pthread_setspecific(psa_workspace_key, ws) != 0) {
      free(ws);
      return NULL;
    }
  }

  if (tsta_psa_workspace_resize(ws, lane_len, block_width) != 0)
    return NULL;
  return ws;
}

/* ── Block alignment task (dispatched by the threadpool) ─────────────── */

static void
tsta_psa_block_alignment(void* p)
{
  psa_block_task* task = (psa_block_task*)p;
  tsta_psa_state* state = task->state;
  int block_i = task->diagonal_index;
  int block_l = task->diagonal_block_offset;
  tsta_psa_workspace* ws = psa_workspace_acquire(state->L, block);
  tsta_psa_workspace local_ws;
  int use_local = 0;

  if (!ws) {
    memset(&local_ws, 0, sizeof(local_ws));
    if (tsta_psa_workspace_resize(&local_ws, state->L, block) != 0)
      return;
    ws = &local_ws;
    use_local = 1;
  }

  tsta_psa_block(state, ws, block_i, block_l);

  if (use_local)
    tsta_psa_workspace_reset(&local_ws);
}

/* ── Aligner create / destroy ────────────────────────────────────────── */

static void
tsta_psa_apply_config(tsta_psa_aligner* aligner, const tsta_config* config)
{
  aligner->config = config ? *config : tsta_config_make_default();
  tsta_init_psa_state(&aligner->state, &aligner->config, block);
}

tsta_psa_aligner*
tsta_psa_aligner_create(const tsta_config* config)
{
  tsta_psa_aligner* aligner = (tsta_psa_aligner*)calloc(1, sizeof(*aligner));
  if (!aligner)
    return NULL;
  tsta_psa_apply_config(aligner, config);
  return aligner;
}

void
tsta_psa_aligner_destroy(tsta_psa_aligner* aligner)
{
  if (!aligner)
    return;
  if (aligner->pool) {
    tsta_threadpool_destroy(aligner->pool);
    aligner->pool = NULL;
  }
  tsta_psa_runtime_buffers_release(&aligner->buffers, &aligner->state, 1);
  tsta_packed_sequence_array_reset_malloc(&aligner->packed_sequences);
  free(aligner);
}

/* ── Core align ──────────────────────────────────────────────────────── */

static int
tsta_psa_aligner_align_internal(tsta_psa_aligner* aligner,
                                const char* sequence1,
                                int sequence1_length,
                                const char* sequence2,
                                int sequence2_length,
                                tsta_psa_result_t* result)
{
  unsigned int tsl;
  int j = 0;
  tsta_psa_state* state;

  if (!aligner)
    return -1;
  state = &aligner->state;
  if (tsta_psa_align_setup(
          state, &aligner->config, sequence1, sequence1_length, sequence2,
          sequence2_length, &aligner->packed_sequences,
          aligner->packed_sequence_slots, &aligner->buffers, result)
      != 0)
    return -1;

  if (!aligner->pool) {
    aligner->pool = tsta_threadpool_create(
        aligner->config.threads > 0 ? aligner->config.threads : 10, 100,
        sizeof(psa_block_task));
    if (!aligner->pool)
      goto cleanup;
  }

  tsl = (unsigned int)((state->length[0] + state->length[1]) / state->L - 1);
  tsta_psa_blockmatrix_init(state);
  for (unsigned int i = 0; i < tsl; i++) {
    if ((int)i <= state->fmaxtag)
      j++;
    else if ((int)i > state->lmaxtag)
      j--;
    for (int l = 0; l < j; l++) {
      psa_block_task task = {
        .diagonal_index = (int)i,
        .diagonal_block_offset = l,
        .state = state,
      };
      if (tsta_threadpool_submit(aligner->pool, tsta_psa_block_alignment,
                                 &task, sizeof(task))
          != 0)
        tsta_psa_block_alignment(&task);
    }
    tsta_threadpool_wait(aligner->pool);
  }

  return tsta_psa_align_finish(&aligner->buffers, state, result);

cleanup:
  tsta_psa_runtime_buffers_release(&aligner->buffers, state, 0);
  return -1;
}

/* ── Public API ──────────────────────────────────────────────────────── */

int
tsta_psa_aligner_align(tsta_psa_aligner* aligner,
                       const char* sequence1,
                       int sequence1_length,
                       const char* sequence2,
                       int sequence2_length,
                       tsta_psa_result_t* result)
{
  if (!result)
    return -1;
  return tsta_psa_aligner_align_internal(aligner, sequence1, sequence1_length,
                                         sequence2, sequence2_length, result);
}

int
tsta_psa_align(const char* sequence1,
               int sequence1_length,
               const char* sequence2,
               int sequence2_length,
               const tsta_config* config,
               tsta_psa_result_t* result)
{
  tsta_config effective_config;
  size_t longest_sequence_length;
  tsta_psa_aligner* aligner;

  effective_config = config ? *config : tsta_config_make_default();
  longest_sequence_length =
      (size_t)(sequence1_length > sequence2_length ? sequence1_length
                                                   : sequence2_length);
  tsta_config_apply_length_fallback(&effective_config, longest_sequence_length,
                                    2);

  aligner = tsta_psa_aligner_create(&effective_config);
  if (!aligner)
    return -1;

  int status = tsta_psa_aligner_align(aligner, sequence1, sequence1_length,
                                      sequence2, sequence2_length, result);
  tsta_psa_aligner_destroy(aligner);
  return status;
}
