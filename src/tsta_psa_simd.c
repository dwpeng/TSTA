/*
 * tsta_psa_simd.c — Standalone SIMD-only pairwise sequence alignment.
 *
 * Single-threaded: no threadpool, no pthread calls, no TLS. Reuses the
 * shared SIMD core (tsta_psa_core.h) and the common helpers
 * (tsta_common.h) instead of reimplementing them; only the serial
 * anti-diagonal driver and the aligner wrapper are module-specific.
 *
 * Output is bit-identical to the threaded tsta_psa_align for the same
 * configuration and inputs.
 */

#include "tsta_psa_simd.h"

#include "tsta_common.h"
#include "tsta_psa_core.h"

#include <stdlib.h>
#include <string.h>

struct tsta_psa_simd_aligner {
  tsta_config config;
  tsta_psa_state state;
  packed_sequence_array_t packed_sequences;
  char* packed_sequence_slots[2];
  tsta_psa_runtime_buffers buffers;
  tsta_psa_workspace workspace;
};

/* ── Aligner lifecycle ───────────────────────────────────────────────── */

static void
tsta_psa_simd_apply_config(tsta_psa_simd_aligner* aligner,
                           const tsta_config* config)
{
  aligner->config = config ? *config : tsta_config_make_default();
  tsta_init_psa_state(&aligner->state, &aligner->config, block);
}

tsta_psa_simd_aligner*
tsta_psa_simd_aligner_create(const tsta_config* config)
{
  tsta_psa_simd_aligner* aligner =
      (tsta_psa_simd_aligner*)calloc(1, sizeof(*aligner));
  if (!aligner)
    return NULL;
  tsta_psa_simd_apply_config(aligner, config);
  return aligner;
}

void
tsta_psa_simd_aligner_destroy(tsta_psa_simd_aligner* aligner)
{
  if (!aligner)
    return;
  tsta_psa_runtime_buffers_release(&aligner->buffers, &aligner->state, 1);
  tsta_packed_sequence_array_reset_malloc(&aligner->packed_sequences);
  tsta_psa_workspace_reset(&aligner->workspace);
  free(aligner);
}

/* ── Serial anti-diagonal driver (threadpool replaced) ───────────────── */

static int
tsta_psa_simd_aligner_align_internal(tsta_psa_simd_aligner* aligner,
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
  if (tsta_psa_align_setup(state, &aligner->config, sequence1,
                           sequence1_length, sequence2, sequence2_length,
                           &aligner->packed_sequences,
                           aligner->packed_sequence_slots, &aligner->buffers)
      != 0)
    return -1;

  if (tsta_psa_workspace_resize(&aligner->workspace, state->L, block) != 0)
    goto cleanup;

  tsl = (unsigned int)((state->length[0] + state->length[1]) / state->L - 1);
  tsta_psa_blockmatrix_init(state);
  for (unsigned int i = 0; i < tsl; i++) {
    if ((int)i <= state->fmaxtag)
      j++;
    else if ((int)i > state->lmaxtag)
      j--;
    for (int l = 0; l < j; l++)
      tsta_psa_block(state, &aligner->workspace, (int)i, l);
  }

  return tsta_psa_align_finish(&aligner->buffers, state, result);

cleanup:
  tsta_psa_runtime_buffers_release(&aligner->buffers, state, 0);
  return -1;
}

/* ── Public API ──────────────────────────────────────────────────────── */

int
tsta_psa_simd_aligner_align(tsta_psa_simd_aligner* aligner,
                            const char* sequence1,
                            int sequence1_length,
                            const char* sequence2,
                            int sequence2_length,
                            tsta_psa_result_t* result)
{
  if (!result)
    return -1;
  return tsta_psa_simd_aligner_align_internal(aligner, sequence1,
                                              sequence1_length, sequence2,
                                              sequence2_length, result);
}

int
tsta_psa_simd_align(const char* sequence1,
                    int sequence1_length,
                    const char* sequence2,
                    int sequence2_length,
                    const tsta_config* config,
                    tsta_psa_result_t* result)
{
  tsta_config effective_config;
  tsta_psa_simd_aligner* aligner;

  if (!result)
    return -1;
  if (config)
    effective_config = *config;
  else
    tsta_config_default(&effective_config);
  tsta_psa_adjust_config(&effective_config);

  aligner = tsta_psa_simd_aligner_create(&effective_config);
  if (!aligner)
    return -1;

  int status =
      tsta_psa_simd_aligner_align(aligner, sequence1, sequence1_length,
                                  sequence2, sequence2_length, result);
  tsta_psa_simd_aligner_destroy(aligner);
  return status;
}
