#ifndef TSTA_PSA_SIMD_H
#define TSTA_PSA_SIMD_H

#include "tsta.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(__CYGWIN__)
#define TSTA_PSA_SIMD_API
#else
#define TSTA_PSA_SIMD_API __attribute__((visibility("default")))
#endif

/* ── Standalone SIMD-only pairwise alignment module ─────────────────────
 *
 * Single-threaded pairwise global alignment accelerated by SIMD only.
 * No threads, no threadpool, no pthread, no TLS. Depends on nothing beyond
 * libc and the SIMD instruction set (same as the main library).
 *
 * Reuses the public types from <tsta.h>:
 *   - configuration: tsta_config (the `threads` field is ignored)
 *   - result:        tsta_psa_result_t
 *
 * `config == NULL` selects built-in defaults
 * (match=2, mismatch=-5, gap_extend=-2, gap_open=-4, block_size=10).
 * Results are bit-identical to the threaded tsta_psa_align for the same
 * configuration and inputs.
 * ─────────────────────────────────────────────────────────────────────── */

typedef struct tsta_psa_simd_aligner tsta_psa_simd_aligner;
/* ── Reusable aligner ───────────────────────────────────────────────── */

TSTA_PSA_SIMD_API tsta_psa_simd_aligner*
tsta_psa_simd_aligner_create(const tsta_config* config);
TSTA_PSA_SIMD_API void
tsta_psa_simd_aligner_destroy(tsta_psa_simd_aligner* aligner);
TSTA_PSA_SIMD_API int
tsta_psa_simd_aligner_align(tsta_psa_simd_aligner* aligner,
                            const char* sequence1,
                            int sequence1_length,
                            const char* sequence2,
                            int sequence2_length,
                            tsta_psa_result_t* result);

/* ── One-shot alignment ─────────────────────────────────────────────── */

TSTA_PSA_SIMD_API int tsta_psa_simd_align(const char* sequence1,
                                          int sequence1_length,
                                          const char* sequence2,
                                          int sequence2_length,
                                          const tsta_config* config,
                                          tsta_psa_result_t* result);

#undef TSTA_PSA_SIMD_API

#ifdef __cplusplus
}
#endif

#endif /* TSTA_PSA_SIMD_H */
