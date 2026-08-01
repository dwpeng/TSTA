#ifndef TSTA_H
#define TSTA_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) || defined(__CYGWIN__)
#ifdef TSTA_BUILD_DLL
#define TSTA_API __declspec(dllexport)
#else
#define TSTA_API __declspec(dllimport)
#endif
#else
#define TSTA_API __attribute__((visibility("default")))
#endif

/* ── Configuration ─────────────────────────────────────────────────── */

/* Traceback is always produced (tsta_psa_result_t / tsta_msa_result_t always
 * include the aligned sequences), so there is no trace_enabled switch. */
typedef struct tsta_config {
  int match;
  int mismatch;
  int gap_extend;
  int gap_open;
  int block_size;
  int threads;
} tsta_config;

/* ── PSA result ────────────────────────────────────────────────────── */

typedef struct tsta_psa_result_t {
  int score;
  int del_count;
  int ins_count;
  int match_count;
  int mismatch_count;
  char* aln[2];
  size_t aln_length;
  char* cigar;
} tsta_psa_result_t;

/* ── MSA result ────────────────────────────────────────────────────── */

typedef struct tsta_msa_result_t {
  int score;
  char** aln;
  size_t sequence_count;
  size_t aln_length;
} tsta_msa_result_t;

/* ── Opaque aligner types ──────────────────────────────────────────── */

typedef struct tsta_psa_aligner tsta_psa_aligner;
typedef struct tsta_msa_aligner tsta_msa_aligner;

/* ── Configuration ─────────────────────────────────────────────────── */

TSTA_API void tsta_config_default(tsta_config* config);
TSTA_API tsta_config tsta_config_make_default(void);

/* ── PSA result lifecycle ──────────────────────────────────────────── */

TSTA_API void tsta_psa_result_init(tsta_psa_result_t* result);
TSTA_API tsta_psa_result_t tsta_psa_result_make(void);
TSTA_API void tsta_psa_result_free(tsta_psa_result_t* result);

/* ── MSA result lifecycle ──────────────────────────────────────────── */

TSTA_API void tsta_msa_result_init(tsta_msa_result_t* result);
TSTA_API tsta_msa_result_t tsta_msa_result_make(void);
TSTA_API void tsta_msa_result_free(tsta_msa_result_t* result);

/* ── PSA aligner ───────────────────────────────────────────────────── */

TSTA_API tsta_psa_aligner* tsta_psa_aligner_create(const tsta_config* config);
TSTA_API void tsta_psa_aligner_destroy(tsta_psa_aligner* aligner);
TSTA_API int tsta_psa_aligner_align(tsta_psa_aligner* aligner,
                                    const char* sequence1,
                                    int sequence1_length,
                                    const char* sequence2,
                                    int sequence2_length,
                                    tsta_psa_result_t* result);

/* ── MSA aligner ───────────────────────────────────────────────────── */

TSTA_API tsta_msa_aligner* tsta_msa_aligner_create(const tsta_config* config);
TSTA_API void tsta_msa_aligner_destroy(tsta_msa_aligner* aligner);
TSTA_API int tsta_msa_aligner_begin(tsta_msa_aligner* aligner,
                                    const char* initial_sequence,
                                    int initial_sequence_length);
TSTA_API int tsta_msa_aligner_add(tsta_msa_aligner* aligner,
                                  const char* sequence,
                                  int sequence_length);
TSTA_API int tsta_msa_aligner_get_result(tsta_msa_aligner* aligner,
                                         tsta_msa_result_t* result);
TSTA_API int tsta_msa_aligner_align(tsta_msa_aligner* aligner,
                                    const char* const* sequences,
                                    const int* sequence_lengths,
                                    size_t sequence_count,
                                    tsta_msa_result_t* result);

/* ── One-shot alignment ────────────────────────────────────────────── */

TSTA_API int tsta_psa_align(const char* sequence1,
                            int sequence1_length,
                            const char* sequence2,
                            int sequence2_length,
                            const tsta_config* config,
                            tsta_psa_result_t* result);

TSTA_API int tsta_msa_align(const char* const* sequences,
                            const int* sequence_lengths,
                            size_t sequence_count,
                            const tsta_config* config,
                            tsta_msa_result_t* result);

#undef TSTA_API

#ifdef __cplusplus
}
#endif

#endif /* TSTA_H */
