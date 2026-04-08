#ifndef TSTA_API_H
#define TSTA_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Configuration for sequence alignment */
typedef struct tsta_config {
  int match;
  int mismatch;
  int gap_extend;
  int gap_open;
  int block_size;
  int threads;
  int trace_enabled;
  size_t msa_trace_dump_threshold;
  int msa_trace_dump_compress;
} tsta_config;

/* Result of pairwise sequence alignment */
typedef struct tsta_psa_result_t {
  int score;
  int del_count;
  int ins_count;
  int match_count;
  int mismatch_count;
  char *aln[2];
  size_t aln_length;
  char *cigar;
} tsta_psa_result_t;

/* Result of multiple sequence alignment */
typedef struct tsta_msa_result_t {
  int score;
  char **aln;
  size_t sequence_count;
  size_t aln_length;
} tsta_msa_result_t;

/* Opaque pointer types for aligner objects */
typedef struct tsta_psa_aligner tsta_psa_aligner;
typedef struct tsta_msa_aligner tsta_msa_aligner;

/* Configuration functions */
void tsta_config_default(tsta_config *config);
tsta_config tsta_config_make_default(void);

/* PSA result management */
void tsta_psa_result_init(tsta_psa_result_t *result);
tsta_psa_result_t tsta_psa_result_make(void);
void tsta_psa_result_free(tsta_psa_result_t *result);

/* MSA result management */
void tsta_msa_result_init(tsta_msa_result_t *result);
tsta_msa_result_t tsta_msa_result_make(void);
void tsta_msa_result_free(tsta_msa_result_t *result);

/* PSA aligner API */
tsta_psa_aligner *tsta_psa_aligner_create(const tsta_config *config);
void tsta_psa_aligner_destroy(tsta_psa_aligner *aligner);
int tsta_psa_aligner_align(tsta_psa_aligner *aligner, const char *sequence1,
                           int sequence1_length, const char *sequence2,
                           int sequence2_length, tsta_psa_result_t *result);

/* MSA aligner API */
tsta_msa_aligner *tsta_msa_aligner_create(const tsta_config *config);
void tsta_msa_aligner_destroy(tsta_msa_aligner *aligner);
int tsta_msa_aligner_begin(tsta_msa_aligner *aligner,
                           const char *initial_sequence,
                           int initial_sequence_length);
int tsta_msa_aligner_add(tsta_msa_aligner *aligner, const char *sequence,
                         int sequence_length);
int tsta_msa_aligner_get_result(tsta_msa_aligner *aligner,
                                tsta_msa_result_t *result);
int tsta_msa_aligner_align(tsta_msa_aligner *aligner,
                           const char *const *sequences,
                           const int *sequence_lengths, size_t sequence_count,
                           tsta_msa_result_t *result);

/* One-shot alignment functions */
int tsta_psa_align(const char *sequence1, int sequence1_length,
                   const char *sequence2, int sequence2_length,
                   const tsta_config *config, tsta_psa_result_t *result);

int tsta_msa_align(const char *const *sequences, const int *sequence_lengths,
                   size_t sequence_count, const tsta_config *config,
                   tsta_msa_result_t *result);

#ifdef __cplusplus
}
#endif

#endif /* TSTA_API_H */