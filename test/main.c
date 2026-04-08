#include "tsta.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

static void assert_msa_result_equal(const tsta_msa_result_t *expected,
                                    const tsta_msa_result_t *actual) {
  assert(expected != NULL);
  assert(actual != NULL);
  assert(expected->score == actual->score);
  assert(expected->sequence_count == actual->sequence_count);
  assert(expected->aln_length == actual->aln_length);
  for (size_t i = 0; i < expected->sequence_count; i++) {
    assert(expected->aln[i] != NULL);
    assert(actual->aln[i] != NULL);
    assert(strcmp(expected->aln[i], actual->aln[i]) == 0);
  }
}

static void fill_repeated_sequence(char *buffer, size_t repeat_count,
                                   const char *motif) {
  size_t motif_length = strlen(motif);

  for (size_t i = 0; i < repeat_count; i++) {
    memcpy(buffer + i * motif_length, motif, motif_length);
  }
  buffer[repeat_count * motif_length] = '\0';
}

int test_psa() {
  const char *seq1 = "ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
  const char *seq2 =
      "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();
  int status;

  status = tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                          &config, &result);
  if (status == 0) {
    printf("Alignment score: %d\n", result.score);
    printf("Aligned sequence 1: %s\n", result.aln[0]);
    printf("Aligned sequence 2: %s\n", result.aln[1]);
    assert(result.aln[0] != NULL);
    assert(result.aln[1] != NULL);
  }
  tsta_psa_result_free(&result);
  return status;
}

int test_msa() {
  const char *seqs[] = {"ACGTAGCTAGCT"};
  int lengths[] = {(int)strlen(seqs[0])};
  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();
  int status = tsta_msa_align(NULL, lengths, 1, &config, &result);
  int status2 = tsta_msa_align(seqs, lengths, 0, &config, &result);
  int status3 = tsta_msa_align(NULL, lengths, 1, &config, &result);
  int status4 = tsta_msa_align(seqs, lengths, 0, &config, &result);

  assert(status == -1);
  assert(status2 == -1);
  assert(status3 == -1);
  assert(status4 == -1);
  tsta_msa_result_free(&result);
  return 0;
}

int test_msa_complex() {
  const char *seqs[] = {
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAAAAAGCTAGGCTAACGTAGCTAGGTAGGCTA",
  };
  int lengths[] = {
      (int)strlen(seqs[0]),
      (int)strlen(seqs[1]),
      (int)strlen(seqs[2]),
      (int)strlen(seqs[3]),
  };
  const size_t seq_count = sizeof(seqs) / sizeof(seqs[0]);

  for (size_t i = 0; i < seq_count; i++) {
    assert(seqs[i] != NULL);
    assert(lengths[i] > 40);
  }

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();
  int status = tsta_msa_align(seqs, lengths, seq_count, &config, &result);
  if (status == 0) {
    printf("MSA score: %d\n", result.score);
    for (size_t i = 0; i < result.sequence_count; i++) {
      printf("Aligned sequence %zu: %s\n", i + 1, result.aln[i]);
      assert(result.aln[i] != NULL);
    }
  }
  tsta_msa_result_free(&result);
  return 0;
}

int test_msa_trace_dump() {
  const char *seqs[] = {
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAAAAAGCTAGGCTAACGTAGCTAGGTAGGCTA",
  };
  int lengths[] = {
      (int)strlen(seqs[0]),
      (int)strlen(seqs[1]),
      (int)strlen(seqs[2]),
      (int)strlen(seqs[3]),
  };
  const size_t seq_count = sizeof(seqs) / sizeof(seqs[0]);
  tsta_config baseline_config = tsta_config_make_default();
  tsta_config spill_config = tsta_config_make_default();
  tsta_config no_compress_config = tsta_config_make_default();
  tsta_msa_result_t baseline = tsta_msa_result_make();
  tsta_msa_result_t spilled = tsta_msa_result_make();
  tsta_msa_result_t spilled_no_compress = tsta_msa_result_make();

  spill_config.msa_trace_dump_threshold = 1;
  spill_config.msa_trace_dump_compress = 1;
  no_compress_config.msa_trace_dump_threshold = 1;
  no_compress_config.msa_trace_dump_compress = 0;

  assert(tsta_msa_align(seqs, lengths, seq_count, &baseline_config,
                        &baseline) == 0);
  assert(tsta_msa_align(seqs, lengths, seq_count, &spill_config, &spilled) ==
         0);
  assert(tsta_msa_align(seqs, lengths, seq_count, &no_compress_config,
                        &spilled_no_compress) == 0);
  assert_msa_result_equal(&baseline, &spilled);
  assert_msa_result_equal(&baseline, &spilled_no_compress);

  tsta_msa_result_free(&baseline);
  tsta_msa_result_free(&spilled);
  tsta_msa_result_free(&spilled_no_compress);
  return 0;
}

int test_msa_trace_dump_chunked() {
  enum {
    repeat_count = 12,
    motif_length = 16,
    sequence_length = repeat_count * motif_length + 1,
  };
  char seq1[sequence_length];
  char seq2[sequence_length];
  char seq3[sequence_length];
  char seq4[sequence_length];
  const char *seqs[] = {seq1, seq2, seq3, seq4};
  int lengths[] = {0, 0, 0, 0};
  const size_t seq_count = sizeof(seqs) / sizeof(seqs[0]);
  tsta_config baseline_config = tsta_config_make_default();
  tsta_config chunked_config = tsta_config_make_default();
  tsta_msa_result_t baseline = tsta_msa_result_make();
  tsta_msa_result_t chunked = tsta_msa_result_make();

  fill_repeated_sequence(seq1, repeat_count, "ACGTAGCTTAGGCTAC");
  fill_repeated_sequence(seq2, repeat_count, "ACGTAGCTTAGGCTAC");
  fill_repeated_sequence(seq3, repeat_count, "ACGTAGCTTAGGCTAC");
  fill_repeated_sequence(seq4, repeat_count, "ACGTAGCTTAGGCTAC");

  seq2[17] = 'T';
  seq2[63] = 'A';
  seq3[31] = 'G';
  seq3[95] = 'C';
  seq4[47] = 'A';
  seq4[111] = 'T';

  lengths[0] = (int)strlen(seqs[0]);
  lengths[1] = (int)strlen(seqs[1]);
  lengths[2] = (int)strlen(seqs[2]);
  lengths[3] = (int)strlen(seqs[3]);

  for (size_t i = 0; i < seq_count; i++) {
    assert(lengths[i] > 160);
  }

  chunked_config.msa_trace_dump_threshold = 1;
  chunked_config.msa_trace_dump_compress = 0;

  assert(tsta_msa_align(seqs, lengths, seq_count, &baseline_config,
                        &baseline) == 0);
  assert(tsta_msa_align(seqs, lengths, seq_count, &chunked_config, &chunked) ==
         0);
  assert_msa_result_equal(&baseline, &chunked);

  tsta_msa_result_free(&baseline);
  tsta_msa_result_free(&chunked);
  return 0;
}

int test_msa_many_sequences() {
  const char *seqs[] = {
      "ACGTACGTACGT", "ACGTACGTACGA", "ACGTACGTACGG",
      "ACGTACGTACGC", "ACGTACGTACGT", "ACGTACGTTCGT",
      "ACGTACGTAAGT", "ACGTACGTACCT", "ACGTACGTACGT",
  };
  int lengths[] = {
      (int)strlen(seqs[0]), (int)strlen(seqs[1]), (int)strlen(seqs[2]),
      (int)strlen(seqs[3]), (int)strlen(seqs[4]), (int)strlen(seqs[5]),
      (int)strlen(seqs[6]), (int)strlen(seqs[7]), (int)strlen(seqs[8]),
  };
  const size_t seq_count = sizeof(seqs) / sizeof(seqs[0]);
  tsta_config baseline_config = tsta_config_make_default();
  tsta_config spill_config = tsta_config_make_default();
  tsta_msa_result_t baseline = tsta_msa_result_make();
  tsta_msa_result_t spilled = tsta_msa_result_make();

  spill_config.msa_trace_dump_threshold = 1;
  spill_config.msa_trace_dump_compress = 0;

  assert(tsta_msa_align(seqs, lengths, seq_count, &baseline_config,
                        &baseline) == 0);
  assert(tsta_msa_align(seqs, lengths, seq_count, &spill_config, &spilled) ==
         0);
  assert_msa_result_equal(&baseline, &spilled);

  tsta_msa_result_free(&baseline);
  tsta_msa_result_free(&spilled);
  return 0;
}

int test_msa_aligner_reuse() {
  const char *seqs[] = {
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA",
  };
  int lengths[] = {
      (int)strlen(seqs[0]),
      (int)strlen(seqs[1]),
      (int)strlen(seqs[2]),
  };
  const size_t seq_count = sizeof(seqs) / sizeof(seqs[0]);
  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t baseline = tsta_msa_result_make();
  tsta_msa_result_t manual = tsta_msa_result_make();
  tsta_msa_result_t reused = tsta_msa_result_make();
  tsta_msa_aligner *aligner;
  int status;

  assert(tsta_msa_align(seqs, lengths, seq_count, &config, &baseline) == 0);

  aligner = tsta_msa_aligner_create(&config);
  assert(aligner != NULL);

  status = tsta_msa_aligner_begin(aligner, seqs[0], lengths[0]);
  assert(status == 0);
  for (size_t i = 1; i < seq_count; i++) {
    status = tsta_msa_aligner_add(aligner, seqs[i], lengths[i]);
    assert(status == 0);
  }
  assert(tsta_msa_aligner_get_result(aligner, &manual) == 0);
  assert_msa_result_equal(&baseline, &manual);

  tsta_msa_result_free(&manual);
  assert(tsta_msa_aligner_align(aligner, seqs, lengths, seq_count, &reused) ==
         0);
  assert_msa_result_equal(&baseline, &reused);

  tsta_msa_aligner_destroy(aligner);
  tsta_msa_result_free(&baseline);
  tsta_msa_result_free(&reused);
  return 0;
}

int test_msa_fallback_clamp() {
  const char *seqs[] = {
      "ACGTAC",
      "ACGTC",
      "ACGTT",
  };
  int lengths[] = {
      (int)strlen(seqs[0]),
      (int)strlen(seqs[1]),
      (int)strlen(seqs[2]),
  };
  const size_t seq_count = sizeof(seqs) / sizeof(seqs[0]);
  tsta_config baseline_config = tsta_config_make_default();
  tsta_config oversized_config = tsta_config_make_default();
  tsta_msa_result_t baseline = tsta_msa_result_make();
  tsta_msa_result_t oversized = tsta_msa_result_make();

  oversized_config.block_size = 4096;
  oversized_config.threads = 4;
  oversized_config.msa_trace_dump_threshold = 0;

  assert(tsta_msa_align(seqs, lengths, seq_count, &baseline_config,
                        &baseline) == 0);
  assert(tsta_msa_align(seqs, lengths, seq_count, &oversized_config,
                        &oversized) == 0);
  assert_msa_result_equal(&baseline, &oversized);

  tsta_msa_result_free(&baseline);
  tsta_msa_result_free(&oversized);
  return 0;
}

int main(void) {
  assert(test_psa() == 0);
  assert(test_msa() == 0);
  assert(test_msa_complex() == 0);
  assert(test_msa_trace_dump() == 0);
  assert(test_msa_trace_dump_chunked() == 0);
  assert(test_msa_many_sequences() == 0);
  assert(test_msa_aligner_reuse() == 0);
  assert(test_msa_fallback_clamp() == 0);
  return 0;
}
