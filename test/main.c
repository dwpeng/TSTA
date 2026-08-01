#include "tsta.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Helpers ─────────────────────────────────────────────────────────── */

static void
assert_msa_result_equal(const tsta_msa_result_t* expected,
                        const tsta_msa_result_t* actual)
{
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

static void
fill_repeated_sequence(char* buffer, size_t repeat_count, const char* motif)
{
  size_t motif_length = strlen(motif);
  for (size_t i = 0; i < repeat_count; i++)
    memcpy(buffer + i * motif_length, motif, motif_length);
  buffer[repeat_count * motif_length] = '\0';
}

/* ── PSA: basic ──────────────────────────────────────────────────────── */

static int
test_psa(void)
{
  const char* seq1 = "ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
  const char* seq2 =
      "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  int status = tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                              &config, &result);
  assert(status == 0);
  assert(result.score != 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.aln_length > 0);
  assert(result.cigar != NULL);
  assert(strlen(result.cigar) > 0);
  /* aligned length >= longest input */
  assert(result.aln_length >= strlen(seq1));
  assert(result.aln_length >= strlen(seq2));
  tsta_psa_result_free(&result);
  return 0;
}

/* ── PSA: very short sequences ───────────────────────────────────────── */

static int
test_psa_short(void)
{
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result;

  /* single character, match */
  result = tsta_psa_result_make();
  assert(tsta_psa_align("A", 1, "A", 1, &config, &result) == 0);
  assert(result.score > 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.match_count > 0);
  tsta_psa_result_free(&result);

  /* single character, mismatch */
  result = tsta_psa_result_make();
  assert(tsta_psa_align("A", 1, "T", 1, &config, &result) == 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.mismatch_count > 0);
  tsta_psa_result_free(&result);

  /* 3 vs 5 characters */
  result = tsta_psa_result_make();
  assert(tsta_psa_align("ACG", 3, "ACGTA", 5, &config, &result) == 0);
  assert(result.aln_length >= 5);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  tsta_psa_result_free(&result);

  return 0;
}

/* ── PSA: identical sequences ────────────────────────────────────────── */

static int
test_psa_identical(void)
{
  const char* seq = "ACGTACGTACGT";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_align(seq, (int)strlen(seq), seq, (int)strlen(seq), &config,
                        &result)
         == 0);
  assert(result.score > 0);
  assert(result.mismatch_count == 0);
  assert(result.del_count == 0);
  assert(result.ins_count == 0);
  /* aligned sequences should match input (no gaps for identical) */
  assert(result.aln_length == strlen(seq));
  tsta_psa_result_free(&result);
  return 0;
}

/* ── PSA: completely different (poly-A vs poly-T) ────────────────────── */

static int
test_psa_divergent(void)
{
  const char* seq1 = "AAAAAAAAAA";
  const char* seq2 = "TTTTTTTTTT";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                        &config, &result)
         == 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  /* For entirely divergent sequences, aligner may prefer gaps over
   * mismatches. Either way, match count should be zero or near-zero. */
  assert(result.match_count == 0);
  tsta_psa_result_free(&result);
  return 0;
}

/* ── PSA: sequences with N (ambiguous base) ──────────────────────────── */

static int
test_psa_N_bases(void)
{
  const char* seq1 = "ACGTNNNACGT";
  const char* seq2 = "ACGTAAACGT";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                        &config, &result)
         == 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  tsta_psa_result_free(&result);
  return 0;
}

/* ── PSA: asymmetric lengths ─────────────────────────────────────────── */

static int
test_psa_asymmetric(void)
{
  char short_seq[256];
  char long_seq[2048];

  fill_repeated_sequence(short_seq, 10, "ACGT"); /* 40 bp */
  fill_repeated_sequence(long_seq, 500, "ACGT"); /* 2000 bp */

  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_align(long_seq, (int)strlen(long_seq), short_seq,
                        (int)strlen(short_seq), &config, &result)
         == 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.aln_length >= strlen(long_seq));
  tsta_psa_result_free(&result);
  return 0;
}

/* ── PSA: long sequences (multi-block) ───────────────────────────────── */

static int
test_psa_long(void)
{
  char seq1[1024];
  char seq2[1024];

  fill_repeated_sequence(seq1, 100, "ACGTACGT"); /* 800 bp */
  /* seq2: same motif with a few SNPs every 80bp */
  {
    const char* pat = "ACGTANGT";
    size_t pos = 0;
    for (int i = 0; i < 100; i++, pos += 8)
      memcpy(seq2 + pos, pat, 8);
    seq2[pos] = '\0';
  }

  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                        &config, &result)
         == 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.aln_length >= strlen(seq1));
  assert(result.cigar != NULL);
  tsta_psa_result_free(&result);
  return 0;
}

/* ── PSA: very long sequences ────────────────────────────────────────── */

static int
test_psa_very_long(void)
{
  char* seq1 = (char*)malloc(4001);
  char* seq2 = (char*)malloc(4001);
  assert(seq1 && seq2);

  fill_repeated_sequence(seq1, 500, "ACGTACGT"); /* 4000 bp */
  fill_repeated_sequence(seq2, 500, "ACGTACNT"); /* 4000 bp, SNP every 8 */

  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                        &config, &result)
         == 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.aln_length >= strlen(seq1));
  tsta_psa_result_free(&result);
  free(seq1);
  free(seq2);
  return 0;
}

/* ── PSA: custom config ──────────────────────────────────────────────── */

static int
test_psa_custom_config(void)
{
  const char* seq1 = "ACGTACGTACGTACGT";
  const char* seq2 = "ACGTANGTACGTACGT";

  tsta_config config;
  tsta_config_default(&config);
  config.match = 5;
  config.mismatch = -4;
  config.gap_open = -10;
  config.gap_extend = -1;
  config.block_size = 8;
  config.threads = 2;

  tsta_psa_result_t result = tsta_psa_result_make();
  assert(tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                        &config, &result)
         == 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  /* match=5 > default=2, should give positive score */
  assert(result.score > 0);
  tsta_psa_result_free(&result);
  return 0;
}

/* ── PSA: aligner reuse ──────────────────────────────────────────────── */

static int
test_psa_aligner_reuse(void)
{
  const char* pairs[][2] = {
    { "ACGTACGT", "ACGTANGT" },
    { "GGGGCCCCAAAATTTT", "GGGGCCCCTTTTAAAA" },
    { "ACGTACGTACGTACGTACGT", "ACGTACNTACGTACGTACGT" },
  };
  const int n_pairs = (int)(sizeof(pairs) / sizeof(pairs[0]));

  tsta_config config = tsta_config_make_default();
  tsta_psa_aligner* aligner = tsta_psa_aligner_create(&config);
  assert(aligner != NULL);

  for (int i = 0; i < n_pairs; i++) {
    tsta_psa_result_t result = tsta_psa_result_make();
    assert(tsta_psa_aligner_align(aligner, pairs[i][0],
                                  (int)strlen(pairs[i][0]), pairs[i][1],
                                  (int)strlen(pairs[i][1]), &result)
           == 0);
    assert(result.aln[0] != NULL);
    assert(result.aln[1] != NULL);
    tsta_psa_result_free(&result);
  }

  tsta_psa_aligner_destroy(aligner);
  return 0;
}

/* ── PSA: error handling ─────────────────────────────────────────────── */

static int
test_psa_errors(void)
{
  tsta_config config = tsta_config_make_default();

  /* null result returns -1 */
  assert(tsta_psa_align("A", 1, "A", 1, &config, NULL) == -1);

  /* null result via aligner */
  tsta_psa_aligner* aligner = tsta_psa_aligner_create(&config);
  assert(aligner != NULL);
  assert(tsta_psa_aligner_align(aligner, "A", 1, "A", 1, NULL) == -1);
  tsta_psa_aligner_destroy(aligner);

  /* null config -> uses default */
  tsta_psa_result_t result = tsta_psa_result_make();
  assert(tsta_psa_align("A", 1, "T", 1, NULL, &result) == 0);
  tsta_psa_result_free(&result);

  /* null aligner destroy (should not crash) */
  tsta_psa_aligner_destroy(NULL);

  /* null aligner create output */
  assert(tsta_psa_aligner_align(NULL, "A", 1, "A", 1, &result) == -1);

  return 0;
}

/* ── PSA: CIGAR validation ───────────────────────────────────────────── */

static int
test_psa_cigar(void)
{
  const char* seq1 = "ACGTACGT";
  const char* seq2 = "ACGTACGT";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                        &config, &result)
         == 0);
  /* For identical sequences, CIGAR should be all M */
  assert(result.cigar != NULL);
  assert(strstr(result.cigar, "M") != NULL);
  assert(result.del_count == 0);
  assert(result.ins_count == 0);
  tsta_psa_result_free(&result);
  return 0;
}

/* ── MSA: basic (single sequence = error) ────────────────────────────── */

static int
test_msa(void)
{
  const char* seqs[] = { "ACGTAGCTAGCT" };
  int lengths[] = { (int)strlen(seqs[0]) };
  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align(NULL, lengths, 1, &config, &result) == -1);
  assert(tsta_msa_align(seqs, lengths, 0, &config, &result) == -1);
  tsta_msa_result_free(&result);
  return 0;
}

/* ── MSA: 2 sequences ────────────────────────────────────────────────── */

static int
test_msa_two(void)
{
  const char* seqs[] = {
    "ACGTACGTACGTACGT",
    "ACGTANGTACGTACGT",
  };
  int lengths[] = { (int)strlen(seqs[0]), (int)strlen(seqs[1]) };

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align(seqs, lengths, 2, &config, &result) == 0);
  assert(result.sequence_count == 2);
  assert(result.aln_length > 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  /* aligned length >= longest input */
  size_t max_input =
      strlen(seqs[0]) > strlen(seqs[1]) ? strlen(seqs[0]) : strlen(seqs[1]);
  assert(result.aln_length >= max_input);
  tsta_msa_result_free(&result);
  return 0;
}

/* ── MSA: complex (original test) ────────────────────────────────────── */

static int
test_msa_complex(void)
{
  const char* seqs[] = {
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

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();
  int status = tsta_msa_align(seqs, lengths, seq_count, &config, &result);
  assert(status == 0);
  assert(result.score != 0);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);
  tsta_msa_result_free(&result);
  return 0;
}

/* ── MSA: many short sequences ───────────────────────────────────────── */

static int
test_msa_many_short(void)
{
  const char* motifs[] = {
    "ACGTACGTACGT", "ACGTANGTACGT", "ACGTACNTACGT", "ACNTACGTACGT",
    "ANGTACGTACGT", "ACGTACGTANGT", "ACGTACGTACNT", "ACGTANGTACNT",
    "ANGTACNTACGT", "ACGTACGTACGT",
  };
  const int n_seqs = (int)(sizeof(motifs) / sizeof(motifs[0])); /* 10 */
  int lengths[10];

  for (int i = 0; i < n_seqs; i++)
    lengths[i] = (int)strlen(motifs[i]);

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align(motifs, lengths, (size_t)n_seqs, &config, &result)
         == 0);
  assert(result.sequence_count == (size_t)n_seqs);
  assert(result.aln_length > 0);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);
  tsta_msa_result_free(&result);
  return 0;
}

/* ── MSA: many sequences (20) ────────────────────────────────────────── */

static int
test_msa_many(void)
{
  const char* seqs[20];
  char seq_bufs[20][64];
  int lengths[20];

  for (int i = 0; i < 20; i++) {
    fill_repeated_sequence(seq_bufs[i], 4, "ACGTACGT"); /* 32 bp */
    /* introduce a SNP at position i */
    if (i < 32)
      seq_bufs[i][i] = (seq_bufs[i][i] == 'A') ? 'T' : 'A';
    seqs[i] = seq_bufs[i];
    lengths[i] = (int)strlen(seqs[i]);
  }

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align(seqs, lengths, 20, &config, &result) == 0);
  assert(result.sequence_count == 20);
  assert(result.aln_length > 0);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);
  tsta_msa_result_free(&result);
  return 0;
}

/* ── MSA: long sequences (two seqs, 800bp) ───────────────────────────── */

static int
test_msa_long_two(void)
{
  char* seqs[2];
  int lengths[2];

  seqs[0] = (char*)malloc(801);
  seqs[1] = (char*)malloc(801);
  assert(seqs[0] && seqs[1]);

  fill_repeated_sequence(seqs[0], 100, "ACGTACGT"); /* 800 bp */
  fill_repeated_sequence(seqs[1], 100, "ACGTANGT"); /* 800 bp */
  lengths[0] = (int)strlen(seqs[0]);
  lengths[1] = (int)strlen(seqs[1]);

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align((const char* const*)seqs, lengths, 2, &config, &result)
         == 0);
  assert(result.sequence_count == 2);
  assert(result.aln_length > 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);

  tsta_msa_result_free(&result);
  free(seqs[0]);
  free(seqs[1]);
  return 0;
}

/* ── MSA: many long sequences (5 sequences, 500bp each) ──────────────── */

static int
test_msa_many_long(void)
{
  const int n_seqs = 5;
  const int repeat = 62; /* 62 * 8 = 496 bp */
  char* seqs[5];
  int lengths[5];
  const char* patterns[] = {
    "ACGTACGT", "ACGTANGT", "ANGTACGT", "ACNTACGT", "ACGTACNT",
  };

  for (int i = 0; i < n_seqs; i++) {
    seqs[i] = (char*)malloc((size_t)(repeat * 8 + 1));
    assert(seqs[i]);
    fill_repeated_sequence(seqs[i], repeat, patterns[i]);
    lengths[i] = (int)strlen(seqs[i]);
  }

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align((const char* const*)seqs, lengths, (size_t)n_seqs,
                        &config, &result)
         == 0);
  assert(result.sequence_count == (size_t)n_seqs);
  assert(result.aln_length >= (size_t)lengths[0]);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);

  tsta_msa_result_free(&result);
  for (int i = 0; i < n_seqs; i++)
    free(seqs[i]);
  return 0;
}

/* ── MSA: long + many (10 sequences, 300bp each) ─────────────────────── */

static int
test_msa_long_many(void)
{
  const int n_seqs = 10;
  const int repeat = 37; /* 37 * 8 = 296 bp */
  char* seqs[10];
  int lengths[10];

  for (int i = 0; i < n_seqs; i++) {
    seqs[i] = (char*)malloc((size_t)(repeat * 8 + 1));
    assert(seqs[i]);
    fill_repeated_sequence(seqs[i], repeat, "ACGTACGT");
    /* introduce unique SNPs */
    if (i < (int)strlen(seqs[i]))
      seqs[i][i * 3] = (seqs[i][i * 3] == 'A') ? 'T' : 'A';
    if (i * 5 + 2 < (int)strlen(seqs[i]))
      seqs[i][i * 5 + 2] = (seqs[i][i * 5 + 2] == 'C') ? 'G' : 'C';
    lengths[i] = (int)strlen(seqs[i]);
  }

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align((const char* const*)seqs, lengths, (size_t)n_seqs,
                        &config, &result)
         == 0);
  assert(result.sequence_count == (size_t)n_seqs);
  assert(result.aln_length > 0);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);

  tsta_msa_result_free(&result);
  for (int i = 0; i < n_seqs; i++)
    free(seqs[i]);
  return 0;
}

/* ── MSA: asymmetric lengths ─────────────────────────────────────────── */

static int
test_msa_asymmetric(void)
{
  const char* seqs[] = {
    "ACGT",                 /* 4 bp */
    "ACGTACGTACGTACGTACGT", /* 20 bp */
    "ACGTACGT",             /* 8 bp */
  };
  int lengths[] = {
    (int)strlen(seqs[0]),
    (int)strlen(seqs[1]),
    (int)strlen(seqs[2]),
  };

  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result = tsta_msa_result_make();

  assert(tsta_msa_align(seqs, lengths, 3, &config, &result) == 0);
  assert(result.sequence_count == 3);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.aln[2] != NULL);
  tsta_msa_result_free(&result);
  return 0;
}

/* ── MSA: aligner reuse ──────────────────────────────────────────────── */

static int
test_msa_aligner_reuse(void)
{
  const char* seqs[] = {
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
  tsta_msa_aligner* aligner;

  assert(tsta_msa_align(seqs, lengths, seq_count, &config, &baseline) == 0);

  aligner = tsta_msa_aligner_create(&config);
  assert(aligner != NULL);

  assert(tsta_msa_aligner_begin(aligner, seqs[0], lengths[0]) == 0);
  for (size_t i = 1; i < seq_count; i++)
    assert(tsta_msa_aligner_add(aligner, seqs[i], lengths[i]) == 0);
  assert(tsta_msa_aligner_get_result(aligner, &manual) == 0);
  assert_msa_result_equal(&baseline, &manual);

  tsta_msa_result_free(&manual);
  assert(tsta_msa_aligner_align(aligner, seqs, lengths, seq_count, &reused)
         == 0);
  assert_msa_result_equal(&baseline, &reused);

  tsta_msa_aligner_destroy(aligner);
  tsta_msa_result_free(&baseline);
  tsta_msa_result_free(&reused);
  return 0;
}

/* ── MSA: fallback clamp ─────────────────────────────────────────────── */

static int
test_msa_fallback_clamp(void)
{
  const char* seqs[] = { "ACGTAC", "ACGTC", "ACGTT" };
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

  assert(tsta_msa_align(seqs, lengths, seq_count, &baseline_config, &baseline)
         == 0);
  assert(
      tsta_msa_align(seqs, lengths, seq_count, &oversized_config, &oversized)
      == 0);
  assert_msa_result_equal(&baseline, &oversized);

  tsta_msa_result_free(&baseline);
  tsta_msa_result_free(&oversized);
  return 0;
}

/* ── MSA: custom config ──────────────────────────────────────────────── */

static int
test_msa_custom_config(void)
{
  const char* seqs[] = {
    "ACGTACGTACGTACGT",
    "ACGTANGTACGTACGT",
    "ACGTACNTACGTACGT",
  };
  int lengths[] = {
    (int)strlen(seqs[0]),
    (int)strlen(seqs[1]),
    (int)strlen(seqs[2]),
  };

  tsta_config config;
  tsta_config_default(&config);
  config.match = 3;
  config.mismatch = -3;
  config.gap_open = -8;
  config.gap_extend = -2;
  config.block_size = 5;
  config.threads = 2;

  tsta_msa_result_t result = tsta_msa_result_make();
  assert(tsta_msa_align(seqs, lengths, 3, &config, &result) == 0);
  assert(result.sequence_count == 3);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);
  tsta_msa_result_free(&result);
  return 0;
}

/* ── MSA: error handling ─────────────────────────────────────────────── */

static int
test_msa_errors(void)
{
  const char* seqs[] = { "ACGT", "ACGT" };
  int lengths[] = { 4, 4 };
  tsta_config config = tsta_config_make_default();
  tsta_msa_result_t result;

  /* null sequences */
  assert(tsta_msa_align(NULL, lengths, 2, &config, &result) == -1);

  /* zero sequence count */
  assert(tsta_msa_align(seqs, lengths, 0, &config, &result) == -1);

  /* null result */
  assert(tsta_msa_align(seqs, lengths, 2, &config, NULL) == -1);

  /* null aligner create -> destroyed safely */
  tsta_msa_aligner_destroy(NULL);

  /* null initial sequence for begin */
  tsta_msa_aligner* aligner = tsta_msa_aligner_create(&config);
  assert(aligner != NULL);
  assert(tsta_msa_aligner_begin(aligner, NULL, 4) == -1);

  /* add without begin */
  tsta_msa_aligner* aligner2 = tsta_msa_aligner_create(&config);
  assert(aligner2 != NULL);
  assert(tsta_msa_aligner_add(aligner2, "ACGT", 4) == -1);
  tsta_msa_aligner_destroy(aligner2);

  /* get_result without begin */
  tsta_msa_result_t empty = tsta_msa_result_make();
  assert(tsta_msa_aligner_get_result(aligner, &empty) == -1);
  tsta_msa_result_free(&empty);

  tsta_msa_aligner_destroy(aligner);
  return 0;
}

/* ── MSA: incremental build with many sequences ──────────────────────── */

static int
test_msa_incremental_many(void)
{
  const int n_seqs = 8;
  const char* motifs[] = {
    "ACGTACGTACGTACGT", "ACGTANGTACGTACGT", "ACGTACNTACGTACGT",
    "ANGTACGTACGTACGT", "ACNTACGTACGTACGT", "ACGTACGTANGTACGT",
    "ACGTACGTACNTACGT", "ACGTACGTACGTANGT",
  };

  tsta_config config = tsta_config_make_default();
  tsta_msa_aligner* aligner = tsta_msa_aligner_create(&config);
  assert(aligner != NULL);

  assert(tsta_msa_aligner_begin(aligner, motifs[0], (int)strlen(motifs[0]))
         == 0);
  for (int i = 1; i < n_seqs; i++)
    assert(tsta_msa_aligner_add(aligner, motifs[i], (int)strlen(motifs[i]))
           == 0);

  tsta_msa_result_t result = tsta_msa_result_make();
  assert(tsta_msa_aligner_get_result(aligner, &result) == 0);
  assert(result.sequence_count == (size_t)n_seqs);
  assert(result.aln_length > 0);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);

  tsta_msa_result_free(&result);
  tsta_msa_aligner_destroy(aligner);
  return 0;
}

/* ── MSA: incremental long sequences ─────────────────────────────────── */

static int
test_msa_incremental_long(void)
{
  const int n_seqs = 4;
  char* seqs[4];
  const char* patterns[] = {
    "ACGTACGT",
    "ACGTANGT",
    "ANGTACGT",
    "ACNTACGT",
  };

  for (int i = 0; i < n_seqs; i++) {
    seqs[i] = (char*)malloc(801);
    assert(seqs[i]);
    fill_repeated_sequence(seqs[i], 100, patterns[i]); /* 800 bp */
  }

  tsta_config config = tsta_config_make_default();
  tsta_msa_aligner* aligner = tsta_msa_aligner_create(&config);
  assert(aligner != NULL);

  assert(tsta_msa_aligner_begin(aligner, seqs[0], (int)strlen(seqs[0])) == 0);
  for (int i = 1; i < n_seqs; i++)
    assert(tsta_msa_aligner_add(aligner, seqs[i], (int)strlen(seqs[i])) == 0);

  tsta_msa_result_t result = tsta_msa_result_make();
  assert(tsta_msa_aligner_get_result(aligner, &result) == 0);
  assert(result.sequence_count == (size_t)n_seqs);
  assert(result.aln_length > 0);
  for (size_t i = 0; i < result.sequence_count; i++)
    assert(result.aln[i] != NULL);

  tsta_msa_result_free(&result);
  tsta_msa_aligner_destroy(aligner);
  for (int i = 0; i < n_seqs; i++)
    free(seqs[i]);
  return 0;
}

/* ── MSA: result lifecycle ───────────────────────────────────────────── */

static int
test_msa_result_lifecycle(void)
{
  /* init */
  tsta_msa_result_t r1;
  tsta_msa_result_init(&r1);
  assert(r1.score == 0);
  assert(r1.aln == NULL);
  assert(r1.sequence_count == 0);
  assert(r1.aln_length == 0);

  /* make */
  tsta_msa_result_t r2 = tsta_msa_result_make();
  assert(r2.score == 0);
  assert(r2.aln == NULL);

  /* free on empty result (safe) */
  tsta_msa_result_free(&r1);
  tsta_msa_result_free(&r2);

  /* free on NULL (safe) */
  tsta_msa_result_free(NULL);

  return 0;
}

/* ── PSA: result lifecycle ───────────────────────────────────────────── */

static int
test_psa_result_lifecycle(void)
{
  tsta_psa_result_t r1;
  tsta_psa_result_init(&r1);
  assert(r1.score == 0);
  assert(r1.aln[0] == NULL);
  assert(r1.aln[1] == NULL);
  assert(r1.cigar == NULL);

  tsta_psa_result_t r2 = tsta_psa_result_make();
  assert(r2.score == 0);

  tsta_psa_result_free(&r1);
  tsta_psa_result_free(&r2);
  tsta_psa_result_free(NULL);

  return 0;
}

/* ── Config lifecycle ────────────────────────────────────────────────── */

static int
test_config_lifecycle(void)
{
  tsta_config c1;
  tsta_config_default(&c1);
  assert(c1.match == 2);
  assert(c1.mismatch == -5);
  assert(c1.gap_extend == -2);
  assert(c1.gap_open == -4);
  assert(c1.block_size == 10);
  assert(c1.threads == 10);

  tsta_config c2 = tsta_config_make_default();
  assert(c2.match == 2);

  /* null config default does not crash */
  tsta_config_default(NULL);

  return 0;
}

/* ── Main ────────────────────────────────────────────────────────────── */

int
main(void)
{
  printf("=== PSA tests ===\n");

  printf("  test_psa ... ");
  assert(test_psa() == 0);
  printf("ok\n");

  printf("  test_psa_short ... ");
  assert(test_psa_short() == 0);
  printf("ok\n");

  printf("  test_psa_identical ... ");
  assert(test_psa_identical() == 0);
  printf("ok\n");

  printf("  test_psa_divergent ... ");
  assert(test_psa_divergent() == 0);
  printf("ok\n");

  printf("  test_psa_N_bases ... ");
  assert(test_psa_N_bases() == 0);
  printf("ok\n");

  printf("  test_psa_asymmetric ... ");
  assert(test_psa_asymmetric() == 0);
  printf("ok\n");

  printf("  test_psa_long ... ");
  assert(test_psa_long() == 0);
  printf("ok\n");

  printf("  test_psa_very_long ... ");
  assert(test_psa_very_long() == 0);
  printf("ok\n");

  printf("  test_psa_custom_config ... ");
  assert(test_psa_custom_config() == 0);
  printf("ok\n");

  printf("  test_psa_aligner_reuse ... ");
  assert(test_psa_aligner_reuse() == 0);
  printf("ok\n");

  printf("  test_psa_errors ... ");
  assert(test_psa_errors() == 0);
  printf("ok\n");

  printf("  test_psa_cigar ... ");
  assert(test_psa_cigar() == 0);
  printf("ok\n");

  printf("  test_psa_result_lifecycle ... ");
  assert(test_psa_result_lifecycle() == 0);
  printf("ok\n");

  printf("=== MSA tests ===\n");

  printf("  test_msa ... ");
  assert(test_msa() == 0);
  printf("ok\n");

  printf("  test_msa_two ... ");
  assert(test_msa_two() == 0);
  printf("ok\n");

  printf("  test_msa_complex ... ");
  assert(test_msa_complex() == 0);
  printf("ok\n");

  printf("  test_msa_many_short ... ");
  assert(test_msa_many_short() == 0);
  printf("ok\n");

  printf("  test_msa_many ... ");
  assert(test_msa_many() == 0);
  printf("ok\n");

  printf("  test_msa_long_two ... ");
  assert(test_msa_long_two() == 0);
  printf("ok\n");

  printf("  test_msa_many_long ... ");
  assert(test_msa_many_long() == 0);
  printf("ok\n");

  printf("  test_msa_long_many ... ");
  assert(test_msa_long_many() == 0);
  printf("ok\n");

  printf("  test_msa_asymmetric ... ");
  assert(test_msa_asymmetric() == 0);
  printf("ok\n");

  printf("  test_msa_aligner_reuse ... ");
  assert(test_msa_aligner_reuse() == 0);
  printf("ok\n");

  printf("  test_msa_fallback_clamp ... ");
  assert(test_msa_fallback_clamp() == 0);
  printf("ok\n");

  printf("  test_msa_custom_config ... ");
  assert(test_msa_custom_config() == 0);
  printf("ok\n");

  printf("  test_msa_errors ... ");
  assert(test_msa_errors() == 0);
  printf("ok\n");

  printf("  test_msa_incremental_many ... ");
  assert(test_msa_incremental_many() == 0);
  printf("ok\n");

  printf("  test_msa_incremental_long ... ");
  assert(test_msa_incremental_long() == 0);
  printf("ok\n");

  printf("  test_msa_result_lifecycle ... ");
  assert(test_msa_result_lifecycle() == 0);
  printf("ok\n");

  printf("  test_config_lifecycle ... ");
  assert(test_config_lifecycle() == 0);
  printf("ok\n");

  printf("\nAll tests passed.\n");
  return 0;
}
