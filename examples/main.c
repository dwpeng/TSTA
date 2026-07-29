#include "tsta.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
main(void)
{
  /* ── Pairwise alignment ──────────────────────────────────────────── */
  {
    const char* s1 = "ACGTACGTACGTACGT";
    const char* s2 = "ACGTANGTACGTACGT";

    tsta_config cfg = tsta_config_make_default();
    tsta_psa_result_t r = tsta_psa_result_make();

    if (tsta_psa_align(s1, (int)strlen(s1), s2, (int)strlen(s2), &cfg, &r)
        != 0) {
      fprintf(stderr, "PSA alignment failed\n");
      return 1;
    }

    printf("=== Pairwise alignment ===\n");
    printf("score: %d  cigar: %s\n", r.score, r.cigar);
    printf("matches=%d  mismatches=%d  insertions=%d  deletions=%d\n",
           r.match_count, r.mismatch_count, r.ins_count, r.del_count);
    printf(">seq1\n%s\n>seq2\n%s\n\n", r.aln[0], r.aln[1]);

    tsta_psa_result_free(&r);
  }

  /* ── PSA aligner reuse ───────────────────────────────────────────── */
  {
    tsta_config cfg = tsta_config_make_default();
    tsta_psa_aligner* a = tsta_psa_aligner_create(&cfg);
    if (!a) {
      fprintf(stderr, "Failed to create PSA aligner\n");
      return 1;
    }

    const char* pairs[][2] = {
      { "ACGTACGTACGT", "ACGTANGTACGT" },
      { "GGGGCCCCAAAA", "GGGGCCCCTTTT" },
    };

    printf("=== PSA aligner reuse ===\n");
    for (int i = 0; i < 2; i++) {
      tsta_psa_result_t r = tsta_psa_result_make();
      if (tsta_psa_aligner_align(a, pairs[i][0], (int)strlen(pairs[i][0]),
                                 pairs[i][1], (int)strlen(pairs[i][1]), &r)
          == 0) {
        printf("pair %d: score=%d cigar=%s\n", i + 1, r.score, r.cigar);
      }
      tsta_psa_result_free(&r);
    }
    printf("\n");

    tsta_psa_aligner_destroy(a);
  }

  /* ── Multiple sequence alignment (one-shot) ──────────────────────── */
  {
    const char* seqs[] = {
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAAAAAGCTAGGCTAACGTAGCTAGGTAGGCTA",
    };
    int lens[] = { (int)strlen(seqs[0]), (int)strlen(seqs[1]),
                   (int)strlen(seqs[2]), (int)strlen(seqs[3]) };
    const size_t count = sizeof(seqs) / sizeof(seqs[0]);

    tsta_config cfg = tsta_config_make_default();
    tsta_msa_result_t r = tsta_msa_result_make();

    if (tsta_msa_align(seqs, lens, count, &cfg, &r) != 0) {
      fprintf(stderr, "MSA alignment failed\n");
      return 1;
    }

    printf("=== Multiple sequence alignment (one-shot) ===\n");
    printf("score: %d  sequences: %zu  length: %zu\n", r.score,
           r.sequence_count, r.aln_length);
    for (size_t i = 0; i < r.sequence_count; i++)
      printf(">seq%zu\n%s\n", i + 1, r.aln[i]);
    printf("\n");

    tsta_msa_result_free(&r);
  }

  /* ── MSA incremental ─────────────────────────────────────────────── */
  {
    tsta_config cfg = tsta_config_make_default();
    tsta_msa_aligner* a = tsta_msa_aligner_create(&cfg);
    if (!a) {
      fprintf(stderr, "Failed to create MSA aligner\n");
      return 1;
    }

    const char* s1 = "ACGTACGTACGTACGT";
    const char* s2 = "ACGTANGTACGTACGT";
    const char* s3 = "ACGTACNTACGTACGT";

    tsta_msa_aligner_begin(a, s1, (int)strlen(s1));
    tsta_msa_aligner_add(a, s2, (int)strlen(s2));
    tsta_msa_aligner_add(a, s3, (int)strlen(s3));

    tsta_msa_result_t r = tsta_msa_result_make();
    tsta_msa_aligner_get_result(a, &r);

    printf("=== MSA incremental ===\n");
    printf("score: %d  sequences: %zu\n", r.score, r.sequence_count);
    for (size_t i = 0; i < r.sequence_count; i++)
      printf(">seq%zu\n%s\n", i + 1, r.aln[i]);

    tsta_msa_result_free(&r);
    tsta_msa_aligner_destroy(a);
  }

  return 0;
}
