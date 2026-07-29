#include "tsta.hpp"

#include <cstdio>
#include <string>
#include <vector>

int
main()
{
  /* ── Pairwise alignment ──────────────────────────────────────────── */
  {
    auto r = tsta::psa_align("ACGTACGTACGTACGT", "ACGTANGTACGTACGT");

    printf("=== Pairwise alignment ===\n");
    printf("score: %d  cigar: %s\n", r.score(), r.cigar_cstr());
    printf("matches=%d  mismatches=%d  insertions=%d  deletions=%d\n",
           r.match_count(), r.mismatch_count(), r.ins_count(), r.del_count());
    printf(">seq1\n%s\n>seq2\n%s\n\n", r.seq1_cstr(), r.seq2_cstr());
  }

  /* ── PSA aligner reuse ───────────────────────────────────────────── */
  {
    tsta::PsaAligner aln;

    auto r1 = aln.align("ACGTACGTACGT", "ACGTANGTACGT");
    printf("=== PSA aligner reuse ===\n");
    printf("pair 1: score=%d cigar=%s\n", r1.score(), r1.cigar_cstr());

    auto r2 = aln.align("GGGGCCCCAAAA", "GGGGCCCCTTTT");
    printf("pair 2: score=%d cigar=%s\n\n", r2.score(), r2.cigar_cstr());
  }

  /* ── Custom config ───────────────────────────────────────────────── */
  {
    tsta::Config cfg;
    cfg.match(5).mismatch(-4).gap_open(-10).gap_extend(-1);

    auto r = tsta::psa_align("ACGTACGTACGTACGT", "ACGTANGTACGTACGT", cfg);

    printf("=== Custom config ===\n");
    printf("match=5 mismatch=-4 gap_open=-10 gap_extend=-1\n");
    printf("score: %d  cigar: %s\n\n", r.score(), r.cigar_cstr());
  }

  /* ── Multiple sequence alignment (one-shot) ──────────────────────── */
  {
    std::vector<std::string> seqs = {
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA",
      "ACGTAGCTAGCTAGGCTAACGTAGCTAAAAAGCTAGGCTAACGTAGCTAGGTAGGCTA",
    };

    auto r = tsta::msa_align(seqs);

    printf("=== Multiple sequence alignment (one-shot) ===\n");
    printf("score: %d  sequences: %zu  length: %zu\n", r.score(),
           r.sequence_count(), r.aln_length());
    for (size_t i = 0; i < r.sequence_count(); i++)
      printf(">seq%zu\n%s\n", i + 1, r.sequence_cstr(i));
    printf("\n");
  }

  /* ── MSA incremental ─────────────────────────────────────────────── */
  {
    tsta::MsaAligner aln;

    aln.begin("ACGTACGTACGTACGT");
    aln.add("ACGTANGTACGTACGT");
    aln.add("ACGTACNTACGTACGT");

    auto r = aln.get_result();

    printf("=== MSA incremental ===\n");
    printf("score: %d  sequences: %zu\n", r.score(), r.sequence_count());
    for (size_t i = 0; i < r.sequence_count(); i++)
      printf(">seq%zu\n%s\n", i + 1, r.sequence_cstr(i));
  }

  return 0;
}
