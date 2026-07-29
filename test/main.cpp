#include "tsta.hpp"

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

static void
test_config_chain()
{
  tsta::Config cfg;
  cfg.match(5).mismatch(-4).gap_open(-10).gap_extend(-1).block_size(8).threads(
      2);
  assert(cfg.match() == 5);
  assert(cfg.mismatch() == -4);
  assert(cfg.gap_open() == -10);
  assert(cfg.gap_extend() == -1);
  assert(cfg.block_size() == 8);
  assert(cfg.threads() == 2);
}

static void
test_psa_one_shot()
{
  auto r = tsta::psa_align("ACGTACGTACGT", "ACGTANGTACGT");
  assert(r.score() != 0);
  assert(!r.seq1().empty());
  assert(!r.seq2().empty());
  assert(!r.cigar().empty());
  assert(r.aln_length() > 0);

  /* string_view overload via std::string */
  std::string s1 = "ACGTACGT";
  std::string s2 = "ACGTANGT";
  auto r2 = tsta::psa_align(s1, s2);
  assert(r2.score() != 0);
}

static void
test_psa_custom_config()
{
  tsta::Config cfg;
  cfg.match(5).mismatch(-4);
  auto r = tsta::psa_align("ACGTACGTACGT", "ACGTANGTACGT", cfg);
  assert(r.score() > 0);
}

static void
test_psa_aligner_reuse()
{
  tsta::PsaAligner aln;
  auto r1 = aln.align("ACGTACGT", "ACGTANGT");
  assert(r1.score() != 0);

  auto r2 = aln.align("GGGGCCCCAAAATTTT", "GGGGCCCCTTTTAAAA");
  assert(r2.score() != 0);

  /* move */
  tsta::PsaAligner aln2 = std::move(aln);
  auto r3 = aln2.align("ACGT", "ACGT");
  assert(r3.score() > 0);
}

static void
test_psa_result_move()
{
  auto r1 = tsta::psa_align("ACGTACGT", "ACGTACGT");
  int score = r1.score();
  tsta::PsaResult r2 = std::move(r1);
  assert(r2.score() == score);
}

static void
test_msa_one_shot()
{
  /* vector<string> overload */
  std::vector<std::string> seqs = {
    "ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA",
    "ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA",
    "ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA",
    "ACGTAGCTAGCTAGGCTAACGTAGCTAAAAAGCTAGGCTAACGTAGCTAGGTAGGCTA",
  };
  auto r = tsta::msa_align(seqs);
  assert(r.sequence_count() == 4);
  assert(r.aln_length() > 0);
  for (size_t i = 0; i < r.sequence_count(); i++)
    assert(!r.sequence(i).empty());
}

static void
test_msa_incremental()
{
  tsta::MsaAligner aln;
  aln.begin("ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA");
  aln.add("ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA");
  aln.add("ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA");
  aln.add("ACGTAGCTAGCTAGGCTAACGTAGCTAAAAAGCTAGGCTAACGTAGCTAGGTAGGCTA");

  auto r = aln.get_result();
  assert(r.sequence_count() == 4);
  assert(r.aln_length() > 0);
  for (size_t i = 0; i < r.sequence_count(); i++)
    assert(!r.sequence(i).empty());
}

static void
test_msa_aligner_reuse()
{
  tsta::MsaAligner aln;

  const char* s1[] = { "ACGT", "ACGT" };
  int l1[] = { 4, 4 };
  auto r1 = aln.align(s1, l1, 2);
  assert(r1.sequence_count() == 2);

  const char* s2[] = { "GGGGCCCC", "GGGGCCCC" };
  int l2[] = { 8, 8 };
  auto r2 = aln.align(s2, l2, 2);
  assert(r2.sequence_count() == 2);
}

static void
test_msa_result_sequences()
{
  auto r = tsta::msa_align(std::vector<std::string>{ "ACGTACGT", "ACGTANGT" });
  auto all = r.sequences();
  assert(all.size() == 2);
  assert(!all[0].empty());
  assert(!all[1].empty());
}

static void
test_msa_result_move()
{
  std::vector<std::string> seqs = { "ACGTACGT", "ACGTANGT", "ACGTACNT" };
  auto r1 = tsta::msa_align(seqs);
  int score = r1.score();
  tsta::MsaResult r2 = std::move(r1);
  assert(r2.score() == score);
  assert(r2.sequence_count() == 3);
}

static void
test_asymmetric_long()
{
  /* long sequence test via C++ wrapper */
  std::string long_seq(2000, 'A');
  for (size_t i = 0; i < long_seq.size(); i += 4)
    long_seq[i] = "ACGT"[i % 4];
  std::string short_seq(40, 'A');
  for (size_t i = 0; i < short_seq.size(); i += 4)
    short_seq[i] = "ACGT"[i % 4];

  auto r = tsta::psa_align(long_seq, short_seq);
  assert(r.score() != 0);
  assert(r.aln_length() >= long_seq.size());
}

int
main()
{
  printf("test_config_chain ... ");
  test_config_chain();
  printf("ok\n");

  printf("test_psa_one_shot ... ");
  test_psa_one_shot();
  printf("ok\n");

  printf("test_psa_custom_config ... ");
  test_psa_custom_config();
  printf("ok\n");

  printf("test_psa_aligner_reuse ... ");
  test_psa_aligner_reuse();
  printf("ok\n");

  printf("test_psa_result_move ... ");
  test_psa_result_move();
  printf("ok\n");

  printf("test_msa_one_shot ... ");
  test_msa_one_shot();
  printf("ok\n");

  printf("test_msa_incremental ... ");
  test_msa_incremental();
  printf("ok\n");

  printf("test_msa_aligner_reuse ... ");
  test_msa_aligner_reuse();
  printf("ok\n");

  printf("test_msa_result_sequences ... ");
  test_msa_result_sequences();
  printf("ok\n");

  printf("test_msa_result_move ... ");
  test_msa_result_move();
  printf("ok\n");

  printf("test_asymmetric_long ... ");
  test_asymmetric_long();
  printf("ok\n");

  printf("\nAll C++ tests passed.\n");
  return 0;
}
