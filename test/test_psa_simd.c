/*
 * test_psa_simd.c — Tests for the standalone SIMD-only pairwise alignment
 * module (tsta_psa_simd.h).
 *
 * Two layers:
 *   1. Standalone API tests mirroring the PSA cases in test/main.c.
 *   2. Cross-check: for the same configuration and inputs, the standalone
 *      module must produce results bit-identical to the threaded library
 *      (tsta_psa_align). This links both static libraries.
 */

#include "tsta.h"
#include "tsta_psa_simd.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_failures = 0;

static void
fill_repeated_sequence(char* buffer, size_t repeat_count, const char* motif)
{
  size_t motif_length = strlen(motif);
  for (size_t i = 0; i < repeat_count; i++)
    memcpy(buffer + i * motif_length, motif, motif_length);
  buffer[repeat_count * motif_length] = '\0';
}

/* ── Deterministic PRNG (xorshift64*) ────────────────────────────────── */

static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;

static uint64_t
rng_next(void)
{
  uint64_t x = rng_state;
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  rng_state = x;
  return x * 0x2545F4914F6CDD1DULL;
}

static int
rng_int(int n)
{
  return (int)(rng_next() % (uint64_t)n);
}

static const char DNA[] = "ACGT";
static const char DNA_N[] = "ACGTN";

static void
gen_random(char* buffer, size_t len, int with_n)
{
  const char* alphabet = with_n ? DNA_N : DNA;
  for (size_t i = 0; i < len; i++)
    buffer[i] = alphabet[rng_int(with_n ? 5 : 4)];
  buffer[len] = '\0';
}

static void
gen_homopolymer(char* buffer, size_t len)
{
  for (size_t i = 0; i < len; i++)
    buffer[i] = (rng_next() & 1) ? 'A' : 'T';
  buffer[len] = '\0';
}

static void
gen_derived(const char* src, char* dst, size_t dst_cap, size_t* out_len)
{
  size_t n = strlen(src);
  size_t d = 0;
  int edits = 1 + rng_int(5);

  for (size_t i = 0; i < n && d + 2 < dst_cap; i++) {
    if (edits > 0 && rng_int(10) == 0) {
      edits--;
      int op = rng_int(3);
      if (op == 0) { /* substitution */
        dst[d++] = DNA[rng_int(4)];
      } else if (op == 1) { /* insertion */
        dst[d++] = src[i];
        dst[d++] = DNA[rng_int(4)];
      } else { /* deletion: skip src[i] */
      }
      continue;
    }
    dst[d++] = src[i];
  }
  dst[d] = '\0';
  *out_len = d;
}

/* ── Cross-check: standalone vs threaded ─────────────────────────────── */

static int
check_equal(const char* s1,
            int l1,
            const char* s2,
            int l2,
            const tsta_config* cfg,
            const char* label)
{
  tsta_psa_result_t rt = tsta_psa_result_make();
  tsta_psa_result_t rs = tsta_psa_result_make();
  int st = tsta_psa_align(s1, l1, s2, l2, cfg, &rt);
  int ss = tsta_psa_simd_align(s1, l1, s2, l2, cfg, &rs);

  if (st != ss) {
    printf("  MISMATCH [%s]: status threaded=%d standalone=%d\n", label, st,
           ss);
    g_failures++;
    goto done;
  }
  if (st != 0) { /* both rejected identically */
    goto done;
  }

  if (rt.score != rs.score || rt.aln_length != rs.aln_length
      || rt.match_count != rs.match_count
      || rt.mismatch_count != rs.mismatch_count || rt.ins_count != rs.ins_count
      || rt.del_count != rs.del_count) {
    printf("  MISMATCH [%s]: counters differ (score %d/%d len %zu/%zu "
           "M %d/%d X %d/%d I %d/%d D %d/%d)\n",
           label, rt.score, rs.score, rt.aln_length, rs.aln_length,
           rt.match_count, rs.match_count, rt.mismatch_count,
           rs.mismatch_count, rt.ins_count, rs.ins_count, rt.del_count,
           rs.del_count);
    g_failures++;
    goto done;
  }
  if (strcmp(rt.cigar, rs.cigar) != 0) {
    printf("  MISMATCH [%s]: cigar \"%s\" vs \"%s\"\n", label, rt.cigar,
           rs.cigar);
    g_failures++;
    goto done;
  }
  if (strcmp(rt.aln[0], rs.aln[0]) != 0 || strcmp(rt.aln[1], rs.aln[1]) != 0) {
    printf("  MISMATCH [%s]: aligned sequences differ\n", label);
    printf("    threaded:  %s / %s\n", rt.aln[0], rt.aln[1]);
    printf("    standalone: %s / %s\n", rs.aln[0], rs.aln[1]);
    g_failures++;
    goto done;
  }

done:
  tsta_psa_result_free(&rt);
  tsta_psa_result_free(&rs);
  return g_failures ? -1 : 0;
}

/* ── Standalone API tests (mirror test/main.c PSA cases) ─────────────── */

static int
test_standalone_basic(void)
{
  const char* seq1 = "ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
  const char* seq2 =
      "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  int status = tsta_psa_simd_align(seq1, (int)strlen(seq1), seq2,
                                   (int)strlen(seq2), &config, &result);
  assert(status == 0);
  assert(result.score != 0);
  assert(result.aln[0] != NULL);
  assert(result.aln[1] != NULL);
  assert(result.aln_length > 0);
  assert(result.cigar != NULL);
  assert(strlen(result.cigar) > 0);
  assert(result.aln_length >= strlen(seq1));
  assert(result.aln_length >= strlen(seq2));
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_short(void)
{
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result;

  result = tsta_psa_result_make();
  assert(tsta_psa_simd_align("A", 1, "A", 1, &config, &result) == 0);
  assert(result.score > 0);
  assert(result.match_count > 0);
  tsta_psa_result_free(&result);

  result = tsta_psa_result_make();
  assert(tsta_psa_simd_align("A", 1, "T", 1, &config, &result) == 0);
  assert(result.mismatch_count > 0);
  tsta_psa_result_free(&result);

  result = tsta_psa_result_make();
  assert(tsta_psa_simd_align("ACG", 3, "ACGTA", 5, &config, &result) == 0);
  assert(result.aln_length >= 5);
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_identical(void)
{
  const char* seq = "ACGTACGTACGT";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_simd_align(seq, (int)strlen(seq), seq, (int)strlen(seq),
                             &config, &result)
         == 0);
  assert(result.score > 0);
  assert(result.mismatch_count == 0);
  assert(result.del_count == 0);
  assert(result.ins_count == 0);
  assert(result.aln_length == strlen(seq));
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_divergent(void)
{
  const char* seq1 = "AAAAAAAAAA";
  const char* seq2 = "TTTTTTTTTT";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_simd_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                             &config, &result)
         == 0);
  assert(result.match_count == 0);
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_cigar(void)
{
  const char* seq1 = "ACGTACGT";
  const char* seq2 = "ACGTACGT";
  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_simd_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                             &config, &result)
         == 0);
  assert(strstr(result.cigar, "M") != NULL);
  assert(result.del_count == 0);
  assert(result.ins_count == 0);
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_asymmetric(void)
{
  char short_seq[256];
  char long_seq[2048];

  fill_repeated_sequence(short_seq, 10, "ACGT");
  fill_repeated_sequence(long_seq, 500, "ACGT");

  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_simd_align(long_seq, (int)strlen(long_seq), short_seq,
                             (int)strlen(short_seq), &config, &result)
         == 0);
  assert(result.aln_length >= strlen(long_seq));
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_long(void)
{
  char seq1[1024];
  char seq2[1024];

  fill_repeated_sequence(seq1, 100, "ACGTACGT");
  {
    const char* pat = "ACGTANGT";
    size_t pos = 0;
    for (int i = 0; i < 100; i++, pos += 8)
      memcpy(seq2 + pos, pat, 8);
    seq2[pos] = '\0';
  }

  tsta_config config = tsta_config_make_default();
  tsta_psa_result_t result = tsta_psa_result_make();

  assert(tsta_psa_simd_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                             &config, &result)
         == 0);
  assert(result.aln_length >= strlen(seq1));
  assert(result.cigar != NULL);
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_custom_config(void)
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
  assert(tsta_psa_simd_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
                             &config, &result)
         == 0);
  assert(result.score > 0);
  tsta_psa_result_free(&result);
  return 0;
}

static int
test_standalone_errors(void)
{
  tsta_config config = tsta_config_make_default();

  assert(tsta_psa_simd_align("A", 1, "A", 1, &config, NULL) == -1);

  tsta_psa_simd_aligner* aligner = tsta_psa_simd_aligner_create(&config);
  assert(aligner != NULL);
  assert(tsta_psa_simd_aligner_align(aligner, "A", 1, "A", 1, NULL) == -1);
  assert(tsta_psa_simd_aligner_align(aligner, NULL, 1, "A", 1, NULL) == -1);
  tsta_psa_simd_aligner_destroy(aligner);

  tsta_psa_simd_aligner_destroy(NULL);

  assert(tsta_psa_simd_aligner_align(NULL, "A", 1, "A", 1, NULL) == -1);

  /* null config -> defaults */
  tsta_psa_result_t result = tsta_psa_result_make();
  assert(tsta_psa_simd_align("A", 1, "T", 1, NULL, &result) == 0);
  tsta_psa_result_free(&result);

  /* empty sequence -> rejected identically to threaded path */
  result = tsta_psa_result_make();
  assert(tsta_psa_simd_align("A", 1, "", 0, &config, &result) == -1);
  tsta_psa_result_free(&result);

  return 0;
}

static int
test_standalone_result_lifecycle(void)
{
  tsta_psa_result_t r1 = tsta_psa_result_make();
  assert(r1.score == 0);
  assert(r1.aln[0] == NULL);
  assert(r1.aln[1] == NULL);
  assert(r1.cigar == NULL);
  tsta_psa_result_free(&r1);
  tsta_psa_result_free(NULL);
  return 0;
}

static int
test_standalone_aligner_reuse(void)
{
  const char* pairs[][2] = {
    { "ACGTACGT", "ACGTANGT" },
    { "GGGGCCCCAAAATTTT", "GGGGCCCCTTTTAAAA" },
    { "ACGTACGTACGTACGTACGT", "ACGTACNTACGTACGTACGT" },
  };
  const int n_pairs = (int)(sizeof(pairs) / sizeof(pairs[0]));

  tsta_config config = tsta_config_make_default();
  tsta_psa_simd_aligner* aligner = tsta_psa_simd_aligner_create(&config);
  assert(aligner != NULL);

  for (int i = 0; i < n_pairs; i++) {
    tsta_psa_result_t result = tsta_psa_result_make();
    assert(tsta_psa_simd_aligner_align(aligner, pairs[i][0],
                                       (int)strlen(pairs[i][0]), pairs[i][1],
                                       (int)strlen(pairs[i][1]), &result)
           == 0);
    assert(result.aln[0] != NULL);
    assert(result.aln[1] != NULL);
    tsta_psa_result_free(&result);
  }

  tsta_psa_simd_aligner_destroy(aligner);
  return 0;
}

static int
test_fallback_clamp(void)
{
  const char* seqs[][2] = {
    { "ACGTACGTACGTACGT", "ACGTANGTACGTACGT" },
    { "ACGT", "ACGTACGTACGT" },
  };
  const int n_pairs = (int)(sizeof(seqs) / sizeof(seqs[0]));

  for (int p = 0; p < n_pairs; p++) {
    const char* s1 = seqs[p][0];
    const char* s2 = seqs[p][1];

    tsta_config baseline_cfg = tsta_config_make_default();
    tsta_config oversized_cfg = tsta_config_make_default();
    oversized_cfg.block_size = 4096;
    oversized_cfg.threads = 4;

    tsta_psa_result_t tb = tsta_psa_result_make();
    tsta_psa_result_t to = tsta_psa_result_make();
    tsta_psa_result_t sb = tsta_psa_result_make();
    tsta_psa_result_t so = tsta_psa_result_make();

    assert(tsta_psa_align(s1, (int)strlen(s1), s2, (int)strlen(s2),
                          &baseline_cfg, &tb)
           == 0);
    assert(tsta_psa_align(s1, (int)strlen(s1), s2, (int)strlen(s2),
                          &oversized_cfg, &to)
           == 0);
    assert(tsta_psa_simd_align(s1, (int)strlen(s1), s2, (int)strlen(s2),
                               &baseline_cfg, &sb)
           == 0);
    assert(tsta_psa_simd_align(s1, (int)strlen(s1), s2, (int)strlen(s2),
                               &oversized_cfg, &so)
           == 0);

    /* oversized block_size must clamp to the same value in both libraries */
    assert(tb.score == to.score);
    assert(tb.score == sb.score);
    assert(tb.score == so.score);
    assert(strcmp(tb.cigar, to.cigar) == 0);
    assert(strcmp(tb.cigar, sb.cigar) == 0);
    assert(strcmp(tb.cigar, so.cigar) == 0);

    tsta_psa_result_free(&tb);
    tsta_psa_result_free(&to);
    tsta_psa_result_free(&sb);
    tsta_psa_result_free(&so);
  }
  return 0;
}

/* ── Cross-check suites ──────────────────────────────────────────────── */

typedef struct config_case {
  tsta_config cfg;
  const char* name;
} config_case;

static void
init_config_cases(config_case cases[6])
{
  tsta_config_default(&cases[0].cfg);
  cases[0].name = "default";

  tsta_config_default(&cases[1].cfg);
  cases[1].cfg.match = 5;
  cases[1].cfg.mismatch = -4;
  cases[1].cfg.gap_open = -10;
  cases[1].cfg.gap_extend = -1;
  cases[1].cfg.block_size = 8;
  cases[1].cfg.threads = 2;
  cases[1].name = "custom scores";

  tsta_config_default(&cases[2].cfg);
  cases[2].cfg.block_size = 1;
  cases[2].cfg.threads = 1;
  cases[2].name = "block_size=1";

  tsta_config_default(&cases[3].cfg);
  cases[3].cfg.block_size = 2;
  cases[3].cfg.threads = 4;
  cases[3].name = "block_size=2";

  tsta_config_default(&cases[4].cfg);
  cases[4].cfg.block_size = 64;
  cases[4].cfg.threads = 8;
  cases[4].name = "block_size=64";

  tsta_config_default(&cases[5].cfg);
  cases[5].cfg.block_size = 4096;
  cases[5].cfg.threads = 4;
  cases[5].name = "block_size=4096 (clamped)";
}

static void
check_random_sweep(void)
{
  config_case cases[6];
  init_config_cases(cases);

  char s1[4096];
  char s2[4096];
  char derived[4096];
  char base[1024];

  for (int c = 0; c < 6; c++) {
    char label[64];
    for (int i = 0; i < 30; i++) {
      int len1, len2;
      int style = rng_int(4);
      switch (style) {
      case 0:
        len1 = 1 + rng_int(200);
        len2 = 1 + rng_int(200);
        gen_random(s1, (size_t)len1, 0);
        gen_random(s2, (size_t)len2, 0);
        break;
      case 1:
        len1 = 1 + rng_int(300);
        len2 = 1 + rng_int(300);
        gen_random(s1, (size_t)len1, 1);
        gen_random(s2, (size_t)len2, 1);
        break;
      case 2:
        len1 = 4 + rng_int(200);
        len2 = 4 + rng_int(200);
        gen_homopolymer(s1, (size_t)len1);
        gen_homopolymer(s2, (size_t)len2);
        break;
      default: {
        size_t b = 8 + rng_int(200);
        gen_random(base, b, 0);
        size_t d1, d2;
        gen_derived(base, s1, sizeof(s1), &d1);
        gen_derived(base, s2, sizeof(s2), &d2);
        len1 = (int)d1;
        len2 = (int)d2;
        break;
      }
      }
      snprintf(label, sizeof(label), "%s/rand%d", cases[c].name, i);
      check_equal(s1, len1, s2, len2, &cases[c].cfg, label);
    }
  }
}

static void
check_long_pairs(void)
{
  config_case cases[6];
  init_config_cases(cases);

  char s1[4096];
  char s2[4096];
  char base[4096];

  for (int c = 0; c < 2; c++) {
    for (int i = 0; i < 6; i++) {
      size_t b = 500 + (size_t)rng_int(3000);
      gen_random(base, b, 0);
      size_t d1 = b, d2 = b;
      gen_derived(base, s1, sizeof(s1), &d1);
      gen_derived(base, s2, sizeof(s2), &d2);
      char label[64];
      snprintf(label, sizeof(label), "%s/long%d (b=%zu)", cases[c].name, i, b);
      check_equal(s1, (int)d1, s2, (int)d2, &cases[c].cfg, label);
    }
  }
}

static void
check_aligner_reuse_cross(void)
{
  tsta_config cfg = tsta_config_make_default();
  cfg.block_size = 10;
  cfg.threads = 4;

  tsta_psa_aligner* ta = tsta_psa_aligner_create(&cfg);
  tsta_psa_simd_aligner* sa = tsta_psa_simd_aligner_create(&cfg);
  assert(ta != NULL);
  assert(sa != NULL);

  const int lens[] = { 10, 100, 1000, 3000 };
  char s1[4096];
  char s2[4096];
  char base[4096];

  for (int i = 0; i < 4; i++) {
    size_t b = (size_t)lens[i];
    gen_random(base, b, 0);
    size_t d1 = b, d2 = b;
    gen_derived(base, s1, sizeof(s1), &d1);
    gen_derived(base, s2, sizeof(s2), &d2);

    tsta_psa_result_t rt = tsta_psa_result_make();
    tsta_psa_result_t rs = tsta_psa_result_make();
    assert(tsta_psa_aligner_align(ta, s1, (int)d1, s2, (int)d2, &rt) == 0);
    assert(tsta_psa_simd_aligner_align(sa, s1, (int)d1, s2, (int)d2, &rs)
           == 0);

    if (rt.score != rs.score || rt.aln_length != rs.aln_length
        || strcmp(rt.cigar, rs.cigar) != 0 || strcmp(rt.aln[0], rs.aln[0]) != 0
        || strcmp(rt.aln[1], rs.aln[1]) != 0) {
      printf("  MISMATCH [aligner-reuse len=%d]: threaded vs standalone\n",
             lens[i]);
      g_failures++;
    }

    tsta_psa_result_free(&rt);
    tsta_psa_result_free(&rs);
  }

  tsta_psa_aligner_destroy(ta);
  tsta_psa_simd_aligner_destroy(sa);
}

/* ── Main ────────────────────────────────────────────────────────────── */

int
main(void)
{
  printf("=== Standalone PSA-SIMD tests ===\n");

  printf("  test_standalone_basic ... ");
  assert(test_standalone_basic() == 0);
  printf("ok\n");

  printf("  test_standalone_short ... ");
  assert(test_standalone_short() == 0);
  printf("ok\n");

  printf("  test_standalone_identical ... ");
  assert(test_standalone_identical() == 0);
  printf("ok\n");

  printf("  test_standalone_divergent ... ");
  assert(test_standalone_divergent() == 0);
  printf("ok\n");

  printf("  test_standalone_cigar ... ");
  assert(test_standalone_cigar() == 0);
  printf("ok\n");

  printf("  test_standalone_asymmetric ... ");
  assert(test_standalone_asymmetric() == 0);
  printf("ok\n");

  printf("  test_standalone_long ... ");
  assert(test_standalone_long() == 0);
  printf("ok\n");

  printf("  test_standalone_custom_config ... ");
  assert(test_standalone_custom_config() == 0);
  printf("ok\n");

  printf("  test_standalone_errors ... ");
  assert(test_standalone_errors() == 0);
  printf("ok\n");

  printf("  test_standalone_result_lifecycle ... ");
  assert(test_standalone_result_lifecycle() == 0);
  printf("ok\n");

  printf("  test_standalone_aligner_reuse ... ");
  assert(test_standalone_aligner_reuse() == 0);
  printf("ok\n");

  printf("  test_fallback_clamp ... ");
  assert(test_fallback_clamp() == 0);
  printf("ok\n");

  printf("=== Cross-check vs threaded tsta_psa_align ===\n");
  check_random_sweep();
  check_long_pairs();
  check_aligner_reuse_cross();

  if (g_failures > 0) {
    printf("\nFAILED: %d cross-check mismatches\n", g_failures);
    return 1;
  }

  printf("\nAll tests passed.\n");
  return 0;
}
