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
    for (size_t i = 0; i < repeat_count; i++)
        memcpy(buffer + i * motif_length, motif, motif_length);
    buffer[repeat_count * motif_length] = '\0';
}

static int test_psa(void) {
    const char *seq1 = "ACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC";
    const char *seq2 =
        "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
    tsta_config config = tsta_config_make_default();
    tsta_psa_result_t result = tsta_psa_result_make();

    int status = tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2),
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

static int test_msa(void) {
    const char *seqs[] = {"ACGTAGCTAGCT"};
    int lengths[] = {(int)strlen(seqs[0])};
    tsta_config config = tsta_config_make_default();
    tsta_msa_result_t result = tsta_msa_result_make();

    assert(tsta_msa_align(NULL, lengths, 1, &config, &result) == -1);
    assert(tsta_msa_align(seqs, lengths, 0, &config, &result) == -1);
    tsta_msa_result_free(&result);
    return 0;
}

static int test_msa_complex(void) {
    const char *seqs[] = {
        "ACGTAGCTAGCTAGGCTAACGTAGCTAGCTTTTTAGGCTAACGTAGCTAGCTAGGCTA",
        "ACGTAGCTAGCTCGGCTAACGTAGAAAAACTAGCTAGGCTAACGTAGCTAGCTAGACTA",
        "ACGTAGCTAGCTAGGCTAACGTAGCTAGAAAAACTTGGCTAACGTAGCTAGCTAGGCTA",
        "ACGTAGCTAGCTAGGCTAACGTAGCTAAAAAGCTAGGCTAACGTAGCTAGGTAGGCTA",
    };
    int lengths[] = {
        (int)strlen(seqs[0]), (int)strlen(seqs[1]),
        (int)strlen(seqs[2]), (int)strlen(seqs[3]),
    };
    const size_t seq_count = sizeof(seqs) / sizeof(seqs[0]);

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

static int test_msa_aligner_reuse(void) {
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

    assert(tsta_msa_align(seqs, lengths, seq_count, &config, &baseline) == 0);

    aligner = tsta_msa_aligner_create(&config);
    assert(aligner != NULL);

    assert(tsta_msa_aligner_begin(aligner, seqs[0], lengths[0]) == 0);
    for (size_t i = 1; i < seq_count; i++)
        assert(tsta_msa_aligner_add(aligner, seqs[i], lengths[i]) == 0);
    assert(tsta_msa_aligner_get_result(aligner, &manual) == 0);
    assert_msa_result_equal(&baseline, &manual);

    tsta_msa_result_free(&manual);
    assert(tsta_msa_aligner_align(aligner, seqs, lengths, seq_count, &reused) == 0);
    assert_msa_result_equal(&baseline, &reused);

    tsta_msa_aligner_destroy(aligner);
    tsta_msa_result_free(&baseline);
    tsta_msa_result_free(&reused);
    return 0;
}

static int test_msa_fallback_clamp(void) {
    const char *seqs[] = {"ACGTAC", "ACGTC", "ACGTT"};
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
    assert(test_msa_aligner_reuse() == 0);
    assert(test_msa_fallback_clamp() == 0);
    return 0;
}
