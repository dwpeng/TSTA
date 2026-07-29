# TSTA

Pairwise and multiple sequence alignment library accelerated by SIMD and threads.

## Introduction

SIMD-accelerated global sequence alignment using striped anti-diagonal DP with difference recurrence. MSA uses Partial Order Alignment (POA) DAG. All data in memory — no disk I/O, no external dependencies beyond libc and pthreads.

## Installation

```bash
git clone https://github.com/dwpeng/TSTA
cd TSTA
make
```

SIMD ISA auto-detected via `-march=native`. Override: `make CFLAGS="-O3 -g -mavx2"`.

Install system-wide:

```bash
sudo make install PREFIX=/usr/local
```

## API

Single header: `#include <tsta.h>`.

### Configuration

```c
tsta_config cfg = tsta_config_make_default();
cfg.match    = 2;     // default
cfg.mismatch = -5;    // default
cfg.gap_open = -4;    // default
cfg.gap_extend = -2;  // default

// block_size / threads auto-clamped per sequence length
```

### Pairwise alignment

```c
#include "tsta.h"
#include <string.h>

const char *s1 = "ACGTAGCTAGCTAGCTAGCTA";
const char *s2 = "AGCTAGCTAGCTAGCTAGC";
tsta_config cfg = tsta_config_make_default();
tsta_psa_result_t r = tsta_psa_result_make();

if (tsta_psa_align(s1, (int)strlen(s1), s2, (int)strlen(s2), &cfg, &r) == 0) {
    printf("score=%d  cigar=%s\n", r.score, r.cigar);
    printf(">1\n%s\n>2\n%s\n", r.aln[0], r.aln[1]);
}
tsta_psa_result_free(&r);
```

### Multiple sequence alignment (one-shot)

```c
const char *seqs[] = {s1, s2, s3};
int lens[] = {(int)strlen(s1), (int)strlen(s2), (int)strlen(s3)};
tsta_config cfg = tsta_config_make_default();
tsta_msa_result_t r = tsta_msa_result_make();

if (tsta_msa_align(seqs, lens, 3, &cfg, &r) == 0) {
    for (size_t i = 0; i < r.sequence_count; i++)
        printf("seq%zu: %s\n", i + 1, r.aln[i]);
}
tsta_msa_result_free(&r);
```

### Multiple sequence alignment (incremental)

```c
tsta_msa_aligner *a = tsta_msa_aligner_create(&cfg);

tsta_msa_aligner_begin(a, s1, (int)strlen(s1));
tsta_msa_aligner_add(a, s2, (int)strlen(s2));
tsta_msa_aligner_add(a, s3, (int)strlen(s3));

tsta_msa_result_t r = tsta_msa_result_make();
tsta_msa_aligner_get_result(a, &r);
// use r.aln[], r.score, r.aln_length ...
tsta_msa_result_free(&r);

tsta_msa_aligner_destroy(a);
```

### Aligner reuse

```c
tsta_msa_aligner *a = tsta_msa_aligner_create(&cfg);

// first MSA
tsta_msa_aligner_align(a, seqs_a, lens_a, count_a, &r1);
tsta_msa_result_free(&r1);

// second MSA — same aligner, same config
tsta_msa_aligner_align(a, seqs_b, lens_b, count_b, &r2);
tsta_msa_result_free(&r2);

tsta_msa_aligner_destroy(a);
```

## Build integration

pkg-config:

```makefile
CFLAGS  += $(shell pkg-config --cflags tsta)
LDFLAGS += $(shell pkg-config --libs tsta)
```

Manual link:

```bash
cc -o myapp myapp.c -ltsta -lpthread
```

## Contact

- Peiyu Zong peiyuzong8@gmail.com (core algorithms)
- Wenpeng Deng 1732889554@qq.com (optimization)
