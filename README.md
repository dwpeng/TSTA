# TSTA

> **Note:** This project is derived from [bxskdh/TSTA](https://github.com/bxskdh/TSTA).
> The core alignment algorithms remain unchanged. Modifications are limited to
> engineering improvements: build system, public API design, C++ wrapper, test
> coverage, code formatting, memory management and threading.

SIMD-accelerated pairwise and multiple sequence alignment library. No external
dependencies beyond libc and pthreads; all data in memory.

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

Useful targets: `make check` (tests), `make check_psa_simd` (standalone vs
threaded cross-check), `make format` / `make check-format` (clang-format),
`make install-hooks` (pre-commit format check).

## Pairwise alignment

```c
#include "tsta.h"
#include <string.h>

tsta_config cfg = tsta_config_make_default();   // match=2, mismatch=-5, gap_open=-4, gap_extend=-2

const char *s1 = "ACGTAGCTAGCTAGCTAGCTA";
const char *s2 = "AGCTAGCTAGCTAGCTAGC";
tsta_psa_result_t r = tsta_psa_result_make();

if (tsta_psa_align(s1, (int)strlen(s1), s2, (int)strlen(s2), &cfg, &r) == 0) {
    printf("score=%d  cigar=%s\n", r.score, r.cigar);
    printf(">1\n%s\n>2\n%s\n", r.aln[0], r.aln[1]);
}
tsta_psa_result_free(&r);
```

`block_size` defaults to 0 (auto: matched to the sequence length, keeping
padding minimal); `threads` defaults to 10. For best throughput on repeated
alignments, reuse an aligner (the one-shot path creates its threadpool per
call and so uses 2 threads). Traceback (aligned sequences + CIGAR) is always
produced.

## Multiple sequence alignment

One-shot:

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

Incremental:

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

Aligners are reusable — runtime buffers are cached across calls, so repeated
alignments on the same aligner are cheap.

## Standalone SIMD module

Single-threaded pairwise alignment with SIMD only — no pthread. Shares the
same core as the threaded library, results are bit-identical. Builds as
`libtsta_psa_simd` (no `-lpthread` needed).

```c
#include "tsta_psa_simd.h"   // reuses tsta_config / tsta_psa_result_t

tsta_config cfg = tsta_config_make_default();
tsta_psa_result_t r = tsta_psa_result_make();
tsta_psa_simd_align(s1, (int)strlen(s1), s2, (int)strlen(s2), &cfg, &r);
tsta_psa_result_free(&r);
```

## Performance

Measured on AVX-512, 8 cores (median of 3 runs; "allocated" = cumulative
malloc/calloc/realloc bytes).

### Pairwise alignment

| Workload | Time | Allocations | Allocated |
|---|---|---|---|
| 1000 × 100 bp, one-shot | 0.26 s | 132k | 29 MB |
| 500 × 500 bp, one-shot | 0.22 s | 266k | 159 MB |
| 200 × 1000 bp, one-shot | 0.17 s | 206k | 228 MB |
| 1000 × 100 bp, aligner reuse | 0.08 s | 5k | 2 MB |
| 500 × 500 bp, aligner reuse | 0.11 s | 3k | 4 MB |
| 200 × 1000 bp, aligner reuse | 0.10 s | 2k | 4 MB |
| 10000 × 10000 bp, single | 0.08 s | 10k | 102 MB |

### Multiple sequence alignment

| Workload | Time | Allocations |
|---|---|---|
| 10 × 500 bp | 0.02 s | 31k |
| 20 × 1000 bp | 0.50 s | 464k |
| 8 × 3000 bp | 0.38 s | 399k |

## Build integration

pkg-config:

```makefile
CFLAGS  += $(shell pkg-config --cflags tsta)
LDFLAGS += $(shell pkg-config --libs tsta)
```

Manual link:

```bash
cc -o myapp myapp.c -ltsta -lpthread        # threaded library
cc -o myapp myapp.c -ltsta_psa_simd         # standalone module (no pthread)
```

## Contact

- Peiyu Zong peiyuzong8@gmail.com (core algorithms)
- Wenpeng Deng 1732889554@qq.com (optimization)
