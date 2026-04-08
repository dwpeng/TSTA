## TSTA
pairwise and multiple sequence alignment accelerated by simd and threads
## Introduction
TSTA is a sequence alignment algorithm that combines SIMD instruction sets and thread acceleration. It uses the difference method, stripe method, and anti-diagonal method, and uses partial order alignment (POA) in multi sequence alignment.<br>

TSTA is committed to achieving faster comparison speed. TSTA only supports global alignment and allows linear and affine gap penalties. It supports SSE4.2/AVX2/AVX512 vectorization, with different choices based on hardware support.
## Installation
```
git clone https://github.com/bxskdh/TSTA.git
cd TSTA
make
```

The build now also produces `libtsta.a`, which exposes the core alignment
API in `include/tsta.h`. The library accepts in-memory sequences directly and
does not depend on the sequence file reader used by the legacy CLI tools.

### Library usage

Structured PSA result (counts, alignments, CIGAR):
```c
#include "tsta.h"
#include <string.h>

tsta_config config = tsta_config_make_default();
tsta_psa_result_t result = tsta_psa_result_make();

if (tsta_psa_align(seq1, (int)strlen(seq1), seq2, (int)strlen(seq2), &config,
                   &result) == 0) {
	printf("score=%d match=%d mismatch=%d del=%d ins=%d\n", result.score,
	       result.match_count, result.mismatch_count, result.del_count,
	       result.ins_count);
	printf("cigar=%s\n", result.cigar);
	printf("aln1=%s\n", result.aln[0]);
	printf("aln2=%s\n", result.aln[1]);
}

tsta_psa_result_free(&result);
```

Incremental MSA usage:
```c
#include "tsta.h"

tsta_config config = tsta_config_make_default();
tsta_msa_result_t result = tsta_msa_result_make();
tsta_msa_aligner *aligner;
int seed_len = (int)strlen(seq0);
int s1_len = (int)strlen(seq1);
int s2_len = (int)strlen(seq2);

aligner = tsta_msa_aligner_create(&config);

if (aligner && tsta_msa_aligner_begin(aligner, seq0, seed_len) == 0) {
	if (seq1) {
		tsta_msa_aligner_add(aligner, seq1, s1_len);
	}
	if (seq2) {
		tsta_msa_aligner_add(aligner, seq2, s2_len);
	}
	if (tsta_msa_aligner_get_result(aligner, &result) == 0) {
		for (size_t i = 0; i < result.sequence_count; i++) {
			printf("seq%zu: %s\n", i + 1, result.aln[i]);
		}
	}
}

tsta_msa_aligner_destroy(aligner);
tsta_msa_result_free(&result);
```

One-shot MSA with raw C strings:
```c
const char *seqs[] = {seq0, seq1, seq2};
int lengths[] = {(int)strlen(seq0), (int)strlen(seq1), (int)strlen(seq2)};
tsta_config config = tsta_config_make_default();
tsta_msa_result_t result = tsta_msa_result_make();

if (tsta_msa_align(seqs, lengths, 3, &config, &result) == 0) {
	for (size_t i = 0; i < result.sequence_count; i++) {
		printf("seq%zu: %s\n", i + 1, result.aln[i]);
	}
}
tsta_msa_result_free(&result);
```

Structured MSA result (each sequence alignment in final MSA):
```c
#include "tsta.h"
#include <string.h>

const char *seqs[] = {seq0, seq1, seq2};
int lengths[] = {(int)strlen(seq0), (int)strlen(seq1), (int)strlen(seq2)};
tsta_config config = tsta_config_make_default();
tsta_msa_result_t result = tsta_msa_result_make();

if (tsta_msa_align(seqs, lengths, 3, &config, &result) == 0) {
	for (size_t i = 0; i < result.sequence_count; i++) {
		printf("seq%zu: %s\n", i + 1, result.aln[i]);
	}
}

tsta_msa_result_free(&result);
```
## Usage and run
### Pairwise sequence alignment
```bash
./TSTA_psa -1 ./example/psa/seq/seqa1.fa -2 ./example/psa/seq/seqb1.fa
# maxsorce=-5
```
### no-backtracing PSA
```bash
./TSTA_psa_notrace -1 ./example/psa/seq/seqa1.fa -2 ./example/psa/seq/seqb1.fa
# maxsorce=-5
```
### Multiple sequence alignment
```bash
./TSTA_msa -i ./example/msa/seq/seq1.fa
# seq:5000, seqN:5120, poa:5000
# poa_len=5000 ,seq_len=5000 ,trace_sub:[num1]=4999 [num2]=4999 ,lastsorce=-5451
# poa_add_len:1885
# seq:5000, seqN:5120, poa:6885
# poa_len=6885 ,seq_len=5000 ,trace_sub:[num1]=6884 [num2]=4999 ,lastsorce=-3101
# poa_add_len:1714
# seq:5000, seqN:5120, poa:8599
# poa_len=8599 ,seq_len=5000 ,trace_sub:[num1]=8598 [num2]=4999 ,lastsorce=-1776
# poa_add_len:1531
# seq:5000, seqN:5120, poa:10130
# poa_len=10130 ,seq_len=5000 ,trace_sub:[num1]=10129 [num2]=4999 ,lastsorce=-870
# poa_add_len:1338
```

> [!Tip]
> While you are running the program in msa mode, we recommend that the sequences used for alignment should be less than 200,000bp and more than 50,000bp. If the sequences are too long, the program may run out of memory.


## Operating system
The Linux platform is supported by TSTA.It run and tested on CentOS Linux 7.
## Contact
- Peiyu Zong peiyuzong8@gmail.com (Designed and implemented core algorithms)
- Wenpeng Deng 1732889554@qq.com (Optimized some of the code) 
