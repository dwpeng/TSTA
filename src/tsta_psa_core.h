#ifndef TSTA_PSA_CORE_H
#define TSTA_PSA_CORE_H

/* ── Shared PSA SIMD core ───────────────────────────────────────────────
 *
 * The reusable pairwise-alignment machinery shared by the threaded driver
 * (tsta_psa.c) and the standalone single-threaded module
 * (tsta_psa_simd.c). Header-only (static inline) so each translation unit
 * compiles exactly what it uses; there is no separate object.
 *
 * These functions are the original tsta_psa.c internals, unchanged except:
 *   - the block computation takes the workspace explicitly (no TLS) and
 *     does not update state->ms (dead; the score is recomputed from the
 *     trace in tsta_psa_capture_result),
 *   - the sequence-cache helpers take the storage container as parameters.
 */

#include "tsta_common.h"

#include <stdio.h>

#define MIN -100
#define I_MIN -2000000000

/* Forward declarations (used by the runtime-buffer helpers below). */
static inline void tsta_psa_trace_context_reset(tsta_psa_trace_context* trace);
static inline void
tsta_psa_trace_context_detach(tsta_psa_trace_context* trace);

/* ── Block workspace ─────────────────────────────────────────────────── */

typedef struct tsta_psa_workspace {
  int lane_len;
  int block_width;
  int* maxsorce;
  char* h_s;
  char* rf;
  char* t_temp;
  char* e_temp;
  char* q_temp;
  int* r_temp;
  char* source[3];
} tsta_psa_workspace;

static inline void
tsta_psa_workspace_reset(tsta_psa_workspace* ws)
{
  if (!ws)
    return;
  free(ws->maxsorce);
  ws->maxsorce = NULL;
  mm_free(ws->h_s);
  ws->h_s = NULL;
  mm_free(ws->rf);
  ws->rf = NULL;
  mm_free(ws->t_temp);
  ws->t_temp = NULL;
  mm_free(ws->e_temp);
  ws->e_temp = NULL;
  mm_free(ws->q_temp);
  ws->q_temp = NULL;
  mm_free(ws->r_temp);
  ws->r_temp = NULL;
  for (int i = 0; i < 3; i++) {
    mm_free(ws->source[i]);
    ws->source[i] = NULL;
  }
  ws->lane_len = 0;
  ws->block_width = 0;
}

static inline int
tsta_psa_workspace_resize(tsta_psa_workspace* ws,
                          int lane_len,
                          int block_width)
{
  if (!ws || lane_len <= 0 || block_width <= 0)
    return -1;

  if (ws->lane_len == lane_len && ws->block_width == block_width
      && ws->maxsorce && ws->h_s && ws->rf && ws->t_temp && ws->e_temp
      && ws->q_temp && ws->r_temp && ws->source[0] && ws->source[1]
      && ws->source[2])
    return 0;

  tsta_psa_workspace_reset(ws);
  ws->maxsorce = (int*)malloc((size_t)lane_len * sizeof(int));
  ws->h_s = (char*)mm_malloc((size_t)lane_len * sizeof(char));
  ws->rf = (char*)mm_malloc((size_t)block_width * sizeof(char));
  ws->t_temp = (char*)mm_malloc((size_t)lane_len * sizeof(char));
  ws->e_temp = (char*)mm_malloc((size_t)lane_len * sizeof(char));
  ws->q_temp = (char*)mm_malloc((size_t)lane_len * sizeof(char));
  ws->r_temp = (int*)mm_malloc((size_t)lane_len * sizeof(int));
  for (int i = 0; i < 3; i++)
    ws->source[i] = (char*)mm_malloc((size_t)lane_len * sizeof(char));

  if (!ws->maxsorce || !ws->h_s || !ws->rf || !ws->t_temp || !ws->e_temp
      || !ws->q_temp || !ws->r_temp || !ws->source[0] || !ws->source[1]
      || !ws->source[2]) {
    tsta_psa_workspace_reset(ws);
    return -1;
  }

  ws->lane_len = lane_len;
  ws->block_width = block_width;
  return 0;
}

/* ── Runtime DP buffers (cached across aligns) ───────────────────────── */

typedef struct tsta_psa_runtime_buffers {
  char* cached_sorce;
  char* cached_esorce;
  int* cached_real;
  char* cached_V;
  char* cached_F;
  size_t cached_main_len;
  size_t cached_aux_len;
  /* Trace matrix cache (grow-only across aligns). Owned by the cache; the
   * aligner's state->trace.matrices borrows these pointers during an align. */
  char** cached_trace_back;
  size_t cached_trace_rows;
  size_t cached_trace_cols;
} tsta_psa_runtime_buffers;

static inline void
tsta_psa_runtime_buffers_release(tsta_psa_runtime_buffers* buffers,
                                 tsta_psa_state* state,
                                 int release_cache)
{
  if (!buffers)
    return;

  state->seq = NULL;
  state->sorce = NULL;
  state->esorce = NULL;
  state->real = NULL;
  state->F = NULL;
  state->V = NULL;
  /* Detach (do not free) the trace matrices — the cache owns them. */
  tsta_psa_trace_context_detach(&state->trace);

  if (!release_cache)
    return;

  free(buffers->cached_sorce);
  free(buffers->cached_esorce);
  free(buffers->cached_real);
  free(buffers->cached_F);
  free(buffers->cached_V);
  buffers->cached_sorce = NULL;
  buffers->cached_esorce = NULL;
  buffers->cached_real = NULL;
  buffers->cached_F = NULL;
  buffers->cached_V = NULL;
  buffers->cached_main_len = 0;
  buffers->cached_aux_len = 0;

  if (buffers->cached_trace_back) {
    for (size_t i = 0; i < buffers->cached_trace_rows; i++)
      mm_free(buffers->cached_trace_back[i]);
    free(buffers->cached_trace_back);
    buffers->cached_trace_back = NULL;
  }
  buffers->cached_trace_rows = 0;
  buffers->cached_trace_cols = 0;
}

static inline int
tsta_psa_runtime_buffers_ensure(tsta_psa_runtime_buffers* buffers,
                                tsta_psa_state* state,
                                size_t main_len,
                                size_t aux_len)
{
  char *new_sorce, *new_esorce, *new_F, *new_V;
  int* new_real;
  char *old_sorce, *old_esorce, *old_F, *old_V;
  int* old_real;

  if (!buffers || main_len == 0 || aux_len == 0)
    return -1;

  if (buffers->cached_main_len >= main_len
      && buffers->cached_aux_len >= aux_len && buffers->cached_sorce
      && buffers->cached_esorce && buffers->cached_real && buffers->cached_F
      && buffers->cached_V) {
    state->sorce = buffers->cached_sorce;
    state->esorce = buffers->cached_esorce;
    state->real = buffers->cached_real;
    state->F = buffers->cached_F;
    state->V = buffers->cached_V;
    return 0;
  }

  new_sorce = (char*)malloc(main_len);
  new_esorce = (char*)malloc(main_len);
  new_real = (int*)malloc(main_len * sizeof(int));
  new_F = (char*)malloc(aux_len);
  new_V = (char*)malloc(aux_len);
  if (!new_sorce || !new_esorce || !new_real || !new_F || !new_V) {
    free(new_sorce);
    free(new_esorce);
    free(new_real);
    free(new_F);
    free(new_V);
    return -1;
  }

  old_sorce = buffers->cached_sorce;
  old_esorce = buffers->cached_esorce;
  old_real = buffers->cached_real;
  old_F = buffers->cached_F;
  old_V = buffers->cached_V;

  buffers->cached_sorce = new_sorce;
  buffers->cached_esorce = new_esorce;
  buffers->cached_real = new_real;
  buffers->cached_F = new_F;
  buffers->cached_V = new_V;
  buffers->cached_main_len = main_len;
  buffers->cached_aux_len = aux_len;

  free(old_sorce);
  free(old_esorce);
  free(old_real);
  free(old_F);
  free(old_V);

  state->sorce = buffers->cached_sorce;
  state->esorce = buffers->cached_esorce;
  state->real = buffers->cached_real;
  state->F = buffers->cached_F;
  state->V = buffers->cached_V;
  return 0;
}

/* ── Trace packing (SIMD byte-shuffle lookup) ────────────────────────── */

/* Pack back/eback/fback (one byte per cell) into a single byte:
 *   bits 0-1: back  (0/1/2)
 *   bits 2-4: eback + 2
 *   bits 5-7: fback + 2
 * vpshufb tables indexed by (value+2)&7 give the pre-shifted field bits. */
static inline void
tsta_psa_pack_trace_row(
    char* dst, const char* s0, const char* s1, const char* s2, int len)
{
  static const unsigned char T2[64] __attribute__((aligned(64))) = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
  };
  static const unsigned char T5[64] __attribute__((aligned(64))) = {
    0, 32, 64, 96, 128, 160, 192, 224, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 32, 64, 96, 128, 160, 192, 224, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 32, 64, 96, 128, 160, 192, 224, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 32, 64, 96, 128, 160, 192, 224, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  __mxxxi vt2 = mm_load((__mxxxi*)T2);
  __mxxxi vt5 = mm_load((__mxxxi*)T5);
  __mxxxi three = mm_set1_epi8(3);
  __mxxxi two = mm_set1_epi8(2);
  __mxxxi seven = mm_set1_epi8(7);
  int k;
  for (k = 0; k + block <= len; k += block) {
    __mxxxi a = mm_load((__mxxxi*)(s0 + k));
    __mxxxi b = mm_load((__mxxxi*)(s1 + k));
    __mxxxi c = mm_load((__mxxxi*)(s2 + k));
    __mxxxi back = mm_and_si(a, three);
    __mxxxi ebm = mm_and_si(mm_add_epi8(c, two), seven);
    __mxxxi fbm = mm_and_si(mm_add_epi8(b, two), seven);
    __mxxxi out = mm_or_epi8(back, mm_or_epi8(mm_shuffle_epi8(vt2, ebm),
                                              mm_shuffle_epi8(vt5, fbm)));
    mm_store((__mxxxi*)(dst + k), out);
  }
  for (; k < len; k++) {
    unsigned char b = (unsigned char)s0[k] & 3u;
    unsigned char eb = ((unsigned char)s2[k] + 2u) & 7u;
    unsigned char fb = ((unsigned char)s1[k] + 2u) & 7u;
    dst[k] = (char)(b | (eb << 2) | (fb << 5));
  }
}

/* ── Core SIMD PSA inner loop (UNCHANGED) ────────────────────────────── */

static inline void
tsta_psa_row(int* maxsorce,
             int y,
             int block_i,
             int block_l,
             int pc1,
             int pc2,
             int pc4,
             char* h_s,
             char* t_temp,
             char* e_temp,
             char* q_temp,
             char* rf,
             int* r_temp,
             char** source,
             tsta_psa_state* state)
{
  __mxxxi zero, s1, s2, Smin, h, b2, e, ev, f, fv, t, s, v, mat, mis, egap,
      ogap, v1, h1, trace, etrace, ftrace, temp1;
  __mask mask, mask1;
  int j = 0;
  zero = mm_setzero();
  egap = mm_set1_epi8(state->E);
  ogap = mm_set1_epi8(state->O + state->E);
  mis = mm_set1_epi8(state->X);
  mat = mm_set1_epi8(state->M);
  Smin = mm_set1_epi8(MIN);
  f = Smin;
  f = mm_insert_epi8(f, state->F[pc2], 0);
  b2 = mm_set1_epi8(state->seq[1][pc2]);

  for (int x = 0; x < state->W; x++) {
    h = mm_load(((__mxxxi*)q_temp) + x);
    mask = mm_cmpeq_epi8(h, b2);
    h = mm_blendv_epi8(mis, mat, mask);
    mm_store(((__mxxxi*)h_s) + x, h);
    t = mm_load(((__mxxxi*)t_temp) + x);
    e = mm_load(((__mxxxi*)e_temp) + x);
    s = mm_max_epi8(h, e);
    s = mm_max_epi8(s, f);
    f = mm_add_epi8(f, egap);
    h1 = mm_add_epi8(s, ogap);
    f = mm_max_epi8(f, h1);
    f = mm_sub_epi8(f, t);
  }

  mm_store((__mxxxi*)rf, f);
  for (int x = 1; x < state->B; x++)
    if (rf[x - 1] + state->W * state->E
            - (r_temp[state->L - state->B + x]
               - r_temp[state->L - state->B + x - 1])
        > rf[x])
      rf[x] = rf[x - 1] + state->W * state->E
              - (r_temp[state->L - state->B + x]
                 - r_temp[state->L - state->B + x - 1]);

  f = mm_load((__mxxxi*)rf);
  temp1 = mm_sub_epi8(f, egap);
  f = mm_slli(f);
  f = mm_insert_epi8(f, state->F[pc2], 0);

  v = mm_sub_epi8(s, t);
  v = mm_max_epi8(temp1, v);
  v = mm_slli(v);
  v = mm_insert_epi8(v, state->V[pc2], 0);

  s1 = mm_set1_epi8(1);
  s2 = mm_add_epi8(s1, s1);
  b2 = mm_set1_epi32(I_MIN);
  for (int x = 0; x < state->W; x++) {
    h1 = mm_load(((__mxxxi*)h_s) + x);
    t = mm_load(((__mxxxi*)t_temp) + x);
    e = mm_load(((__mxxxi*)e_temp) + x);
    s = mm_max_epi8(e, f);
    s = mm_max_epi8(s, h1);
    h = mm_sub_epi8(s, v);
    mm_store(((__mxxxi*)t_temp) + x, h);

    trace = s2;
    mask = mm_cmpeq_epi8(s, f);
    trace = mm_blendv_epi8(trace, zero, mask);
    mask = mm_cmpeq_epi8(s, h1);
    trace = mm_blendv_epi8(trace, s1, mask);
    mm_store(((__mxxxi*)source[0]) + x, trace);

    h1 = mm_add_epi8(s, ogap);
    fv = mm_sub_epi8(f, v);
    f = mm_add_epi8(f, egap);
    mask1 = mm_cmpeq_epi8(f, h1);
    f = mm_max_epi8(f, h1);
    f = mm_sub_epi8(f, t);

    mask = mm_cmpeq_epi8(fv, ogap);
    ftrace = mm_blendv_epi8(s1, s2, mask);
    temp1 = mm_sub_epi8(zero, ftrace);
    mask = mm_and_epi8(mask, mask1);
    ftrace = mm_blendv_epi8(ftrace, temp1, mask);
    mm_store(((__mxxxi*)source[1]) + x, ftrace);

    ev = mm_sub_epi8(e, t);
    e = mm_add_epi8(e, egap);
    mask1 = mm_cmpeq_epi8(e, h1);
    e = mm_max_epi8(e, h1);
    e = mm_sub_epi8(e, v);
    mm_store(((__mxxxi*)e_temp) + x, e);

    mask = mm_cmpeq_epi8(ev, ogap);
    etrace = mm_blendv_epi8(s1, s2, mask);
    temp1 = mm_sub_epi8(zero, etrace);
    mask = mm_and_epi8(mask, mask1);
    etrace = mm_blendv_epi8(etrace, temp1, mask);
    mm_store(((__mxxxi*)source[2]) + x, etrace);

    v = mm_sub_epi8(s, t);
    v1 = mm0_epi8cvt32(v);
    h1 = mm_load(((__mxxxi*)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi*)r_temp) + j, h1);
    j++;
    v1 = mm_epi8cvt32(v, 1);
    h1 = mm_load(((__mxxxi*)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi*)r_temp) + j, h1);
    j++;
    v1 = mm_epi8cvt32(v, 2);
    h1 = mm_load(((__mxxxi*)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi*)r_temp) + j, h1);
    j++;
    v1 = mm_epi8cvt32(v, 3);
    h1 = mm_load(((__mxxxi*)r_temp) + j);
    h1 = mm_add_epi32(v1, h1);
    b2 = mm_max_epi32(b2, h1);
    mm_store(((__mxxxi*)r_temp) + j, h1);
    j++;
  }
  state->F[pc2] = mm_extract_epi8(f);
  state->V[pc2] = mm_extract_epi8(v);
  maxsorce[y] = mm_reduce_max_epi32(b2);
  if (state->trace.matrices.back) {
    /* Pack back (bits 0-1, values 0/1/2), eback (bits 2-4, value+2) and
     * fback (bits 5-7, value+2) into one byte per cell. The traceback
     * matrix is one third of the size of three separate matrices. */
    tsta_psa_pack_trace_row(state->trace.matrices.back[pc2] + pc4, source[0],
                            source[1], source[2], state->L);
  }
}

/* ── DP matrix init (UNCHANGED) ──────────────────────────────────────── */

static inline void
tsta_psa_blockmatrix_init(tsta_psa_state* state)
{
  memset(state->sorce, state->E, state->length[0]);
  memset(state->esorce, state->E + state->E + state->O, state->length[0]);
  state->sorce[0] = state->E + state->O;
  state->esorce[0] = 2 * (state->E + state->O);
  for (int i = 0; i < state->length[0]; i++)
    state->real[i] =
        state->O
        + ((i / state->L * state->L + ((i % state->L) % state->B) * state->W
            + ((i % state->L) / state->B))
           + 1)
              * state->E;

  memset(state->V, state->E, state->length[3]);
  memset(state->F, state->E + state->E + state->O, state->length[3]);
  state->V[0] = state->E + state->O;
  state->F[0] = 2 * (state->E + state->O);
}

/* ── Block computation (UNCHANGED; workspace passed explicitly) ──────── */

static inline void
tsta_psa_block(tsta_psa_state* state,
               tsta_psa_workspace* ws,
               int block_i,
               int block_l)
{
  int* maxsorce = ws->maxsorce;
  int pc0, pc1, pc2, pc4;
  char* h_s = ws->h_s;
  char* rf = ws->rf;
  char* t_temp = ws->t_temp;
  char* e_temp = ws->e_temp;
  char* q_temp = ws->q_temp;
  int* r_temp = ws->r_temp;
  char* source[3] = { ws->source[0], ws->source[1], ws->source[2] };

  if (block_i <= state->lmaxtag)
    pc0 = block_i - block_l;
  else
    pc0 = state->lmaxtag - block_l;
  pc1 = pc0 * state->W;
  pc4 = pc0 * state->L;

  memcpy(t_temp, state->sorce + pc4, state->L);
  memcpy(e_temp, state->esorce + pc4, state->L);
  memcpy(r_temp, state->real + pc4, state->L * sizeof(int));
  for (int i = 0; i < state->L; i++)
    q_temp[i] =
        state->seq[0][pc4 + (i % state->B) * state->W + (i / state->B)];

  for (int i = 0; i < state->L; i++) {
    if (block_i <= state->lmaxtag)
      pc2 = block_l * state->L + i;
    else
      pc2 = (block_l + block_i - state->lmaxtag) * state->L + i;
    if (pc2 >= state->length[3]) {
      for (int s = i; s < state->L; s++)
        maxsorce[s] = I_MIN;
      break;
    }
    tsta_psa_row(maxsorce, i, block_i, block_l, pc1, pc2, pc4, h_s, t_temp,
                 e_temp, q_temp, rf, r_temp, source, state);
  }
  memcpy(state->sorce + pc4, t_temp, state->L);
  memcpy(state->esorce + pc4, e_temp, state->L);
  memcpy(state->real + pc4, r_temp, state->L * sizeof(int));

  /* state->ms is dead (score recomputed from the trace); dropped. */
}

/* ── Trace helpers (UNCHANGED) ───────────────────────────────────────── */

static inline int
tsta_psa_trace_index(const tsta_psa_state* state, int j)
{
  int lane_len = state->length[0];
  return (j / lane_len) * lane_len
         + ((((j % lane_len) % state->W) * state->B)
            + ((j % lane_len) / state->W));
}

static inline void
tsta_psa_trace_result_reset(tsta_psa_trace_result* result)
{
  if (!result)
    return;
  free(result->aligned_seq1);
  free(result->aligned_seq2);
  free(result->formatted_text);
  result->aligned_seq1 = NULL;
  result->aligned_seq2 = NULL;
  result->formatted_text = NULL;
  result->aligned_length = 0;
  result->formatted_length = 0;
}

static inline void
tsta_psa_trace_context_reset(tsta_psa_trace_context* trace)
{
  if (!trace)
    return;

  if (trace->matrices.back) {
    for (int i = 0; i < trace->matrices.rows; i++)
      mm_free(trace->matrices.back[i]);
  }

  free(trace->matrices.back);
  trace->matrices.back = NULL;
  trace->matrices.rows = 0;
  trace->matrices.cols = 0;

  tsta_psa_trace_result_reset(&trace->result);
}

/* Detach the matrices from the trace context without freeing them. The
 * caller (tsta_psa_runtime_buffers) keeps ownership in its cache. */
static inline void
tsta_psa_trace_context_detach(tsta_psa_trace_context* trace)
{
  if (!trace)
    return;
  trace->matrices.back = NULL;
  trace->matrices.rows = 0;
  trace->matrices.cols = 0;
  tsta_psa_trace_result_reset(&trace->result);
}

/* Grow-only trace matrix cache: reuse the aligner's matrices across aligns,
 * reallocating only when the required size grows. The DP writes every cell
 * each align, so reuse introduces no stale reads.
 *
 * Only ONE matrix is allocated: the packed byte per cell stores back (2 bits),
 * eback (3 bits) and fback (3 bits), so the traceback matrix is one third of
 * the size of three separate matrices. */
static inline int
tsta_psa_trace_context_ensure(tsta_psa_runtime_buffers* buffers,
                              tsta_psa_state* state)
{
  tsta_psa_trace_matrices* m = &state->trace.matrices;
  size_t rows = (size_t)state->length[3]; /* real shorter length */
  size_t cols = (size_t)state->length[0]; /* padded longer length */
  char** back = NULL;
  size_t i;

  if (!buffers || rows == 0 || cols == 0)
    return -1;

  if (buffers->cached_trace_rows >= rows && buffers->cached_trace_cols >= cols
      && buffers->cached_trace_back) {
    m->back = buffers->cached_trace_back;
    m->rows = (int)rows;
    m->cols = (int)cols;
    return 0;
  }

  back = (char**)calloc(rows, sizeof(char*));
  if (!back)
    goto fail;
  for (i = 0; i < rows; i++) {
    back[i] = (char*)mm_malloc(cols);
    if (!back[i]) {
      for (size_t k = 0; k < i; k++)
        mm_free(back[k]);
      free(back);
      back = NULL;
      goto fail;
    }
  }

  /* Drop the old cache and adopt the new one. */
  if (buffers->cached_trace_back) {
    for (i = 0; i < buffers->cached_trace_rows; i++)
      mm_free(buffers->cached_trace_back[i]);
    free(buffers->cached_trace_back);
  }
  buffers->cached_trace_back = back;
  buffers->cached_trace_rows = rows;
  buffers->cached_trace_cols = cols;

  m->back = back;
  m->rows = (int)rows;
  m->cols = (int)cols;
  return 0;

fail:
  free(back);
  m->back = NULL;
  m->rows = 0;
  m->cols = 0;
  return -1;
}

static inline int
tsta_psa_trace_build_alignment(tsta_psa_state* state)
{
  tsta_psa_trace_result* result = &state->trace.result;
  int i = state->length[3] - 1;
  int j = state->length[2] - 1;
  int insertion_count = 0;
  size_t write_pos;

  tsta_psa_trace_result_reset(result);

  while (i >= 0 && j >= 0) {
    int current = tsta_psa_trace_index(state, j);
    unsigned char p = (unsigned char)state->trace.matrices.back[i][current];
    int backv = (int)(p & 3);
    int fbv = (int)((p >> 5) & 7) - 2;
    int ebv = (int)((p >> 2) & 7) - 2;
    if (backv == 1) {
      i--;
      j--;
    } else if (backv == 0) {
      if (j - 1 >= 0
          && ((fbv == 1 || fbv == -1)
              || ((fbv == 2 || fbv == -2)
                  && (int)(((unsigned char)state->trace.matrices
                                .back[i][tsta_psa_trace_index(state, j - 1)]
                            >> 5)
                           & 7)
                             - 2
                         < 0))) {
        unsigned char* prevp =
            (unsigned char*)&state->trace.matrices
                .back[i][tsta_psa_trace_index(state, j - 1)];
        *prevp = (unsigned char)(*prevp & ~3u);
      }
      j--;
    } else {
      if (i - 1 >= 0
          && ((ebv == 1 || ebv == -1)
              || ((ebv == 2 || ebv == -2)
                  && (int)(((unsigned char)
                                state->trace.matrices.back[i - 1][current]
                            >> 2)
                           & 7)
                             - 2
                         < 0))) {
        unsigned char* cp =
            (unsigned char*)&state->trace.matrices.back[i - 1][current];
        *cp = (unsigned char)((*cp & ~3u) | 2u);
      }
      i--;
      insertion_count++;
    }
  }
  if (i >= 0)
    insertion_count += i + 1;

  result->aligned_length = (size_t)state->length[2] + (size_t)insertion_count;
  result->aligned_seq1 = (char*)malloc(result->aligned_length + 1);
  result->aligned_seq2 = (char*)malloc(result->aligned_length + 1);
  if (!result->aligned_seq1 || !result->aligned_seq2) {
    tsta_psa_trace_result_reset(result);
    return -1;
  }

  i = state->length[3] - 1;
  j = state->length[2] - 1;
  write_pos = result->aligned_length;
  while (i >= 0 && j >= 0) {
    write_pos--;
    int current = tsta_psa_trace_index(state, j);
    int backv =
        (int)((unsigned char)state->trace.matrices.back[i][current] & 3);
    if (backv == 1) {
      result->aligned_seq1[write_pos] = state->seq[0][j];
      result->aligned_seq2[write_pos] = state->seq[1][i];
      i--;
      j--;
    } else if (backv == 0) {
      result->aligned_seq1[write_pos] = state->seq[0][j];
      result->aligned_seq2[write_pos] = '-';
      j--;
    } else {
      result->aligned_seq1[write_pos] = '-';
      result->aligned_seq2[write_pos] = state->seq[1][i];
      i--;
    }
  }
  while (j >= 0) {
    write_pos--;
    result->aligned_seq1[write_pos] = state->seq[0][j];
    result->aligned_seq2[write_pos] = '-';
    j--;
  }
  while (i >= 0) {
    write_pos--;
    result->aligned_seq2[write_pos] = state->seq[1][i];
    result->aligned_seq1[write_pos] = '-';
    i--;
  }

  result->aligned_seq1[result->aligned_length] = '\0';
  result->aligned_seq2[result->aligned_length] = '\0';
  return 0;
}

static inline int
tsta_psa_build_cigar(const tsta_psa_trace_result* trace,
                     tsta_psa_result_t* result)
{
  size_t i;
  size_t write_pos = 0;
  char prev_op = '\0';
  unsigned int run_len = 0;
  size_t capacity;
  char* cigar;

  if (!trace || !result || !trace->aligned_seq1 || !trace->aligned_seq2)
    return -1;

  capacity = trace->aligned_length > 0 ? trace->aligned_length * 12 + 1 : 2;
  cigar = (char*)malloc(capacity);
  if (!cigar)
    return -1;

  for (i = 0; i < trace->aligned_length; i++) {
    char a = trace->aligned_seq1[i];
    char b = trace->aligned_seq2[i];
    char op;

    if (a == '-' && b != '-') {
      result->ins_count++;
      op = 'I';
    } else if (a != '-' && b == '-') {
      result->del_count++;
      op = 'D';
    } else if (a == b) {
      result->match_count++;
      op = 'M';
    } else {
      result->mismatch_count++;
      op = 'X';
    }

    if (prev_op == '\0') {
      prev_op = op;
      run_len = 1;
      continue;
    }
    if (op == prev_op) {
      run_len++;
      continue;
    }

    write_pos += (size_t)snprintf(cigar + write_pos, capacity - write_pos,
                                  "%u%c", run_len, prev_op);
    prev_op = op;
    run_len = 1;
  }

  if (prev_op != '\0') {
    write_pos += (size_t)snprintf(cigar + write_pos, capacity - write_pos,
                                  "%u%c", run_len, prev_op);
  } else {
    cigar[0] = '0';
    cigar[1] = 'M';
    write_pos = 2;
  }

  cigar[write_pos] = '\0';
  result->cigar = cigar;
  return 0;
}

static inline int
tsta_psa_capture_result(tsta_psa_result_t* result, tsta_psa_state* state)
{
  tsta_psa_trace_result* trace;

  if (!result || !state)
    return -1;
  if (tsta_psa_trace_build_alignment(state) != 0)
    return -1;

  trace = &state->trace.result;
  tsta_psa_result_free(result);
  result->aln_length = trace->aligned_length;
  result->aln[0] = (char*)malloc(trace->aligned_length + 1);
  result->aln[1] = (char*)malloc(trace->aligned_length + 1);
  if (!result->aln[0] || !result->aln[1]) {
    tsta_psa_result_free(result);
    return -1;
  }

  memcpy(result->aln[0], trace->aligned_seq1, trace->aligned_length + 1);
  memcpy(result->aln[1], trace->aligned_seq2, trace->aligned_length + 1);

  if (tsta_psa_build_cigar(trace, result) != 0) {
    tsta_psa_result_free(result);
    return -1;
  }

  /* Compute score from aligned sequences (state->ms tracks the global DP
   * max, not the end-cell score required for global alignment). */
  {
    int score = 0, in_gap1 = 0, in_gap2 = 0;
    for (size_t i = 0; i < trace->aligned_length; i++) {
      char a = trace->aligned_seq1[i];
      char b = trace->aligned_seq2[i];
      if (a == '-') {
        score += state->E;
        if (!in_gap1) {
          score += state->O;
          in_gap1 = 1;
        }
        in_gap2 = 0;
      } else if (b == '-') {
        score += state->E;
        if (!in_gap2) {
          score += state->O;
          in_gap2 = 1;
        }
        in_gap1 = 0;
      } else {
        score += (a == b) ? state->M : state->X;
        in_gap1 = in_gap2 = 0;
      }
    }
    result->score = score;
  }
  return 0;
}

/* ── Sequence preparation (UNCHANGED; cache container passed in) ─────── */

static inline char*
tsta_psa_make_padded_copy(const char* sequence_data,
                          size_t sequence_length,
                          int padded_length)
{
  char* copy;
  size_t length = sequence_length;

  copy = (char*)malloc((size_t)padded_length + 1);
  if (!copy)
    return NULL;

  memset(copy, 'N', (size_t)padded_length);
  if (length > 0)
    memcpy(copy, sequence_data, length);
  copy[padded_length] = '\0';
  return copy;
}

static inline int
tsta_psa_ensure_cached_sequence(packed_sequence_array_t* cache,
                                int slot,
                                int padded_length,
                                const char* sequence_data,
                                size_t sequence_length)
{
  packed_sequence_t* slot_entry;
  char* buffer;
  int required_capacity;

  if (!cache || slot < 0 || slot >= 2 || !sequence_data || padded_length < 0)
    return -1;

  if (tsta_packed_sequence_array_reserve_slot(cache, (size_t)slot, 2) != 0)
    return -1;
  slot_entry = &cache->data[slot];

  required_capacity = padded_length + 1;
  if (slot_entry->len < required_capacity) {
    buffer = tsta_psa_make_padded_copy(sequence_data, sequence_length,
                                       padded_length);
    if (!buffer)
      return -1;
    free(slot_entry->packed_seq);
    slot_entry->packed_seq = buffer;
    slot_entry->len = required_capacity;
    return 0;
  }

  buffer = slot_entry->packed_seq;
  memset(buffer, 'N', (size_t)padded_length);
  if (sequence_length > 0)
    memcpy(buffer, sequence_data, sequence_length);
  buffer[padded_length] = '\0';
  return 0;
}

static inline int
tsta_psa_prepare_sequences(const char* sequence1_data,
                           size_t sequence1_length,
                           const char* sequence2_data,
                           size_t sequence2_length,
                           packed_sequence_array_t* cache,
                           char** slots,
                           tsta_psa_state* state)
{
  const char* longer_data = sequence1_data;
  const char* shorter_data = sequence2_data;
  size_t longer_length = sequence1_length;
  size_t shorter_length = sequence2_length;

  if (sequence2_length > sequence1_length) {
    longer_data = sequence2_data;
    shorter_data = sequence1_data;
    longer_length = sequence2_length;
    shorter_length = sequence1_length;
  }

  state->length[2] = (int)longer_length;
  state->length[3] = (int)shorter_length;
  state->length[0] = (int)longer_length;
  state->length[1] = (int)shorter_length;
  if (state->length[0] % state->L != 0)
    state->length[0] += state->L - state->length[0] % state->L;
  if (state->length[1] % state->L != 0)
    state->length[1] += state->L - state->length[1] % state->L;

  if (tsta_psa_ensure_cached_sequence(cache, 0, state->length[0], longer_data,
                                      longer_length)
      != 0)
    return -1;
  if (tsta_psa_ensure_cached_sequence(cache, 1, state->length[1], shorter_data,
                                      shorter_length)
      != 0)
    return -1;

  slots[0] = cache->data[0].packed_seq;
  slots[1] = cache->data[1].packed_seq;
  state->seq = slots;
  return 0;
}

/* ── Shared align preamble / epilogue ──────────────────────────────────
 * Used by both the threaded driver (tsta_psa.c) and the standalone module
 * (tsta_psa_simd.c); only the DP diagonal loop differs between them. */

static inline int
tsta_psa_align_setup(tsta_psa_state* state,
                     const tsta_config* config,
                     const char* sequence1,
                     int sequence1_length,
                     const char* sequence2,
                     int sequence2_length,
                     packed_sequence_array_t* cache,
                     char** slots,
                     tsta_psa_runtime_buffers* buffers)
{
  tsta_sequence_view_t left, right;
  tsta_config effective;
  size_t longest;

  if (tsta_sequence_view_make(&left, sequence1, sequence1_length) != 0
      || tsta_sequence_view_make(&right, sequence2, sequence2_length) != 0)
    return -1;

  /* Pick the lane length to match the longer sequence: L >= max(len1, len2)
   * rounded up to a SIMD width. This keeps padding minimal for the longer
   * sequence while a large L keeps the SIMD batch efficient (measured:
   * large L beats small L even with more padding). block_size acts as a
   * minimum; the sequence length only raises it. */
  effective = config ? *config : tsta_config_make_default();
  longest = (size_t)(left.length > right.length ? left.length : right.length);
  if (effective.block_size <= 0)
    effective.block_size = 1; /* auto: let the sequence length decide */
  {
    size_t safe = (longest + (size_t)block - 1) / (size_t)block;
    if (safe == 0)
      safe = 1;
    if ((size_t)effective.block_size < safe)
      effective.block_size = (int)safe;
  }
  tsta_init_psa_state(state, &effective, block);

  if (tsta_psa_prepare_sequences(left.sequence, (size_t)left.length,
                                 right.sequence, (size_t)right.length, cache,
                                 slots, state)
      != 0)
    goto fail;

  state->fmaxtag = state->length[1] / state->L - 1;
  state->lmaxtag = state->length[0] / state->L - 1;

  if (tsta_psa_runtime_buffers_ensure(buffers, state, (size_t)state->length[0],
                                      (size_t)state->length[3])
      != 0)
    goto fail;
  /* Traceback is always produced, so the trace matrix is always ensured. */
  if (tsta_psa_trace_context_ensure(buffers, state) != 0)
    goto fail;
  return 0;

fail:
  tsta_psa_runtime_buffers_release(buffers, state, 0);
  return -1;
}

/* One-shot helper: the one-shot path creates its threadpool once per call,
 * so a small thread count wins (spawn/join overhead exceeds the DP speedup —
 * measured: 2 threads fastest across 500..10000 bp). block_size is left to
 * tsta_psa_align_setup, which auto-matches it to the sequence length. Use a
 * reusable aligner to exploit more cores. */
static inline void
tsta_psa_adjust_config(tsta_config* config)
{
  if (config->threads <= 0)
    config->threads = 10;
  if (config->threads > 2)
    config->threads = 2;
  if (config->threads < 1)
    config->threads = 1;
}

static inline int
tsta_psa_align_finish(tsta_psa_runtime_buffers* buffers,
                      tsta_psa_state* state,
                      tsta_psa_result_t* result)
{
  int status = tsta_psa_capture_result(result, state) == 0 ? 0 : -1;
  tsta_psa_runtime_buffers_release(buffers, state, 0);
  return status;
}

#endif /* TSTA_PSA_CORE_H */
