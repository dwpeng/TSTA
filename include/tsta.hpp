#ifndef TSTA_HPP
#define TSTA_HPP

#include "tsta.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace tsta {

/* ── Exception ───────────────────────────────────────────────────────── */

class error : public std::runtime_error {
public:
  explicit error(const std::string& what) : std::runtime_error(what) {}
};

/* ── Configuration ───────────────────────────────────────────────────── */

class Config {
public:
  Config() { tsta_config_default(&cfg_); }

  /* chainable setters */
  Config&
  match(int v)
  {
    cfg_.match = v;
    return *this;
  }
  Config&
  mismatch(int v)
  {
    cfg_.mismatch = v;
    return *this;
  }
  Config&
  gap_open(int v)
  {
    cfg_.gap_open = v;
    return *this;
  }
  Config&
  gap_extend(int v)
  {
    cfg_.gap_extend = v;
    return *this;
  }
  Config&
  block_size(int v)
  {
    cfg_.block_size = v;
    return *this;
  }
  Config&
  threads(int v)
  {
    cfg_.threads = v;
    return *this;
  }

  /* accessors */
  int
  match() const
  {
    return cfg_.match;
  }
  int
  mismatch() const
  {
    return cfg_.mismatch;
  }
  int
  gap_open() const
  {
    return cfg_.gap_open;
  }
  int
  gap_extend() const
  {
    return cfg_.gap_extend;
  }
  int
  block_size() const
  {
    return cfg_.block_size;
  }
  int
  threads() const
  {
    return cfg_.threads;
  }

  const tsta_config*
  raw() const
  {
    return &cfg_;
  }
  tsta_config*
  raw()
  {
    return &cfg_;
  }

private:
  tsta_config cfg_;
};

/* ── PSA result ──────────────────────────────────────────────────────── */

class PsaResult {
public:
  PsaResult() { tsta_psa_result_init(&r_); }

  ~PsaResult() { tsta_psa_result_free(&r_); }

  PsaResult(const PsaResult&) = delete;
  PsaResult& operator=(const PsaResult&) = delete;

  PsaResult(PsaResult&& other) noexcept : r_(other.r_)
  {
    tsta_psa_result_init(&other.r_);
  }
  PsaResult&
  operator=(PsaResult&& other) noexcept
  {
    if (this != &other) {
      tsta_psa_result_free(&r_);
      r_ = other.r_;
      tsta_psa_result_init(&other.r_);
    }
    return *this;
  }

  int
  score() const
  {
    return r_.score;
  }
  int
  del_count() const
  {
    return r_.del_count;
  }
  int
  ins_count() const
  {
    return r_.ins_count;
  }
  int
  match_count() const
  {
    return r_.match_count;
  }
  int
  mismatch_count() const
  {
    return r_.mismatch_count;
  }
  size_t
  aln_length() const
  {
    return r_.aln_length;
  }

  std::string
  seq1() const
  {
    return r_.aln[0] ? std::string(r_.aln[0]) : "";
  }
  std::string
  seq2() const
  {
    return r_.aln[1] ? std::string(r_.aln[1]) : "";
  }
  std::string
  cigar() const
  {
    return r_.cigar ? std::string(r_.cigar) : "";
  }

  const char*
  seq1_cstr() const
  {
    return r_.aln[0];
  }
  const char*
  seq2_cstr() const
  {
    return r_.aln[1];
  }
  const char*
  cigar_cstr() const
  {
    return r_.cigar;
  }

  tsta_psa_result_t*
  raw()
  {
    return &r_;
  }

private:
  tsta_psa_result_t r_;
};

/* ── MSA result ──────────────────────────────────────────────────────── */

class MsaResult {
public:
  MsaResult() { tsta_msa_result_init(&r_); }

  ~MsaResult() { tsta_msa_result_free(&r_); }

  MsaResult(const MsaResult&) = delete;
  MsaResult& operator=(const MsaResult&) = delete;

  MsaResult(MsaResult&& other) noexcept : r_(other.r_)
  {
    tsta_msa_result_init(&other.r_);
  }
  MsaResult&
  operator=(MsaResult&& other) noexcept
  {
    if (this != &other) {
      tsta_msa_result_free(&r_);
      r_ = other.r_;
      tsta_msa_result_init(&other.r_);
    }
    return *this;
  }

  int
  score() const
  {
    return r_.score;
  }
  size_t
  sequence_count() const
  {
    return r_.sequence_count;
  }
  size_t
  aln_length() const
  {
    return r_.aln_length;
  }

  std::string
  sequence(size_t i) const
  {
    if (i >= r_.sequence_count)
      throw error("sequence index out of range");
    return r_.aln[i] ? std::string(r_.aln[i]) : "";
  }

  const char*
  sequence_cstr(size_t i) const
  {
    if (i >= r_.sequence_count)
      throw error("sequence index out of range");
    return r_.aln[i];
  }

  std::vector<std::string>
  sequences() const
  {
    std::vector<std::string> out;
    out.reserve(r_.sequence_count);
    for (size_t i = 0; i < r_.sequence_count; i++)
      out.emplace_back(r_.aln[i] ? r_.aln[i] : "");
    return out;
  }

  tsta_msa_result_t*
  raw()
  {
    return &r_;
  }

private:
  tsta_msa_result_t r_;
};

/* ── PSA aligner ─────────────────────────────────────────────────────── */

class PsaAligner {
public:
  explicit PsaAligner(const Config& config = Config{})
  {
    a_ = tsta_psa_aligner_create(config.raw());
    if (!a_)
      throw error("tsta_psa_aligner_create failed");
  }

  ~PsaAligner()
  {
    if (a_)
      tsta_psa_aligner_destroy(a_);
  }

  PsaAligner(const PsaAligner&) = delete;
  PsaAligner& operator=(const PsaAligner&) = delete;

  PsaAligner(PsaAligner&& other) noexcept : a_(other.a_)
  {
    other.a_ = nullptr;
  }
  PsaAligner&
  operator=(PsaAligner&& other) noexcept
  {
    if (this != &other) {
      tsta_psa_aligner_destroy(a_);
      a_ = other.a_;
      other.a_ = nullptr;
    }
    return *this;
  }

  PsaResult
  align(const char* seq1, int len1, const char* seq2, int len2)
  {
    PsaResult r;
    if (tsta_psa_aligner_align(a_, seq1, len1, seq2, len2, r.raw()) != 0)
      throw error("PSA alignment failed");
    return r;
  }

  PsaResult
  align(std::string_view seq1, std::string_view seq2)
  {
    return align(seq1.data(), (int)seq1.size(), seq2.data(), (int)seq2.size());
  }

private:
  tsta_psa_aligner* a_;
};

/* ── MSA aligner ─────────────────────────────────────────────────────── */

class MsaAligner {
public:
  explicit MsaAligner(const Config& config = Config{})
  {
    a_ = tsta_msa_aligner_create(config.raw());
    if (!a_)
      throw error("tsta_msa_aligner_create failed");
  }

  ~MsaAligner()
  {
    if (a_)
      tsta_msa_aligner_destroy(a_);
  }

  MsaAligner(const MsaAligner&) = delete;
  MsaAligner& operator=(const MsaAligner&) = delete;

  MsaAligner(MsaAligner&& other) noexcept : a_(other.a_)
  {
    other.a_ = nullptr;
  }
  MsaAligner&
  operator=(MsaAligner&& other) noexcept
  {
    if (this != &other) {
      tsta_msa_aligner_destroy(a_);
      a_ = other.a_;
      other.a_ = nullptr;
    }
    return *this;
  }

  /* one-shot */
  MsaResult
  align(const char* const* seqs, const int* lengths, size_t count)
  {
    MsaResult r;
    if (tsta_msa_aligner_align(a_, seqs, lengths, count, r.raw()) != 0)
      throw error("MSA alignment failed");
    return r;
  }

  MsaResult
  align(const std::vector<std::string_view>& seqs)
  {
    std::vector<const char*> ptrs(seqs.size());
    std::vector<int> lens(seqs.size());
    for (size_t i = 0; i < seqs.size(); i++) {
      ptrs[i] = seqs[i].data();
      lens[i] = (int)seqs[i].size();
    }
    return align(ptrs.data(), lens.data(), seqs.size());
  }

  /* incremental */
  MsaAligner&
  begin(const char* seq, int len)
  {
    if (tsta_msa_aligner_begin(a_, seq, len) != 0)
      throw error("tsta_msa_aligner_begin failed");
    return *this;
  }
  MsaAligner&
  begin(std::string_view seq)
  {
    return begin(seq.data(), (int)seq.size());
  }

  MsaAligner&
  add(const char* seq, int len)
  {
    if (tsta_msa_aligner_add(a_, seq, len) != 0)
      throw error("tsta_msa_aligner_add failed");
    return *this;
  }
  MsaAligner&
  add(std::string_view seq)
  {
    return add(seq.data(), (int)seq.size());
  }

  MsaResult
  get_result()
  {
    MsaResult r;
    if (tsta_msa_aligner_get_result(a_, r.raw()) != 0)
      throw error("tsta_msa_aligner_get_result failed");
    return r;
  }

private:
  tsta_msa_aligner* a_;
};

/* ── Free functions: one-shot alignment ──────────────────────────────── */

inline PsaResult
psa_align(const char* seq1,
          int len1,
          const char* seq2,
          int len2,
          const Config& config = Config{})
{
  PsaResult r;
  if (tsta_psa_align(seq1, len1, seq2, len2, config.raw(), r.raw()) != 0)
    throw error("PSA alignment failed");
  return r;
}

inline PsaResult
psa_align(std::string_view seq1,
          std::string_view seq2,
          const Config& config = Config{})
{
  return psa_align(seq1.data(), (int)seq1.size(), seq2.data(),
                   (int)seq2.size(), config);
}

inline MsaResult
msa_align(const char* const* seqs,
          const int* lengths,
          size_t count,
          const Config& config = Config{})
{
  MsaResult r;
  if (tsta_msa_align(seqs, lengths, count, config.raw(), r.raw()) != 0)
    throw error("MSA alignment failed");
  return r;
}

inline MsaResult
msa_align(const std::vector<std::string_view>& seqs,
          const Config& config = Config{})
{
  std::vector<const char*> ptrs(seqs.size());
  std::vector<int> lens(seqs.size());
  for (size_t i = 0; i < seqs.size(); i++) {
    ptrs[i] = seqs[i].data();
    lens[i] = (int)seqs[i].size();
  }
  return msa_align(ptrs.data(), lens.data(), seqs.size(), config);
}

inline MsaResult
msa_align(const std::vector<std::string>& seqs,
          const Config& config = Config{})
{
  std::vector<std::string_view> views;
  views.reserve(seqs.size());
  for (auto& s : seqs)
    views.emplace_back(s);
  return msa_align(views, config);
}

} // namespace tsta

#endif /* TSTA_HPP */
