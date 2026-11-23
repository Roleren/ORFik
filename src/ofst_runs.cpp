#include <Rcpp.h>
#include <vector>
#include <cmath>

// Runs of an integer vector (e.g. OFST seqnames factor codes): one pass, no copy of the input.
// Returns the run values and their 1-based end rows.
// [[Rcpp::export]]
Rcpp::List ofst_code_runs(Rcpp::IntegerVector x) {
  const R_xlen_t n = x.size();
  std::vector<int> values;
  std::vector<double> ends;
  if (n > 0) {
    const int* p = x.begin();
    int current = p[0];
    for (R_xlen_t i = 1; i < n; ++i) {
      if (p[i] != current) {
        values.push_back(current);
        ends.push_back(static_cast<double>(i));
        current = p[i];
      }
    }
    values.push_back(current);
    ends.push_back(static_cast<double>(n));
  }
  return Rcpp::List::create(Rcpp::Named("value") = Rcpp::wrap(values),
                            Rcpp::Named("end") = Rcpp::wrap(ends));
}

// Point coverage sums for P-site style OFST rows in one pass.
// Rows must have each chromosome as one contiguous block with nondecreasing start (OFST
// merge/collapse order); returns NULL otherwise. strand: 1 "+", 2 "-", 3 "*" ("*" counts on both
// strands); size_index: 1..n_size. Returns a list of 2 * (n_size + 1) groups, ordered
// (size 1 forward, size 1 reverse, ..., all sizes forward, all sizes reverse), each
// list(chr, start, score) of positions summed in row order. NA scores count 0.
// [[Rcpp::export]]
SEXP ofst_point_sums(Rcpp::IntegerVector chr, Rcpp::IntegerVector start, Rcpp::IntegerVector strand,
                     Rcpp::IntegerVector size_index, int n_size, SEXP score) {
  const R_xlen_t n = chr.size();
  if (start.size() != n || strand.size() != n || size_index.size() != n || XLENGTH(score) != n)
    Rcpp::stop("inputs must have equal length");
  const bool scoreInt = TYPEOF(score) == INTSXP;
  if (!scoreInt && TYPEOF(score) != REALSXP) Rcpp::stop("score must be integer or double");
  const int* si = scoreInt ? INTEGER(score) : nullptr;
  const double* sd = scoreInt ? nullptr : REAL(score);
  const int* c = chr.begin(); const int* p = start.begin(); const int* st = strand.begin(); const int* sz = size_index.begin();
  // contiguity and order check
  std::vector<char> seen;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (c[i] == NA_INTEGER || p[i] == NA_INTEGER || sz[i] < 1 || sz[i] > n_size || st[i] < 1 || st[i] > 3)
      return R_NilValue;
    if (i == 0 || c[i] != c[i - 1]) {
      if (c[i] < 0) return R_NilValue;
      if (static_cast<size_t>(c[i]) >= seen.size()) seen.resize(c[i] + 1, 0);
      if (seen[c[i]]) return R_NilValue;
      seen[c[i]] = 1;
    } else if (p[i] < p[i - 1]) return R_NilValue;
  }
  const int groups = 2 * (n_size + 1);
  std::vector<std::vector<int> > gc(groups), gp(groups);
  std::vector<std::vector<double> > gs(groups);
  auto add = [&](int g, int ch, int pos, double v) {
    if (!gp[g].empty() && gp[g].back() == pos && gc[g].back() == ch) { gs[g].back() += v; return; }
    gc[g].push_back(ch); gp[g].push_back(pos); gs[g].push_back(v);
  };
  for (R_xlen_t i = 0; i < n; ++i) {
    const double v = scoreInt ? (si[i] == NA_INTEGER ? 0.0 : static_cast<double>(si[i]))
                              : (std::isnan(sd[i]) ? 0.0 : sd[i]);
    const int k = sz[i] - 1;
    if (st[i] != 2) { add(2 * k, c[i], p[i], v); add(2 * n_size, c[i], p[i], v); }
    if (st[i] != 1) { add(2 * k + 1, c[i], p[i], v); add(2 * n_size + 1, c[i], p[i], v); }
  }
  Rcpp::List out(groups);
  for (int g = 0; g < groups; ++g)
    out[g] = Rcpp::List::create(Rcpp::Named("seqnames") = Rcpp::wrap(gc[g]), Rcpp::Named("start") = Rcpp::wrap(gp[g]),
                                Rcpp::Named("score") = Rcpp::wrap(gs[g]));
  return out;
}

namespace {
// Rle runs of one coverage track, built from positions in nondecreasing order. Adjacent equal
// values are merged, so the runs are normalised as S4Vectors::Rle() would store them.
struct RunBuilder {
  std::vector<double> v;
  std::vector<int> l;
  int end = 0;            // bases emitted so far
  int cur_p = -1;         // pending position (1-based) and its summed value
  double cur_v = 0;
  inline void push(double val, int len) {
    if (len <= 0) return;
    if (!v.empty() && v.back() == val) { l.back() += len; return; }
    v.push_back(val); l.push_back(len);
  }
  inline void flush() {
    if (cur_p < 0) return;
    push(0.0, cur_p - end - 1); push(cur_v, 1); end = cur_p; cur_p = -1;
  }
  inline void add(int p, double val) {
    if (p == cur_p) { cur_v += val; return; }
    flush(); cur_p = p; cur_v = val;
  }
  inline void finish(int chr_len) { flush(); push(0.0, chr_len - end); }
};
}

// Coverage runs of one chromosome block of a point OFST table, for every (read length, strand)
// and all read lengths merged, in one pass (no intermediate position sums).
// Rows from + 1 .. to (0-based from, 1-based to) must be one chromosome with nondecreasing
// start in 1..chr_len; returns NULL otherwise. strand: 1 "+", 2 "-", 3 "*" ("*" counts on both
// strands); size_index: 1..n_size. NA scores count 0.
// Returns 2 * (n_size + 1) groups in the order of ofst_point_sums() (size 1 forward, size 1
// reverse, ..., all sizes forward, all sizes reverse), each list(values, lengths). Values are
// integer when score is integer and every summed value fits, else double (as the R path).
// [[Rcpp::export]]
SEXP ofst_chr_rle(Rcpp::IntegerVector start, Rcpp::IntegerVector strand, Rcpp::IntegerVector size_index,
                  int n_size, SEXP score, int chr_len, double from, double to) {
  const R_xlen_t n = start.size();
  if (strand.size() != n || size_index.size() != n || XLENGTH(score) != n)
    Rcpp::stop("inputs must have equal length");
  if (from < 0 || to > n || from > to || n_size < 1 || chr_len < 0) Rcpp::stop("invalid row range or sizes");
  const bool scoreInt = TYPEOF(score) == INTSXP;
  if (!scoreInt && TYPEOF(score) != REALSXP) Rcpp::stop("score must be integer or double");
  const int* si = scoreInt ? INTEGER(score) : nullptr;
  const double* sd = scoreInt ? nullptr : REAL(score);
  const int* p = start.begin(); const int* st = strand.begin(); const int* sz = size_index.begin();
  const R_xlen_t a = static_cast<R_xlen_t>(from), b = static_cast<R_xlen_t>(to);
  const int G = 2 * (n_size + 1);
  // pass 1: validate and count distinct positions per group, so every vector is reserved once
  std::vector<R_xlen_t> cnt(G, 0);
  std::vector<int> last(G, -1);
  for (R_xlen_t i = a; i < b; ++i) {
    if (p[i] == NA_INTEGER || p[i] < 1 || p[i] > chr_len || sz[i] < 1 || sz[i] > n_size ||
        st[i] < 1 || st[i] > 3 || (i > a && p[i] < p[i - 1])) return R_NilValue;
    const int k = sz[i] - 1;
    if (st[i] != 2) {
      if (last[2 * k] != p[i]) { ++cnt[2 * k]; last[2 * k] = p[i]; }
      if (last[2 * n_size] != p[i]) { ++cnt[2 * n_size]; last[2 * n_size] = p[i]; }
    }
    if (st[i] != 1) {
      if (last[2 * k + 1] != p[i]) { ++cnt[2 * k + 1]; last[2 * k + 1] = p[i]; }
      if (last[2 * n_size + 1] != p[i]) { ++cnt[2 * n_size + 1]; last[2 * n_size + 1] = p[i]; }
    }
  }
  std::vector<RunBuilder> bld(G);
  for (int g = 0; g < G; ++g) { bld[g].v.reserve(2 * cnt[g] + 1); bld[g].l.reserve(2 * cnt[g] + 1); }
  for (R_xlen_t i = a; i < b; ++i) {
    const double val = scoreInt ? (si[i] == NA_INTEGER ? 0.0 : static_cast<double>(si[i]))
                                : (std::isnan(sd[i]) ? 0.0 : sd[i]);
    const int k = sz[i] - 1;
    if (st[i] != 2) { bld[2 * k].add(p[i], val); bld[2 * n_size].add(p[i], val); }
    if (st[i] != 1) { bld[2 * k + 1].add(p[i], val); bld[2 * n_size + 1].add(p[i], val); }
  }
  Rcpp::List out(G);
  for (int g = 0; g < G; ++g) {
    RunBuilder& r = bld[g];
    r.finish(chr_len);
    const size_t m = r.v.size();
    bool asInt = scoreInt;
    if (asInt) for (size_t j = 0; j < m; ++j) if (r.v[j] > 2147483647.0) { asInt = false; break; }
    Rcpp::IntegerVector len(r.l.begin(), r.l.end());
    if (asInt) {
      Rcpp::IntegerVector val(m);
      for (size_t j = 0; j < m; ++j) val[j] = static_cast<int>(r.v[j]);
      out[g] = Rcpp::List::create(Rcpp::Named("values") = val, Rcpp::Named("lengths") = len);
    } else {
      out[g] = Rcpp::List::create(Rcpp::Named("values") = Rcpp::NumericVector(r.v.begin(), r.v.end()),
                                  Rcpp::Named("lengths") = len);
    }
    std::vector<double>().swap(r.v); std::vector<int>().swap(r.l);
  }
  return out;
}

// Chromosome blocks of an OFST table for the sidecar index (.ofst.idx): runs of the seqnames
// codes with their 1-based end rows, the score sum per run (0 if score is NULL) and whether
// start is nondecreasing inside the run (FALSE if start is NULL).
// [[Rcpp::export]]
Rcpp::List ofst_block_index(Rcpp::IntegerVector codes, SEXP start, SEXP score) {
  const R_xlen_t n = codes.size();
  const int* c = codes.begin();
  const int* p = (TYPEOF(start) == INTSXP && XLENGTH(start) == n) ? INTEGER(start) : nullptr;
  const bool hasScore = (TYPEOF(score) == INTSXP || TYPEOF(score) == REALSXP) && XLENGTH(score) == n;
  const int* si = hasScore && TYPEOF(score) == INTSXP ? INTEGER(score) : nullptr;
  const double* sd = hasScore && TYPEOF(score) == REALSXP ? REAL(score) : nullptr;
  std::vector<int> values; std::vector<double> ends, sums; std::vector<int> sorted;
  double s = 0; bool srt = p != nullptr;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i > 0 && c[i] != c[i - 1]) {
      values.push_back(c[i - 1]); ends.push_back(static_cast<double>(i)); sums.push_back(s); sorted.push_back(srt);
      s = 0; srt = p != nullptr;
    } else if (i > 0 && srt && (p[i] == NA_INTEGER || p[i] < p[i - 1])) srt = false;
    if (p && p[i] == NA_INTEGER) srt = false;
    if (si) { if (si[i] != NA_INTEGER) s += si[i]; } else if (sd) { if (!std::isnan(sd[i])) s += sd[i]; }
  }
  if (n > 0) { values.push_back(c[n - 1]); ends.push_back(static_cast<double>(n)); sums.push_back(s); sorted.push_back(srt); }
  Rcpp::LogicalVector sorted_out(sorted.size());
  for (size_t j = 0; j < sorted.size(); ++j) sorted_out[j] = sorted[j];
  return Rcpp::List::create(Rcpp::Named("value") = Rcpp::wrap(values), Rcpp::Named("end") = Rcpp::wrap(ends),
                            Rcpp::Named("score") = Rcpp::wrap(sums), Rcpp::Named("start_sorted") = sorted_out);
}
