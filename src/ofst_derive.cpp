// Kernels of ofst_derive(): coverage runs of many chromosome blocks of a sorted point OFST table
// in parallel (OpenMP), and range sums over coverage runs (count tables).
//
// ofst_rle_blocks() is the multi-block, multi-threaded form of ofst_chr_rle() (ofst_runs.cpp):
// every block is split into tasks of about task_rows rows at position borders, each task builds
// its runs independently, and the pieces are joined (equal values at a border are merged), so
// the result equals ofst_chr_rle() on each block. No R API calls inside parallel regions.
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// Runs of one coverage track over positions (pos_begin, pos_end], built from positions in
// nondecreasing order; adjacent equal values are merged (normalised as S4Vectors::Rle()).
struct Runs {
  std::vector<double> v;
  std::vector<int> l;
  int end = 0;            // last position emitted
  int cur_p = -1;         // pending position and its summed value
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
  inline void finish(int pos_end) { flush(); push(0.0, pos_end - end); }
};

struct Task { int block; R_xlen_t a, e; int pos_begin, pos_end; };

inline int n_threads(int threads) { return threads < 1 ? 1 : threads; }

} // namespace

// Coverage runs of several chromosome blocks of a point OFST table, for every (read length,
// strand) group and all read lengths merged.
// start: positions; strand: factor codes of the strand column, mapped by strand_map (code ->
// 1 "+", 2 "-", 3 "*", anything else invalid); size: NULL (all rows size 1) or integer read
// lengths in 1..1e6; score: integer or double, no NA. Rows block_from[b] + 1 .. block_to[b]
// (0-based from, 1-based to) are one chromosome of length chr_len[b] with nondecreasing start.
// Returns NULL if any of this does not hold (the caller then falls back), else
// list(sizes = present read lengths, blocks = per block a list of 2 * (length(sizes) + 1)
// groups in the order of ofst_chr_rle() (size 1 forward, size 1 reverse, ..., all forward, all
// reverse), each list(values, lengths), or NULL for a group without any row in the block).
// Values are integer when score is integer and every summed value fits, else double.
// [[Rcpp::export]]
SEXP ofst_rle_blocks(Rcpp::IntegerVector start, Rcpp::IntegerVector strand, Rcpp::IntegerVector strand_map,
                     SEXP size, SEXP score, Rcpp::NumericVector block_from, Rcpp::NumericVector block_to,
                     Rcpp::IntegerVector chr_len, int threads = 1, double task_rows = 2e6) {
  const R_xlen_t n = start.size();
  const int nb = block_from.size();
  const bool hasSize = !Rf_isNull(size);
  if (strand.size() != n || XLENGTH(score) != n || (hasSize && (TYPEOF(size) != INTSXP || XLENGTH(size) != n)))
    Rcpp::stop("inputs must have equal length (size integer or NULL)");
  if (block_to.size() != nb || chr_len.size() != nb) Rcpp::stop("block_from, block_to and chr_len must have equal length");
  const bool scoreInt = TYPEOF(score) == INTSXP;
  if (!scoreInt && TYPEOF(score) != REALSXP) return R_NilValue;
  const int th = n_threads(threads);
  const int* p = start.begin();
  const int* st = strand.begin();
  const int* sz = hasSize ? INTEGER(size) : nullptr;
  const int* si = scoreInt ? INTEGER(score) : nullptr;
  const double* sd = scoreInt ? nullptr : REAL(score);
  const int nmap = strand_map.size();
  std::vector<int> smap(nmap + 1, 0);
  for (int k = 0; k < nmap; ++k) {
    const int m = strand_map[k];
    smap[k + 1] = (m == 1 || m == 2 || m == 3) ? m : 0;
  }
  for (int b = 0; b < nb; ++b)
    if (!(block_from[b] >= 0 && block_to[b] <= n && block_from[b] < block_to[b]) || chr_len[b] == NA_INTEGER || chr_len[b] < 1)
      return R_NilValue;

  // read lengths present (code k = rank of the size), validity of sizes, strands and scores
  const int max_size = 1000000;
  std::vector<int> gmap;          // size -> 0-based size group, -1 absent
  int S = 1;
  bool bad = false;
  if (hasSize) {
    int mx = 0;
    #pragma omp parallel for num_threads(th) reduction(max:mx) schedule(static)
    for (R_xlen_t i = 0; i < n; ++i) { const int s = sz[i]; if (s > mx) mx = s; }
    if (mx < 1 || mx > max_size) return R_NilValue;
    std::vector<unsigned char> present(mx + 1, 0);
    #pragma omp parallel num_threads(th)
    {
      std::vector<unsigned char> loc(mx + 1, 0);
      bool lbad = false;
      #pragma omp for schedule(static) nowait
      for (R_xlen_t i = 0; i < n; ++i) { const int s = sz[i]; if (s < 1) lbad = true; else loc[s] = 1; }
      #pragma omp critical
      {
        if (lbad) bad = true;
        for (int s = 1; s <= mx; ++s) if (loc[s]) present[s] = 1;
      }
    }
    if (bad) return R_NilValue;
    gmap.assign(mx + 1, -1);
    S = 0;
    for (int s = 1; s <= mx; ++s) if (present[s]) gmap[s] = S++;
  }
  const int G = 2 * (S + 1);

  // tasks: blocks split at position borders, about task_rows rows each
  std::vector<Task> tasks;
  std::vector<int> first_task(nb + 1, 0);
  const R_xlen_t step = std::max<R_xlen_t>(10000, static_cast<R_xlen_t>(task_rows));
  for (int b = 0; b < nb; ++b) {
    first_task[b] = static_cast<int>(tasks.size());
    const R_xlen_t from = static_cast<R_xlen_t>(block_from[b]), to = static_cast<R_xlen_t>(block_to[b]);
    R_xlen_t a = from;
    int pos_begin = 0;
    while (a < to) {
      R_xlen_t e = (to - a > step + step / 2) ? a + step : to;
      while (e < to && p[e] == p[e - 1]) ++e;         // all rows of one position in one task
      if (e < to && p[e] < p[e - 1]) return R_NilValue;  // unsorted at a task border
      const int pos_end = e < to ? p[e - 1] : chr_len[b];
      if (e < to && (p[e - 1] < 1 || p[e - 1] > chr_len[b])) return R_NilValue;
      tasks.push_back({b, a, e, pos_begin, pos_end});
      pos_begin = pos_end;
      a = e;
    }
  }
  first_task[nb] = static_cast<int>(tasks.size());
  const int nt = static_cast<int>(tasks.size());

  // per task: validate and count positions per group (reserve once), then build the runs
  std::vector<std::vector<Runs> > out(nt);
  std::vector<std::vector<R_xlen_t> > counts(nt);
  std::vector<std::vector<char> > big(nt);   // per task and group: a value above the integer range
  #pragma omp parallel for schedule(dynamic, 1) num_threads(th)
  for (int t = 0; t < nt; ++t) {
    if (bad) continue;
    const Task& T = tasks[t];
    const int len = chr_len[T.block];
    std::vector<R_xlen_t> cnt(G, 0);
    std::vector<int> last(G, -1);
    bool tbad = false;
    for (R_xlen_t i = T.a; i < T.e; ++i) {
      const int pi = p[i];
      const int sc = (st[i] >= 1 && st[i] <= nmap) ? smap[st[i]] : 0;
      const int k = hasSize ? ((sz[i] >= 1 && sz[i] < static_cast<int>(gmap.size())) ? gmap[sz[i]] : -1) : 0;
      const bool na_score = scoreInt ? si[i] == NA_INTEGER : std::isnan(sd[i]);
      if (pi == NA_INTEGER || pi <= T.pos_begin || pi > T.pos_end || pi > len || sc == 0 || k < 0 || na_score ||
          (i > T.a && pi < p[i - 1])) { tbad = true; break; }
      if (sc != 2) {
        if (last[2 * k] != pi) { ++cnt[2 * k]; last[2 * k] = pi; }
        if (last[2 * S] != pi) { ++cnt[2 * S]; last[2 * S] = pi; }
      }
      if (sc != 1) {
        if (last[2 * k + 1] != pi) { ++cnt[2 * k + 1]; last[2 * k + 1] = pi; }
        if (last[2 * S + 1] != pi) { ++cnt[2 * S + 1]; last[2 * S + 1] = pi; }
      }
    }
    if (tbad) {
      #pragma omp critical
      bad = true;
      continue;
    }
    std::vector<Runs> bld(G);
    for (int g = 0; g < G; ++g) {
      bld[g].end = T.pos_begin;
      bld[g].v.reserve(2 * cnt[g] + 1); bld[g].l.reserve(2 * cnt[g] + 1);
    }
    for (R_xlen_t i = T.a; i < T.e; ++i) {
      const double val = scoreInt ? static_cast<double>(si[i]) : sd[i];
      const int sc = smap[st[i]];
      const int k = hasSize ? gmap[sz[i]] : 0;
      if (sc != 2) { bld[2 * k].add(p[i], val); bld[2 * S].add(p[i], val); }
      if (sc != 1) { bld[2 * k + 1].add(p[i], val); bld[2 * S + 1].add(p[i], val); }
    }
    std::vector<char> bg(G, 0);
    for (int g = 0; g < G; ++g) {
      bld[g].finish(T.pos_end);
      if (scoreInt) for (size_t j = 0; j < bld[g].v.size(); ++j) if (bld[g].v[j] > 2147483647.0) { bg[g] = 1; break; }
    }
    out[t].swap(bld);
    big[t].swap(bg);
    counts[t].swap(cnt);
  }
  if (bad) return R_NilValue;

  // join the task pieces per (block, group): sizes and types first (serial, allocates R
  // vectors), then the copies in parallel
  struct Job { int b, g; R_xlen_t m; bool asInt; int* lp; int* vi; double* vd; };
  std::vector<Job> jobs;
  Rcpp::List blocks(nb);
  for (int b = 0; b < nb; ++b) {
    Rcpp::List groups(G);
    for (int g = 0; g < G; ++g) {
      R_xlen_t positions = 0, m = 0;
      bool asInt = scoreInt;
      double prev_last = 0; bool has_prev = false;
      for (int t = first_task[b]; t < first_task[b + 1]; ++t) {
        positions += counts[t][g];
        const Runs& r = out[t][g];
        if (r.v.empty()) continue;
        m += static_cast<R_xlen_t>(r.v.size());
        if (has_prev && prev_last == r.v.front()) --m;
        prev_last = r.v.back(); has_prev = true;
        if (big[t][g]) asInt = false;
      }
      if (positions == 0) { groups[g] = R_NilValue; continue; }
      Rcpp::IntegerVector len(Rcpp::no_init(m));   // filled below, in parallel
      Job J = {b, g, m, asInt, len.begin(), nullptr, nullptr};
      if (asInt) {
        Rcpp::IntegerVector val(Rcpp::no_init(m)); J.vi = val.begin();
        groups[g] = Rcpp::List::create(Rcpp::Named("values") = val, Rcpp::Named("lengths") = len);
      } else {
        Rcpp::NumericVector val(Rcpp::no_init(m)); J.vd = val.begin();
        groups[g] = Rcpp::List::create(Rcpp::Named("values") = val, Rcpp::Named("lengths") = len);
      }
      jobs.push_back(J);
    }
    blocks[b] = groups;
  }
  const int nj = static_cast<int>(jobs.size());
  #pragma omp parallel for schedule(dynamic, 1) num_threads(th)
  for (int j = 0; j < nj; ++j) {
    const Job& J = jobs[j];
    R_xlen_t w = 0;
    for (int t = first_task[J.b]; t < first_task[J.b + 1]; ++t) {
      Runs& r = out[t][J.g];
      for (size_t q = 0; q < r.v.size(); ++q) {
        const double x = r.v[q];
        const bool merge = w > 0 && (J.asInt ? static_cast<double>(J.vi[w - 1]) == x : J.vd[w - 1] == x);
        if (merge) { J.lp[w - 1] += r.l[q]; continue; }
        if (J.asInt) J.vi[w] = static_cast<int>(x); else J.vd[w] = x;
        J.lp[w] = r.l[q];
        ++w;
      }
      std::vector<double>().swap(r.v); std::vector<int>().swap(r.l);
    }
  }
  Rcpp::IntegerVector sizes(S);
  if (hasSize) { for (int s = 1; s < static_cast<int>(gmap.size()); ++s) if (gmap[s] >= 0) sizes[gmap[s]] = s; }
  else sizes = Rcpp::IntegerVector::create(1);
  return Rcpp::List::create(Rcpp::Named("sizes") = sizes, Rcpp::Named("blocks") = blocks);
}

// Sums of coverage given as runs (values, lengths; one chromosome strand) over the 1-based closed
// ranges [starts, ends]. Exact for integer coverage (sums below 2^53). One sweep over the runs
// with the sorted range borders (no arrays of the run count). Ranges must lie inside the
// coverage (1 <= start, end <= total length, end >= start - 1); returns NULL otherwise.
// [[Rcpp::export]]
SEXP rle_range_sums(SEXP values, Rcpp::IntegerVector lengths, Rcpp::IntegerVector starts, Rcpp::IntegerVector ends) {
  const R_xlen_t m = lengths.size(), q = starts.size();
  if (XLENGTH(values) != m || ends.size() != q) Rcpp::stop("values/lengths and starts/ends must have equal length");
  const bool isInt = TYPEOF(values) == INTSXP;
  if (!isInt && TYPEOF(values) != REALSXP) Rcpp::stop("values must be integer or double");
  const int* vi = isInt ? INTEGER(values) : nullptr;
  const double* vd = isInt ? nullptr : REAL(values);
  const int* lp = lengths.begin();
  double total = 0;
  for (R_xlen_t j = 0; j < m; ++j) total += lp[j];
  // borders: coverage sum over 1..x is needed at x = start - 1 and x = end
  std::vector<std::pair<double, R_xlen_t> > pts(2 * q);
  for (R_xlen_t i = 0; i < q; ++i) {
    const int s = starts[i], e = ends[i];
    if (s == NA_INTEGER || e == NA_INTEGER || s < 1 || e > total || e < s - 1) return R_NilValue;
    pts[2 * i] = std::make_pair(static_cast<double>(s) - 1, 2 * i);
    pts[2 * i + 1] = std::make_pair(static_cast<double>(e), 2 * i + 1);
  }
  std::sort(pts.begin(), pts.end());
  std::vector<double> upto(2 * q);
  double cend = 0, csum = 0;   // bases and coverage sum before run j
  R_xlen_t j = 0;
  for (const auto& pt : pts) {
    const double x = pt.first;
    while (j < m && cend + lp[j] < x) {
      csum += (isInt ? static_cast<double>(vi[j]) : vd[j]) * lp[j];
      cend += lp[j]; ++j;
    }
    upto[pt.second] = x <= cend ? csum : csum + (isInt ? static_cast<double>(vi[j]) : vd[j]) * (x - cend);
  }
  Rcpp::NumericVector res(q);
  for (R_xlen_t i = 0; i < q; ++i) res[i] = upto[2 * i + 1] - upto[2 * i];
  return res;
}
