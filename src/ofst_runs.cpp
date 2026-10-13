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
