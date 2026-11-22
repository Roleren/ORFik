#include <Rcpp.h>
#include <vector>

// Linear-time kernels for incremental OFST merges and the OFST store (R/utils_ofst_store.R).
// Keys are k integer columns compared lexicographically (e.g. seqnames code, start, strand
// sort key, size). Inputs given as sorted must be strictly increasing (sorted and unique);
// the kernels check this while they scan and stop with an error otherwise.

namespace {

struct KeyCols {
  std::vector<const int*> p;
  R_xlen_t n;
  explicit KeyCols(const Rcpp::List& cols) : p(cols.size()), n(0) {
    for (int c = 0; c < cols.size(); ++c) {
      Rcpp::IntegerVector v = cols[c];
      if (c == 0) n = v.size(); else if (v.size() != n) Rcpp::stop("key columns differ in length");
      p[c] = v.begin();
    }
  }
};

// <0, 0, >0 for row i of a against row j of b
inline int cmp_rows(const KeyCols& a, R_xlen_t i, const KeyCols& b, R_xlen_t j) {
  const size_t k = a.p.size();
  for (size_t c = 0; c < k; ++c) {
    const int x = a.p[c][i], y = b.p[c][j];
    if (x != y) return x < y ? -1 : 1;
  }
  return 0;
}

void check_same_width(const KeyCols& a, const KeyCols& b) {
  if (a.p.size() != b.p.size() || a.p.empty()) Rcpp::stop("key column counts differ or are zero");
}

}  // namespace

// Sorted union of two keyed score tables with summed scores. a and b must each be strictly
// increasing. With drop_zero, keys whose summed score is 0 are left out (removal support);
// a negative sum is an error when forbid_negative is TRUE. Two passes: count, then fill, so
// the output is allocated once at its final size.
// [[Rcpp::export]]
Rcpp::List ofst_sorted_sum_merge(Rcpp::List a_keys, Rcpp::NumericVector a_score,
                                 Rcpp::List b_keys, Rcpp::NumericVector b_score,
                                 bool drop_zero, bool forbid_negative) {
  KeyCols a(a_keys), b(b_keys);
  check_same_width(a, b);
  const R_xlen_t n = a.n, m = b.n;
  if (a_score.size() != n || b_score.size() != m) Rcpp::stop("score length differs from key length");
  const double* sa = a_score.begin(); const double* sb = b_score.begin();
  for (R_xlen_t i = 1; i < n; ++i) if (cmp_rows(a, i - 1, a, i) >= 0) Rcpp::stop("first input is not sorted and unique (row %d)", (double)(i + 1));
  for (R_xlen_t j = 1; j < m; ++j) if (cmp_rows(b, j - 1, b, j) >= 0) Rcpp::stop("second input is not sorted and unique (row %d)", (double)(j + 1));
  // pass 1: count
  R_xlen_t out_n = 0;
  {
    R_xlen_t i = 0, j = 0;
    while (i < n || j < m) {
      double v;
      int c = (i >= n) ? 1 : (j >= m) ? -1 : cmp_rows(a, i, b, j);
      if (c < 0) { v = sa[i]; ++i; } else if (c > 0) { v = sb[j]; ++j; } else { v = sa[i] + sb[j]; ++i; ++j; }
      if (forbid_negative && v < 0) Rcpp::stop("a summed score is negative: removed input is not contained in the base");
      if (!(drop_zero && v == 0)) ++out_n;
    }
  }
  const size_t k = a.p.size();
  Rcpp::List out_keys(k);
  std::vector<int*> op(k);
  for (size_t c = 0; c < k; ++c) { Rcpp::IntegerVector v(out_n); op[c] = v.begin(); out_keys[c] = v; }
  Rcpp::NumericVector out_score(out_n);
  double* os = out_score.begin();
  R_xlen_t i = 0, j = 0, o = 0;
  while (i < n || j < m) {
    double v; const KeyCols* src; R_xlen_t r;
    int c = (i >= n) ? 1 : (j >= m) ? -1 : cmp_rows(a, i, b, j);
    if (c < 0) { v = sa[i]; src = &a; r = i; ++i; }
    else if (c > 0) { v = sb[j]; src = &b; r = j; ++j; }
    else { v = sa[i] + sb[j]; src = &a; r = i; ++i; ++j; }
    if (drop_zero && v == 0) continue;
    for (size_t cc = 0; cc < k; ++cc) op[cc][o] = src->p[cc][r];
    os[o++] = v;
  }
  out_keys.attr("names") = a_keys.attr("names");
  return Rcpp::List::create(Rcpp::Named("keys") = out_keys, Rcpp::Named("score") = out_score);
}

// total[row] += sign * score, in place on a copy (rows 1-based, may repeat). Returns the new
// total; stops if any touched total turns negative when forbid_negative is TRUE.
// [[Rcpp::export]]
Rcpp::NumericVector ofst_store_scatter_add(Rcpp::NumericVector total, Rcpp::IntegerVector row,
                                           Rcpp::NumericVector score, double sign, bool forbid_negative) {
  Rcpp::NumericVector out = Rcpp::clone(total);
  double* t = out.begin();
  const R_xlen_t n = out.size(), m = row.size();
  if (score.size() != m) Rcpp::stop("row and score lengths differ");
  const int* r = row.begin(); const double* s = score.begin();
  for (R_xlen_t j = 0; j < m; ++j) {
    const int x = r[j];
    if (x < 1 || x > n) Rcpp::stop("row id %d outside 1..%d", (double)x, (double)n);
    t[x - 1] += sign * s[j];
  }
  if (forbid_negative) for (R_xlen_t j = 0; j < m; ++j)
    if (t[r[j] - 1] < 0) Rcpp::stop("total score turned negative at row %d (ledger mismatch)", (double)r[j]);
  return out;
}

// TRUE if the key rows are strictly increasing (sorted and unique), one linear pass.
// Lets callers skip a group-by / sort of inputs that are already in merge order.
// [[Rcpp::export]]
bool ofst_keys_sorted_unique(Rcpp::List keys) {
  KeyCols k(keys);
  for (R_xlen_t i = 1; i < k.n; ++i) if (cmp_rows(k, i - 1, k, i) >= 0) return false;
  return true;
}

// Position (1-based) of each query key in a sorted, unique key table (contiguous columns, no
// permutation: cache-friendly), 0 when absent. Query must be strictly increasing.
// [[Rcpp::export]]
Rcpp::IntegerVector ofst_sorted_match(Rcpp::List table, Rcpp::List query) {
  KeyCols t(table), q(query);
  check_same_width(t, q);
  const R_xlen_t n = t.n, m = q.n;
  Rcpp::IntegerVector out(m);
  R_xlen_t i = 0;
  for (R_xlen_t j = 0; j < m; ++j) {
    if (j > 0 && cmp_rows(q, j - 1, q, j) >= 0) Rcpp::stop("query keys are not sorted and unique (row %d)", (double)(j + 1));
    R_xlen_t lo = i, step = 1, hi = i;
    while (hi < n && cmp_rows(t, hi, q, j) < 0) { lo = hi + 1; hi += step; step <<= 1; }
    if (hi > n) hi = n;
    while (lo < hi) {
      R_xlen_t mid = lo + (hi - lo) / 2;
      if (cmp_rows(t, mid, q, j) < 0) lo = mid + 1; else hi = mid;
    }
    i = lo;
    out[j] = (i < n && cmp_rows(t, i, q, j) == 0) ? (int)(i + 1) : 0;
  }
  return out;
}

// Merge a sorted key table (rows ord_a, contiguous sorted keys a) with sorted new keys b
// (rows rows_b) into the merged row order. Keys must not be shared.
// [[Rcpp::export]]
Rcpp::IntegerVector ofst_sorted_merge_order(Rcpp::List a, Rcpp::IntegerVector ord_a,
                                            Rcpp::List b, Rcpp::IntegerVector rows_b) {
  KeyCols ka(a), kb(b);
  check_same_width(ka, kb);
  const R_xlen_t n = ka.n, m = kb.n;
  if (ord_a.size() != n || rows_b.size() != m) Rcpp::stop("order lengths differ from key lengths");
  Rcpp::IntegerVector out(n + m);
  int* o = out.begin(); const int* oa = ord_a.begin(); const int* rb = rows_b.begin();
  R_xlen_t i = 0, j = 0, k = 0;
  while (i < n && j < m) {
    int c = cmp_rows(ka, i, kb, j);
    if (c == 0) Rcpp::stop("a new key is already in the table (row %d)", (double)oa[i]);
    o[k++] = c < 0 ? oa[i++] : rb[j++];
  }
  while (i < n) o[k++] = oa[i++];
  while (j < m) o[k++] = rb[j++];
  return out;
}
