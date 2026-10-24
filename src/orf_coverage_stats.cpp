// Per-group coverage statistics straight from covRle runs (used by detect_ribo_orfs via coveragePerORFStatistics).
// Groups are ORFs/windows given as exons in 5'->3' order; for each group: sum, median, mean over positions, frame sums
// (frame = (position - 1) %% 3) and ORFScore as orfScore(stop3 = FALSE) defines it (only for width > 2, else 0;
// NaN -> 0). Same values as coveragePerTiling(as.data.table, withFrames) + data.table aggregation + orfScore().
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix orfStatsCovRleCPP(List run_values, List run_lengths, IntegerVector rle_id, IntegerVector ex_start,
                                IntegerVector ex_end, LogicalVector ex_minus, IntegerVector group_ends) {
  const int nr = run_values.size();
  std::vector<std::vector<double>> rv(nr);
  std::vector<std::vector<double>> re(nr);  // run ends, 1-based inclusive (double: genome positions can exceed int)
  for (int r = 0; r < nr; ++r) {
    NumericVector v = run_values[r]; NumericVector l = run_lengths[r];
    rv[r].assign(v.begin(), v.end()); re[r].resize(l.size());
    double c = 0; for (R_xlen_t i = 0; i < l.size(); ++i) { c += l[i]; re[r][i] = c; }
  }
  const R_xlen_t ng = group_ends.size();
  NumericMatrix out(ng, 7);
  std::vector<double> v;
  R_xlen_t ex = 0;
  for (R_xlen_t g = 0; g < ng; ++g) {
    v.clear();
    for (; ex < group_ends[g]; ++ex) {
      const int s = ex_start[ex], e = ex_end[ex], w = e - s + 1;
      if (w <= 0) continue;
      const size_t base = v.size(); v.resize(base + w, 0.0);
      const int r = rle_id[ex];
      if (r >= 0) {
        const std::vector<double>& E = re[r]; const std::vector<double>& V = rv[r];
        size_t k = std::lower_bound(E.begin(), E.end(), (double)s) - E.begin();
        for (double p = s; p <= e && k < E.size(); ) {
          const double runEnd = std::min(E[k], (double)e);
          for (double q = p; q <= runEnd; ++q) {
            const int pos = (int)(q - s);
            v[base + (ex_minus[ex] ? (w - 1 - pos) : pos)] = V[k];
          }
          p = runEnd + 1; ++k;
        }
      }
    }
    // minus-strand exons are already in 5'->3' order per group, each read right-to-left
    const int w = (int)v.size();
    double sum = 0, med = 0, mean = 0, f[3] = {0, 0, 0}, score = 0;
    if (w > 0) {
      for (int j = 0; j < w; ++j) { sum += v[j]; f[j % 3] += v[j]; }
      mean = sum / w;
      std::vector<double> t(v); const int h = w / 2;
      std::nth_element(t.begin(), t.begin() + h, t.end());
      if (w % 2) med = t[h]; else med = (*std::max_element(t.begin(), t.begin() + h) + t[h]) / 2.0;
      if (w > 2) {
        const double Ft = (f[0] + f[1] + f[2]) / 3.0;
        double sc = std::log2((f[0] - Ft) * (f[0] - Ft) / Ft + (f[1] - Ft) * (f[1] - Ft) / Ft + (f[2] - Ft) * (f[2] - Ft) / Ft + 1);
        if (f[0] < f[1] || f[0] < f[2]) sc = -sc;
        score = std::isnan(sc) ? 0 : sc;
      } else { f[0] = f[1] = f[2] = 0; }
    }
    out(g, 0) = sum; out(g, 1) = med; out(g, 2) = mean; out(g, 3) = f[0]; out(g, 4) = f[1]; out(g, 5) = f[2]; out(g, 6) = score;
  }
  colnames(out) = CharacterVector::create("sum", "median", "mean", "F1", "F2", "F3", "ORFScores");
  return out;
}
