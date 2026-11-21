#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Profiled-gene Poisson (conditional multinomial within gene) evaluation for codon_dwell().
//
// Model: E[y_gi] = mu_g * exp(eta_gi), eta_gi = sum_k beta[nf * k + c_ik - 1] (c_ik = 0 is the reference codon of
// offset k and adds nothing). mu_g is profiled out exactly: mu_g = n_g / sum_i exp(eta_gi), so the expected count is
// m_gi = n_g * p_gi with p_gi the within-gene softmax of eta. The design is implicit one-hot, stored as codes.
//
// codes  : raw matrix (all design positions x offsets), values 0..nf (0 = reference codon)
// rows   : 0-based design rows used in this fit, grouped by gene
// y      : counts, one per entry of rows
// starts : length G + 1, gene g owns rows[starts[g] .. starts[g + 1] - 1]
// beta   : length nf * offsets
// Returns the gradient (score summed over genes), the Hessian part sum_i m_gi x_i x_i' (lower triangle filled,
// symmetric), the per-gene expected feature sums scaled by sqrt(n_g) (P x G; the profiled Hessian is
// hessian_raw - xp %*% t(xp)), optionally the per-gene scores (P x G, for the gene-clustered sandwich), and the
// deviance and null (uniform within gene) deviance.
// threads: OpenMP threads (genes in parallel for the gradient/scores, chunks of positions for the Hessian; each thread
// keeps its own accumulators, summed at the end, so results do not depend on the thread count beyond rounding).
// [[Rcpp::export]]
List codon_dwell_eval_cpp(RawMatrix codes, IntegerVector rows, NumericVector y,
                          IntegerVector starts, NumericVector beta, int nf,
                          bool hessian = true, bool scores = false, int threads = 1) {
  const int N = rows.size(), K = codes.ncol(), G = starts.size() - 1;
  const int P = nf * K;
  const size_t NR = codes.nrow();
  if (y.size() != N || beta.size() != P || G < 1 || starts[G] != N)
    stop("codon_dwell_eval_cpp: incompatible inputs");
  const unsigned char* cd = RAW(codes);
  const int* rw = INTEGER(rows);
  const double* yy = REAL(y);
  const double* bb = REAL(beta);
  for (int i = 0; i < N; i++) if (rw[i] < 0 || (size_t)rw[i] >= NR)
    stop("codon_dwell_eval_cpp: row index out of bounds");

  for (int g = 0; g < G; g++) if (starts[g + 1] < starts[g]) stop("codon_dwell_eval_cpp: starts must be increasing");
  // codes <= nf is checked once per design on the R side (codon_dwell_check_design)
#ifdef _OPENMP
  const int nt = std::max(1, std::min(threads, omp_get_num_procs()));
#else
  const int nt = 1;
#endif
  std::vector<unsigned char> cc((size_t)N * K);
  std::vector<double> mu(N);
  NumericMatrix xp(hessian ? P : 0, hessian ? G : 0);
  NumericMatrix score(scores ? P : 0, scores ? G : 0);
  double* xp_ptr = hessian ? REAL(xp) : nullptr;
  double* score_ptr = scores ? REAL(score) : nullptr;
  std::vector<std::vector<double>> grt(nt, std::vector<double>(P, 0.0));
  std::vector<std::vector<double>> sgt(nt, std::vector<double>(P, 0.0));
  double dev = 0, dev0 = 0;
#ifdef _OPENMP
  #pragma omp parallel for num_threads(nt) schedule(dynamic, 8) reduction(+:dev,dev0)
#endif
  for (int g = 0; g < G; g++) {
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    const int a = starts[g], b = starts[g + 1], m = b - a;
    if (m <= 0) continue;
    double n = 0, mx = -INFINITY;
    std::vector<double> ex(m);
    for (int i = a; i < b; i++) {
      double v = 0; n += yy[i];
      const size_t r = rw[i];
      for (int k = 0; k < K; k++) {
        unsigned char c = cd[r + NR * k];
        cc[(size_t)i * K + k] = c;
        if (c) v += bb[nf * k + c - 1];
      }
      ex[i - a] = v; mx = std::max(mx, v);
    }
    double z = 0;
    for (int j = 0; j < m; j++) { ex[j] = std::exp(ex[j] - mx); z += ex[j]; }
    double* xg = hessian ? xp_ptr + (size_t)g * P : nullptr;
    double* sg = scores ? score_ptr + (size_t)g * P : sgt[tid].data();
    if (!scores) std::fill(sgt[tid].begin(), sgt[tid].end(), 0.0);
    for (int i = a; i < b; i++) {
      const double p = ex[i - a] / z, yi = yy[i];
      mu[i] = n * p;
      if (yi > 0) {
        dev += 2 * yi * (std::log(yi / n) - std::log(p));
        dev0 += 2 * yi * std::log(yi * m / n);
      }
      const double s = yi - mu[i];
      for (int k = 0; k < K; k++) {
        unsigned char c = cc[(size_t)i * K + k];
        if (c) { int id = nf * k + c - 1; if (hessian) xg[id] += p; sg[id] += s; }
      }
    }
    double* gr = grt[tid].data();
    for (int j = 0; j < P; j++) gr[j] += sg[j];
    if (hessian) { const double sn = std::sqrt(n); for (int j = 0; j < P; j++) xg[j] *= sn; }
  }
  NumericVector grad(P);
  for (int t = 0; t < nt; t++) for (int j = 0; j < P; j++) grad[j] += grt[t][j];
  if (!hessian) {
    return List::create(_["gradient"] = grad, _["score"] = score,
                        _["deviance"] = dev, _["null_deviance"] = dev0);
  }
  // Hessian sum_i mu_i x_i x_i': accumulated per chunk of positions and per offset pair, into contiguous
  // (nf + 1) x (nf + 1) blocks (slot 0 = reference, discarded), so updates stay in cache.
  const int B1 = nf + 1, BB = B1 * B1, NP = K * (K + 1) / 2, CH = 16384;
  const int nch = (N + CH - 1) / CH;
  std::vector<std::vector<double>> ht(nt);
#ifdef _OPENMP
  #pragma omp parallel num_threads(nt)
#endif
  {
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    std::vector<double>& h = ht[tid];
    h.assign((size_t)NP * BB, 0.0);
    std::vector<unsigned char> ct((size_t)CH * K);
#ifdef _OPENMP
    #pragma omp for schedule(dynamic, 1)
#endif
    for (int ch = 0; ch < nch; ch++) {
      const int a = ch * CH, b = std::min(N, a + CH), m = b - a;
      for (int i = 0; i < m; i++) for (int k = 0; k < K; k++) ct[(size_t)k * CH + i] = cc[(size_t)(a + i) * K + k];
      const double* mc = mu.data() + a;
      int pr = 0;
      for (int k = 0; k < K; k++) for (int l = 0; l <= k; l++, pr++) {
        double* blk = h.data() + (size_t)pr * BB;
        const unsigned char* ck = ct.data() + (size_t)k * CH;
        const unsigned char* cl = ct.data() + (size_t)l * CH;
        for (int i = 0; i < m; i++) blk[cl[i] * B1 + ck[i]] += mc[i];
      }
    }
  }
  std::vector<double> h((size_t)NP * BB, 0.0);
  for (int t = 0; t < nt; t++) for (size_t j = 0; j < h.size(); j++) h[j] += ht[t][j];
  NumericMatrix H(P, P);
  int pr = 0;
  for (int k = 0; k < K; k++) for (int l = 0; l <= k; l++, pr++) {
    const double* blk = h.data() + (size_t)pr * BB;
    for (int c2 = 1; c2 <= nf; c2++) for (int c1 = 1; c1 <= nf; c1++) {
      const double v = blk[c2 * B1 + c1];
      if (v == 0) continue;
      const int r = nf * k + c1 - 1, c = nf * l + c2 - 1;
      H(r, c) += v;
      if (k != l) H(c, r) += v;
    }
  }
  return List::create(_["gradient"] = grad, _["hessian_raw"] = H, _["xp"] = xp, _["score"] = score,
                      _["deviance"] = dev, _["null_deviance"] = dev0);
}

// Per offset and code (0..nf): number of positions and summed counts over the selected rows (feature support for
// codon_dwell(): codons without positions or reads at an offset are not estimated).
// [[Rcpp::export]]
List codon_dwell_feature_counts_cpp(RawMatrix codes, IntegerVector rows, NumericVector y, int nf) {
  const int N = rows.size(), K = codes.ncol();
  const size_t NR = codes.nrow();
  if (y.size() != N) stop("codon_dwell_feature_counts_cpp: incompatible inputs");
  const unsigned char* cd = RAW(codes);
  NumericMatrix n(nf + 1, K), ys(nf + 1, K);
  for (int k = 0; k < K; k++) {
    const unsigned char* ck = cd + NR * k;
    double* nk = &n(0, k); double* yk = &ys(0, k);
    for (int i = 0; i < N; i++) {
      const size_t r = rows[i];
      if (r >= NR) stop("codon_dwell_feature_counts_cpp: row index out of bounds");
      const unsigned char c = ck[r];
      if (c > nf) stop("codon_dwell_feature_counts_cpp: code larger than number of features");
      nk[c] += 1; yk[c] += y[i];
    }
  }
  return List::create(_["n"] = n, _["y"] = ys);
}
