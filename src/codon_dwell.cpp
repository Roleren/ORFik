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
// codon_dwell(): codons without positions or reads at an offset are not estimated). y2 (optional, same length as y)
// gives a second summed count (used for the start values: expected counts under uniform density within gene).
// [[Rcpp::export]]
List codon_dwell_feature_counts_cpp(RawMatrix codes, IntegerVector rows, NumericVector y, int nf,
                                    Nullable<NumericVector> y2_ = R_NilValue, int threads = 1) {
  const int N = rows.size(), K = codes.ncol();
  const size_t NR = codes.nrow();
  if (y.size() != N) stop("codon_dwell_feature_counts_cpp: incompatible inputs");
  NumericVector y2v;
  const bool has2 = y2_.isNotNull();
  if (has2) { y2v = NumericVector(y2_); if (y2v.size() != N) stop("codon_dwell_feature_counts_cpp: incompatible y2"); }
  const unsigned char* cd = RAW(codes);
  const int* rw = INTEGER(rows);
  const double* yy = REAL(y);
  const double* y2 = has2 ? REAL(y2v) : nullptr;
  for (int i = 0; i < N; i++) if (rw[i] < 0 || (size_t)rw[i] >= NR)
    stop("codon_dwell_feature_counts_cpp: row index out of bounds");
  for (size_t j = 0; j < NR * (size_t)K; j++) if (cd[j] > nf)
    stop("codon_dwell_feature_counts_cpp: code larger than number of features");
  NumericMatrix n(nf + 1, K), ys(nf + 1, K), ys2(has2 ? nf + 1 : 0, has2 ? K : 0);
  double* np = REAL(n); double* yp = REAL(ys); double* y2p = has2 ? REAL(ys2) : nullptr;
#ifdef _OPENMP
  const int nt = std::max(1, std::min(threads, omp_get_num_procs()));
  #pragma omp parallel for num_threads(nt) schedule(dynamic, 1)
#endif
  for (int k = 0; k < K; k++) {
    const unsigned char* ck = cd + NR * k;
    double* nk = np + (size_t)(nf + 1) * k; double* yk = yp + (size_t)(nf + 1) * k;
    double* y2k = has2 ? y2p + (size_t)(nf + 1) * k : nullptr;
    for (int i = 0; i < N; i++) {
      const unsigned char c = ck[rw[i]];
      nk[c] += 1; yk[c] += yy[i];
      if (has2) y2k[c] += y2[i];
    }
  }
  return List::create(_["n"] = n, _["y"] = ys, _["y2"] = ys2);
}

// Weighted tabulation: out[bin[i]] += w[i] for bins 1..n (bins outside 1..n are ignored).
// [[Rcpp::export]]
NumericVector codon_dwell_tabulate_cpp(IntegerVector bin, NumericVector w, int n) {
  const R_xlen_t N = bin.size();
  if (w.size() != N) stop("codon_dwell_tabulate_cpp: incompatible inputs");
  NumericVector out(n);
  const int* b = INTEGER(bin); const double* ww = REAL(w); double* o = REAL(out);
  for (R_xlen_t i = 0; i < N; i++) { const int k = b[i]; if (k >= 1 && k <= n) o[k - 1] += ww[i]; }
  return out;
}

// P-site -> CDS lookup for codon_dwell_map(): per read (chromosome index into the read seqlevels, strand 1 "+",
// 2 "-", 3 "*" = both, position = start, or end on "-"), find the CDS exon containing it in every layer of
// non-overlapping exons (keys chromosome-strand id * 2^32 + position, sorted within layer; layer_end = cumulative
// layer sizes). A per layer and strand hint makes position-sorted reads mostly O(1). Returns, per hit, the read
// index (1-based), CDS index (gene) and the 1-based position in the CDS (tx).
// [[Rcpp::export]]
List codon_dwell_map_cpp(IntegerVector chr, IntegerVector strand, IntegerVector start, IntegerVector end,
                         IntegerVector lvmap, NumericVector key_start, NumericVector key_end,
                         IntegerVector layer_end, IntegerVector ex_start, IntegerVector ex_end,
                         LogicalVector ex_minus, IntegerVector ex_gene, IntegerVector ex_tx0) {
  const R_xlen_t n = chr.size();
  const int NL = layer_end.size();
  const double* ks = REAL(key_start); const double* ke = REAL(key_end);
  if (strand.size() != n || start.size() != n || end.size() != n) stop("codon_dwell_map_cpp: incompatible inputs");
  std::vector<int> qo, go, to;
  qo.reserve(n); go.reserve(n); to.reserve(n);
  std::vector<int> hint(2 * NL, -1);
  for (R_xlen_t i = 0; i < n; i++) {
    const int c = chr[i];
    if (c == NA_INTEGER || c < 1 || c > lvmap.size()) continue;
    const int d = lvmap[c - 1];
    if (d == NA_INTEGER) continue;
    const int st = strand[i];
    const int pos = st == 2 ? end[i] : start[i];
    for (int s = 0; s < 2; s++) {
      if (st == 1 && s == 1) break;
      if (st == 2 && s == 0) continue;
      const double key = ((double)(d - 1) * 2 + s) * 4294967296.0 + pos;
      int a = 0;
      for (int l = 0; l < NL; l++) {
        const int b = layer_end[l];
        int& h = hint[2 * l + s];
        int j;
        if (h >= a && h < b && ks[h] <= key && (h + 1 == b || key < ks[h + 1])) {
          j = h;
        } else {
          const double* it = std::upper_bound(ks + a, ks + b, key);
          j = (int)(it - ks) - 1;
          if (j < a) { a = b; continue; }
          h = j;
        }
        if (key <= ke[j]) {
          qo.push_back((int)i + 1);
          go.push_back(ex_gene[j]);
          to.push_back(ex_tx0[j] + (ex_minus[j] ? ex_end[j] - pos : pos - ex_start[j]) + 1);
        }
        a = b;
      }
    }
  }
  return List::create(_["idx"] = wrap(qo), _["gene"] = wrap(go), _["tx"] = wrap(to));
}

// Counts of codon_dwell_counts(): gene_total[g] += w (all codons and frames), y[row] += wy for P-site codons inside
// the trimmed window (row = row_start[g] + codon - trim_start, 1-based).
// [[Rcpp::export]]
List codon_dwell_counts_cpp(IntegerVector gene, IntegerVector tx, NumericVector w, NumericVector wy,
                            IntegerVector ncodons, NumericVector row_start, int trim_start, int trim_end,
                            int nrow) {
  const R_xlen_t n = gene.size();
  const int G = ncodons.size();
  if (tx.size() != n || w.size() != n || wy.size() != n || row_start.size() != G)
    stop("codon_dwell_counts_cpp: incompatible inputs");
  NumericVector gt(G), y(nrow);
  for (R_xlen_t i = 0; i < n; i++) {
    const int g = gene[i] - 1;
    if (g < 0 || g >= G) stop("codon_dwell_counts_cpp: gene index out of bounds");
    gt[g] += w[i];
    const int cod = (tx[i] - 1) / 3 + 1;
    if (cod > trim_start && cod <= ncodons[g] - trim_end) {
      const R_xlen_t r = (R_xlen_t)row_start[g] + cod - trim_start - 1;
      if (r < 0 || r >= nrow) stop("codon_dwell_counts_cpp: row out of bounds");
      y[r] += wy[i];
    }
  }
  return List::create(_["y"] = y, _["gene_total"] = gt);
}
