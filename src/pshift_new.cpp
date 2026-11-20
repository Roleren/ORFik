// pshiftR core (backclaude, 2026-11-17): one streaming pass over an unshifted ofst for ORFik-equivalent P-site detection.
// Per read: 5' end (strand aware, from start + reference span of the CIGAR), size (aligned length without N gaps, as
// ORFik readWidths(along.reference = TRUE, N.regions.removed = TRUE)), chromosome gate (seqlevels shared with the CDS set),
// CDS-overlap gate (reads overlapping any CDS exon, ORFik countOverlapsW(cds, footprints, "score") > 0), and for accepted
// sizes the counts at transcript positions -30..+149 relative to the CDS start of every eligible transcript.
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <set>
#include <cstring>
using namespace Rcpp;

static const int WLO = -30, WHI = 149, WLEN = WHI - WLO + 1;  // 180 window positions (ORFik start window + first 150 CDS nt)
// extended storage (opt-in): start -45..-31, stop window -60..+30 relative to the last CDS nt, CDS-body 5'-end frame counts
static const int XS_LO = -45, XS_N = 15, ST_LO = -60, ST_N = 91, EXT = XS_N + ST_N + 3;

// ---- disjoint segments of intervals in key space (key = ((chr * 2 + minus) << 30) | pos) ----
// [[Rcpp::export]]
List psr_segments_cpp(NumericVector ks, NumericVector ke) {
  const R_xlen_t n = ks.size();
  std::vector<std::pair<int64_t, int>> ev; ev.reserve(2 * n);
  std::vector<int64_t> b; b.reserve(2 * n);
  for (R_xlen_t i = 0; i < n; ++i) {
    int64_t s = (int64_t) ks[i], e = (int64_t) ke[i];
    if (e < s) stop("interval with end < start");
    ev.push_back({s, (int)(i + 1)}); ev.push_back({e + 1, -(int)(i + 1)});
    b.push_back(s); b.push_back(e + 1);
  }
  std::sort(b.begin(), b.end()); b.erase(std::unique(b.begin(), b.end()), b.end());
  std::sort(ev.begin(), ev.end());
  std::set<int> active;
  std::vector<double> ss, se; std::vector<int> ptr(1, 0), ids;
  size_t q = 0;
  for (size_t i = 0; i + 1 < b.size(); ++i) {
    while (q < ev.size() && ev[q].first == b[i]) {
      int v = ev[q].second;
      if (v > 0) active.insert(v - 1); else active.erase(-v - 1);
      ++q;
    }
    if (active.empty()) continue;
    ss.push_back((double) b[i]); se.push_back((double)(b[i + 1] - 1));
    for (int a : active) ids.push_back(a);
    ptr.push_back((int) ids.size());
  }
  return List::create(_["seg_start"] = wrap(ss), _["seg_end"] = wrap(se), _["seg_ptr"] = wrap(ptr), _["seg_exon"] = wrap(ids));
}

// ---- CIGAR table: span along reference (with N) and aligned size (M, D, =, X; N removed) for factor levels ----
// [[Rcpp::export]]
DataFrame psr_cigar_cpp(CharacterVector cig) {
  const R_xlen_t n = cig.size();
  IntegerVector span(n), size(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const char *p = CHAR(STRING_ELT(cig, i));
    long num = 0; int sp = 0, sz = 0;
    for (; *p; ++p) {
      if (*p >= '0' && *p <= '9') { num = num * 10 + (*p - '0'); continue; }
      switch (*p) {
        case 'M': case 'D': case '=': case 'X': sp += num; sz += num; break;
        case 'N': sp += num; break;
        default: break;  // S, H, I, P
      }
      num = 0;
    }
    span[i] = sp; size[i] = sz;
  }
  return DataFrame::create(_["span"] = span, _["size"] = size);
}

struct Seg {
  std::vector<int64_t> s, e; const int *ptr = nullptr, *ex = nullptr; bool ok = false;
};
static Seg make_seg(const List &idx) {
  Seg S;
  NumericVector s = idx["seg_start"], e = idx["seg_end"];
  S.s.assign(s.begin(), s.end()); S.e.assign(e.begin(), e.end());
  IntegerVector ptr = idx["seg_ptr"], ex = idx["seg_exon"];
  S.ptr = ptr.begin(); S.ex = ex.begin(); S.ok = true; return S;
}
static inline long find_seg(const Seg &I, int64_t key) {
  auto it = std::upper_bound(I.s.begin(), I.s.end(), key);
  if (it == I.s.begin()) return -1;
  long j = (long)(it - I.s.begin()) - 1;
  return key <= I.e[j] ? j : -1;
}

struct Acc {
  Seg TX, CDS;  // TX: segments of eligible transcript exons (exon -> tx, offset); CDS: segments of CDS exons
  IntegerVector tx_ex_tx, tx_ex_start, tx_ex_end, tx_ex_off, tx_ex_minus, cds_ex_tx, cds_start;  // keep alive
  List txidx, cdsidx;
  int ntx = 0, nsz = 0; bool ext = false; int blen = WLEN; std::vector<int> cds_end;
  std::vector<int> size_slot;         // size -> accepted size index or -1
  std::vector<int> slot;              // tx * nsz + k -> block index or -1
  std::vector<double> store;          // blocks of WLEN
  std::vector<int> blk_tx, blk_sz;
  std::vector<double> tot;            // per size (0..255) all reads passing the chromosome gate
  std::vector<char> gate;             // per tx: any read overlaps a CDS exon
  double reads_total = 0, reads_dropped_chr = 0, reads_skipped = 0;
};

// [[Rcpp::export]]
SEXP psr_new_cpp(List txidx, IntegerVector ex_tx, IntegerVector ex_start, IntegerVector ex_end, IntegerVector ex_off,
                 IntegerVector ex_minus, List cdsidx, IntegerVector cds_ex_tx, IntegerVector cds_start,
                 IntegerVector sizes, IntegerVector cds_end, bool extended) {
  Acc *a = new Acc();
  a->ext = extended; a->blen = WLEN + (extended ? EXT : 0); a->cds_end.assign(cds_end.begin(), cds_end.end());
  a->txidx = txidx; a->cdsidx = cdsidx;
  a->tx_ex_tx = ex_tx; a->tx_ex_start = ex_start; a->tx_ex_end = ex_end; a->tx_ex_off = ex_off; a->tx_ex_minus = ex_minus;
  a->cds_ex_tx = cds_ex_tx; a->cds_start = cds_start;
  a->TX = make_seg(txidx); a->CDS = make_seg(cdsidx);
  a->ntx = cds_start.size(); a->nsz = sizes.size();
  a->size_slot.assign(256, -1);
  for (int k = 0; k < sizes.size(); ++k) if (sizes[k] >= 0 && sizes[k] < 256) a->size_slot[sizes[k]] = k;
  a->slot.assign((size_t) a->ntx * a->nsz, -1);
  a->tot.assign(256, 0.0); a->gate.assign(a->ntx, 0);
  XPtr<Acc> p(a, true); return p;
}

// one chunk of ofst rows. chr: chromosome id per row (-1 = not in the CDS seqlevels), strand: 1 '+', 2 '-', 3 '*'
// [[Rcpp::export]]
void psr_add_cpp(SEXP ptr, IntegerVector chr, IntegerVector start, IntegerVector strand, IntegerVector cig,
                 NumericVector score, IntegerVector cig_span, IntegerVector cig_size, CharacterVector cig_levels) {
  XPtr<Acc> A(ptr); Acc &a = *A;
  const R_xlen_t n = chr.size();
  const int *txs = a.tx_ex_start.begin(), *txe = a.tx_ex_end.begin(), *txo = a.tx_ex_off.begin(), *txm = a.tx_ex_minus.begin(),
            *txt = a.tx_ex_tx.begin(), *cs = a.cds_start.begin();
  (void) txs;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double sc = score[i];
    a.reads_total += sc;
    int ch = chr[i], st = start[i], cc = cig[i];
    if (ch < 0) { a.reads_dropped_chr += sc; continue; }
    if (cc == NA_INTEGER || st == NA_INTEGER) { a.reads_skipped += sc; continue; }
    const int sp = cig_span[cc - 1], sz = cig_size[cc - 1];
    if (sz < 0 || sz >= 256) { a.reads_skipped += sc; continue; }
    a.tot[sz] += sc;
    const bool minus = strand[i] == 2;
    const int pos5 = minus ? st + sp - 1 : st, last = st + sp - 1;
    const int64_t base = ((int64_t)(ch * 2 + (minus ? 1 : 0))) << 30;
    // CDS overlap gate: any CDS exon overlapping an aligned block (5' end only)
    {
      auto gate_block = [&](int bs, int be) {
        const int64_t k0 = base | (int64_t) bs, k1 = base | (int64_t) be;
        auto it = std::upper_bound(a.CDS.s.begin(), a.CDS.s.end(), k0);
        long j = (it == a.CDS.s.begin()) ? 0 : (long)(it - a.CDS.s.begin()) - 1;
        for (; j < (long) a.CDS.s.size() && a.CDS.s[j] <= k1; ++j) {
          if (a.CDS.e[j] < k0) continue;
          for (int q = a.CDS.ptr[j]; q < a.CDS.ptr[j + 1]; ++q) a.gate[a.cds_ex_tx[a.CDS.ex[q]]] = 1;
        }
      };
      gate_block(pos5, pos5);  // ORFik gates on the one-based 5' end (convertToOneBasedRanges runs before countOverlapsW)
    }
    const int k = a.size_slot[sz];
    if (k < 0) continue;
    const int64_t key = base | (int64_t) pos5;
    long j = find_seg(a.TX, key);
    if (j < 0) continue;
    for (int q = a.TX.ptr[j]; q < a.TX.ptr[j + 1]; ++q) {
      const int x = a.TX.ex[q], t = txt[x];
      const int tp = txo[x] + (txm[x] ? (txe[x] - pos5) : (pos5 - txs[x]));   // 1-based transcript coordinate
      const int rel = tp - cs[t];                                                       // 0 at the CDS start
      int col = -1;
      if (rel >= WLO && rel <= WHI) col = rel - WLO;
      int xcol = -1;
      if (a.ext) {
        const int ce = a.cds_end[t], rs = tp - ce;
        if (rel >= XS_LO && rel < WLO) xcol = WLEN + (rel - XS_LO);
        else if (rs >= ST_LO && rs < ST_LO + ST_N) xcol = WLEN + XS_N + (rs - ST_LO);
        if (rel >= 30 && tp <= ce - 30) {  // CDS body (outside the first 10 and last 10 codons): 5'-end frame counts
          int &sl2 = a.slot[(size_t) t * a.nsz + k];
          if (sl2 < 0) { sl2 = (int) a.blk_tx.size(); a.blk_tx.push_back(t); a.blk_sz.push_back(k); a.store.resize(a.store.size() + a.blen, 0.0); }
          a.store[(size_t) sl2 * a.blen + WLEN + XS_N + ST_N + (rel % 3)] += sc;
        }
      }
      if (col < 0 && xcol < 0) continue;
      int &sl = a.slot[(size_t) t * a.nsz + k];
      if (sl < 0) { sl = (int) a.blk_tx.size(); a.blk_tx.push_back(t); a.blk_sz.push_back(k); a.store.resize(a.store.size() + a.blen, 0.0); }
      if (col >= 0) a.store[(size_t) sl * a.blen + col] += sc;
      if (xcol >= 0) a.store[(size_t) sl * a.blen + xcol] += sc;
    }
  }
}

// [[Rcpp::export]]
List psr_result_cpp(SEXP ptr) {
  XPtr<Acc> A(ptr); Acc &a = *A;
  const int nb = (int) a.blk_tx.size();
  NumericMatrix W(nb, a.blen);
  for (int b = 0; b < nb; ++b) for (int c = 0; c < a.blen; ++c) W(b, c) = a.store[(size_t) b * a.blen + c];
  LogicalVector gate(a.ntx); for (int t = 0; t < a.ntx; ++t) gate[t] = a.gate[t] != 0;
  return List::create(_["tx"] = wrap(a.blk_tx), _["size_idx"] = wrap(a.blk_sz), _["win"] = W, _["gate"] = gate,
                      _["size_total"] = wrap(a.tot), _["reads_total"] = a.reads_total, _["reads_dropped_chr"] = a.reads_dropped_chr,
                      _["reads_skipped"] = a.reads_skipped);
}
