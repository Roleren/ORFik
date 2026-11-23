#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <memory>
#include <thread>
#include <vector>

// Linear-time kernels for incremental OFST merges and the OFST store (R/utils_ofst_store.R).
// Keys are k integer columns compared lexicographically (e.g. seqnames code, start, strand
// sort key, size). Inputs given as sorted must be strictly increasing (sorted and unique);
// the kernels check this while they scan and stop with an error otherwise.
//
// Packed keys: before matching, the k columns are packed into W 64-bit words per row (one
// comparison per word instead of one per column). Each column takes ceil(log2(max - min + 1))
// bits after subtracting its minimum over all tables of the call, in column order from the
// most significant bit of word 0, so that word-wise unsigned comparison gives the same order as
// the column-wise comparison. A column never straddles two words; every column fits in 32 bits,
// so W <= ceil(k / 2). Typical OFST keys (seq, start, strand, size) take one word.
// Packing is per call and in memory only: the store files keep the plain columns.
// Kernels that use several threads (std::thread, fork-safe, no R API inside threads) split
// the key range at sampled splitter keys, so equal keys always land in the same thread.

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

// ---- packing ---------------------------------------------------------------------------------

struct Packer {
  size_t k = 0;
  int W = 1;
  std::vector<int64_t> lo;
  std::vector<int> bits, word, shift;
};

// Layout from the column ranges of all tables (empty tables are ignored). Missing values
// (NA_integer_) are not valid keys.
// With na_out, a missing value sets *na_out instead of stopping (the layout is then unusable).
Packer make_packer(const std::vector<const KeyCols*>& sets, bool* na_out = nullptr) {
  Packer pk;
  pk.k = sets[0]->p.size();
  for (const KeyCols* s : sets) if (s->p.size() != pk.k) Rcpp::stop("key column counts differ");
  if (pk.k == 0) Rcpp::stop("no key columns");
  std::vector<int64_t> hi(pk.k, INT64_MIN);
  pk.lo.assign(pk.k, INT64_MAX);
  for (const KeyCols* s : sets) for (size_t c = 0; c < pk.k; ++c) {
    const int* x = s->p[c];
    int mn = INT_MAX, mx = INT_MIN;
    bool na = false;
    for (R_xlen_t i = 0; i < s->n; ++i) {
      const int v = x[i];
      na |= (v == NA_INTEGER);
      mn = v < mn ? v : mn; mx = v > mx ? v : mx;
    }
    if (na) {
      if (na_out) { *na_out = true; return pk; }
      Rcpp::stop("missing value (NA) in key column %d", (int)(c + 1));
    }
    if (s->n) { pk.lo[c] = std::min<int64_t>(pk.lo[c], mn); hi[c] = std::max<int64_t>(hi[c], mx); }
  }
  pk.bits.assign(pk.k, 0); pk.word.assign(pk.k, 0); pk.shift.assign(pk.k, 0);
  int w = 0, free_bits = 64;
  for (size_t c = 0; c < pk.k; ++c) {
    int b = 0;
    if (hi[c] > pk.lo[c]) { uint64_t r = (uint64_t)(hi[c] - pk.lo[c]); while (r) { ++b; r >>= 1; } }   // <= 32
    if (b > free_bits) { ++w; free_bits = 64; }
    pk.bits[c] = b; pk.word[c] = w; pk.shift[c] = free_bits - b;
    free_bits -= b;
  }
  pk.W = w + 1;
  return pk;
}

template <int W> using Key = std::array<uint64_t, W>;

template <int W>
inline Key<W> pack_row(const KeyCols& kc, R_xlen_t i, const Packer& pk) {
  Key<W> key{};
  for (size_t c = 0; c < pk.k; ++c) {
    if (!pk.bits[c]) continue;
    key[pk.word[c]] |= (uint64_t)((int64_t)kc.p[c][i] - pk.lo[c]) << pk.shift[c];
  }
  return key;
}

// Packed rows of a table; returns the first row (1-based) that is not strictly greater than
// its predecessor, or 0 when the table is sorted and unique (check_sorted).
template <int W>
R_xlen_t pack_table(const KeyCols& kc, const Packer& pk, Key<W>* out, bool check_sorted) {
  R_xlen_t bad = 0;
  for (R_xlen_t i = 0; i < kc.n; ++i) {
    out[i] = pack_row<W>(kc, i, pk);
    if (check_sorted && !bad && i > 0 && !(out[i - 1] < out[i])) bad = i + 1;
  }
  return bad;
}

int n_threads(int threads, double work) {
  int t = std::max(1, threads);
  if (work < 2e5) t = 1;                      // thread start-up is not worth it
  const unsigned hw = std::thread::hardware_concurrency();
  if (hw > 0 && t > (int)hw) t = (int)hw;
  return t;
}

template <class F>
void run_threads(int t, F f) {
  if (t <= 1) { f(0); return; }
  std::vector<std::thread> th;
  th.reserve(t - 1);
  for (int i = 1; i < t; ++i) th.emplace_back(f, i);
  f(0);
  for (auto& x : th) x.join();
}

// ---- k-way match of sorted tables (the store's key union) ------------------------------------

template <int W>
struct KwayResult {
  std::vector<int> ord;           // merged key order: row ids (provisional -(i + 1) for new key i)
  std::vector<Key<W>> newkeys;    // packed new keys, in key order
  std::vector<Key<W>> uni;        // the key union of the range (only when the sorted keys are wanted)
};

template <int W>
inline uint64_t unpack_col(const Key<W>& key, const Packer& pk, size_t c) {
  if (!pk.bits[c]) return 0;
  const uint64_t mask = pk.bits[c] == 64 ? ~0ull : ((1ull << pk.bits[c]) - 1);
  return (key[pk.word[c]] >> pk.shift[c]) & mask;
}

// Position of each key of sorted run q[0..m) in the sorted superset U (galloping from the
// previous hit; every key of q is in U).
template <int W>
void locate_in(const Key<W>* U, size_t nu, const Key<W>* q, R_xlen_t m, std::vector<uint32_t>& out) {
  out.resize(m);
  size_t i = 0;
  for (R_xlen_t j = 0; j < m; ++j) {
    const Key<W>& key = q[j];
    if (U[i] < key) {
      size_t lo = i + 1, step = 1, hi = i + 1;
      while (hi < nu && U[hi] < key) { lo = hi + 1; hi += step; step <<= 1; }
      if (hi > nu) hi = nu;
      i = std::lower_bound(U + lo, U + hi, key) - U;
    }
    out[j] = (uint32_t)i;
  }
}

// One key range [lo, hi) of every stream. Streams: s = 0 the stored keys (sorted; ids from
// old_order, may be empty), s >= 1 the new tables. (1) Union U of the range: set_union
// tree (branch-light linear merges, duplicates collapse level by level). (2) Each stream is
// located in U by galloping. Stored keys keep their ids; the others get provisional ids in U
// order; rows[s - 1][p] receives the id of key p of stream s.
template <int W>
void kway_range(const std::vector<const Key<W>*>& S, const std::vector<R_xlen_t>& lo,
                const std::vector<R_xlen_t>& hi, const int* old_order, const std::vector<int*>& rows,
                KwayResult<W>& res, bool want_union) {
  const int ns = (int)S.size();
  // Union tree: always merge the two smallest spans (Huffman order), so the large stored span
  // is copied once at the end instead of once per tree level.
  struct Span { const Key<W>* p; size_t n; std::shared_ptr<std::vector<Key<W>>> own; };
  std::vector<Span> cur;
  for (int s = 0; s < ns; ++s) if (hi[s] > lo[s]) cur.push_back({S[s] + lo[s], (size_t)(hi[s] - lo[s]), nullptr});
  auto bigger = [](const Span& a, const Span& b) { return a.n > b.n; };   // min-heap on size
  std::make_heap(cur.begin(), cur.end(), bigger);
  while (cur.size() > 1) {
    std::pop_heap(cur.begin(), cur.end(), bigger); Span a = cur.back(); cur.pop_back();
    std::pop_heap(cur.begin(), cur.end(), bigger); Span b = cur.back(); cur.pop_back();
    auto u = std::make_shared<std::vector<Key<W>>>(a.n + b.n);
    auto e = std::set_union(a.p, a.p + a.n, b.p, b.p + b.n, u->begin());
    u->resize(e - u->begin());
    a.own.reset(); b.own.reset();   // frees merged intermediates
    cur.push_back({u->data(), u->size(), u});
    std::push_heap(cur.begin(), cur.end(), bigger);
  }
  if (cur.empty()) return;
  const Key<W>* U = cur[0].p;
  const size_t nu = cur[0].n;
  if (want_union) res.uni.assign(U, U + nu);
  std::vector<int> id(nu, 0);
  std::vector<uint32_t> loc;
  if (hi[0] > lo[0]) {
    locate_in<W>(U, nu, S[0] + lo[0], hi[0] - lo[0], loc);
    for (size_t j = 0; j < loc.size(); ++j) id[loc[j]] = old_order[lo[0] + (R_xlen_t)j];
  }
  int n_new = 0;
  for (size_t u = 0; u < nu; ++u) if (id[u] == 0) { id[u] = -(++n_new); res.newkeys.push_back(U[u]); }
  for (int s = 1; s < ns; ++s) {
    if (hi[s] <= lo[s]) continue;
    locate_in<W>(U, nu, S[s] + lo[s], hi[s] - lo[s], loc);
    int* r = rows[s - 1] + lo[s];
    for (size_t j = 0; j < loc.size(); ++j) r[j] = id[loc[j]];
  }
  res.ord.swap(id);
}

template <int W>
Rcpp::List kway_impl(const KeyCols& stored, const int* old_order, const std::vector<KeyCols>& inputs,
                     const Packer& pk, double n_old, int threads, bool stored_sorted, bool want_union) {
  auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  std::vector<double> tm{now()};
  const int ns = (int)inputs.size() + 1;
  std::vector<std::vector<Key<W>>> packed(ns);
  std::vector<const KeyCols*> src(ns);
  src[0] = &stored;
  for (int s = 1; s < ns; ++s) src[s] = &inputs[s - 1];
  double total = 0;
  for (int s = 0; s < ns; ++s) { packed[s].resize(src[s]->n); total += (double)src[s]->n; }
  const int T = n_threads(threads, total);
  // pack (and check the order) in parallel over streams
  std::vector<R_xlen_t> bad(ns, 0);
  {
    std::atomic<int> next(0);
    run_threads(std::min(T, ns), [&](int) {
      for (int s; (s = next++) < ns;) {
        if (s > 0) { bad[s] = pack_table<W>(*src[s], pk, packed[s].data(), true); continue; }
        // stored keys: gathered into key order through old_order while packing (or, when they
        // are given in key order already, read sequentially)
        Key<W>* out = packed[0].data();
        for (R_xlen_t i = 0; i < stored.n; ++i) {
          out[i] = pack_row<W>(stored, stored_sorted ? i : (R_xlen_t)old_order[i] - 1, pk);
          if (!bad[0] && i > 0 && !(out[i - 1] < out[i])) bad[0] = i + 1;
        }
      }
    });
  }
  for (int s = 0; s < ns; ++s) if (bad[s]) {
    if (s == 0) Rcpp::stop("stored keys are not sorted and unique in the given order (position %d)", (double)bad[s]);
    Rcpp::stop("input %d: keys are not sorted and unique (row %d)", s, (double)bad[s]);
  }
  std::vector<const Key<W>*> S(ns);
  for (int s = 0; s < ns; ++s) S[s] = packed[s].data();
  tm.push_back(now());
  // outputs allocated here (R API only in this thread)
  Rcpp::List rows_out(ns - 1);
  std::vector<int*> rows(ns - 1);
  for (int s = 1; s < ns; ++s) { Rcpp::IntegerVector v(src[s]->n); rows[s - 1] = v.begin(); rows_out[s - 1] = v; }
  // splitters: sampled keys, then per-thread [lo, hi) positions in every stream
  std::vector<Key<W>> split;
  if (T > 1) {
    std::vector<Key<W>> sample;
    const double want = 256.0 * T;
    for (int s = 0; s < ns; ++s) {
      const R_xlen_t n = src[s]->n;
      if (!n) continue;
      const R_xlen_t m = std::max<R_xlen_t>(1, (R_xlen_t)(want * (double)n / total));
      for (R_xlen_t i = 0; i < m; ++i) sample.push_back(S[s][(R_xlen_t)((double)i * (double)n / (double)m)]);
    }
    std::sort(sample.begin(), sample.end());
    for (int t = 1; t < T; ++t) split.push_back(sample[(size_t)((double)t * sample.size() / T)]);
  }
  const int nt = (int)split.size() + 1;
  std::vector<std::vector<R_xlen_t>> bound(nt + 1, std::vector<R_xlen_t>(ns));
  for (int s = 0; s < ns; ++s) {
    bound[0][s] = 0; bound[nt][s] = src[s]->n;
    for (int t = 1; t < nt; ++t)
      bound[t][s] = std::lower_bound(S[s], S[s] + src[s]->n, split[t - 1]) - S[s];
  }
  std::vector<KwayResult<W>> res(nt);
  run_threads(nt, [&](int t) {   // thread-local results (no false sharing), moved out at the end
    KwayResult<W> local;
    kway_range<W>(S, bound[t], bound[t + 1], old_order, rows, local, want_union);
    res[t] = std::move(local);
  });
  tm.push_back(now());
  std::vector<double> new_off(nt + 1, 0), ord_off(nt + 1, 0);
  for (int t = 0; t < nt; ++t) {
    new_off[t + 1] = new_off[t] + (double)res[t].newkeys.size();
    ord_off[t + 1] = ord_off[t] + (double)res[t].ord.size();
  }
  const double n_new = new_off[nt];
  if (n_old + n_new > (double)INT_MAX) Rcpp::stop("more than 2^31 - 1 keys in one partition");
  packed.clear(); packed.shrink_to_fit();
  Rcpp::IntegerVector ord((R_xlen_t)ord_off[nt]);
  int* o = ord.begin();
  const size_t k = pk.k;
  Rcpp::List new_keys(k);
  std::vector<int*> nk(k);
  for (size_t c = 0; c < k; ++c) { Rcpp::IntegerVector v((R_xlen_t)n_new); nk[c] = v.begin(); new_keys[c] = v; }
  // all keys in key order (the union), when wanted
  Rcpp::List all_keys(want_union ? k : 0);
  std::vector<int*> ak(k, nullptr);
  if (want_union)
    for (size_t c = 0; c < k; ++c) { Rcpp::IntegerVector v((R_xlen_t)ord_off[nt]); ak[c] = v.begin(); all_keys[c] = v; }
  // pass 2: final ids of new keys (n_old + global new index), merged order, new key columns
  run_threads(nt, [&](int t) {
    const int base = (int)(n_old + new_off[t]);
    for (int s = 1; s < ns; ++s) {
      int* r = rows[s - 1];
      for (R_xlen_t p = bound[t][s]; p < bound[t + 1][s]; ++p) if (r[p] < 0) r[p] = base - r[p];
    }
    int* oo = o + (R_xlen_t)ord_off[t];
    for (size_t i = 0; i < res[t].ord.size(); ++i) { const int v = res[t].ord[i]; oo[i] = v < 0 ? base - v : v; }
    const R_xlen_t first = (R_xlen_t)new_off[t];
    for (size_t i = 0; i < res[t].newkeys.size(); ++i)
      for (size_t c = 0; c < k; ++c)
        nk[c][first + (R_xlen_t)i] = (int)((int64_t)unpack_col<W>(res[t].newkeys[i], pk, c) + pk.lo[c]);
    if (want_union) {
      const R_xlen_t off = (R_xlen_t)ord_off[t];
      for (size_t c = 0; c < k; ++c) {
        int* a = ak[c] + off;
        const int64_t lo_c = pk.lo[c];
        for (size_t i = 0; i < res[t].uni.size(); ++i) a[i] = (int)((int64_t)unpack_col<W>(res[t].uni[i], pk, c) + lo_c);
      }
      std::vector<Key<W>>().swap(res[t].uni);
    }
  });
  tm.push_back(now());
  Rcpp::NumericVector secs = Rcpp::NumericVector::create(
    Rcpp::Named("pack") = tm[1] - tm[0], Rcpp::Named("match") = tm[2] - tm[1], Rcpp::Named("finish") = tm[3] - tm[2]);
  return Rcpp::List::create(Rcpp::Named("rows") = rows_out, Rcpp::Named("order") = ord,
                            Rcpp::Named("new_keys") = new_keys, Rcpp::Named("keys") = all_keys,
                            Rcpp::Named("words") = pk.W,
                            Rcpp::Named("threads") = nt, Rcpp::Named("seconds") = secs);
}

// ---- packed two-pointer kernels ---------------------------------------------------------------

template <int W>
Rcpp::IntegerVector match_impl(const KeyCols& t, const KeyCols& q, const Packer& pk) {
  std::vector<Key<W>> pt(t.n), pq(q.n);
  pack_table<W>(t, pk, pt.data(), false);
  const R_xlen_t bad = pack_table<W>(q, pk, pq.data(), true);
  if (bad) Rcpp::stop("query keys are not sorted and unique (row %d)", (double)bad);
  const R_xlen_t n = t.n, m = q.n;
  Rcpp::IntegerVector out(m);
  R_xlen_t i = 0;
  for (R_xlen_t j = 0; j < m; ++j) {
    const Key<W>& key = pq[j];
    R_xlen_t lo = i, step = 1, hi = i;   // galloping search from the previous hit
    while (hi < n && pt[hi] < key) { lo = hi + 1; hi += step; step <<= 1; }
    if (hi > n) hi = n;
    i = std::lower_bound(pt.begin() + lo, pt.begin() + hi, key) - pt.begin();
    out[j] = (i < n && pt[i] == key) ? (int)(i + 1) : 0;
  }
  return out;
}

template <int W>
Rcpp::IntegerVector merge_order_impl(const KeyCols& ka, const int* oa, const KeyCols& kb, const int* rb,
                                     const Packer& pk) {
  std::vector<Key<W>> pa(ka.n), pb(kb.n);
  pack_table<W>(ka, pk, pa.data(), false);
  pack_table<W>(kb, pk, pb.data(), false);
  const R_xlen_t n = ka.n, m = kb.n;
  Rcpp::IntegerVector out(n + m);
  int* o = out.begin();
  R_xlen_t i = 0, j = 0, k = 0;
  while (i < n && j < m) {
    if (pa[i] == pb[j]) Rcpp::stop("a new key is already in the table (row %d)", (double)oa[i]);
    o[k++] = pa[i] < pb[j] ? oa[i++] : rb[j++];
  }
  while (i < n) o[k++] = oa[i++];
  while (j < m) o[k++] = rb[j++];
  return out;
}

// ---- delta updates: locate new tables in the stored key streams -------------------------------

// Pack a large table over t threads (contiguous chunks); returns the first row (1-based) that is
// not strictly greater than its predecessor, or 0.
template <int W>
R_xlen_t pack_table_par(const KeyCols& kc, const Packer& pk, Key<W>* out, int t) {
  const R_xlen_t n = kc.n;
  t = std::max(1, std::min<int>(t, (int)std::max<R_xlen_t>(1, n / 1000000)));
  std::vector<R_xlen_t> bad(t, 0);
  run_threads(t, [&](int i) {
    const R_xlen_t a = n * i / t, b = n * (i + 1) / t;
    for (R_xlen_t r = a; r < b; ++r) {
      out[r] = pack_row<W>(kc, r, pk);
      if (!bad[i] && r > a && !(out[r - 1] < out[r])) bad[i] = r + 1;
    }
  });
  for (int i = 0; i < t; ++i) {
    if (bad[i]) return bad[i];
    const R_xlen_t a = n * i / t;
    if (i > 0 && a > 0 && a < n && !(out[a - 1] < out[a])) return a + 1;
  }
  return 0;
}

// Stored streams S[s] (sorted, disjoint, ids ids[s]) and new tables Q[q] (sorted, unique).
// Every key of a new table gets its stored id, or a new id n_old + rank in the sorted union of
// all keys found in no stream. Work is split over the new tables (one thread each).
template <int W>
Rcpp::List locate_impl(const std::vector<KeyCols>& streams, const std::vector<const int*>& ids,
                       const std::vector<KeyCols>& inputs, const Packer& pk, double n_old, int threads) {
  auto now = [] { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); };
  std::vector<double> tm{now()};
  const int ns = (int)streams.size(), nq = (int)inputs.size();
  double work = 0;
  for (auto& x : streams) work += (double)x.n;
  for (auto& x : inputs) work += (double)x.n;
  const int T = n_threads(threads, work);
  std::vector<std::vector<Key<W>>> S(ns), Q(nq);
  for (int s = 0; s < ns; ++s) {   // large streams: chunked over all threads
    S[s].resize(streams[s].n);
    const R_xlen_t bad = pack_table_par<W>(streams[s], pk, S[s].data(), T);
    if (bad) Rcpp::stop("stored key stream %d is not sorted and unique (row %d)", s + 1, (double)bad);
  }
  std::vector<R_xlen_t> qbad(nq, 0);
  {
    std::atomic<int> next(0);
    run_threads(std::min(T, std::max(1, nq)), [&](int) {
      for (int q; (q = next++) < nq;) { Q[q].resize(inputs[q].n); qbad[q] = pack_table<W>(inputs[q], pk, Q[q].data(), true); }
    });
  }
  for (int q = 0; q < nq; ++q) if (qbad[q]) Rcpp::stop("input %d: keys are not sorted and unique (row %d)", q + 1, (double)qbad[q]);
  tm.push_back(now());
  Rcpp::List rows_out(nq);
  std::vector<int*> rows(nq);
  for (int q = 0; q < nq; ++q) { Rcpp::IntegerVector v(inputs[q].n); rows[q] = v.begin(); rows_out[q] = v; }
  // pass 1: stored ids by galloping search in every stream; misses collected per table
  std::vector<std::vector<Key<W>>> miss(nq);
  std::vector<std::vector<R_xlen_t>> miss_pos(nq);
  {
    std::atomic<int> next(0);
    run_threads(std::min(T, std::max(1, nq)), [&](int) {
      for (int q; (q = next++) < nq;) {
        std::vector<size_t> cur(ns, 0);
        const Key<W>* qq = Q[q].data();
        int* r = rows[q];
        for (R_xlen_t j = 0; j < inputs[q].n; ++j) {
          const Key<W>& key = qq[j];
          int id = 0;
          for (int s = 0; s < ns && !id; ++s) {
            const Key<W>* U = S[s].data();
            const size_t nu = S[s].size();
            size_t i = cur[s];
            if (i < nu && U[i] < key) {
              size_t lo = i + 1, step = 1, hi = i + 1;
              while (hi < nu && U[hi] < key) { lo = hi + 1; hi += step; step <<= 1; }
              if (hi > nu) hi = nu;
              i = std::lower_bound(U + lo, U + hi, key) - U;
              cur[s] = i;
            }
            if (i < nu && U[i] == key) id = ids[s][i];
          }
          r[j] = id;
          if (!id) { miss[q].push_back(key); miss_pos[q].push_back(j); }
        }
      }
    });
  }
  tm.push_back(now());
  // union of the misses (Huffman order of pairwise set unions)
  std::vector<std::shared_ptr<std::vector<Key<W>>>> cur;
  for (int q = 0; q < nq; ++q) if (!miss[q].empty()) cur.push_back(std::make_shared<std::vector<Key<W>>>(miss[q]));
  auto bigger = [](const std::shared_ptr<std::vector<Key<W>>>& a, const std::shared_ptr<std::vector<Key<W>>>& b) { return a->size() > b->size(); };
  std::make_heap(cur.begin(), cur.end(), bigger);
  while (cur.size() > 1) {
    std::pop_heap(cur.begin(), cur.end(), bigger); auto a = cur.back(); cur.pop_back();
    std::pop_heap(cur.begin(), cur.end(), bigger); auto b = cur.back(); cur.pop_back();
    auto u = std::make_shared<std::vector<Key<W>>>(a->size() + b->size());
    auto e = std::set_union(a->begin(), a->end(), b->begin(), b->end(), u->begin());
    u->resize(e - u->begin());
    cur.push_back(u);
    std::push_heap(cur.begin(), cur.end(), bigger);
  }
  std::vector<Key<W>> U;
  if (!cur.empty()) U.swap(*cur[0]);
  cur.clear();
  const double n_new = (double)U.size();
  if (n_old + n_new > (double)INT_MAX) Rcpp::stop("more than 2^31 - 1 keys in one partition");
  // pass 2: new ids
  {
    std::atomic<int> next(0);
    run_threads(std::min(T, std::max(1, nq)), [&](int) {
      std::vector<uint32_t> loc;
      for (int q; (q = next++) < nq;) {
        if (miss[q].empty()) continue;
        locate_in<W>(U.data(), U.size(), miss[q].data(), (R_xlen_t)miss[q].size(), loc);
        int* r = rows[q];
        for (size_t j = 0; j < loc.size(); ++j) r[miss_pos[q][j]] = (int)(n_old + 1 + loc[j]);
        std::vector<Key<W>>().swap(miss[q]);
      }
    });
  }
  const size_t k = pk.k;
  Rcpp::List new_keys(k);
  for (size_t c = 0; c < k; ++c) {
    Rcpp::IntegerVector v((R_xlen_t)n_new);
    for (size_t i = 0; i < U.size(); ++i) v[i] = (int)((int64_t)unpack_col<W>(U[i], pk, c) + pk.lo[c]);
    new_keys[c] = v;
  }
  tm.push_back(now());
  Rcpp::NumericVector secs = Rcpp::NumericVector::create(
    Rcpp::Named("pack") = tm[1] - tm[0], Rcpp::Named("locate") = tm[2] - tm[1], Rcpp::Named("finish") = tm[3] - tm[2]);
  return Rcpp::List::create(Rcpp::Named("rows") = rows_out, Rcpp::Named("new_keys") = new_keys,
                            Rcpp::Named("words") = pk.W, Rcpp::Named("threads") = T, Rcpp::Named("seconds") = secs);
}

#define OFST_DISPATCH_W(pk, CALL)                                         \
  switch ((pk).W) {                                                       \
    case 1: return CALL(1);                                               \
    case 2: return CALL(2);                                               \
    case 3: return CALL(3);                                               \
    case 4: return CALL(4);                                               \
    default: Rcpp::stop("too many key columns to pack (%d words)", (pk).W); \
  }

}  // namespace

// Sorted union of two keyed score tables with summed scores. a and b must each be strictly
// increasing. With drop_zero, keys whose summed score is 0 are left out (removal support);
// a negative sum is an error when forbid_negative is TRUE. Two passes: count, then fill, so
// the output is allocated once at its final size. Keys that pack into one word are compared as
// one 64-bit integer.
// [[Rcpp::export]]
Rcpp::List ofst_sorted_sum_merge(Rcpp::List a_keys, Rcpp::NumericVector a_score,
                                 Rcpp::List b_keys, Rcpp::NumericVector b_score,
                                 bool drop_zero, bool forbid_negative) {
  KeyCols a(a_keys), b(b_keys);
  check_same_width(a, b);
  const R_xlen_t n = a.n, m = b.n;
  if (a_score.size() != n || b_score.size() != m) Rcpp::stop("score length differs from key length");
  const double* sa = a_score.begin(); const double* sb = b_score.begin();
  Packer pk = make_packer({&a, &b});
  const bool packed = pk.W == 1;
  std::vector<Key<1>> pa, pb;
  if (packed) {
    pa.resize(n); pb.resize(m);
    R_xlen_t bad = pack_table<1>(a, pk, pa.data(), true);
    if (bad) Rcpp::stop("first input is not sorted and unique (row %d)", (double)bad);
    bad = pack_table<1>(b, pk, pb.data(), true);
    if (bad) Rcpp::stop("second input is not sorted and unique (row %d)", (double)bad);
  } else {
    for (R_xlen_t i = 1; i < n; ++i) if (cmp_rows(a, i - 1, a, i) >= 0) Rcpp::stop("first input is not sorted and unique (row %d)", (double)(i + 1));
    for (R_xlen_t j = 1; j < m; ++j) if (cmp_rows(b, j - 1, b, j) >= 0) Rcpp::stop("second input is not sorted and unique (row %d)", (double)(j + 1));
  }
  auto cmp = [&](R_xlen_t i, R_xlen_t j) -> int {
    if (packed) return pa[i][0] < pb[j][0] ? -1 : (pa[i][0] > pb[j][0] ? 1 : 0);
    return cmp_rows(a, i, b, j);
  };
  // pass 1: count
  R_xlen_t out_n = 0;
  {
    R_xlen_t i = 0, j = 0;
    while (i < n || j < m) {
      double v;
      int c = (i >= n) ? 1 : (j >= m) ? -1 : cmp(i, j);
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
    int c = (i >= n) ? 1 : (j >= m) ? -1 : cmp(i, j);
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

// TRUE if the key rows are strictly increasing (sorted and unique), one linear pass on packed
// keys. Lets callers skip a group-by / sort of inputs that are already in merge order.
template <int W>
bool sorted_unique_impl(const KeyCols& k, const Packer& pk) {
  if (k.n < 2) return true;
  Key<W> prev = pack_row<W>(k, 0, pk);
  for (R_xlen_t i = 1; i < k.n; ++i) {
    const Key<W> cur = pack_row<W>(k, i, pk);
    if (!(prev < cur)) return false;
    prev = cur;
  }
  return true;
}

// [[Rcpp::export]]
bool ofst_keys_sorted_unique(Rcpp::List keys) {
  KeyCols k(keys);
  if (k.p.empty()) Rcpp::stop("no key columns");
  bool na = false;
  Packer pk = make_packer({&k}, &na);
  if (na) {   // NA keys (sort first as INT_MIN): column-wise check
    for (R_xlen_t j = 1; j < k.n; ++j) if (cmp_rows(k, j - 1, k, j) >= 0) return false;
    return true;
  }
#define OFST_CALL(W) sorted_unique_impl<W>(k, pk)
  OFST_DISPATCH_W(pk, OFST_CALL)
#undef OFST_CALL
}

// Number of 64-bit words the keys of these tables pack into (see "Packed keys" above).
// [[Rcpp::export]]
int ofst_packed_key_words(Rcpp::List tables) {
  if (!tables.size()) Rcpp::stop("no tables");
  std::vector<KeyCols> kc;
  kc.reserve(tables.size());
  for (int i = 0; i < tables.size(); ++i) kc.emplace_back(Rcpp::as<Rcpp::List>(tables[i]));
  std::vector<const KeyCols*> sets;
  for (auto& x : kc) sets.push_back(&x);
  return make_packer(sets).W;
}

// Position (1-based) of each query key in a sorted, unique key table, 0 when absent. Query
// must be strictly increasing. Packed keys, galloping search from the previous hit.
// [[Rcpp::export]]
Rcpp::IntegerVector ofst_sorted_match(Rcpp::List table, Rcpp::List query) {
  KeyCols t(table), q(query);
  check_same_width(t, q);
  Packer pk = make_packer({&t, &q});
#define OFST_CALL(W) match_impl<W>(t, q, pk)
  OFST_DISPATCH_W(pk, OFST_CALL)
#undef OFST_CALL
}

// Merge a sorted key table (rows ord_a, contiguous sorted keys a) with sorted new keys b
// (rows rows_b) into the merged row order. Keys must not be shared.
// [[Rcpp::export]]
Rcpp::IntegerVector ofst_sorted_merge_order(Rcpp::List a, Rcpp::IntegerVector ord_a,
                                            Rcpp::List b, Rcpp::IntegerVector rows_b) {
  KeyCols ka(a), kb(b);
  check_same_width(ka, kb);
  if (ord_a.size() != ka.n || rows_b.size() != kb.n) Rcpp::stop("order lengths differ from key lengths");
  Packer pk = make_packer({&ka, &kb});
#define OFST_CALL(W) merge_order_impl<W>(ka, ord_a.begin(), kb, rows_b.begin(), pk)
  OFST_DISPATCH_W(pk, OFST_CALL)
#undef OFST_CALL
}

// The store's update kernel: one k-way pass over the stored keys (rows in id order; `order`
// lists their row ids in key order) and the new tables (each sorted and unique). Returns
//   rows:     per new table, the row id of each of its keys (stored id, or a new id);
//   order:    row ids of the union in key order (the new order file);
//   new_keys: the keys not stored yet, in key order; they get ids n_old + 1, n_old + 2, ...;
//   keys:     with want_keys, all keys of the union in key order (the new sorted key file);
//   words:    packed 64-bit words per key; threads: key ranges processed in parallel.
// threads > 1 splits the key range over std::threads. With stored_sorted, `stored` is already
// in key order (the sorted key file), so it is read sequentially instead of through `order`.
// [[Rcpp::export]]
Rcpp::List ofst_store_kway(Rcpp::List stored, Rcpp::IntegerVector order, Rcpp::List inputs,
                           double n_old, int threads, bool stored_sorted = false, bool want_keys = false) {
  KeyCols st(stored);
  if (order.size() != st.n) Rcpp::stop("order length differs from stored key length");
  for (R_xlen_t i = 0; i < st.n; ++i)
    if (order[i] < 1 || order[i] > st.n) Rcpp::stop("order is not a permutation of the stored rows");
  if ((double)st.n != n_old) Rcpp::stop("n_old differs from the stored key count");
  std::vector<KeyCols> in;
  in.reserve(inputs.size());
  for (int i = 0; i < inputs.size(); ++i) {
    in.emplace_back(Rcpp::as<Rcpp::List>(inputs[i]));
    check_same_width(st, in.back());
  }
  std::vector<const KeyCols*> sets{&st};
  for (auto& x : in) sets.push_back(&x);
  Packer pk = make_packer(sets);
#define OFST_CALL(W) kway_impl<W>(st, order.begin(), in, pk, n_old, threads, stored_sorted, want_keys)
  OFST_DISPATCH_W(pk, OFST_CALL)
#undef OFST_CALL
}

// One pass over a partition's scores for the streaming export: non-zero rows per sequence
// code (seq may be empty: all rows count for code 1) and whether every non-zero score is a
// whole number within the int32 range (then the exported score column is integer).
// [[Rcpp::export]]
Rcpp::List ofst_store_score_scan(Rcpp::NumericVector score, Rcpp::IntegerVector seq, int n_codes) {
  const R_xlen_t n = score.size();
  const bool by_seq = seq.size() > 0;
  if (by_seq && seq.size() != n) Rcpp::stop("seq and score differ in length");
  std::vector<double> counts(std::max(1, n_codes), 0.0);
  bool int_ok = true;
  const double* s = score.begin();
  const int* q = by_seq ? seq.begin() : nullptr;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double v = s[i];
    if (v == 0) continue;
    if (by_seq) {
      const int c = q[i];
      if (c < 1 || c > n_codes) Rcpp::stop("sequence code out of range");
      counts[c - 1] += 1;
    } else {
      counts[0] += 1;
    }
    if (int_ok && !(v == std::floor(v) && std::fabs(v) <= 2147483647.0)) int_ok = false;
  }
  return Rcpp::List::create(Rcpp::Named("counts") = Rcpp::wrap(counts), Rcpp::Named("int_ok") = int_ok);
}

// The store's delta-update kernel. `streams`: stored key tables, each in key order (the sorted
// key file of the base generation, then the pending blocks of new keys), disjoint; `ids`: the
// row id of every stored key (per stream). `inputs`: new tables, each sorted and unique.
// Returns rows (per input, the row id of each key: a stored id, or n_old + rank among the keys
// found in no stream) and new_keys (those keys, in key order). Unlike ofst_store_kway, the
// stored keys are only searched (galloping per input), never merged, so the cost is a pass to
// pack them plus O(input rows * log gap).
// [[Rcpp::export]]
Rcpp::List ofst_store_locate(Rcpp::List streams, Rcpp::List ids, Rcpp::List inputs, double n_old, int threads) {
  if (streams.size() != ids.size()) Rcpp::stop("one id vector per stream is required");
  std::vector<KeyCols> st, in;
  std::vector<const int*> id;
  std::vector<Rcpp::IntegerVector> keep;   // (keeps coerced id vectors alive)
  st.reserve(streams.size()); in.reserve(inputs.size());
  for (int s = 0; s < streams.size(); ++s) {
    st.emplace_back(Rcpp::as<Rcpp::List>(streams[s]));
    keep.push_back(Rcpp::as<Rcpp::IntegerVector>(ids[s]));
    if (keep.back().size() != st.back().n) Rcpp::stop("stream %d: id length differs from key length", s + 1);
  }
  for (auto& v : keep) id.push_back(v.begin());
  for (int i = 0; i < inputs.size(); ++i) in.emplace_back(Rcpp::as<Rcpp::List>(inputs[i]));
  std::vector<const KeyCols*> sets;
  for (auto& x : st) sets.push_back(&x);
  for (auto& x : in) sets.push_back(&x);
  if (sets.empty()) Rcpp::stop("no key tables");
  for (auto* x : sets) check_same_width(*sets[0], *x);
  Packer pk = make_packer(sets);
#define OFST_CALL(W) locate_impl<W>(st, id, in, pk, n_old, threads)
  OFST_DISPATCH_W(pk, OFST_CALL)
#undef OFST_CALL
}
