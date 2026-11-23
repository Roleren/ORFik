// bigWig writer for coverage given as Rle runs (one track per chromosome).
//
// Writes the UCSC bigWig format (version 4): header, chromosome B+ tree, bedGraph data blocks
// of up to 1024 items, the R-tree index over the blocks, 7 zoom levels (256 bp * 4^k) with
// their R-tree indexes and the total summary. Data and zoom blocks are deflated in parallel
// (OpenMP, zlib); everything else is written serially. No R API calls inside parallel regions.
// Reference: Kent et al. 2010, Bioinformatics 26:2204, and the UCSC kent source (bbiWrite.c,
// bPlusTree.c, cirTree.c), which rtracklayer embeds.
#include <Rcpp.h>
#include <zlib.h>
#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

// File output that tracks its own offset (no ftello/fseeko, which are not portable for > 2 GB).
struct Out {
  FILE* f = nullptr;
  uint64_t pos = 0;
  bool err = false;
  void write(const void* d, size_t n) {
    if (n && fwrite(d, 1, n, f) != n) err = true;
    pos += n;
  }
  template <class T> void put(T x) { write(&x, sizeof(T)); }
  void zeros(size_t n) { static const char z[256] = {0}; while (n) { size_t k = std::min<size_t>(n, 256); write(z, k); n -= k; } }
};

template <class T> inline void put(std::string& b, T x) { b.append(reinterpret_cast<const char*>(&x), sizeof(T)); }

// zlib deflate with a reusable stream: deflateReset() instead of a new deflateInit() per block
// (compress2() allocates and clears ~270 KB of state per call, which dominates for the many
// small blocks of fragmented assemblies). Output bytes equal compress2() at the same level.
// One Deflater per thread; never calls R.
struct Deflater {
  z_stream zs;
  int level = -1;
  bool ready = false;
  Deflater() = default;
  Deflater(const Deflater&) = delete;
  Deflater& operator=(const Deflater&) = delete;
  ~Deflater() { if (ready) deflateEnd(&zs); }
  bool run(const char* src, size_t n, int lvl, std::string& out) {
    if (!ready || lvl != level) {
      if (ready) deflateEnd(&zs);
      memset(&zs, 0, sizeof(zs));
      ready = deflateInit(&zs, lvl) == Z_OK;
      if (!ready) return false;
      level = lvl;
    } else if (deflateReset(&zs) != Z_OK) return false;
    const uLong bound = deflateBound(&zs, static_cast<uLong>(n));
    out.resize(bound);
    zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(src));
    zs.avail_in = static_cast<uInt>(n);
    zs.next_out = reinterpret_cast<Bytef*>(&out[0]);
    zs.avail_out = static_cast<uInt>(bound);
    if (deflate(&zs, Z_FINISH) != Z_STREAM_END) return false;
    out.resize(zs.total_out);
    return true;
  }
};

inline int thread_id() {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

struct IdxItem { uint32_t c0, s, c1, e; uint64_t off, size; };
struct ZoomRec { uint32_t chrom, start, end, valid; float mn, mx, sum, sumsq; };
static_assert(sizeof(ZoomRec) == 32, "zoom record must be 32 bytes");

// Chromosome B+ tree as UCSC bptFileBulkIndexToOpenFile: full blocks (zero padded), any
// number of levels, so there is no limit on the number of sequences.
void write_chrom_tree(Out& o, const std::vector<std::string>& names, const std::vector<uint32_t>& sizes,
                      uint32_t blockSize, uint32_t keySize) {
  const uint64_t n = names.size();
  const uint32_t valSize = 8;
  o.put<uint32_t>(0x78CA8C91u); o.put<uint32_t>(blockSize); o.put<uint32_t>(keySize); o.put<uint32_t>(valSize);
  o.put<uint64_t>(n); o.put<uint64_t>(0);
  int levels = 1;
  for (uint64_t c = n; c > blockSize; c = (c + blockSize - 1) / blockSize) ++levels;
  uint64_t indexOffset = o.pos;
  const uint64_t bytesIndexBlock = 4 + static_cast<uint64_t>(blockSize) * (keySize + 8);
  const uint64_t bytesLeafBlock = 4 + static_cast<uint64_t>(blockSize) * (keySize + valSize);
  std::string key(keySize, '\0');
  auto putKey = [&](uint64_t i) {
    std::fill(key.begin(), key.end(), '\0');
    memcpy(&key[0], names[i].data(), names[i].size());
    o.write(key.data(), keySize);
  };
  for (int level = levels - 1; level > 0; --level) {
    uint64_t slotSizePer = 1;
    for (int k = 0; k < level; ++k) slotSizePer *= blockSize;
    const uint64_t nodeSizePer = slotSizePer * blockSize;
    const uint64_t nodeCount = (n + nodeSizePer - 1) / nodeSizePer;
    const uint64_t endLevel = indexOffset + nodeCount * bytesIndexBlock;
    uint64_t nextChild = endLevel;
    const uint64_t nextBlock = level == 1 ? bytesLeafBlock : bytesIndexBlock;
    for (uint64_t i = 0; i < n; i += nodeSizePer) {
      uint64_t countOne = (n - i + slotSizePer - 1) / slotSizePer;
      if (countOne > blockSize) countOne = blockSize;
      o.put<uint8_t>(0); o.put<uint8_t>(0); o.put<uint16_t>(static_cast<uint16_t>(countOne));
      for (uint64_t j = 0; j < countOne; ++j) {
        putKey(i + j * slotSizePer);
        o.put<uint64_t>(nextChild);
        nextChild += nextBlock;
      }
      o.zeros(static_cast<size_t>((blockSize - countOne) * (keySize + 8)));
    }
    indexOffset = o.pos;
  }
  for (uint64_t i = 0; i < n; i += blockSize) {
    const uint64_t countOne = std::min<uint64_t>(blockSize, n - i);
    o.put<uint8_t>(1); o.put<uint8_t>(0); o.put<uint16_t>(static_cast<uint16_t>(countOne));
    for (uint64_t j = 0; j < countOne; ++j) {
      putKey(i + j);
      o.put<uint32_t>(static_cast<uint32_t>(i + j)); o.put<uint32_t>(sizes[i + j]);
    }
    o.zeros(static_cast<size_t>((blockSize - countOne) * (keySize + valSize)));
  }
}

// R-tree (cirTree) over blocks sorted by (chrom, start). Nodes are written root first; only
// the last node of a level is partial, so child offsets follow from full node sizes.
void write_rtree(Out& o, const std::vector<IdxItem>& items, uint64_t endOffset, uint32_t itemsPerSlot) {
  const uint32_t B = 256;
  std::vector<std::vector<IdxItem> > lv(1, items);
  while (lv.back().size() > B) {
    const std::vector<IdxItem>& c = lv.back();
    std::vector<IdxItem> up;
    for (size_t j = 0; j < c.size(); j += B) {
      IdxItem x = c[j];
      for (size_t i = j; i < std::min(c.size(), j + B); ++i)
        if (c[i].c1 > x.c1 || (c[i].c1 == x.c1 && c[i].e > x.e)) { x.c1 = c[i].c1; x.e = c[i].e; }
      up.push_back(x);
    }
    lv.push_back(up);
  }
  IdxItem all = {0, 0, 0, 0, 0, 0};
  if (!items.empty()) {
    all = items[0];
    for (const IdxItem& x : lv.back())
      if (x.c1 > all.c1 || (x.c1 == all.c1 && x.e > all.e)) { all.c1 = x.c1; all.e = x.e; }
  }
  o.put<uint32_t>(0x2468ACE0u); o.put<uint32_t>(B); o.put<uint64_t>(items.size());
  o.put(all.c0); o.put(all.s); o.put(all.c1); o.put(all.e); o.put<uint64_t>(endOffset);
  o.put<uint32_t>(itemsPerSlot); o.put<uint32_t>(0);
  const int K = static_cast<int>(lv.size()) - 1;
  std::vector<uint64_t> off(K + 1);
  uint64_t q = o.pos;
  for (int k = K; k >= 0; --k) {
    off[k] = q;
    const uint64_t nodes = (lv[k].size() + B - 1) / B;
    q += std::max<uint64_t>(nodes, 1) * 4 + lv[k].size() * (k == 0 ? 32 : 24);
  }
  for (int k = K; k >= 0; --k) {
    const std::vector<IdxItem>& c = lv[k];
    size_t j = 0;
    do {
      const size_t cnt = std::min<size_t>(B, c.size() - j);
      o.put<uint8_t>(k == 0); o.put<uint8_t>(0); o.put<uint16_t>(static_cast<uint16_t>(cnt));
      for (size_t i = j; i < j + cnt; ++i) {
        o.put(c[i].c0); o.put(c[i].s); o.put(c[i].c1); o.put(c[i].e);
        if (k == 0) { o.put(c[i].off); o.put(c[i].size); }
        else o.put<uint64_t>(off[k - 1] + static_cast<uint64_t>(i) * (4 + static_cast<uint64_t>(B) * (k - 1 == 0 ? 32 : 24)));
      }
      j += cnt;
    } while (j < c.size());
  }
}


// One sequence handed to the writer: raw pointers into R vectors that the R side keeps alive
// (the writer holds a reference) until the file is closed.
struct ChromJob { int c; R_xlen_t m; const int* lp; const int* vi; const double* vd; };

// bigWig writer that takes one sequence at a time, in any order (bw_stream_*), optionally on a
// background thread so that R continues while blocks are deflated. Data blocks are written in
// the order sequences arrive; the R-tree index is sorted by sequence id, and zoom records are
// kept per sequence and written at close in id order. Sequences added in id order give the same
// bytes as one bw_write_runs() call. The background thread never calls R.
class BwWriter {
public:
  BwWriter(const std::string& path_, const std::vector<std::string>& names_, const std::vector<uint32_t>& sizes_,
           int threads_, int level_, bool zero_runs_, uint32_t blockSize_, bool async_)
    : path(path_), names(names_), sizes(sizes_), threads(std::max(1, threads_)), level(level_),
      zero_runs(zero_runs_), async(async_), pool(std::max(1, threads_)), zrec(names_.size()),
      added(names_.size(), 0) {
    batch = std::max<size_t>(64, static_cast<size_t>(threads) * 64);
    size_t keySize = 1;
    for (const std::string& s : names) keySize = std::max(keySize, s.size());
    o.f = fopen(path.c_str(), "wb");
    if (!o.f) throw std::runtime_error("cannot open " + path + " for writing");
    o.zeros(headerSize());
    chromTreeOffset = o.pos;
    write_chrom_tree(o, names, sizes, blockSize_, static_cast<uint32_t>(keySize));
    fullDataOffset = o.pos;
    o.put<uint64_t>(0);   // number of data blocks, patched at close
    if (async) worker = std::thread([this]() { this->run(); });
  }
  ~BwWriter() { abort_(); }

  // Queue (async) or write (sync) one sequence. Called from the R thread only.
  void add(const ChromJob& j) {
    if (closed) throw std::runtime_error("bigWig writer is closed");
    if (j.c < 0 || j.c >= static_cast<int>(names.size()) || added[j.c]) throw std::runtime_error("sequence id invalid or added twice");
    added[j.c] = 1;
    if (!async) { process(j); if (!error.empty()) throw std::runtime_error(error); return; }
    { std::lock_guard<std::mutex> lk(mu); queue.push_back(j); }
    cv.notify_one();
  }
  std::string pending_error() { std::lock_guard<std::mutex> lk(mu); return error; }

  // Finish the file: wait for the background thread, then index, zoom levels and header.
  std::vector<double> close() {
    if (closed) throw std::runtime_error("bigWig writer is already closed");
    join_();
    closed = true;
    if (!error.empty()) { fclose(o.f); o.f = nullptr; throw std::runtime_error(error); }
    if (!flush_data()) { fclose(o.f); o.f = nullptr; throw std::runtime_error("zlib compression failed"); }
    // R-tree items sorted by (sequence id, start): blocks of one sequence are already in order
    std::stable_sort(idx.begin(), idx.end(), [](const IdxItem& a, const IdxItem& b) {
      return a.c0 < b.c0 || (a.c0 == b.c0 && a.s < b.s); });
    const uint32_t per = 1024, zper = 512;
    const int Z = static_cast<int>(red.size());
    // zoom records of all sequences in id order, deflated in blocks of zper in parallel
    std::vector<std::vector<ZoomRec> > all(Z);
    for (int z = 0; z < Z; ++z) {
      size_t tot = 0;
      for (size_t c = 0; c < zrec.size(); ++c) if (!zrec[c].empty()) tot += zrec[c][z].size();
      all[z].reserve(tot);
      for (size_t c = 0; c < zrec.size(); ++c) if (!zrec[c].empty()) {
        all[z].insert(all[z].end(), zrec[c][z].begin(), zrec[c][z].end());
        std::vector<ZoomRec>().swap(zrec[c][z]);
      }
    }
    std::vector<std::vector<IdxItem> > zidx(Z);
    std::vector<std::vector<std::string> > zblk(Z);
    std::vector<std::pair<int, size_t> > tasks;
    for (int z = 0; z < Z; ++z) {
      const size_t nb = (all[z].size() + zper - 1) / zper;
      zblk[z].resize(nb); zidx[z].resize(nb);
      for (size_t s = 0; s < nb; ++s) tasks.push_back(std::make_pair(z, s));
    }
    const R_xlen_t nt = static_cast<R_xlen_t>(tasks.size());
    bool err = false;
    #pragma omp parallel for schedule(dynamic, 16) num_threads(threads) if (nt > 1)
    for (R_xlen_t t = 0; t < nt; ++t) {
      const int z = tasks[t].first; const size_t s = tasks[t].second;
      const std::vector<ZoomRec>& rec = all[z];
      const size_t a = s * zper, e = std::min(rec.size(), a + zper);
      if (!pool[thread_id()].run(reinterpret_cast<const char*>(&rec[a]), (e - a) * sizeof(ZoomRec), level, zblk[z][s])) {
        #pragma omp critical
        err = true;
      }
      zidx[z][s] = {rec[a].chrom, rec[a].start, rec[e - 1].chrom, rec[e - 1].end, 0, 0};
    }
    if (err) { fclose(o.f); o.f = nullptr; throw std::runtime_error("zlib compression failed"); }
    const uint64_t fullIndexOffset = o.pos;
    write_rtree(o, idx, fullIndexOffset, per);
    std::vector<uint64_t> zdata(Z), zindex(Z);
    for (int z = 0; z < Z; ++z) {
      zdata[z] = o.pos;
      o.put<uint32_t>(static_cast<uint32_t>(all[z].size()));
      for (size_t s = 0; s < zblk[z].size(); ++s) {
        zidx[z][s].off = o.pos; zidx[z][s].size = zblk[z][s].size();
        o.write(zblk[z][s].data(), zblk[z][s].size());
        std::string().swap(zblk[z][s]);
      }
      zindex[z] = o.pos;
      write_rtree(o, zidx[z], zindex[z], zper);
    }
    o.put<uint32_t>(0x888FFC26u);   // trailing magic
    const double total = static_cast<double>(o.pos);
    if (fseek(o.f, 0, SEEK_SET) != 0) o.err = true;
    o.put<uint32_t>(0x888FFC26u); o.put<uint16_t>(4); o.put<uint16_t>(static_cast<uint16_t>(Z));
    o.put<uint64_t>(chromTreeOffset); o.put<uint64_t>(fullDataOffset); o.put<uint64_t>(fullIndexOffset);
    o.put<uint16_t>(0); o.put<uint16_t>(0); o.put<uint64_t>(0);
    o.put<uint64_t>(64 + 24 * Z);                                               // total summary offset
    o.put<uint32_t>(std::max<uint32_t>(24 + per * 12, zper * sizeof(ZoomRec)));  // uncompress buffer size
    o.put<uint64_t>(0);
    for (int z = 0; z < Z; ++z) { o.put<uint32_t>(red[z]); o.put<uint32_t>(0); o.put<uint64_t>(zdata[z]); o.put<uint64_t>(zindex[z]); }
    if (!covered) mn = mx = 0;
    o.put<uint64_t>(covered); o.put<double>(mn); o.put<double>(mx); o.put<double>(sum); o.put<double>(sumsq);
    if (fullDataOffset > 2000000000ULL || fseek(o.f, static_cast<long>(fullDataOffset), SEEK_SET) != 0) o.err = true;
    o.put<uint64_t>(nsect);
    const bool closeErr = fclose(o.f) != 0;
    o.f = nullptr;
    if (o.err || closeErr) throw std::runtime_error("error writing " + path + " (disk full?)");
    return {total, static_cast<double>(nitems), static_cast<double>(nsect)};
  }

  std::vector<Rcpp::RObject> keep;   // R vectors referenced by queued jobs (R thread only)

private:
  static size_t headerSize() { return 64 + 24 * 7 + 40; }   // header, zoom headers, total summary
  const std::vector<uint32_t> red = {256, 1024, 4096, 16384, 65536, 262144, 1048576};
  std::string path;
  std::vector<std::string> names;
  std::vector<uint32_t> sizes;
  int threads, level;
  bool zero_runs, async;
  Out o;
  uint64_t chromTreeOffset = 0, fullDataOffset = 0;
  std::vector<IdxItem> idx;
  uint64_t nsect = 0, nitems = 0, covered = 0;
  double mn = INFINITY, mx = -INFINITY, sum = 0, sumsq = 0;
  std::vector<Deflater> pool;                      // one zlib stream per thread
  std::vector<std::vector<std::vector<ZoomRec> > > zrec;   // [sequence][zoom level]
  std::vector<char> added;
  size_t batch = 64;
  std::vector<std::string> qraw; std::vector<IdxItem> qidx;
  bool closed = false;
  std::thread worker;
  std::mutex mu;
  std::condition_variable cv;
  std::deque<ChromJob> queue;
  bool stop = false;
  std::string error;

  void run() {
    for (;;) {
      ChromJob j;
      {
        std::unique_lock<std::mutex> lk(mu);
        cv.wait(lk, [this]() { return stop || !queue.empty(); });
        if (queue.empty()) return;   // stop requested and nothing left
        j = queue.front(); queue.pop_front();
        if (!error.empty()) continue;
      }
      process(j);
    }
  }
  void join_() {
    if (!worker.joinable()) return;
    { std::lock_guard<std::mutex> lk(mu); stop = true; }
    cv.notify_one();
    worker.join();
  }
  void abort_() {
    if (worker.joinable()) {
      { std::lock_guard<std::mutex> lk(mu); queue.clear(); stop = true; }
      cv.notify_one();
      worker.join();
    }
    if (o.f) { fclose(o.f); o.f = nullptr; }
  }
  void set_error(const std::string& e) { std::lock_guard<std::mutex> lk(mu); if (error.empty()) error = e; }

  bool flush_data() {
    const R_xlen_t nq = static_cast<R_xlen_t>(qraw.size());
    std::vector<std::string> comp(nq);
    bool err = false;
    #pragma omp parallel for schedule(dynamic, 4) num_threads(threads) if (nq > 1)
    for (R_xlen_t q = 0; q < nq; ++q) {
      if (!pool[thread_id()].run(qraw[q].data(), qraw[q].size(), level, comp[q])) {
        #pragma omp critical
        err = true;
      }
    }
    for (R_xlen_t q = 0; q < nq; ++q) {
      qidx[q].off = o.pos; qidx[q].size = comp[q].size();
      o.write(comp[q].data(), comp[q].size());
      idx.push_back(qidx[q]);
    }
    nsect += nq;
    qraw.clear(); qidx.clear();
    return !err;
  }

  // Data blocks, total summary and zoom records of one sequence (no R API).
  void process(const ChromJob& j) {
    try {
      const uint32_t per = 1024;
      const int Z = static_cast<int>(red.size());
      const int c = j.c;
      std::vector<uint32_t> is, ie; std::vector<float> iv;
      is.reserve(j.m); ie.reserve(j.m); iv.reserve(j.m);
      uint32_t pos = 0;
      for (R_xlen_t k = 0; k < j.m; ++k) {
        const int L = j.lp[k];
        if (L == 0) continue;
        const double x = j.vi ? static_cast<double>(j.vi[k]) : j.vd[k];
        if (x != 0 || zero_runs) { is.push_back(pos); ie.push_back(pos + L); iv.push_back(static_cast<float>(x)); }
        pos += L;
      }
      const size_t ni = is.size();
      if (!ni) return;
      nitems += ni;
      for (size_t a = 0; a < ni; a += per) {
        const size_t e = std::min(ni, a + per);
        std::string raw; raw.reserve(24 + (e - a) * 12);
        put<uint32_t>(raw, c); put<uint32_t>(raw, is[a]); put<uint32_t>(raw, ie[e - 1]);
        put<uint32_t>(raw, 0); put<uint32_t>(raw, 0);
        put<uint8_t>(raw, 1); put<uint8_t>(raw, 0); put<uint16_t>(raw, static_cast<uint16_t>(e - a));
        for (size_t i = a; i < e; ++i) { put(raw, is[i]); put(raw, ie[i]); put(raw, iv[i]); }
        qraw.push_back(std::move(raw));
        qidx.push_back({static_cast<uint32_t>(c), is[a], static_cast<uint32_t>(c), ie[e - 1], 0, 0});
        if (qraw.size() >= batch && !flush_data()) { set_error("zlib compression failed"); return; }
      }
      for (size_t i = 0; i < ni; ++i) {
        const double w = ie[i] - is[i], x = iv[i];
        covered += static_cast<uint64_t>(w); sum += x * w; sumsq += x * x * w;
        if (x < mn) mn = x;
        if (x > mx) mx = x;
      }
      std::vector<std::vector<ZoomRec> >& zr = zrec[c];
      zr.assign(Z, std::vector<ZoomRec>());
      #pragma omp parallel for schedule(dynamic, 1) num_threads(std::min(threads, Z)) if (ni > 100000)
      for (int z = 0; z < Z; ++z) {
        const uint32_t r = red[z];
        std::vector<ZoomRec>& out = zr[z];
        ZoomRec cur = {static_cast<uint32_t>(c), 0, 0, 0, 0, 0, 0, 0};
        uint32_t bin = UINT32_MAX;
        for (size_t i = 0; i < ni; ++i) {
          uint32_t a = is[i];
          const float x = iv[i];
          while (a < ie[i]) {
            const uint32_t b = a / r;
            const uint32_t e = static_cast<uint32_t>(std::min<uint64_t>(ie[i], static_cast<uint64_t>(b + 1) * r));
            const uint32_t w = e - a;
            if (b != bin) {
              if (bin != UINT32_MAX) out.push_back(cur);
              bin = b;
              cur = {static_cast<uint32_t>(c), a, e, w, x, x, x * w, x * x * w};
            } else {
              cur.end = e; cur.valid += w; cur.mn = std::min(cur.mn, x); cur.mx = std::max(cur.mx, x);
              cur.sum += x * w; cur.sumsq += x * x * w;
            }
            a = e;
          }
        }
        if (bin != UINT32_MAX) out.push_back(cur);
      }
    } catch (const std::exception& e) {
      set_error(std::string("bigWig writer: ") + e.what());
    } catch (...) {
      set_error("bigWig writer: unknown error");
    }
  }
};

// Checks names and sizes; returns the B+ tree block size.
uint32_t check_chroms(Rcpp::CharacterVector chrom, Rcpp::IntegerVector chrom_size, int chrom_block_size,
                      std::vector<std::string>& names, std::vector<uint32_t>& sizes) {
  const int nchr = chrom.size();
  if (chrom_size.size() != nchr) Rcpp::stop("chrom and chrom_size must have equal length");
  if (nchr < 1) Rcpp::stop("no sequences to write");
  if (chrom_block_size < 2 || chrom_block_size > 65535) Rcpp::stop("chrom_block_size must be 2..65535");
  const bool multiLevel = nchr > chrom_block_size;
  names.resize(nchr); sizes.resize(nchr);
  for (int i = 0; i < nchr; ++i) {
    names[i] = Rcpp::as<std::string>(chrom[i]);
    if (names[i].empty()) Rcpp::stop("empty sequence name");
    if (multiLevel && i && names[i] <= names[i - 1])
      Rcpp::stop("more than chrom_block_size sequences: chrom must be unique and sorted bytewise (C locale)");
    if (chrom_size[i] == NA_INTEGER || chrom_size[i] < 0) Rcpp::stop("invalid size of sequence %s", names[i]);
    sizes[i] = static_cast<uint32_t>(chrom_size[i]);
  }
  if (!multiLevel) {
    std::vector<std::string> sorted_names(names);
    std::sort(sorted_names.begin(), sorted_names.end());
    if (std::adjacent_find(sorted_names.begin(), sorted_names.end()) != sorted_names.end())
      Rcpp::stop("chrom must be unique");
  }
  return static_cast<uint32_t>(std::max(1, std::min(nchr, chrom_block_size)));
}

// Validates the runs of sequence c (R thread) and returns the job.
ChromJob make_job(int c, const std::string& name, uint32_t size, SEXP v, SEXP l) {
  if ((TYPEOF(v) != INTSXP && TYPEOF(v) != REALSXP) || TYPEOF(l) != INTSXP || XLENGTH(v) != XLENGTH(l))
    Rcpp::stop("runs of %s: values must be integer or double and lengths integer, of equal length", name);
  const R_xlen_t m = XLENGTH(l);
  const int* lp = INTEGER(l);
  const int* vi = TYPEOF(v) == INTSXP ? INTEGER(v) : nullptr;
  const double* vd = TYPEOF(v) == REALSXP ? REAL(v) : nullptr;
  double total = 0;
  for (R_xlen_t j = 0; j < m; ++j) {
    if (lp[j] == NA_INTEGER || lp[j] < 0) Rcpp::stop("invalid run length in %s", name);
    if (vi ? vi[j] == NA_INTEGER : std::isnan(vd[j])) Rcpp::stop("NA coverage value in %s", name);
    total += lp[j];
  }
  if (total != size) Rcpp::stop("runs of %s do not sum to the sequence length", name);
  return {c, m, lp, vi, vd};
}

} // namespace

// Write a bigWig file from Rle runs.
// chrom: unique sequence names; their order is the order of the chromosome ids (rtracklayer keeps
// the seqlevels order). Above chrom_block_size sequences the chromosome B+ tree has several
// levels, and then the names must be in bytewise (C locale) order for name lookups. chrom_size:
// their lengths; values / lengths: run values (integer or double, no NA) and run lengths per
// chromosome (NULL = no data). zero_runs = TRUE writes zero runs as items, like
// rtracklayer::export.bw(GRanges(<RleList>)), so import.bw() returns the same GRanges.
// level: zlib level (rtracklayer/UCSC use 6). threads: OpenMP threads for deflate.
// chrom_block_size: B+ tree block size; the default (65535, the largest) gives one leaf for all
// sequences, as rtracklayer writes it.
// Returns c(file size, data items, data blocks).
// [[Rcpp::export]]
Rcpp::NumericVector bw_write_runs(std::string path, Rcpp::CharacterVector chrom, Rcpp::IntegerVector chrom_size,
                                  Rcpp::List values, Rcpp::List lengths, int threads = 1, int level = 6,
                                  bool zero_runs = true, int chrom_block_size = 65535) {
  const int nchr = chrom.size();
  if (values.size() != nchr || lengths.size() != nchr)
    Rcpp::stop("chrom, chrom_size, values and lengths must have equal length");
  if (level < 0 || level > 9) Rcpp::stop("level must be 0..9");
  std::vector<std::string> names; std::vector<uint32_t> sizes;
  const uint32_t blockSize = check_chroms(chrom, chrom_size, chrom_block_size, names, sizes);
  std::vector<ChromJob> jobs;
  for (int i = 0; i < nchr; ++i)
    if (!Rf_isNull(values[i])) jobs.push_back(make_job(i, names[i], sizes[i], values[i], lengths[i]));
  std::vector<double> res;
  try {
    BwWriter w(path, names, sizes, threads, level, zero_runs, blockSize, false);
    for (const ChromJob& j : jobs) w.add(j);
    res = w.close();
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
  return Rcpp::NumericVector(res.begin(), res.end());
}

// Streaming form of bw_write_runs(): open the file with all sequence names and lengths, add the
// runs of one sequence at a time in any order (bw_stream_add, chrom_id 1-based), then close.
// async = TRUE deflates and writes on a background thread, so the R session continues while
// earlier sequences are compressed. Sequences never added are written without data.
// [[Rcpp::export]]
SEXP bw_stream_open(std::string path, Rcpp::CharacterVector chrom, Rcpp::IntegerVector chrom_size,
                    int threads = 1, int level = 6, bool zero_runs = true, int chrom_block_size = 65535,
                    bool async = true) {
  if (level < 0 || level > 9) Rcpp::stop("level must be 0..9");
  std::vector<std::string> names; std::vector<uint32_t> sizes;
  const uint32_t blockSize = check_chroms(chrom, chrom_size, chrom_block_size, names, sizes);
  BwWriter* w = nullptr;
  try {
    w = new BwWriter(path, names, sizes, threads, level, zero_runs, blockSize, async);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
  Rcpp::XPtr<BwWriter> ptr(w, true);
  ptr.attr("chrom") = chrom;
  ptr.attr("chrom_size") = chrom_size;
  return ptr;
}

// [[Rcpp::export]]
SEXP bw_stream_add(SEXP writer, int chrom_id, SEXP values, SEXP lengths) {
  Rcpp::XPtr<BwWriter> w(writer);
  if (!w.get()) Rcpp::stop("bigWig writer is closed");
  Rcpp::CharacterVector chrom = w.attr("chrom");
  Rcpp::IntegerVector chrom_size = w.attr("chrom_size");
  if (chrom_id < 1 || chrom_id > chrom.size()) Rcpp::stop("chrom_id out of range");
  const int c = chrom_id - 1;
  ChromJob j = make_job(c, Rcpp::as<std::string>(chrom[c]), static_cast<uint32_t>(chrom_size[c]), values, lengths);
  const std::string pe = w->pending_error();
  if (!pe.empty()) Rcpp::stop(pe);
  w->keep.push_back(Rcpp::RObject(values));
  w->keep.push_back(Rcpp::RObject(lengths));
  try {
    w->add(j);
  } catch (const std::exception& e) {
    Rcpp::stop(e.what());
  }
  return R_NilValue;
}

// [[Rcpp::export]]
Rcpp::NumericVector bw_stream_close(SEXP writer) {
  Rcpp::XPtr<BwWriter> w(writer);
  if (!w.get()) Rcpp::stop("bigWig writer is closed");
  std::vector<double> res;
  std::string err;
  try {
    res = w->close();
  } catch (const std::exception& e) {
    err = e.what();
  }
  w.release();   // deletes the writer (joins nothing: the thread has finished)
  if (!err.empty()) Rcpp::stop(err);
  return Rcpp::NumericVector(res.begin(), res.end());
}

// Stop a stream without finishing the file (background thread stopped, file closed); the R
// side removes the partial file.
// [[Rcpp::export]]
SEXP bw_stream_abort(SEXP writer) {
  Rcpp::XPtr<BwWriter> w(writer);
  if (w.get()) w.release();
  return R_NilValue;
}
