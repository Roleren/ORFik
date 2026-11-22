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

inline double run_value(SEXP v, R_xlen_t j) {
  return TYPEOF(v) == INTSXP ? static_cast<double>(INTEGER(v)[j]) : REAL(v)[j];
}

} // namespace

// Write a bigWig file from Rle runs.
// chrom: unique sequence names; their order is the order of the chromosome ids (rtracklayer keeps
// the seqlevels order). Above chrom_block_size sequences the chromosome B+ tree has several
// levels, and then the names must be in bytewise (C locale) order for name lookups. chrom_size: their lengths; values / lengths: run values (integer or double, no NA) and run
// lengths per chromosome (NULL = no data). zero_runs = TRUE writes zero runs as items, like
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
  if (chrom_size.size() != nchr || values.size() != nchr || lengths.size() != nchr)
    Rcpp::stop("chrom, chrom_size, values and lengths must have equal length");
  if (nchr < 1) Rcpp::stop("no sequences to write");
  if (level < 0 || level > 9) Rcpp::stop("level must be 0..9");
  if (threads < 1) threads = 1;
  if (chrom_block_size < 2 || chrom_block_size > 65535) Rcpp::stop("chrom_block_size must be 2..65535");
  const bool multiLevel = nchr > chrom_block_size;
  const uint32_t per = 1024, zper = 512;
  const std::vector<uint32_t> red = {256, 1024, 4096, 16384, 65536, 262144, 1048576};
  const int Z = static_cast<int>(red.size());

  std::vector<std::string> names(nchr);
  std::vector<uint32_t> sizes(nchr);
  size_t keySize = 1;
  for (int i = 0; i < nchr; ++i) {
    names[i] = Rcpp::as<std::string>(chrom[i]);
    if (names[i].empty()) Rcpp::stop("empty sequence name");
    if (multiLevel && i && names[i] <= names[i - 1])
      Rcpp::stop("more than chrom_block_size sequences: chrom must be unique and sorted bytewise (C locale)");
    if (chrom_size[i] == NA_INTEGER || chrom_size[i] < 0) Rcpp::stop("invalid size of sequence %s", names[i]);
    sizes[i] = static_cast<uint32_t>(chrom_size[i]);
    keySize = std::max(keySize, names[i].size());
    SEXP v = values[i], l = lengths[i];
    if (Rf_isNull(v)) continue;
    if ((TYPEOF(v) != INTSXP && TYPEOF(v) != REALSXP) || TYPEOF(l) != INTSXP || XLENGTH(v) != XLENGTH(l))
      Rcpp::stop("runs of %s: values must be integer or double and lengths integer, of equal length", names[i]);
    double total = 0;
    const int* lp = INTEGER(l);
    for (R_xlen_t j = 0; j < XLENGTH(l); ++j) {
      if (lp[j] == NA_INTEGER || lp[j] < 0) Rcpp::stop("invalid run length in %s", names[i]);
      if (TYPEOF(v) == INTSXP ? INTEGER(v)[j] == NA_INTEGER : std::isnan(REAL(v)[j]))
        Rcpp::stop("NA coverage value in %s", names[i]);
      total += lp[j];
    }
    if (total != sizes[i]) Rcpp::stop("runs of %s do not sum to the sequence length", names[i]);
  }
  if (!multiLevel) {
    std::vector<std::string> sorted_names(names);
    std::sort(sorted_names.begin(), sorted_names.end());
    if (std::adjacent_find(sorted_names.begin(), sorted_names.end()) != sorted_names.end())
      Rcpp::stop("chrom must be unique");
  }
  const uint32_t blockSize = static_cast<uint32_t>(std::max(1, std::min(nchr, chrom_block_size)));

  Out o;
  o.f = fopen(path.c_str(), "wb");
  if (!o.f) Rcpp::stop("cannot open %s for writing", path);
  const size_t headerSize = 64 + 24 * Z + 40;   // header, zoom headers, total summary
  o.zeros(headerSize);
  const uint64_t chromTreeOffset = o.pos;
  write_chrom_tree(o, names, sizes, blockSize, static_cast<uint32_t>(keySize));
  const uint64_t fullDataOffset = o.pos;
  o.put<uint64_t>(0);   // number of data blocks, patched at the end

  std::vector<IdxItem> idx;
  uint64_t nsect = 0, nitems = 0, covered = 0;
  double mn = INFINITY, mx = -INFINITY, sum = 0, sumsq = 0;
  std::vector<std::vector<IdxItem> > zidx(Z);
  std::vector<std::vector<std::string> > zblk(Z);
  std::vector<uint32_t> zcount(Z, 0);
  std::vector<std::vector<ZoomRec> > pend(Z);   // zoom records not yet in a block
  std::vector<Deflater> pool(threads);          // one zlib stream per thread
  // queued raw data blocks, deflated in parallel and written in order by flush_data()
  const size_t batch = std::max<size_t>(64, static_cast<size_t>(threads) * 64);
  std::vector<std::string> qraw; std::vector<IdxItem> qidx;
  auto flush_data = [&]() {
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
  };
  // Deflate the pending zoom records in blocks of zper (all of them when final), in parallel.
  auto compress_zoom = [&](bool final) {
    std::vector<std::pair<int, size_t> > tasks;
    std::vector<size_t> base(Z), nrec(Z);
    for (int z = 0; z < Z; ++z) {
      nrec[z] = final ? pend[z].size() : pend[z].size() / zper * zper;
      const size_t nb = (nrec[z] + zper - 1) / zper;
      base[z] = zblk[z].size();
      zblk[z].resize(base[z] + nb); zidx[z].resize(base[z] + nb);
      for (size_t s = 0; s < nb; ++s) tasks.push_back(std::make_pair(z, s));
      zcount[z] += static_cast<uint32_t>(nrec[z]);
    }
    const R_xlen_t nt = static_cast<R_xlen_t>(tasks.size());
    bool err = false;
    #pragma omp parallel for schedule(dynamic, 16) num_threads(threads) if (nt > 1)
    for (R_xlen_t t = 0; t < nt; ++t) {
      const int z = tasks[t].first; const size_t s = tasks[t].second;
      const std::vector<ZoomRec>& rec = pend[z];
      const size_t a = s * zper, e = std::min(nrec[z], a + zper);
      if (!pool[thread_id()].run(reinterpret_cast<const char*>(&rec[a]), (e - a) * sizeof(ZoomRec), level, zblk[z][base[z] + s])) {
        #pragma omp critical
        err = true;
      }
      zidx[z][base[z] + s] = {rec[a].chrom, rec[a].start, rec[e - 1].chrom, rec[e - 1].end, 0, 0};
    }
    for (int z = 0; z < Z; ++z) pend[z].erase(pend[z].begin(), pend[z].begin() + nrec[z]);
    return !err;
  };

  for (int c = 0; c < nchr; ++c) {
    SEXP v = values[c];
    if (Rf_isNull(v)) continue;
    const int* lp = INTEGER(lengths[c]);
    const R_xlen_t m = XLENGTH(v);
    std::vector<uint32_t> is, ie; std::vector<float> iv;
    is.reserve(m); ie.reserve(m); iv.reserve(m);
    uint32_t pos = 0;
    for (R_xlen_t j = 0; j < m; ++j) {
      if (lp[j] == 0) continue;
      const double x = run_value(v, j);
      if (x != 0 || zero_runs) { is.push_back(pos); ie.push_back(pos + lp[j]); iv.push_back(static_cast<float>(x)); }
      pos += lp[j];
    }
    const size_t ni = is.size();
    if (!ni) continue;
    nitems += ni;
    // data blocks: raw blocks are queued (across sequences) and deflated in parallel batches
    for (size_t a = 0; a < ni; a += per) {
      const size_t e = std::min(ni, a + per);
      std::string raw; raw.reserve(24 + (e - a) * 12);
      put<uint32_t>(raw, c); put<uint32_t>(raw, is[a]); put<uint32_t>(raw, ie[e - 1]);
      put<uint32_t>(raw, 0); put<uint32_t>(raw, 0);
      put<uint8_t>(raw, 1); put<uint8_t>(raw, 0); put<uint16_t>(raw, static_cast<uint16_t>(e - a));
      for (size_t i = a; i < e; ++i) { put(raw, is[i]); put(raw, ie[i]); put(raw, iv[i]); }
      qraw.push_back(std::move(raw));
      qidx.push_back({static_cast<uint32_t>(c), is[a], static_cast<uint32_t>(c), ie[e - 1], 0, 0});
      if (qraw.size() >= batch && !flush_data()) { fclose(o.f); Rcpp::stop("zlib compression failed"); }
    }
    // total summary
    for (size_t i = 0; i < ni; ++i) {
      const double w = ie[i] - is[i], x = iv[i];
      covered += static_cast<uint64_t>(w); sum += x * w; sumsq += x * x * w;
      if (x < mn) mn = x;
      if (x > mx) mx = x;
    }
    // zoom records per level (items split at bin borders), levels in parallel; records are
    // appended to pend[z] and full blocks of zper records (which may span sequences) deflated
    #pragma omp parallel for schedule(dynamic, 1) num_threads(std::min(threads, Z)) if (ni > 100000)
    for (int z = 0; z < Z; ++z) {
      const uint32_t r = red[z];
      std::vector<ZoomRec>& out = pend[z];
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
    if (!compress_zoom(false)) { fclose(o.f); Rcpp::stop("zlib compression failed"); }
  }
  if (!flush_data() || !compress_zoom(true)) { fclose(o.f); Rcpp::stop("zlib compression failed"); }

  const uint64_t fullIndexOffset = o.pos;
  write_rtree(o, idx, fullIndexOffset, per);
  std::vector<uint64_t> zdata(Z), zindex(Z);
  for (int z = 0; z < Z; ++z) {
    zdata[z] = o.pos;
    o.put<uint32_t>(zcount[z]);
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
  // header, zoom headers and total summary
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
  // number of data blocks (fullDataOffset is small: right after the chromosome tree)
  if (fullDataOffset > 2000000000ULL || fseek(o.f, static_cast<long>(fullDataOffset), SEEK_SET) != 0) o.err = true;
  o.put<uint64_t>(nsect);
  const bool closeErr = fclose(o.f) != 0;
  if (o.err || closeErr) Rcpp::stop("error writing %s (disk full?)", path);
  return Rcpp::NumericVector::create(total, static_cast<double>(nitems), static_cast<double>(nsect));
}
