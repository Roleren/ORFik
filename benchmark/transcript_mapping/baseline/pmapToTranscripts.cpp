// Map genomic exon ranges to their paired transcripts.
#include <Rcpp.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <cstdint>
using namespace Rcpp;

static IntegerVector plain_indices(const IntegerVector& indices) {
  if (ATTRIB(indices) == R_NilValue) return indices;
  IntegerVector result(indices.size());
  std::copy(indices.begin(), indices.end(), result.begin());
  return result;
}


// [[Rcpp::export]]
List pmapToTranscriptsCPP(const IntegerVector& xStart, const IntegerVector& xEnd,
                       const IntegerVector& transcriptStart,
                       const IntegerVector& transcriptEnd,
                       const IntegerVector& indices, const char& direction,
                       const IntegerVector& nExonCumSum) {
  const R_xlen_t nq = xStart.size(), ne = transcriptStart.size();
  if (xEnd.size() != nq || nExonCumSum.size() != nq ||
      transcriptEnd.size() != ne || indices.size() != ne)
    stop("Inconsistent mapping vector lengths");
  if (direction != '+' && direction != '-') stop("Invalid mapping direction");
  if (ne > std::numeric_limits<int>::max()) stop("Too many transcript exons");
  IntegerVector starts(nq), ends(nq, -1);
  if (!nq) return List::create(_["index"] = plain_indices(indices),
    _["ranges"] = List::create(starts, ends));
  IntegerVector prefix(ne), groupEnd(ne);
  const int* xStart_ptr = xStart.begin();
  const int* xEnd_ptr = xEnd.begin();
  const int* transcriptStart_ptr = transcriptStart.begin();
  const int* transcriptEnd_ptr = transcriptEnd.begin();
  const int* indices_ptr = indices.begin();
  const int* nExonCumSum_ptr = nExonCumSum.begin();
  int* prefix_ptr = prefix.begin();
  int* groupEnd_ptr = groupEnd.begin();
  int* starts_ptr = starts.begin();
  int* ends_ptr = ends.begin();
  for (R_xlen_t begin = 0; begin < ne;) {
    if (indices_ptr[begin] == NA_INTEGER || indices_ptr[begin] < 1) stop("Invalid transcript indices");
    R_xlen_t last = begin + 1;
    while (last < ne && indices_ptr[last] == indices_ptr[begin]) ++last;
    int64_t total = 0;
    bool ordered = true;
    for (R_xlen_t j = begin; j < last; ++j) {
      if (transcriptStart_ptr[j] == NA_INTEGER || transcriptEnd_ptr[j] == NA_INTEGER ||
          transcriptStart_ptr[j] < 1 || transcriptEnd_ptr[j] < transcriptStart_ptr[j])
        stop("Invalid transcript exon coordinates");
      prefix_ptr[j] = static_cast<int>(total);
      total += static_cast<int64_t>(transcriptEnd_ptr[j]) - transcriptStart_ptr[j] + 1;
      if (total > std::numeric_limits<int>::max()) stop("Transcript width exceeds integer range");
      groupEnd_ptr[j] = last;
      if (j > begin && (direction == '+' ? transcriptStart_ptr[j] <= transcriptEnd_ptr[j - 1]
                                         : transcriptEnd_ptr[j] >= transcriptStart_ptr[j - 1]))
        ordered = false;
    }
    groupEnd_ptr[begin] = ordered ? last : -last;
    begin = last;
  }
  for (R_xlen_t i = 0; i < nq; ++i) {
    if ((i & 65535) == 0) checkUserInterrupt();
    const int base = nExonCumSum_ptr[i];
    if (base < 0 || base >= ne) continue;
    if (base && indices_ptr[base] == indices_ptr[base - 1])
      stop("Transcript offset must identify its first exon");
    if (xStart_ptr[i] == NA_INTEGER || xEnd_ptr[i] == NA_INTEGER ||
        xStart_ptr[i] < 1 || xEnd_ptr[i] < xStart_ptr[i]) continue;
    const int encodedEnd = groupEnd_ptr[base];
    const int last = std::abs(encodedEnd);
    int j = base;
    if (encodedEnd > 0 && last - base > 8) {
      if (direction == '+')
        j = std::lower_bound(transcriptEnd_ptr + base, transcriptEnd_ptr + last,
                             xEnd_ptr[i]) - transcriptEnd_ptr;
      else
        j = std::lower_bound(transcriptStart_ptr + base, transcriptStart_ptr + last,
                             xStart_ptr[i], std::greater<int>()) - transcriptStart_ptr;
    } else {
      while (j < last && !(transcriptStart_ptr[j] <= xStart_ptr[i] && transcriptEnd_ptr[j] >= xEnd_ptr[i])) ++j;
    }
    if (j == last || transcriptStart_ptr[j] > xStart_ptr[i] || transcriptEnd_ptr[j] < xEnd_ptr[i]) continue;
    if (direction == '+') {
      starts_ptr[i] = prefix_ptr[j] + (xStart_ptr[i] - transcriptStart_ptr[j] + 1);
      ends_ptr[i] = prefix_ptr[j] + (xEnd_ptr[i] - transcriptStart_ptr[j] + 1);
    } else {
      starts_ptr[i] = prefix_ptr[j] + (transcriptEnd_ptr[j] - xEnd_ptr[i] + 1);
      ends_ptr[i] = prefix_ptr[j] + (transcriptEnd_ptr[j] - xStart_ptr[i] + 1);
    }
  }
  return List::create(_["index"] = plain_indices(indices), _["ranges"] = List::create(starts, ends));
}
