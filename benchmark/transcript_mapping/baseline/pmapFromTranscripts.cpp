// Transcript coordinate mapping with bounded exon searches.
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
List pmapFromTranscriptsCPP(const IntegerVector& xStart, const IntegerVector& xEnd,
                         const IntegerVector& transcriptStart,
                         const IntegerVector& transcriptEnd,
                         const IntegerVector& indices, const char& direction,
                         const bool removeEmpty) {
  const R_xlen_t nq = xStart.size(), ne = transcriptStart.size();
  if (xEnd.size() != nq || transcriptEnd.size() != ne || indices.size() != ne)
    stop("Inconsistent mapping vector lengths");
  if (direction != '+' && direction != '-') stop("Invalid mapping direction");
  if (ne > std::numeric_limits<int>::max()) stop("Too many transcript exons");
  if (!nq) {
    if (ne) stop("Query and transcript group counts differ");
    return List::create(_["index"] = IntegerVector(),
      _["ranges"] = List::create(IntegerVector(), IntegerVector()));
  }
  IntegerVector cumulative(ne), first(nq), mappedFirst(nq), mappedLast(nq);
  const int* xStart_ptr = xStart.begin();
  const int* xEnd_ptr = xEnd.begin();
  const int* transcriptStart_ptr = transcriptStart.begin();
  const int* transcriptEnd_ptr = transcriptEnd.begin();
  const int* indices_ptr = indices.begin();
  int* cumulative_ptr = cumulative.begin();
  int* first_ptr = first.begin();
  int* mappedFirst_ptr = mappedFirst.begin();
  int* mappedLast_ptr = mappedLast.begin();
  R_xlen_t q = 0, outputCount = removeEmpty ? 0 : ne;
  for (R_xlen_t begin = 0; begin < ne;) {
    if (q == nq) stop("Query and transcript group counts differ");
    if (indices_ptr[begin] == NA_INTEGER || indices_ptr[begin] < 1) stop("Invalid transcript indices");
    R_xlen_t finish = begin + 1;
    while (finish < ne && indices_ptr[finish] == indices_ptr[begin]) ++finish;
    int64_t total = 0;
    for (R_xlen_t j = begin; j < finish; ++j) {
      if (transcriptStart_ptr[j] == NA_INTEGER || transcriptEnd_ptr[j] == NA_INTEGER ||
          transcriptStart_ptr[j] < 1 || transcriptEnd_ptr[j] < transcriptStart_ptr[j])
        stop("Invalid transcript exon coordinates");
      total += static_cast<int64_t>(transcriptEnd_ptr[j]) - transcriptStart_ptr[j] + 1;
      if (total > std::numeric_limits<int>::max()) stop("Transcript width exceeds integer range");
      cumulative_ptr[j] = total;
    }
    if (xStart_ptr[q] == NA_INTEGER || xEnd_ptr[q] == NA_INTEGER ||
        xStart_ptr[q] < 1 || xEnd_ptr[q] < xStart_ptr[q] || xEnd_ptr[q] > total)
      stop("Transcript query coordinates are out of bounds");
    first_ptr[q] = begin;
    mappedFirst_ptr[q] = std::lower_bound(cumulative_ptr + begin,
      cumulative_ptr + finish, xStart_ptr[q]) - cumulative_ptr;
    mappedLast_ptr[q] = std::lower_bound(cumulative_ptr + begin,
      cumulative_ptr + finish, xEnd_ptr[q]) - cumulative_ptr;
    if (removeEmpty) outputCount += mappedLast_ptr[q] - mappedFirst_ptr[q] + 1;
    begin = finish; ++q;
  }
  if (q != nq) stop("Query and transcript group counts differ");
  IntegerVector starts(outputCount), ends(outputCount);
  IntegerVector resultIndex = removeEmpty ? IntegerVector(outputCount) : plain_indices(indices);
  int* starts_ptr = starts.begin();
  int* ends_ptr = ends.begin();
  int* resultIndex_ptr = resultIndex.begin();
  R_xlen_t out = 0;
  for (R_xlen_t i = 0; i < nq; ++i) {
    if ((i & 65535) == 0) checkUserInterrupt();
    for (int j = mappedFirst_ptr[i]; j <= mappedLast_ptr[i]; ++j) {
      const int previous = j == first_ptr[i] ? 0 : cumulative_ptr[j - 1];
      const int s = std::max(xStart_ptr[i], previous + 1) - previous - 1;
      const int e = std::min(xEnd_ptr[i], cumulative_ptr[j]) - previous - 1;
      const R_xlen_t row = removeEmpty ? out++ : j;
      if (direction == '+') {
        starts_ptr[row] = transcriptStart_ptr[j] + s;
        ends_ptr[row] = transcriptStart_ptr[j] + e;
      } else {
        starts_ptr[row] = transcriptEnd_ptr[j] - e;
        ends_ptr[row] = transcriptEnd_ptr[j] - s;
      }
      if (removeEmpty) resultIndex_ptr[row] = indices_ptr[j];
    }
  }
  return List::create(_["index"] = resultIndex, _["ranges"] = List::create(starts, ends));
}
