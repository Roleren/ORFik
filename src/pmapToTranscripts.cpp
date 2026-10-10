#include "transcript_mapping_R_helpers.h"
using namespace Rcpp;
namespace mapping=orfik::mapping;
namespace binding=orfik::r_mapping;
// [[Rcpp::export]]
List pmapToTranscriptsCPP(const IntegerVector& xStart,const IntegerVector& xEnd,
                         const IntegerVector& transcriptStart,const IntegerVector& transcriptEnd,
                         const IntegerVector& indices,const char& direction,const IntegerVector& nExonCumSum) {
  IntegerVector starts(xStart.size()),ends(xStart.size());
  mapping::map_to(binding::view(xStart),binding::view(xEnd),binding::view(transcriptStart),binding::view(transcriptEnd),
                  binding::view(indices),direction,binding::view(nExonCumSum),starts.begin(),ends.begin(),checkUserInterrupt);
  return binding::result(binding::plain_indices(indices),starts,ends);
}
