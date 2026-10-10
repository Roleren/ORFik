#include "transcript_mapping_R_helpers.h"
using namespace Rcpp;
namespace mapping=orfik::mapping;
namespace binding=orfik::r_mapping;
// [[Rcpp::export]]
List pmapFromTranscriptsCPP(const IntegerVector& xStart,const IntegerVector& xEnd,
                           const IntegerVector& transcriptStart,const IntegerVector& transcriptEnd,
                           const IntegerVector& indices,const char& direction,const bool removeEmpty) {
  mapping::require(xStart.size()==xEnd.size() && transcriptStart.size()==transcriptEnd.size() && transcriptStart.size()==indices.size(),"Inconsistent mapping vector lengths");
  mapping::require(direction=='+' || direction=='-',"Invalid mapping direction");
  auto partitions=mapping::partitions_from_indices(binding::view(indices));
  mapping::require(partitions.size()==static_cast<std::size_t>(xStart.size()),"Query and transcript group counts differ");
  std::vector<int> strands(partitions.size(),direction=='+'?1:2),groups(partitions.size());
  for(std::size_t i=0;i<groups.size();++i)groups[i]=i;
  mapping::Catalogue tx(binding::view(transcriptStart),binding::view(transcriptEnd),{partitions.data(),partitions.size()},{strands.data(),strands.size()});
  auto plan=mapping::plan_from(binding::view(xStart),binding::view(xEnd),{groups.data(),groups.size()},tx,removeEmpty);
  IntegerVector starts(plan.count),ends(plan.count),queryGroups(plan.count);
  mapping::fill_from(binding::view(xStart),binding::view(xEnd),tx,plan,removeEmpty,starts.begin(),ends.begin(),queryGroups.begin(),nullptr,checkUserInterrupt);
  IntegerVector resultIndex=removeEmpty?IntegerVector(plan.count):binding::plain_indices(indices);
  if(removeEmpty)for(std::size_t i=0;i<plan.count;++i)resultIndex[i]=indices[tx.boundaries[queryGroups[i]-1]];
  return binding::result(resultIndex,starts,ends);
}
