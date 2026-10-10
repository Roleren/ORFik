#ifndef ORFIK_TRANSCRIPT_MAPPING_R_HELPERS_H
#define ORFIK_TRANSCRIPT_MAPPING_R_HELPERS_H
#include <Rcpp.h>
#include "../inst/include/ORFik/transcript_mapping_core.h"
namespace orfik { namespace r_mapping {
inline mapping::IntView view(const Rcpp::IntegerVector& x) {return {x.begin(),static_cast<std::size_t>(x.size())};}
inline Rcpp::IntegerVector plain_indices(const Rcpp::IntegerVector& indices) {
  if(ATTRIB(indices)==R_NilValue)return indices;
  return Rcpp::IntegerVector(indices.begin(),indices.end());
}
inline Rcpp::List result(const Rcpp::IntegerVector& index,const Rcpp::IntegerVector& starts,const Rcpp::IntegerVector& ends) {
  return Rcpp::List::create(Rcpp::_["index"]=index,Rcpp::_["ranges"]=Rcpp::List::create(starts,ends));
}
}}
#endif
