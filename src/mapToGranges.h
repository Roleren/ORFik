
#ifndef PKG_mapToGRanges_H
#define PKG_mapToGRanges_H

Rcpp::S4 map_to_GRangesC(Rcpp::IntegerMatrix ORFdef,
      Rcpp::IntegerMatrix grangesObj,
      Rcpp::CharacterMatrix grangesObjStrings,
      std::string transcriptName);

#endif
