// Find many orfs faster by putting all the data into Cpp directly
//NB! Slower for few, because of memory overhead
// Cpp format: Webkit: 80 character max line

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>

#include <Rcpp.h>
#include "findORFs.h"

using vi = std::vector<int>;
using string = std::string;

using namespace Rcpp;

// Input is fasta sequences and fiveUTRs as grangesList
// Output is GRanges object with ranges mapped to genomic coordinates
// Instead of passing back and forth from r, do it all in C

// [[Rcpp::export]]
List get_all_ORFs_as_List(
    CharacterVector fastaSeqs,
    std::string startCodon,
    std::string stopCodon,
    bool longestORF,
    int minimumLength)
{
  size_t startSize = 2000000;
  vi result_index(startSize); // index of input that output belongs to
  size_t index = 0;
  std::vector<vi> result_value(2, vi(startSize)); // output values
  // std::cout << "size is: "<< fastaSeqs.size() << std::endl;
  // create result in C++
  for (size_t i = 0; i < fastaSeqs.size(); i++) {

    std::string fastaSeq = static_cast<std::string>(fastaSeqs[i]);
    std::vector<int> ORFdef = get_all_orfs_as_vector(fastaSeq, startCodon,
                                                     stopCodon,longestORF,
                                                     minimumLength);
    //check if resize is needed
    if((index + (ORFdef.size()/2)) > startSize){
      startSize = startSize*2;
      result_index.resize(startSize);
      result_value[0].resize(startSize);
      result_value[1].resize(startSize);
    }
    for (int j = 0; j < (ORFdef.size() / 2); j++) {
      result_value[0][index + j] = ORFdef[j * 2]; // set starts
      result_value[1][index + j] = ORFdef[(j * 2) + 1]; // set ends
      result_index[index + j] = i + 1;
    }
    index += (ORFdef.size() / 2); // Increase counter
  }

  result_index.resize(index);
  result_value[0].resize(index);
  result_value[1].resize(index);
  // then return to R, e.g., allowing many orf per input
  return List::create(
    Named("index") = wrap(result_index),
    Named("orf") = wrap(result_value));
}
