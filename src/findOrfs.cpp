// Find many orfs faster by putting all the data into Cpp directly
// Cpp format: Webkit: 80 character max line

#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

#include <Rcpp.h>
#include "findORFsHelpers.h"

using vi = std::vector<int>;
using string = std::string;

using namespace Rcpp;

// [[Rcpp::export]]
List orfs_as_List(
    CharacterVector fastaSeqs,
    std::string startCodon,
    std::string stopCodon,
    int minimumLength)
{
  vi result_index; // index of input that output belongs to
  vi all_orfs;

  for (int i = 0; i < fastaSeqs.size(); i++) {
    std::string fastaSeq = static_cast<std::string>(fastaSeqs[i]);
    vi ORFdef = orfs_as_vector(fastaSeq, startCodon, stopCodon, minimumLength);

    all_orfs.insert(all_orfs.end(), ORFdef.begin(), ORFdef.end());
    result_index.insert(result_index.end(), ORFdef.size() / 2, i + 1);
  }
  // all_orfs is an interlaced vector. We de-interlace it into two vectors.
  std::vector<vi> result_value(2);
  result_value[0].resize(all_orfs.size() / 2);
  result_value[1].resize(all_orfs.size() / 2);
  for (size_t i = 0; i < all_orfs.size() / 2; i++) {
    result_value[0][i] = all_orfs[2 * i];
    result_value[1][i] = all_orfs[2 * i + 1];
  }
  // then return to R, e.g., allowing many orf per input
  return List::create(
    Named("index") = wrap(result_index),
    Named("orf") = wrap(result_value));
}
