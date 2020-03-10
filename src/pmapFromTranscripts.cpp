// Mapping from transcript to genomic coordinates
// Cpp format: Webkit: 80 character max line

#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

#include <Rcpp.h>

using vi = std::vector<int>;
using string = std::string;

using namespace Rcpp;

// Integer vector of cumulative sums of exons
// (1,1,1,2,2,3,3,3) -> (0,3,5) [this is now indices of first of every]
vi nExonsCumSums(const int &xSize, const int &indexSize, const vi &indices){
  vi nExons(xSize);
  int ref = indices[0];
  int counter = 0;
  int current = 0;
  for(auto i = 0; i < indexSize; i++){
    if(indices[i] != ref){
      nExons[counter] = i - current;
      ref = indices[i];
      counter++;
      current = i;
    }
  }
  nExons[counter] = indexSize - current;

  vi nExonCumSum(xSize);
  nExonCumSum[0] = 0;
  if(xSize > 1) {
    for(auto i = 1;i < xSize; i++){
      nExonCumSum[i] = nExons[i-1] + nExonCumSum[i - 1];
    }
  }
  return(nExonCumSum);
}

//Forward strand calculation
void pmapPositive(vi &all_orfs, const int &xSize, int &counter,
                  const vi &xWidths, const vi &nExonCumSum,
                  const vi &transcriptStart, const vi &xStart,
                  const vi &txWidthsCumSum, const vi &txWidths,
                  const vi &transcriptEnd){
  int currentWidth = 0;
  int currentStart = 0;
  int dif = 0;
  bool firstExon, firstAssigned = true;
  for (auto i = 0; i < xSize; i++) { // per range to map
    currentWidth = xWidths[i];
    counter = nExonCumSum[i];
    currentStart = transcriptStart[counter];
    firstAssigned = firstExon = true; //reset
    while(true) {
      if (txWidthsCumSum[counter] >= xStart[i]) { // if not jump to next exon
        currentStart = transcriptStart[counter];
        if (firstExon) {
          dif = xStart[i] - 1;
          all_orfs[2 * counter] = currentStart + dif;
        } else if (firstAssigned) { // set starts
          dif = xStart[i] - txWidthsCumSum[counter - 1] - 1;
          all_orfs[2 * counter] = currentStart + dif;
        } else {
          all_orfs[2 * counter] = currentStart;
          dif = 0;
        }
        firstAssigned = false;

        if (currentWidth <= (txWidths[counter] - dif)) { //set ends
          all_orfs[2 * counter + 1] = all_orfs[2 * counter] + currentWidth - 1;
          break;
        } else {
          all_orfs[2 * counter + 1] = transcriptEnd[counter];
          currentWidth = currentWidth -
            (all_orfs[2 * counter + 1] - all_orfs[2 * counter] + 1);
        }
      }
      counter++;
      firstExon = false;
    }
  }
}

//Reverse strand calculation
void pmapNegative(vi &all_orfs, const int &xSize, int &counter,
                  const vi &xWidths, const vi &nExonCumSum,
                  const vi &transcriptStart, const vi &xStart,
                  const vi &txWidthsCumSum, const vi &txWidths,
                  const vi &transcriptEnd){
  int currentWidth = 0;
  int currentStart = 0;
  int dif = 0;
  bool firstExon, firstAssigned = true;
  for (auto i = 0; i < xSize; i++) {
    currentWidth = xWidths[i];
    counter = nExonCumSum[i];
    currentStart = transcriptEnd[counter];
    firstAssigned = firstExon = true; //reset
    while(true) {
      if (txWidthsCumSum[counter] >= xStart[i]) { // if not jump to next exon
        currentStart = transcriptEnd[counter];
        if (firstExon) {
          dif = xStart[i] - 1;
          all_orfs[2 * counter + 1] = currentStart - dif;
        } else if (firstAssigned) { // set starts
          dif = xStart[i] - txWidthsCumSum[counter - 1] - 1;
          all_orfs[2 * counter + 1] = currentStart - dif;
        } else {
          all_orfs[2 * counter + 1] = currentStart;
          dif = 0;
        }
        firstAssigned = false;

        if (currentWidth <= (txWidths[counter] - dif)) { //set ends
          all_orfs[2 * counter] = all_orfs[2 * counter + 1] - currentWidth + 1;
          break;
        } else {
          all_orfs[2 * counter] = transcriptStart[counter];
          currentWidth = currentWidth -
            (all_orfs[2 * counter + 1] - all_orfs[2 * counter] + 1);
        }
      }
      counter++;
      firstExon = false;
    }
  }
}

// xStart, xEnd, the ranges to be mapped coordinates
// transcriptStart, transcriptEnd, the transcripts to map against
// indices are the number of exons per range, (1,1,1,2,2) is 3 exons then 2.
// [[Rcpp::export]]
List pmapFromTranscriptsCPP(const std::vector<int> &xStart,
                            const std::vector<int> &xEnd,
                            const std::vector<int> &transcriptStart,
                            const std::vector<int> &transcriptEnd,
                            const std::vector<int> &indices,
                            const char &direction,
                            const bool removeEmpty)
{
  vi result_index; // index of input that output belongs to
  size_t xSize = xStart.size();
  size_t indexSize = transcriptStart.size();


  vi xWidths(xSize);
  for(auto i = 0;i < xSize; i++){
    xWidths[i] = xEnd[i] - xStart[i] + 1;
  }
  vi txWidths(indexSize);
  for(auto i = 0; i < indexSize; i++){
    txWidths[i] = transcriptEnd[i] - transcriptStart[i] + 1;
  }

  vi txWidthsCumSum(indexSize);
  txWidthsCumSum[0] = txWidths[0];
  for(auto i = 1; i < indexSize; i++){
    if(indices[i] == indices[i-1]){
      txWidthsCumSum[i] = txWidths[i] + txWidthsCumSum[i - 1];
    } else {
      txWidthsCumSum[i] = txWidths[i];
    }
  }
  vi nExonCumSum = nExonsCumSums(xSize, indexSize, indices);

  int counter = 0;
  vi all_orfs(2 * indexSize);
  if(direction == '+'){
    pmapPositive(all_orfs, xSize, counter, xWidths, nExonCumSum,
                 transcriptStart, xStart, txWidthsCumSum,
                 txWidths, transcriptEnd);
  } else {
    pmapNegative(all_orfs, xSize, counter, xWidths, nExonCumSum,
                 transcriptStart, xStart, txWidthsCumSum,
                 txWidths, transcriptEnd);
  }

  // all_orfs is an interlaced vector. We de-interlace it into two vectors.
  std::vector<vi> result_value(2);
  result_value[0].resize(indexSize);
  result_value[1].resize(indexSize);

  if(removeEmpty) {
    counter = 0;
    vi indicesNew(indexSize);
    for (size_t i = 0; i < indexSize; i++) {
      if(all_orfs[2 * i] != 0) {
        result_value[0][counter] = all_orfs[2 * i];
        result_value[1][counter] = all_orfs[2 * i + 1];
        indicesNew[counter] = indices[i];
        counter++;
      }
    }
    result_value[0].resize(counter);
    result_value[1].resize(counter);
    indicesNew.resize(counter);
    // then return to R, e.g., allowing many orf per input
    return List::create(
      Named("index") = wrap(indicesNew),
      Named("ranges") = wrap(result_value));

  } else {
    for (size_t i = 0; i < indexSize; i++) {
      result_value[0][i] = all_orfs[2 * i];
      result_value[1][i] = all_orfs[2 * i + 1];
    }
    // then return to R, e.g., allowing many orf per input
    return List::create(
      Named("index") = wrap(indices),
      Named("ranges") = wrap(result_value));
  }
}

