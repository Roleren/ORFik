// Mapping to transcript from genomic coordinates
// Cpp format: Webkit: 80 character max line

#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

#include <Rcpp.h>

using vi = std::vector<int>;
using string = std::string;

using namespace Rcpp;

//Forward strand calculation
void pmapToPositive(vi &all_orfs, const int &xSize, int &counter,
                  const vi &xWidths, const vi &nExonCumSum,
                  const vi &transcriptStart, const vi &xStart,
                  const vi &txWidthsCumSum, const vi &txWidths,
                  const vi &transcriptEnd, const vi &xEnd){
  int currentWidth = 0;
  int currentStart = 0;
  int dif = 0;
  int txSize = transcriptStart.size();
  bool firstExon, firstAssigned = true;

  for (auto i = 0; i < xSize; i++) { // per range to map
    currentWidth = xWidths[i];
    counter = nExonCumSum[i];
    firstAssigned = firstExon = true; //reset
    while(counter <= txSize) {
      currentStart = transcriptStart[counter];
      if (xStart[i] >= currentStart && transcriptEnd[counter] >= xStart[i]) { // if not jump to next exon
        if (firstExon) {
          dif = xStart[i] + 1;
          all_orfs[2 * i] = dif - currentStart;
        } else if (firstAssigned) { // set starts
          dif = xStart[i] + txWidthsCumSum[counter - 1] + 1;
          all_orfs[2 * i] = dif - currentStart;
        }
        firstAssigned = false;


        if (!firstAssigned && (xEnd[i] <= transcriptEnd[counter])) { //set ends
          if (firstExon) {
            all_orfs[2 * i + 1] = xEnd[i] - currentStart + 1;
          } else {
            all_orfs[2 * i + 1] = xEnd[i] -
              currentStart + txWidthsCumSum[counter - 1] + 1;
          }
          break;
        }
      }
      counter++;
      firstExon = false;
    }
  }
  return;
}

//Reverse strand calculation
void pmapToNegative(vi &all_orfs, const int &xSize, int &counter,
                  const vi &xWidths, const vi &nExonCumSum,
                  const vi &transcriptStart, const vi &xStart,
                  const vi &txWidthsCumSum, const vi &txWidths,
                  const vi &transcriptEnd, const vi &xEnd){
  int currentWidth = 0;
  int currentStart = 0;
  int dif = 0;
  bool firstExon, firstAssigned = true;
  int txSize = transcriptStart.size();
  for (auto i = 0; i < xSize; i++) { // per range to map
    currentWidth = xWidths[i];
    counter = nExonCumSum[i];
    firstAssigned = firstExon = true; //reset
    while(counter <= txSize) {
      currentStart = transcriptEnd[counter];
      if (xStart[i] >= transcriptStart[counter] && currentStart >= xStart[i]) { // if not jump to next exon
        if (firstExon) { //set starts (ends for - strand)
          dif = xStart[i];
          all_orfs[2 * i + 1] = currentStart - dif  + 1;
        } else if (firstAssigned) {
          dif = xStart[i] ;
          all_orfs[2 * i + 1] = currentStart - dif + txWidthsCumSum[counter - 1] + 1;
        }
        firstAssigned = false;


        if (!firstAssigned && (xEnd[i] <= currentStart)) { //set ends (starts for - strand)
          if (firstExon) {
            all_orfs[2 * i] = currentStart - xEnd[i] + 1;
          } else {
            all_orfs[2 * i] = currentStart - xEnd[i]
               + txWidthsCumSum[counter - 1] + 1;
          }
          break;
        }
      }
      counter++;
      firstExon = false;
    }
  }
  return;
}

// xStart, xEnd, the ranges to be mapped coordinates
// transcriptStart, transcriptEnd, the transcripts to map against
// indices are the number of exons per range, (1,1,1,2,2) is 3 exons then 2,
// in transcripts
// [[Rcpp::export]]
List pmapToTranscriptsCPP(const std::vector<int> &xStart,
                            const std::vector<int> &xEnd,
                            const std::vector<int> &transcriptStart,
                            const std::vector<int> &transcriptEnd,
                            const std::vector<int> &indices,
                            const char &direction,
                            const std::vector<int> &nExonCumSum)
{
  vi result_index; // index of input that output belongs to
  size_t xSize = xStart.size();
  size_t indexSize = transcriptStart.size();


  vi xWidths(xSize);
  for(auto i = 0;i < xSize; i++){ // Width per exon in x
    xWidths[i] = xEnd[i] - xStart[i] + 1;
  }
  vi txWidths(indexSize);
  for(auto i = 0; i < indexSize; i++){ // Width per exon in tx
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

  int counter = 0;
  vi all_orfs(2 * xSize);

  if(direction == '+'){
    pmapToPositive(all_orfs, xSize, counter, xWidths, nExonCumSum,
                 transcriptStart, xStart, txWidthsCumSum,
                 txWidths, transcriptEnd, xEnd);
  } else {
    pmapToNegative(all_orfs, xSize, counter, xWidths, nExonCumSum,
                 transcriptStart, xStart, txWidthsCumSum,
                 txWidths, transcriptEnd, xEnd);
  }
  // all_orfs is an interlaced vector. We de-interlace it into two vectors.
  std::vector<vi> result_value(2);
  result_value[0].resize(xSize);
  result_value[1].resize(xSize);

  for (size_t i = 0; i < xSize; i++) {
    result_value[0][i] = all_orfs[2 * i];
    result_value[1][i] = all_orfs[2 * i + 1];
  }
  // then return to R, e.g., allowing many orf per input
  return List::create(
    Named("index") = wrap(indices),
    Named("ranges") = wrap(result_value));
}

