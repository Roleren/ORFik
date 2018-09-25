// Find many orfs faster by putting all the data into Cpp directly
// Cpp format: Webkit: 80 character max line

#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <locale>

#include <Rcpp.h>
#include "findORFsHelpers.h"

using vi = std::vector<int>;
using namespace Rcpp;

Function GRangesC("GRanges", Environment::namespace_env("GenomicRanges"));
Function IRangesC("IRanges", Environment::namespace_env("IRanges"));

// Rule for undefined behavior is to return (space)
char complement(char n)
{
  switch (n) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'G':
    return 'C';
  case 'C':
    return 'G';
  case 'N':
    return 'N';
  case 'a':
    return 't';
  case 't':
    return 'a';
  case 'g':
    return 'c';
  case 'c':
    return 'g';
  }

  return ' ';
}

// all orfs in fasta file
// can do procaryotic circular genomes
// Takes a fasta file and finds all orfs,
// including orfs that span from before stop to over start position of genome
// TODO: We should add the seqinfo from the fasta here, so that especially
//  if genome is circular, we add lengths of chromosomes
// We could also add name of genome from name of fasta file
// And get the real seqnames
// [[Rcpp::export]]
S4 findORFs_fasta(std::string file,
                       std::string startCodon,
                       std::string stopCodon,
                       bool longestORF,
                       int minimumLength,
                       bool isCircular)
{
  std::vector<int> all_orfs;
  std::vector<std::string> Seqnames;
  std::vector<int> strands;
  std::ifstream in(file.c_str());
  in.get(); // remove first '>'
  std::string rec;
  int n = 0;
  int chromoLength;
  while (getline(in, rec, '>')) { // For each chromosome
    int newLineLoc = rec.find('\n');
    std::string header = rec.substr(0, newLineLoc);
    header = header.substr(0, header.find(' '));
    std::string fastaSeq = rec.substr(newLineLoc + 1,
                                      rec.length() - newLineLoc - 2);
    chromoLength = fastaSeq.length() + 1;
    // Check if sequences is lowercase
    if (n == 0) {
      std::locale loc;
      if (std::islower(fastaSeq.at(0), loc)) {
        std::transform(startCodon.begin(), startCodon.end(),
                       startCodon.begin(), ::tolower);
        std::transform(stopCodon.begin(), stopCodon.end(),
                         stopCodon.begin(), ::tolower);
      }
    }

    // get all orfs for start to stop
    std::vector<int> ORFdef = orfs_as_vector(fastaSeq, startCodon,
                                             stopCodon, longestORF,
                                             minimumLength);
    all_orfs.insert(all_orfs.end(), ORFdef.begin(), ORFdef.end());
    Seqnames.insert(Seqnames.end(), ORFdef.size() / 2, header);
    strands.insert(strands.end(), ORFdef.size() / 2, 1);

    //Definitions if isCircular is TRUE
    int rescaleStart = fastaSeq.length() - 4000;
    int rescaleStop = 4000;
    std::string startStopBoundary;
    std::vector<int> ORFdefBoundary;

    if (isCircular) {
      // Now do stop/start boundary, +/-4000, only keep the ones
      // who are overlapping start/stop boundary
      startStopBoundary = fastaSeq.substr(
        rescaleStart, rescaleStop);
      startStopBoundary.append(fastaSeq.substr(0, rescaleStop));
      ORFdefBoundary = orfs_as_vector(startStopBoundary,
                                                       startCodon,
                                                       stopCodon, longestORF,
                                                       minimumLength);
      // now filter out wrong ones, the ones that does not contain point
      // 4000 in range start/stop
      std::vector<int> ORFdefOverlap; // <- vector only for valid ones

      for (size_t i = 0; i < ORFdefBoundary.size() / 2; i++) {
        if (ORFdefBoundary[2 * i] < rescaleStop) { // start
          if (ORFdefBoundary[2 * i + 1] >= rescaleStop) { // stop
            ORFdefOverlap.push_back(ORFdefBoundary[2 * i] + rescaleStart);
            ORFdefOverlap.push_back(ORFdefBoundary[2 * i + 1] + rescaleStart);
          }
        }
      }
      all_orfs.insert(all_orfs.end(), ORFdefOverlap.begin(),
                      ORFdefOverlap.end());
      Seqnames.insert(Seqnames.end(), ORFdefOverlap.size() / 2, header);
      strands.insert(strands.end(), ORFdefOverlap.size() / 2, 1);
    }

    // now do on reverse compliment
    std::reverse(fastaSeq.begin(), fastaSeq.end());
    std::transform(fastaSeq.begin(), fastaSeq.end(), fastaSeq.begin(),
                   complement);

    ORFdef = orfs_as_vector(fastaSeq, startCodon,
                            stopCodon, longestORF,
                            minimumLength); // <-get all orfs for start to stop
    //Following standard, negative strand should have stop <= start
    for(size_t i = 0; i < ORFdef.size(); i++)
      ORFdef[i] = chromoLength - ORFdef[i];
    all_orfs.insert(all_orfs.end(), ORFdef.rbegin(), ORFdef.rend());
    Seqnames.insert(Seqnames.end(), ORFdef.size() / 2, header);
    strands.insert(strands.end(), ORFdef.size() / 2, -1);

    if (isCircular) {
      // Now do stop/start boundary, +/-4000, only keep the ones
      // who are overlapping start/stop boundary
      startStopBoundary = fastaSeq.substr(rescaleStart, rescaleStop);
      startStopBoundary.append(fastaSeq.substr(0, rescaleStop));
      ORFdefBoundary = orfs_as_vector(startStopBoundary, startCodon,
                                      stopCodon, longestORF,
                                      minimumLength);

      // now filter out wrong ones, the ones that does
      // not contain point 4000 in range start/stop
      std::vector<int> ORFdefOverlapMin; // <- vector only for valid ones

      for (size_t i = 0; i < ORFdefBoundary.size() / 2; i++) {
        if (ORFdefBoundary[2 * i] < rescaleStop) { // start
          if (ORFdefBoundary[2 * i + 1] >= rescaleStop) { // stop
            ORFdefOverlapMin.push_back(ORFdefBoundary[2 * i] + rescaleStart);
            ORFdefOverlapMin.push_back(ORFdefBoundary[2 * i + 1] + rescaleStart);
          }
        }
      }
      for(size_t i = 0; i < ORFdefOverlapMin.size(); i++)
        ORFdefOverlapMin[i] = chromoLength - ORFdefOverlapMin[i];
      all_orfs.insert(all_orfs.end(),
                      ORFdefOverlapMin.rbegin(), ORFdefOverlapMin.rend());
      Seqnames.insert(Seqnames.end(), ORFdefOverlapMin.size() / 2, header);
      strands.insert(strands.end(), ORFdefOverlapMin.size() / 2, -1);
    }
    n++;
  }

  // all_orfs is an interlaced vector. We de-interlace it into two vectors.
  std::vector<vi> result_value(2);
  result_value[0].resize(all_orfs.size() / 2);
  result_value[1].resize(all_orfs.size() / 2);
  for (size_t i = 0; i < all_orfs.size() / 2; i++) {
    result_value[0][i] = all_orfs[2 * i];
    result_value[1][i] = all_orfs[2 * i + 1];
  }

  return (GRangesC(Seqnames,
                   IRangesC(wrap(result_value[0]),
                            wrap(result_value[1])),
                            strands));
}
