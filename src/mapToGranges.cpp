//Map either a matrix or IRanges to GRanges
//Cpp format: Webkit: 80 character max line

#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdlib.h>

#include <Rcpp.h>

using ui = unsigned int;
using vi = std::vector<int>;
using string = std::string;
using vc = std::vector<char>;
using vs = std::vector<string>;

using namespace Rcpp;


Function GRanges("GRanges", Environment::namespace_env("GenomicRanges"));
//Environment envGenomicRanges = Environment("package:GenomicRanges");
//Function GRanges = envGenomicRanges["GRanges"];

Function IRanges("IRanges", Environment::namespace_env("IRanges"));

// Environment envIRanges("package:IRanges");
// Function IRanges = envIRanges["IRanges"];


//Global variables
String chrom;
String strands;

//Get start of IRanges object
ui start(const IntegerMatrix& mat, int index)
{
  return mat(index, 0);
}

//Get end of IRanges object
ui end(const IntegerMatrix& mat, int index)
{
  return mat(index, 1);
}

//Get width of GRanges object
ui width(const IntegerMatrix mat, int index)
{
  return (mat(index, 1) - mat(index, 0) + 1);
}

ui length(string s)
{
  return s.size();
}

//Get length of character vector
ui length(vc s)
{
  return s.size();
}

//Get length of Matrix
ui length(const IntegerMatrix& s)
{
  return s.nrow();
}

//Stop cpp from running, but continue R
void stop(string s)
{
  Rcpp::stop(s);
}

//Resize a matrix to size n, rcpp does not include a standard for this
IntegerMatrix resize(const IntegerMatrix& x, int n)
{
  int oldsize = x.nrow();
  IntegerMatrix y(n, 2);

  for (int i = 0; i < oldsize; i++) {
    y(i, 0) = x(i, 0);
    y(i, 1) = x(i, 1);
  }
  return y;
}

//Expand a matrix by size + 1 and, set the last start and end to the extra row
IntegerMatrix c(IntegerMatrix& one, ui start, ui end, ui last)
{

  if (one(0, 0) != 0)
    one = resize(one, one.nrow() + 1);

  one(last, 0) = start;
  one(last, 1) = end;

  return one;
}

//Combination method for two integer matrices
IntegerMatrix c(IntegerMatrix& one, IntegerMatrix& two)
{
  ui oneOldSize = one.nrow();
  ui twoOldSize = two.nrow();

  one = resize(one, oneOldSize + twoOldSize);

  for (int i = 0; i < twoOldSize; i++) {
    one(oneOldSize + i, 0) = two(i, 0);
    one(oneOldSize + i, 1) = two(i, 1);
  }
  return one;
}

//if non unique elements exists, return '0' as error, else return element [0]
String Unique(const CharacterVector& s)
{

  String temp = s[0];

  for (unsigned int i = 1; i < s.size(); i++) {
    if (temp != s[i]) {
      return ('0');
    }
  }
  return temp;
}

//Get all the correct names for the ranges
CharacterVector getAllNames(int numberOfOrfs, string txName, vi& nexons)
{
  CharacterVector orfNames(numberOfOrfs);
  int index = 0;
  int inc = 0;

  for (int& i : nexons) {
    for (int l = 0; l < i + 1; l++) {
      orfNames[index + l] = txName + "_" + std::to_string(inc + 1);
    }
    index += i + 1;
    inc++;
  }
  return orfNames;
}

//repete string
CharacterVector rep(string s, ui times)
{
  CharacterVector vec(times, s);
  return vec;
}


//ORFdef are the orfs as IRanges
//grangesObj is object that orfs are from
//transcriptName is the name of the grangesObj
//The Functions are function from R used to call constructors

S4 forPositiveStrand(IntegerMatrix ORFdef, IntegerMatrix grangesObj,
                     CharacterMatrix grangesObjStrings, string transcriptName)
{
  //Definitions
  ui startingPos;
  ui endingPos;
  ui j; //exon counter
  ui a; //ranges counter
  S4 returnRanges;
  IntegerMatrix ORFranges;
  vi nexons;
  IntegerMatrix ORFGranges;

  for (ui i = 0; i < length(ORFdef); i++) { //For each orf in list

    ORFranges = IntegerMatrix(1, 2);
    startingPos = ORFdef(i, 0) - 1;
    endingPos = ORFdef(i, 1) - 1;

    j = 0;
    a = 0;
    //While loops to jump to next exon, if IRange exeeds first exon width
    while (startingPos > width(grangesObj, j)) {
      startingPos = startingPos - width(grangesObj, j);
      endingPos = endingPos - width(grangesObj, j);
      j = j + 1;
    }

    while (endingPos > width(grangesObj, j)) {

      ORFranges = c(ORFranges, start(grangesObj, j) + startingPos,
                    end(grangesObj, j), a);
      startingPos = 0;
      endingPos = endingPos - width(grangesObj, j);
      j = j + 1;
      a = a + 1;
    }

    ORFranges = c(ORFranges, start(grangesObj, j) + startingPos,
                  start(grangesObj, j) + endingPos, a);

    if (i == 0) //If first
      ORFGranges = ORFranges;
    else
      ORFGranges = c(ORFGranges, ORFranges);

    nexons.push_back(a); //push number of ranges for number transcript _#
  }

  S4 iranges = IRanges(ORFGranges(_, 0), ORFGranges(_, 1));
  CharacterVector names = getAllNames(length(ORFGranges), transcriptName,
                                      nexons);

  returnRanges = GRanges(rep(chrom, length(ORFGranges)),
                         iranges,
                         rep(strands, length(ORFGranges)),
                         Named("names") = names);
  return returnRanges;
}

S4 forNegativeStrand(IntegerMatrix ORFdef, IntegerMatrix grangesObj,
                     CharacterMatrix grangesObjStrings, string transcriptName)
{
  //Definitions
  ui startingPos;
  ui endingPos;
  ui j; //exon counter
  ui a; //ranges counter
  S4 returnRanges;
  IntegerMatrix ORFranges;

  vi nexons;
  IntegerMatrix ORFGranges;

  for (ui i = 0; i < length(ORFdef); i++) { //For each orf in list

    ORFranges = IntegerMatrix(1, 2);
    startingPos = ORFdef(i, 0) - 1;
    endingPos = ORFdef(i, 1) - 1;
    a = 0;
    j = 0;
    //While loops to jump to next exon, if IRange exeeds first exon width
    while (startingPos > width(grangesObj, j)) {
      startingPos = startingPos - width(grangesObj, j);
      endingPos = endingPos - width(grangesObj, j);
      j = j + 1;
    }

    while (endingPos > width(grangesObj, j)) {
      ORFranges = c(ORFranges, start(grangesObj, j),
                    end(grangesObj, j) - startingPos, a);
      startingPos = 0;
      endingPos = endingPos - width(grangesObj, j);
      j = j + 1;
      a = a + 1;
    }

    ORFranges = c(ORFranges, end(grangesObj, j) - endingPos,
                  end(grangesObj, j) - startingPos, a);

    if (i == 0) //if first orf
      ORFGranges = ORFranges;
    else
      ORFGranges = c(ORFGranges, ORFranges);

    nexons.push_back(a); //push number of ranges for number transcript _#
  }

  S4 iranges = IRanges(ORFGranges(_, 0), ORFGranges(_, 1));
  CharacterVector names = getAllNames(length(ORFGranges), transcriptName,
                                      nexons);

  returnRanges = GRanges(rep(chrom, length(ORFGranges)),
                         iranges,
                         rep(strands, length(ORFGranges)),
                         Named("names") = names);
  return returnRanges;
}

// Map from transcript coordinates to genomic coordinates
// Slow version, implemented in c for speedup, not recomended to use for over 5000 ranges.
// ORFdef: intmatrix of transcriptranges, 2 columns, col1 = start, col2 = end
// grangesObjString: charmatrix of transcript strings,
// 2 columns, col1 = seqnames, col2 = strands
// [[Rcpp::export]]
S4 map_to_GRangesC(IntegerMatrix ORFdef,
                   IntegerMatrix grangesObj, CharacterMatrix grangesObjStrings,
                   std::string transcriptName)
{

  chrom = Unique(grangesObjStrings(_, 0));
  strands = Unique(grangesObjStrings(_, 1));

  if (chrom == "0") {
    stop("Different chromosomes in GRanges object");
  }
  if (strands == "0") {
    stop("Different strands in GRanges object");
  }

  S4 returnRanges("GRanges");
  if (length(ORFdef) == 0) {
    return (returnRanges);
  }

  if (strands == "+")
    returnRanges = forPositiveStrand(ORFdef, grangesObj,
                                     grangesObjStrings, transcriptName);

  else
    returnRanges = forNegativeStrand(ORFdef, grangesObj,
                                     grangesObjStrings, transcriptName);

  return (returnRanges);
}
