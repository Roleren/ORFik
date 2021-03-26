//Uses Knuth–Morris–Pratt algorithm to search global string for substring
//Cpp format: Webkit: 80 character max line

#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

#include <Rcpp.h>

using vi = std::vector<int>;
using string = std::string;

using namespace Rcpp;


Function IRangesA("IRanges", Environment::namespace_env("IRanges"));

////////////--Functions--/////////////
// z - scoring list
// return: indices of all z >= substring length
vi get_index_list(const vi z, const  string& working_string,
                  const string& substring)
{
  size_t counter = 0;
  const size_t subSize = substring.size();
  const int subSizeInt = static_cast<unsigned int>(subSize);
  for (size_t i = subSize; i < working_string.size(); ++i)
    if (z[i] >= subSizeInt)
      counter++;


  vi indeces(counter, 0);
  counter = 0;

  for (size_t i = subSize; i < working_string.size(); ++i)
    if (z[i] >= subSizeInt)
      indeces[counter++] = i - subSize;
  return indeces;
}

// The string-searching algorithm
// s - the string to search, z - the scoring vector
// To understand, look up Knut-morris-pratt and z algorithm.
void calc_z(const string& s, vi& z)
{
  const int len = s.size();
  z.resize(len);

  int l = 0, r = 0;
  // For each position in string length l
  for (int i = 1; i < len; ++i)
    if (z[i - l] + i <= r)
      z[i] = z[i - l];
    else {
      l = i;
      if (i > r)
        r = i;
      for (z[i] = r - i; r < len; ++r, ++z[i])
        if (s[r] != s[z[i]])
          break;
      --r;
    }
}

//returns a list of indeces for searched substring given mainstring
vi return_outputs_of_substring(string& main_string, const string& substring)
{

  string working_string = substring + main_string;
  vi z;
  calc_z(working_string, z);
  z = get_index_list(z, working_string, substring);
  return (z);
}

//Find all orf in either frame 0,1 or 2.
vi find_orfs_in_specific_frame(const vi &frameS,const vi &frameE,
                               const int endSize, const int max_size)
{
  vi res(max_size * 2, -1);
  unsigned int counter = 0;
  //For each valid start u, find valid stop v that is greater than u.

  //binary search for end that is > start
  for (auto& u : frameS) {
    if (u == -1)
      break;
    auto it = std::lower_bound(frameE.begin(),
                               frameE.begin() + endSize, u + 1);
    if (it != frameE.end()) {
      auto& v = *it;
      res[counter] = u + 1;
      res[counter + 1] = v + 1;
      counter = counter + 2;
    }
  }
  res.resize(counter);
  return res;
}
// Combine all three frames
vi find_matched_startends(const vi& starts, const vi& ends, const int max_size)
{

  //the 3 possible frames of orfs
  vi sFrame0(max_size, -1);
  vi sFrame1(max_size, -1);
  vi sFrame2(max_size, -1);

  vi eFrame0(max_size, -1);
  vi eFrame1(max_size, -1);
  vi eFrame2(max_size, -1);

  //Counters for lists
  int zeroC = 0;
  int oneC = 0;
  int twoC = 0;
  int modHit = 0;
  for (auto i : starts) {
    modHit = (i % 3);
    if (modHit == 0)
      sFrame0[zeroC++] = i;
    else if (modHit == 1)
      sFrame1[oneC++] = i;
    else
      sFrame2[twoC++] = i;
  }

  zeroC = 0;
  oneC = 0;
  twoC = 0;
  modHit = 0;
  for (auto i : ends) {
    modHit = (i % 3);
    if (modHit == 0)
      eFrame0[zeroC++] = i + 2;
    else if (modHit == 1)
      eFrame1[oneC++] = i + 2;
    else
      eFrame2[twoC++] = i + 2;
  }

  vi tempRes;
  vi res = find_orfs_in_specific_frame(sFrame0, eFrame0, zeroC, max_size);
  tempRes.insert(tempRes.end(), res.begin(), res.end());
  res = find_orfs_in_specific_frame(sFrame1, eFrame1, oneC, max_size);
  tempRes.insert(tempRes.end(), res.begin(), res.end());
  res = find_orfs_in_specific_frame(sFrame2, eFrame2, twoC, max_size);
  tempRes.insert(tempRes.end(), res.begin(), res.end());

  return tempRes;
}
//Get all hits of either starts or ends
vi get_all_hits(string& main_string, string s)
{
  std::stringstream sStream(s);
  string segment;

  vi tempS;
  while (getline(sStream, segment, '|')) {
    vi starts = return_outputs_of_substring(main_string, segment);
    tempS.insert(tempS.end(), starts.begin(), starts.end());
  }
  sort(tempS.begin(), tempS.end());

  return tempS;
}

// Return ranges as vector, only for internal c++ use!!
vi orfs_as_vector(std::string &main_string, const std::string s,
                                const std::string e, int minimumLength)
{
  minimumLength = 6 + (minimumLength * 3) - 1;

  vi tempStarts = get_all_hits(main_string, s); //Get starts
  vi tempEnds = get_all_hits(main_string, e); //Get ends
  size_t max_size = main_string.length(); //maximun number of orfs

  vi res = find_matched_startends(tempStarts, tempEnds, max_size);
  int nHits = 0; //How many uorfs have current max length
  int cl; //length of current orf

  vi maxUORF(max_size, 0);
  for (size_t i = 0; i < res.size(); i = i + 2) {
    cl = (res[i + 1] - res[i]);
    if (cl >= minimumLength) {
      maxUORF[nHits * 2] = res[i];
      maxUORF[(nHits * 2) + 1] = res[i + 1];
      nHits++;
    }
  }
  //Resize
  maxUORF.resize(nHits*2);

  return maxUORF; //Returns as matrix
}

// Now used in the fast version of orf finding, use _as_IRanges for
// direct use in R.
IntegerMatrix orfs_as_matrix(std::string &main_string, std::string s,
                             std::string e, int minimumLength)
{

  vi maxUORF =  orfs_as_vector(main_string, s, e, minimumLength);
  size_t uorfSize = maxUORF.size();

  if (uorfSize == 0) {
    IntegerMatrix a;
    return a;
  }

  IntegerMatrix mat(uorfSize/2, 2);

  for (int i = 0; i < uorfSize/2; i++) {
    mat(i, 0) = maxUORF[i * 2];
    mat(i, 1) = maxUORF[(i * 2) + 1];
  }

  return mat; //Returns as matrix
}

//main_string is a dna sequence[ATCG] as string
//s is start codons allowed, seperated by "|"
//e is end codons (stop codons) allowed
//Minimum length filters the list to only contain orfs longer...
//or equal to this number of triplets
// [[Rcpp::export]]
S4 orfs_as_IRanges(std::string &main_string, const std::string s,
                           const std::string e, int minimumLength)
{
  size_t minLength = 6 + (minimumLength * 3) - 1;
  if (main_string.length() < 6 ||
      main_string.length() < minLength) {
    S4 I("IRanges");
    return I;
  }
  IntegerMatrix mat = orfs_as_matrix(main_string, s, e, minimumLength);

  if (mat.ncol() == 0) {
    S4 I("IRanges");
    return I;
  }

  return IRangesA(mat(_, 0), mat(_, 1));
}

