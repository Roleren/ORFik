//Uses Knuth–Morris–Pratt algorithm to search global string for substring
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>

#include <Rcpp.h>

using vi = std::vector<int>;
using string = std::string;

using namespace Rcpp;

////////////--Functions--/////////////

vi get_index_list(vi z,string &working_string,string &substring){
    int counter = 0;
    int subSize = substring.size();

    for(int i = subSize; i < working_string.size(); ++i)
           if(z[i] >=subSize)
                counter++;
    vi indeces(counter,0);
    counter = 0;


    for(int i = subSize; i < working_string.size(); ++i)
           if(z[i] >=subSize)
               indeces[counter++] =  i - substring.size();
    return indeces;
}

//The string-searching algorithm

void calc_z (string &s, vi & z)
{
    int len = s.size();
    z.resize (len);

    int l = 0, r = 0;
    for (int i=1; i<len; ++i)
        if (z[i-l]+i <= r)
            z[i] = z[i-l];
        else
        {
            l = i;
            if (i > r) r = i;
            for (z[i] = r-i; r<len; ++r, ++z[i])
                if (s[r] != s[z[i]])
                    break;
            --r;
        }
}
//returns a list of indeces for searched substring given mainstring
vi return_outputs_of_substring(string &main_string,string &substring){


    string working_string = substring + main_string;
    vi z;
    calc_z(working_string, z);
    // print_indeces(z,working_string,substring);
    z = get_index_list(z,working_string,substring);
    return(z);
}
//Find all orf in either frame 0,1 or 2.
vi find_orfs_in_specific_frame(vi frameS,vi frameE,int max_size){
    vi res(max_size*2,-1);
    int counter = 0;
    //bool used[max_size];

    for(auto u : frameS){

        if(u != -1){
            for(auto v : frameE){
                if((v != -1) && v > u){
                    res[counter] = u + 1;
                    res[counter+1] = v + 1;
                    counter = counter + 2;
                    break;
                }
            }
        }
    }
    res.resize(counter);
    return res;
}
// Combine all three frames
vi find_matched_startends(vi& starts,vi& ends, int max_size){

    //the 3 possible frames of orfs
    vi sFrame0(max_size,-1);
    vi sFrame1(max_size,-1);
    vi sFrame2(max_size,-1);

    vi eFrame0(max_size,-1);
    vi eFrame1(max_size,-1);
    vi eFrame2(max_size,-1);

    //Counters for lists
    int zeroC = 0;
    int oneC = 0;
    int twoC = 0;
    int modHit  = 0;
    for(auto i : starts){
        modHit = ( i % 3 );
        if(modHit == 0)
            sFrame0[zeroC++] = i;
        else if(modHit == 1)
            sFrame1[oneC++] = i;
        else
            sFrame2[twoC++] = i;

    }

    zeroC = 0;
    oneC = 0;
    twoC = 0;
    modHit  = 0;
    for(auto i : ends){
        modHit = ( i % 3 );
        if(modHit == 0)
            eFrame0[zeroC++] = i + 2;
        else if(modHit == 1)
            eFrame1[oneC++] = i + 2;
        else
            eFrame2[twoC++] = i + 2;
    }

    vi tempRes;
    vi res = find_orfs_in_specific_frame(sFrame0,eFrame0,max_size);
    tempRes.insert( tempRes.end(), res.begin(), res.end() );
    res = find_orfs_in_specific_frame(sFrame1,eFrame1,max_size);
    tempRes.insert( tempRes.end(), res.begin(), res.end() );
    res = find_orfs_in_specific_frame(sFrame2,eFrame2,max_size);
    tempRes.insert( tempRes.end(), res.begin(), res.end() );

    return tempRes;
}
//Get all hits of either starts or ends
vi get_all_hits(string &main_string,string s){
  std::stringstream sStream(s);
  string segment;

  vi tempS;
  while( getline(sStream, segment, '|') )
  {
    vi starts = return_outputs_of_substring(main_string,segment);
    tempS.insert( tempS.end(), starts.begin(), starts.end() );
  }
  sort(tempS.begin(),tempS.end());

  return tempS;
}

S4 IRangesCC(IntegerMatrix & Mat){
  S4 I("IRanges");

  I.slot("start") = Mat( _, 0);
  I.slot("width") = Mat( _, 1) - Mat(_,0) + 1;

  return I;
}



// [[Rcpp::export]]
IntegerMatrix get_all_orfs_as_matrix(std::string main_string,std::string s,std::string e ,bool longestORF, int minimumLength){

    minimumLength = 6 + (minimumLength*3) -1;

    vi tempStarts = get_all_hits(main_string,s); //Get starts
    vi tempEnds = get_all_hits(main_string,e);  //Get ends


    int max_size =  ( main_string.length() ); //maximun number of orfs are third of total

    vi res = find_matched_startends(tempStarts,tempEnds,max_size);
    //for(auto u : res)
    //  D(u);
    int tempMax = 0;
    int nHits = 0; //How many uorfs have current max length
    int cl; //length of current orf

    vi maxUORF(max_size,0);
    for(int i = 0; i < res.size();i=i+2){
        cl = (res[i+1]-res[i]);
        if( cl  >= minimumLength){
            if(longestORF){
                if( cl > tempMax){
                    tempMax = res[i+1]-res[i];
                    maxUORF[0] = res[i];
                    maxUORF[1] = res[i+1];
                    nHits = 1;
                }
                else if( cl == tempMax){
                    maxUORF[nHits*2] = res[i];
                    maxUORF[(nHits*2) + 1] = res[i+1];
                    nHits++;
                }
            }
            else{ //if longestORF == false
                //D(" hit" );
                maxUORF[nHits*2] = res[i];
                maxUORF[(nHits*2) + 1] = res[i+1];
                nHits++;
            }
        }
    }

    if(tempMax == 0 && longestORF == true){
        IntegerMatrix a;
        return a;
    }

    IntegerMatrix mat(nHits,2);

    for(int i = 0; i < nHits;i++){
        mat(i,0) = maxUORF[i*2];
        mat(i,1) = maxUORF[(i*2)+1];
    }

    return mat; //Returns as matrix
}
//This method does not follow the API of GRanges!
//Main string is fasta sequence
//s is start codons allowed, seperated by "|"
//e is end codons (stop codons) allowed
//longest ORF is a bool returning only the longest or all orfs
//Minimum length filters the list to only contain orf longer or equal to this number of triplets
// [[Rcpp::export]]
S4 get_all_orfs_as_IRanges(std::string main_string,std::string s,std::string e ,bool longestORF, int minimumLength){
  if(main_string.length() < 6 || main_string.length() < (6 + (minimumLength*3) -1) ){
    S4 I("IRanges");
    return I;
  }
  IntegerMatrix mat = get_all_orfs_as_matrix(main_string,s,e ,longestORF, minimumLength);
  S4 ir = IRangesCC(mat);
  return ir;
}

// [[Rcpp::export]]
CharacterVector read_fasta(std::string file) {
  CharacterVector records;
  std::ifstream in(file.c_str());
  in.get(); // remove first '>'
  std::string rec;
  while(getline(in,rec,'>')){
    int newLineLoc = rec.find('\n');
    std::string header = rec.substr(0,newLineLoc);
    std::string sequence = rec.substr(newLineLoc+1, rec.length()-newLineLoc-2);
    sequence.erase(std::remove(sequence.begin(),sequence.end(),'\n'),sequence.end());
    records[header]=sequence;
  }
  return(records);
}



