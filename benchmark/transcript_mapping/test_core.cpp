#include "../../inst/include/ORFik/transcript_mapping_core.h"
#include <cassert>
#include <iostream>
#include <climits>
using namespace orfik::mapping;
int main() {
  int s[]={100,200,400},e[]={109,219,409},p[]={2,3},strand[]={1,2};
  Catalogue tx({s,3},{e,3},{p,2},{strand,2});
  auto a=tx.to_transcript(205,208,0);assert(a.start==16 && a.end==19);
  a=tx.to_transcript(402,405,1);assert(a.start==5 && a.end==8);
  a=tx.to_transcript(402,405,0);assert(a.start==0 && a.end==-1);
  int qs[]={5,2,16},qe[]={15,8,20},g[]={0,1,0};
  auto plan=plan_from({qs,3},{qe,3},{g,3},tx,true);
  assert(plan.count==4);std::vector<int> os(plan.count),oe(plan.count),og(plan.count),op(3);
  fill_from({qs,3},{qe,3},tx,plan,true,os.data(),oe.data(),og.data(),op.data());
  assert(os[0]==104 && oe[0]==109 && os[1]==200 && oe[1]==204);
  assert(os[2]==402 && oe[2]==408 && os[3]==205 && oe[3]==209);
  assert(op[0]==2 && op[1]==3 && op[2]==4);
  bool rejected=false;try {tx.window(0,2,0);}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
  int ms[]={1,INT_MAX},me[]={INT_MAX,INT_MAX},mp[]={2},mt[]={1};
  rejected=false;try {Catalogue overflow({ms,2},{me,2},{mp,1},{mt,1});}catch(const std::invalid_argument&){rejected=true;}assert(rejected);
  for(int direction=1;direction<=2;++direction) {
    std::vector<int> ss(200),ee(200);for(int j=0;j<200;++j){ss[j]=100+j*30;ee[j]=ss[j]+19;}
    if(direction==2){std::reverse(ss.begin(),ss.end());std::reverse(ee.begin(),ee.end());}
    int end=200;Catalogue many({ss.data(),200},{ee.data(),200},{&end,1},{&direction,1});
    for(int j=0;j<200;++j){a=many.to_transcript(ss[j]+1,ss[j]+5,0);assert(a.start==j*20+(direction==1?2:15));}
  }
  std::cout << "Standalone mapping core checks passed without R or Rcpp.\n";
}
