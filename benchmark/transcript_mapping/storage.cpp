#include <Rcpp.h>
#include <unordered_map>
#include "../../inst/include/ORFik/transcript_mapping_core.h"
using namespace Rcpp;
namespace mapping=orfik::mapping;
static mapping::IntView view(const IntegerVector& x) {return {x.begin(),static_cast<std::size_t>(x.size())};}
// Same kernel; vary only access/ownership representation.
// [[Rcpp::export]]
List candidate_storage(IntegerVector qs,IntegerVector qe,IntegerVector ts,IntegerVector te,
                       IntegerVector indices,char direction,IntegerVector offsets,int mode) {
  if(mode==3) {
    std::vector<int> s(qs.begin(),qs.end()),e(qe.begin(),qe.end()),a(ts.begin(),ts.end()),b(te.begin(),te.end()),g(indices.begin(),indices.end()),o(offsets.begin(),offsets.end());
    std::vector<int> os(qs.size()),oe(qs.size());
    mapping::map_to({s.data(),s.size()},{e.data(),e.size()},{a.data(),a.size()},{b.data(),b.size()},{g.data(),g.size()},direction,{o.data(),o.size()},os.data(),oe.data());
    IntegerVector output_s=wrap(os),output_e=wrap(oe);
    return List::create(output_s,output_e);

  }
  IntegerVector starts(qs.size()),ends(qs.size());
  if(mode==0) {
    mapping::map_to(view(qs),view(qe),view(ts),view(te),view(indices),direction,view(offsets),starts.begin(),ends.begin());
  } else if(mode==1) {
    auto partitions=mapping::partitions_from_indices(view(indices));std::vector<int> strands(partitions.size(),direction=='+'?1:2);
    mapping::Catalogue tx(view(ts),view(te),{partitions.data(),partitions.size()},{strands.data(),strands.size()});
    std::vector<int> groups(ts.size(),-1);for(std::size_t g=0;g<partitions.size();++g)groups[tx.boundaries[g]]=g;
    for(R_xlen_t i=0;i<qs.size();++i) {
      mapping::Interval r{0,-1};const int base=offsets[i];if(base>=0 && base<ts.size())r=tx.to_transcript(qs[i],qe[i],groups[base]);starts[i]=r.start;ends[i]=r.end;
    }
  } else if(mode==2) {
    std::vector<int> s(qs.begin(),qs.end()),e(qe.begin(),qe.end()),a(ts.begin(),ts.end()),b(te.begin(),te.end()),g(indices.begin(),indices.end()),o(offsets.begin(),offsets.end());
    mapping::map_to({s.data(),s.size()},{e.data(),e.size()},{a.data(),a.size()},{b.data(),b.size()},{g.data(),g.size()},direction,{o.data(),o.size()},starts.begin(),ends.begin());
  } else stop("Unknown storage mode");
  return List::create(starts,ends);
}
