#include <Rcpp.h>
#include <unordered_map>
#include "../inst/include/ORFik/transcript_mapping_core.h"
#include <cerrno>
#include <cstdlib>
#include <cctype>
using namespace Rcpp;
namespace mapping=orfik::mapping;
static mapping::IntView view(const IntegerVector& x) {return {x.begin(),static_cast<std::size_t>(x.size())};}
static S4 slot_s4(const S4& x,const char* n) {return S4(x.slot(n));}
static IntegerVector expand_rle(const S4& x) {
  IntegerVector values=x.slot("values"),lengths=x.slot("lengths");
  mapping::require(values.size()==lengths.size(),"Invalid Rle slots");
  std::size_t total=0;for(int n:lengths){mapping::require(n>=0,"Invalid Rle length");total+=n;}
  mapping::require(total<=static_cast<std::size_t>(INT_MAX),"Rle is too long");
  IntegerVector out(total);int* p=out.begin();
  for(R_xlen_t i=0;i<values.size();++i) p=std::fill_n(p,lengths[i],values[i]);
  return out;
}
struct RangesView {
  S4 ranges; IntegerVector starts,widths,ends,partitions,strands,chromosomes;
  CharacterVector levels; SEXP names, flat_names; bool grouped,genomic;
  RangesView(SEXP object):ranges(object),names(R_NilValue) {
    grouped=Rf_inherits(object,"CompressedGRangesList")||Rf_inherits(object,"CompressedIRangesList");
    S4 flat=grouped?slot_s4(ranges,"unlistData"):ranges;
    genomic=Rf_inherits(flat,"GRanges");
    if(grouped) {S4 partition=slot_s4(ranges,"partitioning");partitions=partition.slot("end");names=partition.slot("NAMES");}
    S4 ir=genomic?slot_s4(flat,"ranges"):flat;
    mapping::require(Rf_inherits(ir,"IRanges"),"Expected IRanges, GRanges or a compressed ranges list");
    starts=ir.slot("start");widths=ir.slot("width");flat_names=ir.slot("NAMES");
    mapping::require(starts.size()==widths.size() && starts.size()<=INT_MAX,"Invalid range slots");
    ends=IntegerVector(starts.size());for(R_xlen_t j=0;j<starts.size();++j) {
      int64_t e=static_cast<int64_t>(starts[j])+widths[j]-1;
      mapping::require(widths[j]>=0 && e>=INT_MIN && e<=INT_MAX,"Invalid range coordinates");ends[j]=e;
    }
    if(!grouped) {partitions=seq_len(starts.size());names=ir.slot("NAMES");}
    int previous=0;for(int end:partitions) {mapping::require(end>=previous && end<=starts.size(),"Invalid query partition");previous=end;}
    mapping::require(previous==starts.size(),"Incomplete query partition");
    mapping::require(names==R_NilValue || Rf_xlength(names)==partitions.size(),"Invalid range names");
    if(genomic) {
      strands=expand_rle(slot_s4(flat,"strand")); S4 seq=slot_s4(flat,"seqnames");chromosomes=expand_rle(seq);
      IntegerVector vals=seq.slot("values");levels=vals.attr("levels");
      mapping::require(strands.size()==starts.size() && chromosomes.size()==starts.size(),"Invalid genomic slots");
      for(int value:strands) mapping::require(value>=1 && value<=3,"Invalid strand code");
      for(int value:chromosomes) mapping::require(value>=1 && value<=levels.size(),"Invalid sequence code");
    }
  }
  int group_strand(int g) const {return genomic && partitions[g]>0?strands[partitions[g]-1]:1;}
};
static S4 make_iranges(const List& templates,IntegerVector starts,IntegerVector ends,SEXP names=R_NilValue) {
  S4 result=clone(S4(templates["ir"]));IntegerVector widths(starts.size());
  for(R_xlen_t i=0;i<starts.size();++i) widths[i]=mapping::checked_int(static_cast<int64_t>(ends[i])-starts[i]+1);
  result.slot("start")=starts;result.slot("width")=widths;result.slot("NAMES")=names;return result;
}
static S4 make_rle(S4 prototype,const IntegerVector& codes,CharacterVector levels) {
  S4 result=clone(prototype);std::vector<int> values,lengths;
  for(int code:codes) {
    mapping::require(code>=1 && code<=levels.size(),"Invalid factor code");
    if(!values.empty() && values.back()==code) ++lengths.back();else {values.push_back(code);lengths.push_back(1);}
  }
  IntegerVector v=wrap(values);v.attr("levels")=levels;v.attr("class")="factor";
  result.slot("values")=v;result.slot("lengths")=wrap(lengths);return result;
}
static S4 empty_mcols(S4 prototype,int n) {S4 result=clone(prototype);result.slot("nrows")=n;return result;}
static S4 make_granges(const List& templates,S4 ir,IntegerVector seqcodes,CharacterVector seqlevels,
                       IntegerVector strands,S4 seqinfo) {
  S4 result=clone(S4(templates["gr"]));
  result.slot("ranges")=ir;result.slot("seqnames")=make_rle(slot_s4(result,"seqnames"),seqcodes,seqlevels);
  result.slot("strand")=make_rle(slot_s4(result,"strand"),strands,CharacterVector::create("+","-","*"));
  result.slot("seqinfo")=seqinfo;
  result.slot("elementMetadata")=empty_mcols(slot_s4(result,"elementMetadata"),seqcodes.size());return result;
}
static S4 make_list(const List& templates,S4 flat,IntegerVector ends,SEXP names,bool genomic) {
  S4 result=clone(S4(templates[genomic?"grl":"irl"]));S4 partition=slot_s4(result,"partitioning");
  partition.slot("end")=ends;partition.slot("NAMES")=names;result.slot("partitioning")=partition;
  result.slot("unlistData")=flat;
  if(static_cast<SEXP>(result.slot("elementMetadata"))!=R_NilValue) result.slot("elementMetadata")=empty_mcols(slot_s4(result,"elementMetadata"),ends.size());
  return result;
}
static S4 make_seqinfo(const List& templates,CharacterVector levels,IntegerVector lengths) {
  S4 result=clone(slot_s4(S4(templates["gr"]),"seqinfo"));result.slot("seqnames")=levels;
  result.slot("seqlengths")=lengths;result.slot("is_circular")=LogicalVector(levels.size(),NA_LOGICAL);
  result.slot("genome")=CharacterVector(levels.size(),NA_STRING);return result;
}
static SEXP map_to_objects_impl(SEXP query,SEXP reference,bool ignore_strand,bool set_lengths,bool reduce_ranges,List templates) {
  RangesView x(query),tx(reference); const int ng=x.partitions.size();
  mapping::require(ng==tx.partitions.size(),"Paired query and transcript lengths differ");
  IntegerVector txstrands(ng);for(int g=0;g<ng;++g) txstrands[g]=ignore_strand?1:tx.group_strand(g);
  mapping::Catalogue catalogue(view(tx.starts),view(tx.ends),view(tx.partitions),view(txstrands));
  if(x.genomic && tx.genomic) for(SEXP level:x.levels) {
    bool found=false;for(SEXP other:tx.levels) if(level==other){found=true;break;}
    mapping::require(found,"subscript contains out-of-bounds indices");
  }
  IntegerVector starts(x.starts.size()),ends(x.starts.size()),seqcodes(x.starts.size()),strands(x.starts.size());
  CharacterVector group_labels(ng);
  for(int g=0;g<ng;++g) group_labels[g]=tx.names==R_NilValue?Rf_mkChar(std::to_string(g+1).c_str()):STRING_ELT(tx.names,g);
  std::vector<int> levelgroups;std::vector<SEXP> labels;std::unordered_map<SEXP,int> lookup;
  for(int g=0,begin=0;g<ng;++g) {
    checkUserInterrupt();const int stop=x.partitions[g];int64_t query_width=0;
    for(int j=begin;j<stop;++j) query_width+=x.widths[j];
    mapping::require(query_width<=catalogue.widths[g],"Invalid ranges to map, check them. One has width bigger than its reference");
    SEXP label=group_labels[g];
    mapping::require(!x.genomic || label!=NA_STRING,"Transcript names cannot be missing for genomic output");
    int code=0;
    if(stop>begin && x.genomic) {
      auto found=lookup.find(label);
      if(found==lookup.end()) {code=labels.size()+1;lookup.emplace(label,code);labels.push_back(label);levelgroups.push_back(g);}else code=found->second;
    }
    for(int j=begin;j<stop;++j) {
      auto mapped=catalogue.to_transcript(x.starts[j],x.ends[j],g);
      if(x.genomic && tx.genomic && (tx.partitions[g]==0 ||
         STRING_ELT(x.levels,x.chromosomes[j]-1)!=STRING_ELT(tx.levels,tx.chromosomes[tx.partitions[g]-1]-1))) mapped={0,-1};
      starts[j]=mapped.start;ends[j]=mapped.end;seqcodes[j]=code;
      strands[j]=ignore_strand?3:tx.group_strand(g);
      if(mapped.start==0 || (!ignore_strand && x.genomic && strands[j]!=x.strands[j])) strands[j]=3;
    }
    begin=stop;
  }
  IntegerVector partitions=x.partitions;
  bool multiple=false;for(int g=0,begin=0;g<ng;++g){if(partitions[g]-begin>1)multiple=true;begin=partitions[g];}
  if(x.grouped && reduce_ranges && multiple) {
    auto reduced=mapping::reduce_groups(view(starts),view(ends),view(partitions),view(strands),view(seqcodes),x.genomic);
    starts=wrap(reduced.starts);ends=wrap(reduced.ends);strands=wrap(reduced.strands);seqcodes=wrap(reduced.sequences);partitions=wrap(reduced.partitions);
  }
  S4 result=make_iranges(templates,starts,ends,x.grouped?R_NilValue:x.names);
  if(x.genomic) {
    CharacterVector levels(labels.size());IntegerVector lengths(labels.size(),NA_INTEGER);
    for(std::size_t i=0;i<labels.size();++i) {levels[i]=labels[i];if(set_lengths && tx.names==R_NilValue) lengths[i]=catalogue.widths[levelgroups[i]];}
    result=make_granges(templates,result,seqcodes,levels,strands,make_seqinfo(templates,levels,lengths));
  }
  if(x.grouped) result=make_list(templates,result,partitions,x.names,x.genomic);
  return result;
}
static SEXP map_from_objects_impl(SEXP query,SEXP reference,bool drop,List templates,IntegerVector indices) {
  RangesView x(query),tx(reference);mapping::require(tx.grouped && tx.genomic,"transcripts must be a compressed GRangesList");
  mapping::require(indices.size()==x.starts.size(),"length of ranges != indices");
  IntegerVector txstrands(tx.partitions.size());for(int g=0;g<txstrands.size();++g) txstrands[g]=tx.group_strand(g);
  mapping::Catalogue catalogue(view(tx.starts),view(tx.ends),view(tx.partitions),view(txstrands));
  IntegerVector groups(indices.size());for(int i=0;i<groups.size();++i) {mapping::require(indices[i]>=1 && indices[i]<=tx.partitions.size(),"Invalid transcript index");groups[i]=indices[i]-1;}
  auto plan=mapping::plan_from(view(x.starts),view(x.ends),view(groups),catalogue,drop);
  IntegerVector starts(plan.count),ends(plan.count),querygroups(plan.count),partitions(indices.size());
  mapping::fill_from(view(x.starts),view(x.ends),catalogue,plan,drop,starts.begin(),ends.begin(),querygroups.begin(),partitions.begin(),checkUserInterrupt);
  IntegerVector chromosomes(plan.count),strands(plan.count);CharacterVector range_names(plan.count),list_names(indices.size());
  for(int i=0;i<indices.size();++i) list_names[i]=tx.names==R_NilValue?NA_STRING:STRING_ELT(tx.names,indices[i]-1);
  for(std::size_t i=0;i<plan.count;++i) {
    const int g=indices[querygroups[i]-1]-1;
    chromosomes[i]=tx.chromosomes[tx.partitions[g]-1];strands[i]=tx.group_strand(g);range_names[i]=list_names[querygroups[i]-1];
  }
  S4 flat=slot_s4(S4(reference),"unlistData");
  S4 result=make_granges(templates,make_iranges(templates,starts,ends,tx.names==R_NilValue?R_NilValue:static_cast<SEXP>(range_names)),chromosomes,tx.levels,strands,slot_s4(flat,"seqinfo"));
  return make_list(templates,result,partitions,tx.names==R_NilValue?R_NilValue:static_cast<SEXP>(list_names),true);
}

static bool logical_flag(SEXP value,const char* name) {
  if(TYPEOF(value)!=LGLSXP || Rf_xlength(value)!=1 || LOGICAL(value)[0]==NA_LOGICAL)
    stop("%s must be a single non-missing logical",name);
  return LOGICAL(value)[0];
}
static int object_groups(SEXP object) {
  if(Rf_inherits(object,"CompressedGRangesList") || Rf_inherits(object,"CompressedIRangesList"))
    return Rf_xlength(slot_s4(S4(object),"partitioning").slot("end"));
  if(Rf_inherits(object,"GRanges"))return Rf_xlength(slot_s4(S4(object),"ranges").slot("start"));
  if(Rf_inherits(object,"IRanges"))return Rf_xlength(S4(object).slot("start"));
  stop("Expected IRanges, GRanges or a compressed ranges list");return 0;
}
// [[Rcpp::export]]
SEXP pmapToObjectsCPP(SEXP query,SEXP reference,SEXP ignoreStrand,SEXP querySorted,
                     SEXP referenceSorted,SEXP reduceRanges,SEXP setLengths,List templates) {
  const bool ignore=logical_flag(ignoreStrand,"ignore.strand");
  const bool qs=logical_flag(querySorted,"x.is.sorted"),ts=logical_flag(referenceSorted,"tx.is.sorted");
  const bool reduce=logical_flag(reduceRanges,"reduce.ranges"),lengths=logical_flag(setLengths,"set.seqlengths");
  int nq=object_groups(query),nt=object_groups(reference);if(!nq)return query;
  RObject x(query),tx(reference);
  if(nq!=nt) {
    Function rep=Environment::base_env()["rep"];
    if(nq==1)x=rep(x,nt);else if(nt==1)tx=rep(tx,nq);
    else stop("recycling is supported when length(x) == 1 or length(transcripts) == 1; otherwise the lengths must match");
  }
  if((!qs && Rf_inherits(x,"CompressedGRangesList")) || (!ts && Rf_inherits(tx,"CompressedGRangesList"))) {
    Function sort=Environment::namespace_env("ORFik")["sortPerGroup"];
    if(!qs && Rf_inherits(x,"CompressedGRangesList"))x=sort(x,ignore);
    if(!ts && Rf_inherits(tx,"CompressedGRangesList"))tx=sort(tx,ignore);
  }
  return map_to_objects_impl(x,tx,ignore,lengths,reduce,templates);
}
static int parse_index(SEXP string,int limit) {
  if(string==NA_STRING)stop("Invalid transcript index");
  char* end;errno=0;const char* text=CHAR(string);long index=std::strtol(text,&end,0);
  while(*end && std::isspace(static_cast<unsigned char>(*end)))++end;
  if(errno || end==text || *end || index<1 || index>limit)stop("Invalid transcript index");
  return index;
}
// [[Rcpp::export]]
SEXP pmapFromObjectsCPP(SEXP query,SEXP reference,SEXP removeEmpty,List templates) {
  const bool drop=logical_flag(removeEmpty,"removeEmpty");
  int nq=object_groups(query),nt=object_groups(reference);
  if(!nq) {Function subset=Environment::base_env()["["];return subset(reference,IntegerVector());}
  RangesView x(query);if(x.names==R_NilValue && nq!=nt)stop("invalid names of ranges, when length(x) != length(transcripts)");
  if(x.genomic && x.grouped)stop("x must be IRanges, IRangesList or GRanges");
  IntegerVector indices(x.starts.size());
  for(int g=0,begin=0;g<nq;++g) {
    for(int j=begin;j<x.partitions[g];++j) {
      if(x.grouped && x.names!=R_NilValue) {
        if(x.flat_names!=R_NilValue)stop("Invalid transcript index: named inner ranges are not integer transcript identifiers");
        indices[j]=parse_index(STRING_ELT(x.names,g),nt);
      } else if(x.flat_names!=R_NilValue) indices[j]=parse_index(STRING_ELT(x.flat_names,j),nt);
      else indices[j]=j+1;
    }
    begin=x.partitions[g];
  }
  return map_from_objects_impl(query,reference,drop,templates,indices);
}
