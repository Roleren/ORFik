#ifndef ORFIK_TRANSCRIPT_MAPPING_CORE_H
#define ORFIK_TRANSCRIPT_MAPPING_CORE_H
#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace orfik { namespace mapping {
struct IntView {
  const int* data;
  std::size_t size;
  int operator[](std::size_t i) const { return data[i]; }
};
inline void require(bool condition, const char* message) {
  if (!condition) throw std::invalid_argument(message);
}
inline int checked_int(int64_t x) {
  require(x >= 0 && x <= std::numeric_limits<int>::max(), "Coordinate exceeds integer range");
  return static_cast<int>(x);
}
struct Interval { int start, end; };
// Borrow reference coordinates; own only the reusable search index.
class Catalogue {
public:
  IntView starts, ends;
  std::vector<int> boundaries, prefix, cumulative, widths;
  std::vector<unsigned char> reverse, ordered;
  Catalogue(IntView s, IntView e, IntView partitions, IntView strands)
    : starts(s), ends(e), prefix(s.size), cumulative(s.size) {
    require(s.size == e.size, "Inconsistent mapping vector lengths");
    require(s.size <= static_cast<std::size_t>(std::numeric_limits<int>::max()), "Too many transcript exons");
    require(strands.size == partitions.size, "Inconsistent transcript strands");
    boundaries.push_back(0);
    for (std::size_t g = 0; g < partitions.size; ++g) {
      const int begin = boundaries.back(), end = partitions[g];
      require(end >= begin && static_cast<std::size_t>(end) <= s.size, "Invalid transcript partition");
      require(strands[g] == 1 || strands[g] == 2 || strands[g] == 3, "Invalid transcript strand");
      bool sorted = true; const bool rev = strands[g] != 1; int64_t total = 0;
      for (int j = begin; j < end; ++j) {
        require(s[j] >= 1 && e[j] >= s[j], "Invalid transcript exon coordinates");
        prefix[j] = checked_int(total);
        total += static_cast<int64_t>(e[j]) - s[j] + 1;
        cumulative[j] = checked_int(total);
        if (j > begin && (rev ? e[j] >= s[j-1] : s[j] <= e[j-1])) sorted = false;
      }
      widths.push_back(checked_int(total)); reverse.push_back(rev); ordered.push_back(sorted);
      boundaries.push_back(end);
    }
    require(static_cast<std::size_t>(boundaries.back()) == s.size, "Incomplete transcript partition");
  }
  Interval to_transcript(int s, int e, std::size_t group) const {
    require(group < widths.size(), "Invalid transcript index");
    if (s < 1 || e < s) return {0,-1};
    const int begin = boundaries[group], end = boundaries[group+1]; int j = begin;
    if (ordered[group] && end-begin > 8) {
      if (!reverse[group]) j = std::lower_bound(ends.data+begin, ends.data+end, e)-ends.data;
      else j = std::lower_bound(starts.data+begin, starts.data+end, s, std::greater<int>())-starts.data;
    } else while (j < end && !(starts[j] <= s && ends[j] >= e)) ++j;
    if (j == end || starts[j] > s || ends[j] < e) return {0,-1};
    if (!reverse[group]) return {prefix[j]+(s-starts[j]+1),prefix[j]+(e-starts[j]+1)};
    return {prefix[j]+(ends[j]-e+1),prefix[j]+(ends[j]-s+1)};
  }
  Interval window(int s, int e, std::size_t group) const {
    require(group < widths.size(), "Invalid transcript index");
    // An empty query (e == s - 1, e.g. a window clipped to nothing at a transcript end) maps to an empty
    // genomic range at that position, as before the 1.29.16 rewrite (windowPerGroup upstream of an ORF
    // starting at transcript position 1 or 2).
    const bool empty = e == s - 1;
    require(s >= 1 && (e >= s || empty) && e <= widths[group] && s <= widths[group] + (empty ? 1 : 0),
            "Transcript query coordinates are out of bounds");
    const int begin = boundaries[group], end = boundaries[group+1];
    if (empty) {
      const int k = static_cast<int>(std::lower_bound(cumulative.data()+begin,cumulative.data()+end,
                                                      std::min(s, widths[group]))-cumulative.data());
      return {k, k};
    }
    return {static_cast<int>(std::lower_bound(cumulative.data()+begin,cumulative.data()+end,s)-cumulative.data()),
            static_cast<int>(std::lower_bound(cumulative.data()+begin,cumulative.data()+end,e)-cumulative.data())};
  }
  Interval to_genome(int s, int e, int exon) const {
    const int a = std::max(s,prefix[exon]+1)-prefix[exon]-1;
    const int b = std::min(e,cumulative[exon])-prefix[exon]-1;
    // Exon orientation is supplied by the caller's group.
    return {a,b};
  }
};
struct ReducedRanges {
  std::vector<int> starts, ends, strands, sequences, partitions;
};
inline ReducedRanges reduce_groups(IntView starts, IntView ends, IntView partitions,
                                    IntView strands, IntView sequences, bool genomic) {
  ReducedRanges out;
  int begin=0;
  for(std::size_t g=0;g<partitions.size;++g) {
    std::vector<int> rows;
    for(int j=begin;j<partitions[g];++j) rows.push_back(j);
    std::sort(rows.begin(),rows.end(),[&](int a,int b) {
      if(genomic && strands[a]!=strands[b]) return strands[a]<strands[b];
      if(starts[a]!=starts[b]) return starts[a]<starts[b];
      return ends[a]<ends[b];
    });
    const std::size_t group_begin=out.starts.size();
    for(int j:rows) {
      if(out.starts.size()>group_begin && (!genomic || out.strands.back()==strands[j]) &&
         static_cast<int64_t>(starts[j])<=static_cast<int64_t>(out.ends.back())+1) {
        out.ends.back()=std::max(out.ends.back(),ends[j]);
      } else {
        out.starts.push_back(starts[j]);out.ends.push_back(ends[j]);
        out.strands.push_back(strands[j]);out.sequences.push_back(sequences[j]);
      }
    }
    out.partitions.push_back(checked_int(out.starts.size()));begin=partitions[g];
  }
  return out;
}
inline std::vector<int> partitions_from_indices(IntView indices) {
  std::vector<int> result;
  for (std::size_t j=0; j<indices.size;) {
    require(indices[j] >= 1, "Invalid transcript indices");
    std::size_t end=j+1; while(end<indices.size && indices[end]==indices[j]) ++end;
    result.push_back(checked_int(end)); j=end;
  }
  return result;
}
inline void map_to(IntView qs, IntView qe, IntView ts, IntView te, IntView indices,
                   char direction, IntView offsets, int* out_s, int* out_e,
                   void (*interrupt)() = nullptr) {
  require(qs.size == qe.size && qs.size == offsets.size && ts.size == te.size && ts.size == indices.size,
          "Inconsistent mapping vector lengths");
  require(direction == '+' || direction == '-', "Invalid mapping direction");
  if (!qs.size) return;
  auto parts=partitions_from_indices(indices); std::vector<int> strands(parts.size(),direction=='+'?1:2);
  Catalogue tx(ts,te,{parts.data(),parts.size()},{strands.data(),strands.size()});
  std::vector<int> groups(ts.size,-1);
  for(std::size_t g=0;g<parts.size();++g) groups[tx.boundaries[g]]=g;
  for(std::size_t i=0;i<qs.size;++i) {
    if(interrupt && !(i&65535)) interrupt();
    const int base=offsets[i]; Interval result{0,-1};
    if(base>=0 && static_cast<std::size_t>(base)<ts.size) {
      require(groups[base]>=0,"Transcript offset must identify its first exon");
      result=tx.to_transcript(qs[i],qe[i],groups[base]);
    }
    out_s[i]=result.start;out_e[i]=result.end;
  }
}
struct FromPlan {
  std::vector<Interval> windows;
  std::vector<int> groups;
  std::size_t count=0;
};
inline FromPlan plan_from(IntView qs, IntView qe, IntView groups, const Catalogue& tx, bool drop) {
  require(qs.size==qe.size && qs.size==groups.size,"Inconsistent mapping vector lengths");
  FromPlan plan; plan.windows.reserve(qs.size);plan.groups.reserve(qs.size);
  for(std::size_t i=0;i<qs.size;++i) {
    require(groups[i]>=0 && static_cast<std::size_t>(groups[i])<tx.widths.size(),"Invalid transcript index");
    const int g=groups[i]; auto w=tx.window(qs[i],qe[i],g);plan.windows.push_back(w);plan.groups.push_back(g);
    plan.count+=drop?w.end-w.start+1:tx.boundaries[g+1]-tx.boundaries[g];
    require(plan.count<=static_cast<std::size_t>(std::numeric_limits<int>::max()),"Too many mapped ranges");
  }
  return plan;
}
inline void fill_from(IntView qs, IntView qe, const Catalogue& tx, const FromPlan& plan, bool drop,
                      int* out_s,int* out_e,int* out_group,int* out_partition,
                      void (*interrupt)()=nullptr) {
  std::size_t row=0;
  for(std::size_t i=0;i<qs.size;++i) {
    if(interrupt && !(i&65535)) interrupt();
    const int g=plan.groups[i]; auto w=plan.windows[i];
    const int begin=drop?w.start:tx.boundaries[g],end=drop?w.end+1:tx.boundaries[g+1];
    for(int j=begin;j<end;++j,++row) {
      Interval r{0,0};
      if(j>=w.start && j<=w.end) {
        auto offset=tx.to_genome(qs[i],qe[i],j);
        r=tx.reverse[g]?Interval{tx.ends[j]-offset.end,tx.ends[j]-offset.start}:
                           Interval{tx.starts[j]+offset.start,tx.starts[j]+offset.end};
      }
      out_s[row]=r.start;out_e[row]=r.end;out_group[row]=i+1;
    }
    if(out_partition) out_partition[i]=row;
  }
}
}}
#endif
