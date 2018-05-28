#ifndef MUTATIONS_HPP
#define MUTATIONS_HPP

#include "gzstream.hpp"
#include "data.hpp"
#include "anc.hpp"

#include <iostream>
#include <deque>

struct SNPInfo{

  std::string rs_id;
  int snp_id;
  int pos, dist;

  int tree;
  std::deque<int> branch;
  std::vector<int> freq;
  bool flipped = false;
  float age_begin = 0.0, age_end = 0.0;

  std::string upstream_base = "NA", downstream_base = "NA";
  std::string mutation_type = "NA";

};


class Mutations{

  private:

    int N;
    int L;

    int num_flips, num_notmappingmutations;

  public:

    std::string header;
    std::vector<SNPInfo> info;

    Mutations(){};
    Mutations(Data& data);
    void Init(Data& data);

    //////////////////////////////////////////////////////

    void GetAge(AncesTree& anc);
   
    void Read(igzstream& is);
    void Read(const std::string& filename);
    void Dump(const std::string& filename);


    void ReadShortFormat(const std::vector<std::string>& filenames);
    void DumpShortFormat(const std::string& filename);
    void DumpShortFormat(const std::string& filename, const int section_startpos, const int section_endpos);

    int GetNumFlippedMutations(){return num_flips;}
    int GetNumNotMappingMutations(){return num_notmappingmutations;}

};


#endif //MUTATIONS_HPP 
