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

typedef std::vector<SNPInfo> Muts; //I could change this to deque but deque has no splice

class Mutations{

  private:

    int N;
    int L;

    int num_flips, num_notmappingmutations;

  public:

    std::string header;
    Muts info;

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

class AncMutIterators{

  private:

    igzstream is;
    Muts::iterator pit_mut;
    Mutations mut;

    int N, num_trees;
    int tree_index_in_anc, tree_index_in_mut;
    double num_bases_tree_persists;
    std::string line, filename_anc, filename_mut;

  public:

    std::vector<double> sample_ages;
    AncMutIterators(){};
    AncMutIterators(const std::string& filename_anc, const std::string& filename_mut);

    void OpenFiles(const std::string& i_filename_anc, const std::string& i_filename_mut);
    void CloseFiles(){
      if(is.rdbuf() -> is_open()) is.close(); //close if stream is still open
    }
    int NumTips(){
      return(N);
    }
    int NumSnps(){
      return(mut.info.size());
    }
    int NumTrees(){
      return(num_trees);
    }

    Muts::iterator mut_begin(){
      return(mut.info.begin());
    }
    Muts::iterator mut_end(){
      return(mut.info.end());
    }
    int get_treecount(){
      return(tree_index_in_anc);
    }

    double NextTree(MarginalTree& mtr, Muts::iterator& it_mut);
    double FirstSNP(MarginalTree& mtr, Muts::iterator& it_mut);
    double NextSNP(MarginalTree& mtr, Muts::iterator& it_mut);


};


#endif //MUTATIONS_HPP 
