#ifndef ARG_BUILDER_HPP
#define ARG_BUILDER_HPP

////////////////////////////
// Class for building tree sequences.
/////////////////////////////

#include "fast_log.hpp"
#include "collapsed_matrix.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"
#include "fast_painting.hpp"
#include "tree_builder.hpp"

#include <utility> //for std::pair
#include <limits>
#include <utility>
#include <tgmath.h>
#include <random>
#include <cassert>

//////////////////////////////////////////////////////

class Correlation{

  private:

    int N;
    float N_float;

  public:

    Correlation(int N): N(N){
      N_float = (float) N;
    }
    float Pearson(const Leaves& set1, const Leaves& set2);

};

struct PropagateStructLocal{

  int num_carriers = 0;
  int num_flipped_carriers = 0;
  int best_branch = -1;
  int best_flipped_branch = -1;

};

struct PropagateStructGlobal{

  int num_correct_carriers, num_correct_noncarriers;
  int num_incorrect_carriers, num_incorrect_noncarriers;
  int best_branch, best_flipped_branch;
  int min, flipped_min;

};

struct EquivalentNode{

  int node1;
  int node2;
  float corr;

  EquivalentNode(){};
  EquivalentNode(int node1, int node2, float corr): node1(node1), node2(node2), corr(corr){};

  bool operator > (const EquivalentNode& n) const{
    return (corr > n.corr);
  }

};

//////////////////////////////////////////////////////


//This class is used in AncesTreeBuilder to obtain the actual distance matrix
//as a weighted mean of the distance matrices on SNPs around the current
//SNP. (We have to do this because we don't paint every SNP, rather we only
//paint SNPs at which the sequence that is being painted has a mutation.)
class DistanceMeasure{

  private:

    int N, L;
    int section;
    int section_startpos, section_endpos;
    const float scale = -1.0;

    std::vector<CollapsedMatrix<float>> top;
    std::vector<CollapsedMatrix<float>>* topology; //This is inefficient, but not by much I think
    std::vector<std::vector<float>> log;
    std::vector<std::vector<float>>* logscales;
    Data* data;

  public:

    CollapsedMatrix<float> matrix;
    std::vector<int> v_snp_prev; 
    std::vector<double> v_rpos_prev, v_rpos_next;

    DistanceMeasure(Data& idata, int section): section(section){
      data      = &idata; 
      N    = idata.N;
      L    = idata.L;

      v_snp_prev.resize(N,0);
      v_rpos_prev.resize(N);
      v_rpos_next.resize(N);
      matrix.resize(N, N);
      top.resize(N);
      log.resize(N);

      section_startpos = -1;
      section_endpos   = -1;
    }

		void Assign(std::vector<CollapsedMatrix<float>>& itop, std::vector<std::vector<float>>& ilog, const int isection_startpos, const int isection_endpos, const int snp);
    void GetTopologyWithRepaint(const int snp); //call this function whenever data.pos[snp] does not lie within [section_startpos, section_endpos]
    void GetMatrix(const int snp);

};

class AncesTreeBuilder{

  private:

    int N, N_total, root;
    int L;
    int section;

    //for mapping mutations
    int thr; //threshold for mismatches when placing mutations on tree.
    
    //for finding equivalent branches
    float threshold_brancheq; //threshold for R2 value for regarding two branches as equivalent
    std::vector<std::vector<int>> potential_branches;

    //for random flipping
    std::mt19937 rng;

    //////////////////
    //A few helper functions used for building AncesTrees 
    int MapMutation(Tree& tree, Leaves& sequences_carrying_mutations, std::uniform_real_distribution<double>& dist_unif, const int snp, float& min_value);
    int MapMutation(Tree& tree, Leaves& sequences_carrying_mutations, const int snp, float& min_value);
    int ForceMapMutation(Tree& tree, Leaves& sequences_carrying_mutations, const int snp, const bool force = false);
    
    int PropagateMutationExact(Node& node, std::deque<int>& branches, std::deque<int>& branches_flipped, Leaves& sequences_carrying_mutations);
    void PropagateMutationGlobal(Node& node, Leaves& sequences_carrying_mutations, PropagateStructGlobal& report);
    void PropagateMutationLocal(Node& node, std::deque<int>& branches, std::deque<int>& branches_flipped, Leaves& sequences_carrying_mutations, PropagateStructLocal& report);

    void UpdateBranchSNPbegin(Tree& tree, int snp);
    void UpdateBranchSNPend(Tree& tree, int snp);

  public:

    Mutations mutations;
    AncesTreeBuilder(Data& data);

    /////////////////////////////////////////////////////
    //This is using the painting to calculate a distance measure
    //and then the MinMatch hierarchical clustering algorithm to
    //build tree sequences that adapt whenever a recombination
    //event is detected.

    void BuildTopology(const int section, const int section_startpos, const int section_endpos, Data& data, AncesTree& anc, const int seed, const bool ancestral_state);
    void AssociateTrees(std::vector<AncesTree>& v_anc, const std::string& dirname = "./");
		int OptimizeParameters(const int section, const int section_startpos, const int section_endpos, Data& data, const int seed);

    void PreCalcPotentialBranches();
    void BranchAssociation(Tree& ref_tree, Tree& tree, std::vector<int>& equivalent_branches);

    int IsSNPMapping(Tree& tree, Leaves& sequences_carrying_mutations, int snp){
      float min_value;
      if(MapMutation(tree, sequences_carrying_mutations, snp, min_value) > 2){
        ForceMapMutation(tree, sequences_carrying_mutations, snp, true);
        return 2;
      }else{
        return 1;
      }
    }

    void ForceMap(Tree& tree, Leaves& sequences_carrying_mutations, int snp){
      ForceMapMutation(tree, sequences_carrying_mutations, snp, true);
    }

};


#endif //ARG_BUILDER_HPP
