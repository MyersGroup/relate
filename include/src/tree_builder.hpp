#ifndef TREE_BUILDER_HPP
#define TREE_BUILDER_HPP

////////////////////////////
// Class for Hierarchical clustering.
/////////////////////////////

#include "fast_log.hpp"
#include "collapsed_matrix.hpp"
#include "data.hpp"
#include "anc.hpp"

#include <sys/types.h>
#include <math.h>
#include <unistd.h>
#include <algorithm>
#include <ctime>
#include <deque>
#include <set>
#include <random>

struct Candidate{
  int lin1 = -1;
  int lin2 = -1;
  double dist = std::numeric_limits<float>::infinity();
  double dist2 = std::numeric_limits<float>::infinity();

  Candidate(){};
  Candidate(int lin1, int lin2, double dist): lin1(lin1), lin2(lin2), dist(dist){};

};

//typedef std::pair<std::list<int>::iterator, std::list<int>::iterator> PairIterator; 
typedef std::list<Candidate> candidates_type;

class MinMatch{

  private:

    std::mt19937 rng;
    int N, L, N_total, Ne;
    float threshold;

    std::vector<int> convert_index; //this will convert the value in cluster_index to the actual index which is between 0 - (2N-1)
    std::vector<float> cluster_size; //size of cluster, accessed using cluster_index
    std::deque<int> cluster_index; //the cluster index will always stay between 0 - N-1 so that I can access the distance matrix.
  
    std::vector<Candidate> mcandidates, mcandidates_sym;
    std::vector<Candidate> candidates_to_check;
    int candidates_to_check_size = 0;

    CollapsedMatrix<float> sym_d;
    float sym_dist, dist_random;
    Candidate best_candidate, best_sym_candidate;

    std::vector<float> min_values, min_values_sym;  //Stores the min values of each row of the distance matrix
    std::vector<int> min_updated; //0: not updated, 1:updated but unchanged, 2: updated and changed
    float min_value_old;

    //only for initial assignment of candidates and min_values
    std::vector<float>::iterator it_min_values_it; //iterator for min_values[i]
    std::vector<float>::iterator it_min_values_jt; //iterator for min_values[j]

    int i, j; //indices of cluster that are being merged
    int conv_i, conv_j; //indices converted back to [N, 2*N-1]
    float added_cluster_size; 
    float dkj, dki; //old values of d[k][j] and d[k][i], respectively
    float djk, dik;
    std::vector<float>::iterator dk_it, dj_it, di_it, d_it; //iterators pointing to beginning of row k, j, i of matrix d
    std::vector<float>::iterator it_d_it_jt; //iterator for d[i][j]

    std::vector<int> updated_cluster;

		void Initialize(CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif);
    void InitializeSym(CollapsedMatrix<float>& sym_d, CollapsedMatrix<float>& d);
		void Coalesce(const int i, const int j, CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif);
    void CoalesceSym(const int i, const int j, CollapsedMatrix<float>& sym_d); 

  public:

    MinMatch(Data& data);
 
    void QuickBuild(CollapsedMatrix<float>& d, Tree& tree);
    void UPGMA(CollapsedMatrix<float>& d, Tree& tree);

};

class InferBranchLengths{

  private:

    std::mt19937 rng;

    int N, L, N_total, Ne;
    float convergence_threshold;

    int root;
    bool is_avg_increasing;
    float max_diff = 0.0, diff; //diff of coordinates
    int count; //count number of iterations
    bool accept; //use to check if proposal in MCMC is accepted.
    float log_likelihood_ratio; //store log of ratio of likelihoods

    Node n;
    int child_left_label, child_right_label;
    
    double tb, tb_new, tb_child_left, tb_child_left_new, tb_child_right, tb_child_right_new;
    float n_num_events, child_left_num_events, child_right_num_events;

    double delta_tau, tau_old, tau_new; //time while k ancestors
    int num_lineages;
    double k_choose_2;

    std::vector<double> avg_prev, avg;
    std::vector<double> coordinates;
    std::vector<double>::iterator it_avg, it_avg_prev, it_coords;

    std::vector<float> logF;
    std::vector<float> mut_rate;

    std::vector<int> sorted_indices; //node indices in order of coalescent events
    std::vector<int> order; //order of coalescent events
    std::vector<int>::iterator it_sorted_indices;

    //update avg
    std::vector<int> last_update;
    std::vector<int>::iterator it_last_update;
    std::vector<double> last_coordinates;
    std::vector<double>::iterator it_last_coords;
    int update_node1 = -1, update_node2 = -1;

    //EM
    std::vector<double> old_branch_length;
    std::deque<int> spanning_branches;

    void InitializeBranchLengths(Tree& tree);
    void InitializeMCMC(const Data& data, Tree& tree);

    float LikelihoodGivenTimes(Tree& tree);
    void UpdateAvg(Tree& tree);
    void logFactorial(int max);

    float ChangeTimeWhilekAncestors(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);
    float ChangeTimeWhilekAncestorsVP(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif);

    void SwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);
    void RandomSwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);

  public:

    InferBranchLengths(const Data& data);

    //this is a post-processing step
    void GetCoordinates(Node& n, std::vector<double>& coords);
    //constant population size MCMC
    void MCMC(const Data& data, Tree& tree, const int seed = std::time(0) + getpid());
    //variable population size MCMC
    void MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed = std::time(0) + getpid());
    //variable population size MCMC, for Relate (i.e. correct initialization etc)
    void MCMCVariablePopulationSizeForRelate(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed = std::time(0) + getpid());
    //variable population size MCMC resample
    void MCMCVariablePopulationSizeSample(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, int num_proposals, const bool init, const int seed = std::time(0) + getpid());
    void EM(const Data& data, Tree& tree, bool called_as_main = false);

};


#endif //TREE_BUILDER_HPP
