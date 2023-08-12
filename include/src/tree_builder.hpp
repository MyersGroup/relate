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
  double dist3 = std::numeric_limits<float>::infinity();
  bool replace = false;
  bool is_from_tmpl = false;

  Candidate(){};
  Candidate(int lin1, int lin2, double dist): lin1(lin1), lin2(lin2), dist(dist){};

  Candidate& operator =(const Candidate& a);

};

//typedef std::pair<std::list<int>::iterator, std::list<int>::iterator> PairIterator; 
typedef std::list<Candidate> candidates_type;

class MinMatch{

  private:

    std::mt19937 rng;
    int N, L, N_total, Ne;
    float threshold, threshold_CF;

    std::vector<int> convert_index; //this will convert the value in cluster_index to the actual index which is between 0 - (2N-1)
    std::vector<float> cluster_size; //size of cluster, accessed using cluster_index
    std::deque<int> cluster_index; //the cluster index will always stay between 0 - N-1 so that I can access the distance matrix.
    std::vector<int> conv_to_tmpl;

    std::vector<Candidate> mcandidates, mcandidates_sym;
    std::vector<Candidate> candidates_to_check;
    int candidates_to_check_size = 0;

    CollapsedMatrix<float> sym_d;
    float sym_dist, dist_random;
    double age;
    Candidate best_candidate, best_sym_candidate, cand;
    std::vector<double> unique_sample_ages;
    std::vector<int> sample_ages_count;

    std::vector<float> min_values, min_values_sym, min_values_CF;  //Stores the min values of each row of the distance matrix
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
    int num_nodes;
    std::vector<float>::iterator dk_it, dj_it, di_it, d_it; //iterators pointing to beginning of row k, j, i of matrix d
    std::vector<float>::iterator it_d_it_jt; //iterator for d[i][j]

    std::vector<int> updated_cluster;

		void Initialize(CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, Tree *tmpl_tree = NULL);
		void Initialize(CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages, Tree *tmpl_tree = NULL);
    void InitializeSym(CollapsedMatrix<float>& sym_d, CollapsedMatrix<float>& d);
		void Coalesce(const int i, const int j, CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, Tree *tmpl_tree = NULL);
		void Coalesce(const int i, const int j, CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages, Tree *tmpl_tree = NULL);
    void CoalesceSym(const int i, const int j, CollapsedMatrix<float>& sym_d); 


    void Initialize(CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif);
    void Initialize(CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages);
    void Coalesce(const int i, const int j, CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif);
    void Coalesce(const int i, const int j, CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages);

  public:

    MinMatch(Data& data);
 
    void QuickBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& sample_ages, Tree *tmpl_tree = NULL);   
    void QuickBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& sample_ages, const CollapsedMatrix<float>& d_prior);
    
    void TestBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& sample_ages);
    void SlowBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& sample_ages);
    void UPGMA(CollapsedMatrix<float>& d, Tree& tree);

};

#endif //TREE_BUILDER_HPP
