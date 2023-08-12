#ifndef BRANCH_LENGTH_ESTIMATOR_HPP
#define BRANCH_LENGTH_ESTIMATOR_HPP

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

class EstimateBranchLengthsWithSampleAge{

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
    int child_left_label, child_right_label, parent_label;

    double tb, tb_new, tb_child_left, tb_child_left_new, tb_child_right, tb_child_right_new;
    float n_num_events, child_left_num_events, child_right_num_events;

    double delta_tau, tau_old, tau_new; //time while k ancestors
    double k_choose_2;

    double tau_below, tau_above, T, x1, x2, tau_new_below, tau_new_above;

    std::vector<double> avg_prev, avg;
    std::vector<double> coordinates;
    std::vector<double>::iterator it_avg, it_avg_prev, it_coords;

    std::vector<float> logt_pos, logt_neg;
    std::vector<float> mut_rate;
    std::vector<double> sample_age;
  
    std::vector<int> num_lineages, num_lineages_new;
    std::vector<int> sorted_indices, sorted_indices_new; //node indices in order of coalescent events
    std::vector<int> order, order_new; //order of coalescent events
    std::vector<std::vector<int>> remaining, remaining_new;
    std::vector<int>::iterator it_sorted_indices, it_order;

    //update avg
    std::vector<int> last_update;
    std::vector<int>::iterator it_last_update;
    std::vector<double> last_coordinates;
    std::vector<double>::iterator it_last_coords;
    int update_node1 = -1, update_node2 = -1, update_node3 = -1;

    float log_deltat(float t);

    void InitializeBranchLengths(Tree& tree);
    void InitializeOrder(Tree& tree);
    void InitializeMCMC(const Data& data, Tree& tree);

    void UpdateAvg(Tree& tree);

    //calculate coalescent prior TODO: optimize these functions
    double CalculatePrior(std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages);
    double CalculatePrior(int k_start, int k_end, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages);
    double CalculatePrior(const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages);
    double CalculatePrior(int k_start, int k_end, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages);
    double CalculatePrior(const Tree& tree, const std::vector<double>& epoch, const std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages);
    double CalculatePrior(int k_start, int k_end, const Tree& tree, const std::vector<double>& epoch, const std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages);


    //constant Ne
    void UpdateOneEvent(Tree& tree, int node_k, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif);

    //variable Ne
    void UpdateOneEventVP(Tree& tree, int node_k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif);
    void UpdateOneEventVP(Tree& tree, int node_k, const std::vector<double>& epoch, const std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif);


    void SwitchOrder(Tree& tree, int node_k, std::uniform_real_distribution<double>& dist_unif);
    void SwitchTopo(Tree& tree, std::vector<Leaves>& desc, const std::vector<double>& epoch, std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, int node_k, std::uniform_real_distribution<double>& dist_unif);

    void RandomSwitchOrder(Tree& tree, int node_k, std::uniform_real_distribution<double>& dist_unif);

  public:

    EstimateBranchLengthsWithSampleAge(const Data& data);
    EstimateBranchLengthsWithSampleAge(const Data& data, const std::vector<double>& sample_age_input);

    //this is a post-processing step
    void GetCoordinates(Node& n, std::vector<double>& coords);
    //constant population size MCMC
    void MCMC(const Data& data, Tree& tree, const int seed = std::time(0) + getpid());
    //variable Ne MCMC
    void MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed = std::time(0) + getpid());
    void MCMCVariablePopulationSizeForRelate(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed = std::time(0) + getpid());
    void MCMCCoalRatesForRelate(const Data& data, Tree& tree, const std::vector<int>& membership, const std::vector<double>& epoch, std::vector<std::vector<std::vector<double>>>& coal_rate, const int seed = std::time(0) + getpid());

    void MCMCCoalRatesSample(const Data& data, Tree& tree, const std::vector<int>& membership, const std::vector<double>& epoch, std::vector<std::vector<std::vector<double>>>& coal_rate, int num_proposals, const bool init, const int seed = std::time(0) + getpid());
    void MCMCVariablePopulationSizeSample(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, int num_proposals, const bool init, const int seed = std::time(0) + getpid());

};

#endif //BRANCH_LENGTH_ESTIMATOR_HPP
