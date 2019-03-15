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


class EstimateBranchLengths{

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
    int num_lineages;
    double k_choose_2;

    double tau_below, tau_above, T, x1, x2, tau_new_below, tau_new_above;

    std::vector<double> avg_prev, avg;
    std::vector<double> coordinates;
    std::vector<double>::iterator it_avg, it_avg_prev, it_coords;

    std::vector<float> logt_pos, logt_neg;
    std::vector<float> mut_rate;

    std::vector<int> sorted_indices; //node indices in order of coalescent events
    std::vector<int> order; //order of coalescent events
    std::vector<int>::iterator it_sorted_indices;

    //update avg
    std::vector<int> last_update;
    std::vector<int>::iterator it_last_update;
    std::vector<double> last_coordinates;
    std::vector<double>::iterator it_last_coords;
    int update_node1 = -1, update_node2 = -1, update_node3 = -1;
 
    float log_deltat(float t);

    void InitializeBranchLengths(Tree& tree);
    void InitializeMCMC(const Data& data, Tree& tree);

    void UpdateAvg(Tree& tree);

    void UpdateOneEvent(Tree& tree, int k, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif);
    void ChangeTimeWhilekAncestors(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);

    void UpdateOneEventVP(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif);
    void ChangeTimeWhilekAncestorsVP(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif);
    void ChangeTimeWhilekAncestorsVP_new(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif);


    void SwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);
    void RandomSwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);
    void InitialiseEventOrder(Tree& tree, std::uniform_real_distribution<double>& dist_unif);
    void InitialiseEventOrder2(Tree& tree, std::uniform_real_distribution<double>& dist_unif);

  public:

    EstimateBranchLengths(const Data& data);

    //this is a post-processing step
    void GetCoordinates(Node& n, std::vector<double>& coords);
    //constant population size MCMC
    void MCMC(const Data& data, Tree& tree, const int seed = std::time(0) + getpid());
    //variable population size MCMC
    void MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed = std::time(0) + getpid());
    //variable population size MCMC, for Relate (i.e. correct initialization etc)
    void MCMCVariablePopulationSizeForRelate(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed = std::time(0) + getpid());

};

