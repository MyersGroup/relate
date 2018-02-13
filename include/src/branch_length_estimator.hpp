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
    int count, count_accept, count_proposal; //count number of iterations
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

    void InitializeBranchLengths(Tree& tree);
    void InitializeMCMC(const Data& data, Tree& tree);

    void UpdateAvg(Tree& tree);

    void ChangeTimeWhilekAncestors(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);
    void ChangeTimeWhilekAncestorsVP(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif);
    void SwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);
    void RandomSwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif);

  public:

    EstimateBranchLengths(const Data& data);

    //this is a post-processing step
    void GetCoordinates(Node& n, std::vector<double>& coords);
    void MCMC(const Data& data, Tree& tree, const int seed = std::time(0) + getpid());
    void MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed = std::time(0) + getpid());

};

