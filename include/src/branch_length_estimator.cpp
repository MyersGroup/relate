#include "branch_length_estimator.hpp"


//////////////////////////////////////////

EstimateBranchLengthsWithSampleAge::EstimateBranchLengthsWithSampleAge(const Data& data, const std::vector<double>& sample_age_input){
  N       = data.N;
  N_total = 2*N - 1;
  L       = data.L;
  Ne      = data.Ne;

  logt_pos.resize(10000);
  for(int i = 0; i < 10000; i++){
    logt_pos[i] = log(1.0 + 0.0001*((float) i));
  }
  logt_neg.resize(1000);
  for(int i = 0; i < 1000; i++){
    logt_neg[i] = log(1.0 - 0.0001*((float)i));
  }

  assert(sample_age_input.size() == N);
  sample_age.resize(N);
  std::vector<double>::const_iterator it_sample_age_input = sample_age_input.begin();
  for(std::vector<double>::iterator it_sample_age = sample_age.begin(); it_sample_age != sample_age.end(); it_sample_age++){
    *it_sample_age = *it_sample_age_input/data.Ne;
    it_sample_age_input++;
  }

  num_lineages.resize(N_total);
  coordinates.resize(N_total);
  sorted_indices.resize(N_total); //node indices in order of coalescent events
  order.resize(N_total); //order of coalescent events
};

//MCMC
void
EstimateBranchLengthsWithSampleAge::InitializeBranchLengths(Tree& tree){

  int node_i, num_lins;

  num_lins = 0;
  double ages = sample_age[*sorted_indices.begin()];
  std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
  for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    if(*it_sorted_indices >= N){
      for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        num_lineages[*it_sorted_indices_start] = num_lins; 
      }
      num_lins--;
      num_lineages[*it_sorted_indices] = num_lins;
      it_sorted_indices_start++;
    }else if(ages < sample_age[*it_sorted_indices]){      
      for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        num_lineages[*it_sorted_indices_start] = num_lins; 
      }
      ages = sample_age[*it_sorted_indices];
      num_lins++;
    }else{
      num_lins++;
    }
  }
  num_lineages_new = num_lineages;

  //initialize using coalescent prior
  coordinates.resize(N_total);
  std::fill(coordinates.begin(), coordinates.end(), 0.0);
  for(int i = 0; i < N; i++){
    coordinates[i] = sample_age[i];
  }

  //for each node determine upper limit of age
  int j = 1, i = 1;
  double age_upper = coordinates[sorted_indices[0]];
  for(i = 1; i < N_total; i++){
    if(sorted_indices[i] < N){
      age_upper = coordinates[sorted_indices[i]];
      for(; j < i; j++){
        assert(sorted_indices[j] >= N);
        coordinates[sorted_indices[j]] = age_upper;
      }
      j = i+1;
    }
  }

  for(int i = 0; i < N_total; i++){
    node_i = sorted_indices[i];
    if(node_i >= N){
      num_lins = num_lineages[sorted_indices[i-1]];
      assert(num_lins > 0);
      if(coordinates[node_i] > 0){
        double tmp = coordinates[node_i];
        assert(tmp >= coordinates[sorted_indices[i-1]]);
        coordinates[node_i] = (tmp + coordinates[sorted_indices[i-1]])/10.0;
        //coordinates[node_i] = (tmp + coordinates[sorted_indices[i-1]])/sorted_indices.size();
        //coordinates[node_i] = coordinates[sorted_indices[i-1]];
        assert(coordinates[node_i] <= tmp);
      }else{
        coordinates[node_i] = coordinates[sorted_indices[i-1]] + 2.0/( num_lins * (num_lins - 1.0) ); // determined by the prior
      }
      (*tree.nodes[node_i].child_left).branch_length  = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_left).label];
      (*tree.nodes[node_i].child_right).branch_length = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_right).label];
    }
  }

  for(int i = 0; i < N_total-1; i++){
    //std::cerr << coordinates[sorted_indices[i+1]] << " ";
    assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
  }
  //std::cerr << std::endl;

}

void
EstimateBranchLengthsWithSampleAge::InitializeOrder(Tree& tree){

  //initialize
  //1. sort coordinate vector to obtain sorted_indices
  //2. sort sorted_indices to obtain order

  //strategy:
  //assign pseudo coordinates to nodes, by giving them a lower bound on age + epsilon
  //then use code as if coordinates are given
  std::vector<double> pseudo_coords(N_total, 0.0);
  double epsilon = 1.0/log(N);
  //double epsilon = 1/N_total;
  //double epsilon = 0.0;
  for(int i = 0; i < N; i++){
    pseudo_coords[i] = sample_age[i];
    int k1 = i, k2 = i;
    while(k2 < root){
      k1 = k2;
      k2 = (*tree.nodes[k2].parent).label;
      if(pseudo_coords[k2] < pseudo_coords[k1] + epsilon){
        pseudo_coords[k2] = std::nextafter(pseudo_coords[k1]+epsilon, pseudo_coords[k1]+epsilon+1);
      }
    }
  } 

  std::size_t m1(0);
  std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
  std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
      return std::tie(pseudo_coords[i1],i1) < std::tie(pseudo_coords[i2],i2); });

  //obtain order of coalescent events
  std::fill(order.begin(), order.end(), 0);
  std::size_t m2(0);
  std::generate(std::begin(order), std::end(order), [&]{ return m2++; });
  std::sort(std::begin(order), std::end(order), [&](int i1, int i2) { return sorted_indices[i1] < sorted_indices[i2]; } );

  ////////////////////////////////

  if(0){
    int num_lins = 0;
    double age = pseudo_coords[*sorted_indices.begin()];
    std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
    for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      if(pseudo_coords[*it_sorted_indices] > age){
        for(; it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages[*it_sorted_indices_start] = num_lins;          
        }
        age = pseudo_coords[*it_sorted_indices_start];
      }
      if(*it_sorted_indices < N){
        num_lins++;
      }else{
        num_lins--;
      }
      assert(num_lins >= 1);

      if(it_sorted_indices != sorted_indices.begin()){
        assert(pseudo_coords[*it_sorted_indices] >= pseudo_coords[*std::prev(it_sorted_indices,1)]);
      }
    }
  }

  sorted_indices_new = sorted_indices;
  order_new          = order;
  //num_lineages_new   = num_lineages;

  //debug
  for(int i = 0; i < N_total-1; i++){
    assert(order[sorted_indices[i]] == i);
    assert(order[i] < order[(*tree.nodes[i].parent).label]);
  }

}

void
EstimateBranchLengthsWithSampleAge::InitializeMCMC(const Data& data, Tree& tree){

  mut_rate.resize(N_total);
  for(int i = 0; i < N_total; i++){
    int snp_begin = tree.nodes[i].SNP_begin;
    int snp_end   = tree.nodes[i].SNP_end;

    assert(snp_end < data.pos.size());
    //if(snp_end >= data.pos.size()) snp_end = data.pos.size()-1;
    mut_rate[i]            = 0.0;
    for(int snp = snp_begin; snp < snp_end; snp++){
      mut_rate[i]         += data.pos[snp];
    }

    if(snp_begin > 0){
      snp_begin--;
      mut_rate[i]         += 0.5 * data.pos[snp_begin];
    }
    if(snp_end < data.L-1){
      mut_rate[i]         += 0.5 * data.pos[snp_end];
    }

    mut_rate[i]           *= data.Ne * data.mu;
  }

  //initialize
  //1. sort coordinate vector to obtain sorted_indices
  //2. sort sorted_indices to obtain order

  order.resize(N_total);
  sorted_indices.resize(N_total);

}

void
EstimateBranchLengthsWithSampleAge::UpdateAvg(Tree& tree){


  if(update_node1 != -1){

    if(update_node2 != -1){

      it_avg         = std::next(avg.begin(), update_node1);
      it_coords      = std::next(coordinates.begin(), update_node1);
      it_last_update = std::next(last_update.begin(), update_node1);
      it_last_coords = std::next(last_coordinates.begin(), update_node1);
      *it_avg         += ((count - *it_last_update) * (*it_last_coords - *it_avg) + *it_coords - *it_last_coords)/count;
      *it_last_update  = count;
      *it_last_coords  = *it_coords;

      it_avg         = std::next(avg.begin(), update_node2);
      it_coords      = std::next(coordinates.begin(), update_node2);
      it_last_update = std::next(last_update.begin(), update_node2);
      it_last_coords = std::next(last_coordinates.begin(), update_node2);
      *it_avg         += ((count - *it_last_update) * (*it_last_coords - *it_avg) + *it_coords - *it_last_coords)/count;
      *it_last_update  = count;
      *it_last_coords  = *it_coords;

      update_node1 = -1;
      update_node2 = -1;

    }else if(update_node3 != -1){

      it_avg         = std::next(avg.begin(), update_node1);
      it_coords      = std::next(coordinates.begin(), update_node1);
      it_last_update = std::next(last_update.begin(), update_node1);
      it_last_coords = std::next(last_coordinates.begin(), update_node1);
      *it_avg         += ((count - *it_last_update) * (*it_last_coords - *it_avg) + *it_coords - *it_last_coords)/count;
      *it_last_update  = count;
      *it_last_coords  = *it_coords;

      update_node1 = -1;
      update_node3 = -1;

    }else{

      for(std::vector<int>::iterator it_node = std::next(sorted_indices.begin(), update_node1); it_node != sorted_indices.end(); it_node++){
        avg[*it_node]             += ((count - last_update[*it_node]) * (last_coordinates[*it_node] - avg[*it_node]) + coordinates[*it_node] - last_coordinates[*it_node])/count;
        last_update[*it_node]      = count;
        last_coordinates[*it_node] = coordinates[*it_node];
      }
      update_node1 = -1;

    }

  }

}

float
EstimateBranchLengthsWithSampleAge::log_deltat(float t){
  if(t >= 0){
    //return(fast_log(1+t));
    if(t < 1){
      return(logt_pos[(int) (t*10000)]);
    }else{
      return(fast_log(1.0+t));
    }  
  }else{
    //return(fast_log(1+t));
    if(t > -0.1){
      return(logt_neg[(int) (-t*10000)]);
    }else{
      return(fast_log(1.0+t));
    }  
  }
}


//////////////////////////////////////////

//This function switches the order of coalescent events randomly for initialisation.
void
EstimateBranchLengthsWithSampleAge::RandomSwitchOrder(Tree& tree, int node_k, std::uniform_real_distribution<double>& dist_unif){

  int node_swap_k;

  //check order of parent and order of children
  //choose a random number in this range (not including)
  //swap order with that node

  //node_k = sorted_indices[k];

  //cannot swap a tip
  int k = order[node_k];
  if(node_k >= N){

    int parent_order    = order[(*tree.nodes[node_k].parent).label];
    int child_order     = order[(*tree.nodes[node_k].child_left).label];
    int child_order_alt = order[(*tree.nodes[node_k].child_right).label];
    if(child_order < child_order_alt) child_order = child_order_alt;
    //if(child_order < N) child_order = N-1;
    assert(child_order < parent_order);

    if( parent_order - child_order > 2 ){

      std::uniform_int_distribution<int> d_swap(child_order+1, parent_order-1);
      int new_order = d_swap(rng);

      //this can't be a sample
      if(sorted_indices[new_order] >= N){

        //swap
        node_swap_k     = sorted_indices[new_order];
        parent_order    = order[(*tree.nodes[node_swap_k].parent).label];
        child_order     = order[(*tree.nodes[node_swap_k].child_left).label];
        child_order_alt = order[(*tree.nodes[node_swap_k].child_right).label];
        if(child_order < child_order_alt) child_order = child_order_alt;
        if(child_order < N) child_order = N-1;

        if(child_order < k && k < parent_order){ 

          if(new_order != k){
            sorted_indices[k]         = node_swap_k;
            sorted_indices[new_order] = node_k;
            order[node_k]             = new_order;
            order[node_swap_k]        = k;
            //int num_lins = num_lineages[node_k];
            //num_lineages[node_k]      = num_lineages[node_swap_k];
            //num_lineages[node_swap_k] = num_lins;
          }

        }

      }

    }

  }

}

void
EstimateBranchLengthsWithSampleAge::SwitchOrder(Tree& tree, int node_k, std::uniform_real_distribution<double>& dist_unif){

  //Idea:
  //
  //For node_k = sorted_indices[k], look up order of parent and children
  //Then, choose a node with order between children and parent uniformly at random.
  //Check whether chosen node has order of children < k and order of parent > k
  //If so, this is a valid transition - calculate likelihood ratio and accept transition with this probability
  //Otherwise, reject transition immediately
  //
  //Notice that the transition probabilities are symmetric: 
  //Distribution d_swap(child_order+1, parent_order-1); does not change after a transition.
  //

  //This update is constant time.
  accept = true;
  int node_swap_k;
  log_likelihood_ratio = 0.0;

  //check order of parent and order of children
  //choose a random number in this range (not including)
  //swap order with that node

  //node_k = sorted_indices[k];

  int k = order[node_k];
  if(node_k >= N){

    int parent_order    = order[(*tree.nodes[node_k].parent).label];
    int child_order     = order[(*tree.nodes[node_k].child_left).label];
    int child_order_alt = order[(*tree.nodes[node_k].child_right).label];
    if(child_order < child_order_alt) child_order = child_order_alt;
    //if(child_order < N) child_order = N-1;

    //std::cerr << k << " " << node_k << " " << child_order << " " << parent_order << std::endl;
    //std::cerr << coordinates[node_k] << " " << coordinates[sorted_indices[parent_order]] << " " << coordinates[sorted_indices[child_order]] << " " << coordinates[sorted_indices[child_order_alt]] << std::endl;
    assert(child_order < parent_order);

    if( parent_order - child_order > 2 ){

      std::uniform_int_distribution<int> d_swap(child_order+1, parent_order-1);
      int new_order = d_swap(rng);

      if(sorted_indices[new_order] >= N){

        //calculate odds
        node_swap_k     = sorted_indices[new_order];
        parent_order    = order[(*tree.nodes[node_swap_k].parent).label];
        child_order     = order[(*tree.nodes[node_swap_k].child_left).label];
        child_order_alt = order[(*tree.nodes[node_swap_k].child_right).label];
        if(child_order < child_order_alt) child_order = child_order_alt;
        //if(child_order < N) child_order = N-1;

        //std::cerr << node_k << " " << node_swap_k << " " << N << std::endl;
        if(child_order < k && k < parent_order){

          //branch length of node and children change
          delta_tau = coordinates[node_swap_k] - coordinates[node_k];
          //assert(delta_tau >= 0.0);

          child_left_label  = (*tree.nodes[node_k].child_left).label;
          child_right_label = (*tree.nodes[node_k].child_right).label;

          n_num_events           = tree.nodes[node_k].num_events;
          child_left_num_events  = tree.nodes[child_left_label].num_events;
          child_right_num_events = tree.nodes[child_right_label].num_events;

          tb                 = tree.nodes[node_k].branch_length;
          tb_child_left      = tree.nodes[child_left_label].branch_length;
          tb_child_right     = tree.nodes[child_right_label].branch_length;

          //mutation and recombination part
          if(tb == 0.0){
            log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          }else if(tb <= delta_tau){
            log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          }else{

            if(tb_child_left == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
            }else if(tb_child_left <= -delta_tau){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
            }else{

              if(tb_child_right == 0.0){
                log_likelihood_ratio  = std::numeric_limits<float>::infinity();
              }else if(tb_child_right <= -delta_tau){
                log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
              }else{

                log_likelihood_ratio += (mut_rate[node_k] - mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
                if(n_num_events >= 0.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
                if(child_right_num_events >= 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
                if(child_left_num_events >= 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

                delta_tau  *= -1.0;
                child_left_label  = (*tree.nodes[node_swap_k].child_left).label;
                child_right_label = (*tree.nodes[node_swap_k].child_right).label;

                n_num_events           = tree.nodes[node_swap_k].num_events;
                child_left_num_events  = tree.nodes[child_left_label].num_events;
                child_right_num_events = tree.nodes[child_right_label].num_events;

                tb                 = tree.nodes[node_swap_k].branch_length;
                tb_child_left      = tree.nodes[child_left_label].branch_length;
                tb_child_right     = tree.nodes[child_right_label].branch_length;

                if(tb == 0.0){
                  log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                }else if(tb <= delta_tau){
                  log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
                }else{

                  if(tb_child_left == 0.0){
                    log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                  }else if(tb_child_left <= -delta_tau){
                    log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
                  }else{

                    if(tb_child_right == 0.0){
                      log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                    }else if(tb_child_right <= -delta_tau){
                      log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
                    }else{ 
                      log_likelihood_ratio += (mut_rate[node_swap_k] - mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
                      if(n_num_events >= 0.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
                      if(child_right_num_events >= 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
                      if(child_left_num_events >= 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);
                    }

                  }
                }

              }
            }
          } 
          assert(!std::isnan(log_likelihood_ratio));

          accept = true;
          if(log_likelihood_ratio < 0.0){
            //accept with probability exp(log_likelihood_ratio)
            if(dist_unif(rng) > exp(log_likelihood_ratio)){
              accept = false;
            }
          }

          //update coordinates and sorted_indices
          if(accept && new_order != k){
            //order of nodes in k - new_order decreases by one.

            sorted_indices[k]         = node_swap_k;
            sorted_indices[new_order] = node_k;
            order[node_k]             = new_order;
            order[node_swap_k]        = k;
            int num_lins = num_lineages[node_k];
            num_lineages[node_k]      = num_lineages[node_swap_k];
            num_lineages[node_swap_k] = num_lins;

            double tmp_coords                       = coordinates[node_k];
            coordinates[node_k]                     = coordinates[node_swap_k];
            coordinates[node_swap_k]                = tmp_coords;
            update_node1 = node_k;
            update_node2 = node_swap_k;

            //calculate new branch lengths
            tree.nodes[node_k].branch_length                 = coordinates[(*tree.nodes[node_k].parent).label]  - coordinates[node_k];
            if(tree.nodes[node_k].branch_length < 0.0) tree.nodes[node_k].branch_length = 0.0;
            assert(tree.nodes[node_k].branch_length >= 0.0); 
            child_left_label                                 = (*tree.nodes[node_k].child_left).label;
            tree.nodes[child_left_label].branch_length       = coordinates[node_k] - coordinates[child_left_label];
            if(tree.nodes[child_left_label].branch_length < 0.0) tree.nodes[child_left_label].branch_length = 0.0;
            assert(tree.nodes[child_left_label].branch_length >= 0.0);
            child_right_label                                = (*tree.nodes[node_k].child_right).label;
            tree.nodes[child_right_label].branch_length      = coordinates[node_k] - coordinates[child_right_label];
            if(tree.nodes[child_right_label].branch_length < 0.0) tree.nodes[child_right_label].branch_length = 0.0;
            assert(tree.nodes[child_right_label].branch_length >= 0.0);

            tree.nodes[node_swap_k].branch_length            = coordinates[(*tree.nodes[node_swap_k].parent).label]  - coordinates[node_swap_k];
            if(tree.nodes[node_k].branch_length < 0.0) tree.nodes[node_k].branch_length = 0.0;
            assert(tree.nodes[node_swap_k].branch_length >= 0.0); 
            child_left_label                                 = (*tree.nodes[node_swap_k].child_left).label;
            tree.nodes[child_left_label].branch_length       = coordinates[node_swap_k] - coordinates[child_left_label];
            if(tree.nodes[child_left_label].branch_length < 0.0) tree.nodes[child_left_label].branch_length = 0.0;
            assert(tree.nodes[child_left_label].branch_length >= 0.0);
            child_right_label                                = (*tree.nodes[node_swap_k].child_right).label;
            tree.nodes[child_right_label].branch_length      = coordinates[node_swap_k] - coordinates[child_right_label];
            if(tree.nodes[child_right_label].branch_length < 0.0) tree.nodes[child_right_label].branch_length = 0.0;
            assert(tree.nodes[child_right_label].branch_length >= 0.0);

          }
        }

      }
    }

  }

}

//This changes the time while k ancestors remain. Want to only make a constant number of updates for each k
void
EstimateBranchLengthsWithSampleAge::ChangeTimeWhilekAncestors(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif){

  //This step is O(N)
  k_choose_2 = num_lineages[sorted_indices[k]] * (num_lineages[sorted_indices[k]]-1.0)/2.0;
  tau_old    = coordinates[sorted_indices[k]] - coordinates[sorted_indices[k-1]];

  log_likelihood_ratio = 0.0; 
  if(tau_old > 0.0){
    //choose tau_new according to Gamma(alpha, alpha/tau_old)
    tau_new   = -fast_log(dist_unif(rng)) * tau_old;
    delta_tau = tau_new - tau_old;
    //calculate ratio of proposals
    log_likelihood_ratio = fast_log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);
    assert(tau_new > 0.0);
  }else{
    //choose tau_new according to Gamma(alpha, alpha/tau_old) 
    tau_new   = -fast_log(dist_unif(rng))/k_choose_2;
    tau_old   = 0.0;
    delta_tau = tau_new;
    //calculate ratio of proposals
    log_likelihood_ratio = fast_log(1.0/(tau_new*k_choose_2)) + tau_new*k_choose_2;
    assert(tau_new > 0.0);
  }

  //coalescent prior
  log_likelihood_ratio -= k_choose_2 * delta_tau;

  //assert(order[node_k] == k);
  int count_number_of_spanning_branches = 0;
  it_sorted_indices = std::next(sorted_indices.begin(), k);
  for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){

    if(order[(*tree.nodes[*it_sorted_indices].child_left).label] < k){
      count_number_of_spanning_branches++;

      n = *tree.nodes[*it_sorted_indices].child_left;
      assert(order[n.label] < k);
      tb     = n.branch_length;
      tb_new = tb + delta_tau;

      if(n.num_events > 0.0){ 
        if(tb == 0.0){
          log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          break;
        }else if(tb_new <= 0.0){
          log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          break;
        }else{
          log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
          //log_likelihood_ratio += n.num_events * fast_log(tb_new/tb);
          log_likelihood_ratio += n.num_events * log_deltat(delta_tau/tb);
        }
      }else{
        if(tb == 0.0){
          log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          break;
        }else if(tb_new <= 0.0){
          log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          break;
        }else{
          log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
        }
      }
    }
    if(order[(*tree.nodes[*it_sorted_indices].child_right).label] < k){
      count_number_of_spanning_branches++;

      n = *tree.nodes[*it_sorted_indices].child_right;
      assert(order[n.label] < k);
      tb     = n.branch_length;
      tb_new = tb + delta_tau;

      if(n.num_events > 0.0){ 
        if(tb == 0.0){
          log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          break;
        }else if(tb_new <= 0.0){
          log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          break;
        }else{
          log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
          //log_likelihood_ratio += n.num_events * fast_log(tb_new/tb);
          log_likelihood_ratio += n.num_events * log_deltat(delta_tau/tb);
        }
      }else{
        if(tb == 0.0){
          log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          break;
        }else if(tb_new <= 0.0){
          log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          break;
        }else{
          log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
        }
      }
    }

    if(count_number_of_spanning_branches == num_lineages[sorted_indices[k]]) break;
  }

  //decide whether to accept proposal or not.
  accept = true;
  if(log_likelihood_ratio < 0.0){
    //accept with probability exp(log_p)
    if(dist_unif(rng) > exp(log_likelihood_ratio)){
      accept = false;
    }
  }

  //update coordinates and sorted_indices
  if(accept){
    //calculate new branch lengths
    it_sorted_indices = std::next(sorted_indices.begin(), k);
    update_node1 = k;
    for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      coordinates[*it_sorted_indices]            += delta_tau;
      child_left_label                            = (*tree.nodes[*it_sorted_indices].child_left).label;
      tree.nodes[child_left_label].branch_length  = coordinates[*it_sorted_indices] - coordinates[child_left_label];
      child_right_label                           = (*tree.nodes[*it_sorted_indices].child_right).label;
      tree.nodes[child_right_label].branch_length = coordinates[*it_sorted_indices] - coordinates[child_right_label];
    }
  }

}

//propose a new time and change events only up to t+bound
void
EstimateBranchLengthsWithSampleAge::ChangeTimeWhilekAncestorsVP(Tree& tree, int node, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif){

  //This step is O(N)
  int k = order[node];

  assert(node == sorted_indices[k]);
  //coorindates[sorted_indices[k-1]] determines the epoch at the lower end.
  //from there, I have to propose a new time by drawing a time for the first epoch, and if it is exceeding the epoch, for the next epoch etc.

  //for(int i = 0; i < N_total; i++){
  //  std::cerr << num_lineages[i] << " ";
  //}
  //std::cerr << std::endl;

  double age = coordinates[node];
  if(sorted_indices[k] < N){
    //std::cerr << node << " " << age << " " << sample_age[node] << std::endl;
    assert(age == sample_age[node]);
    while(sorted_indices[k] < N){ 
      k++;
      if(sorted_indices[k] < N){
        if(sample_age[sorted_indices[k]] != age) break;
      }
    }
    k--;
  }
  assert(age == coordinates[sorted_indices[k]]);
  node      = sorted_indices[k];
  tau_old   = coordinates[sorted_indices[k+1]] - age;
  k_choose_2 = num_lineages[node] * (num_lineages[node]-1.0)/2.0;
  log_likelihood_ratio = 0.0;

  double min_tip = std::numeric_limits<float>::infinity();
  for(it_order = order.begin(); it_order != std::next(order.begin(), N); it_order++){
    if(*it_order > k){
      double bl = tree.nodes[sorted_indices[*it_order]].branch_length;
      if(bl < min_tip) min_tip = bl;
    }
  }

  if(tau_old > 0.0){

    if(tau_old <= min_tip){
      tau_new   = -log(dist_unif(rng)) * tau_old;
      delta_tau = tau_new - tau_old;
      //calculate ratio of proposals
      log_likelihood_ratio = log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);   
    }else{
      tau_new   = -log(dist_unif(rng)) * (min_tip) + tau_old - min_tip;
      delta_tau = tau_new - tau_old;
      assert(min_tip+delta_tau >= 0.0);
      //calculate ratio of proposals
      log_likelihood_ratio = log(min_tip/(min_tip+delta_tau)) + ((min_tip + delta_tau)/min_tip) - (min_tip/(min_tip+delta_tau));
    }
    assert(tau_new > 0.0);

  }else{
    tau_new   = -log(dist_unif(rng)) * 1.0/k_choose_2;
    tau_old   = 0.0;
    delta_tau = tau_new;
    //calculate ratio of proposals
    log_likelihood_ratio = log(1.0/(tau_new*k_choose_2)) + tau_new*k_choose_2;
    assert(tau_new > 0.0);
  }

  ////////////////////////////////////////////////
  //coalescent prior  
  int ep_begin = 0;
  while(age >= epoch[ep_begin]){
    ep_begin++;
    if(ep_begin == (int)epoch.size()) break;
  }
  ep_begin--;
  assert(ep_begin > -1);
  assert(age >= epoch[ep_begin]);
  if( age >= epoch[ep_begin+1]  ){
    assert(ep_begin == (int) epoch.size() - 1);
  }

  double tmp_tau, delta_tmp_tau, tmp_tau_remaining, lower_coord;
  int ep            = ep_begin;
  tmp_tau           = tau_new;
  tmp_tau_remaining = tmp_tau;
  lower_coord       = age;

  int k_tmp = k;
  int num_lineages_tmp = num_lineages[sorted_indices[k_tmp]];
  float k_choose_2_tmp;
  bool is_sample = false, has_coal = false, had_sample = false;
  int end_ep = -1;
  // For k_tmp == k, I am actually changing the time while num_lineages ancestors remain.
  // For k_tmp > k, I need to check whether these were pushed to older epochs by the update, and adjust coalescent times.
  while(k_tmp < 2*N-2){

    //std::cerr << "node: " << k_tmp << " " << sorted_indices[k_tmp] << " " << num_lineages[180] << " " << num_lineages_tmp << std::endl; 
    k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
    assert(num_lineages_tmp >= 1);

    k_tmp++;
    is_sample = false;
    if(sorted_indices[k_tmp] < N){
      age = sample_age[sorted_indices[k_tmp]];
      while(sorted_indices[k_tmp] < N){ 
        k_tmp++;
        if(sorted_indices[k_tmp] < N){
          if(sample_age[sorted_indices[k_tmp]] != age) break;
        }
      }
      k_tmp--;
      if(sorted_indices[k_tmp] < N) is_sample = true;
    }
    num_lineages_tmp = num_lineages[sorted_indices[k_tmp]];

    if(ep < epoch.size() - 1){

      //std::cerr << "new: " << k_tmp << " " << is_sample << " " << has_coal << " " << ep << " " << epoch[ep] << " " << epoch[ep+1] << " " << coordinates[sorted_indices[k_tmp-1]] << " " << lower_coord << " " << lower_coord + delta_tau << std::endl;
      if(has_coal){
        tmp_tau       = coordinates[sorted_indices[k_tmp]] - lower_coord;
        if(!is_sample) tmp_tau += delta_tau;
        delta_tmp_tau = epoch[ep+1] - lower_coord;
        lower_coord   = coordinates[sorted_indices[k_tmp]];
        if(!is_sample) lower_coord += delta_tau;
        //update k_choose_2
      }else{
        assert(coordinates[sorted_indices[k_tmp-1]] >= epoch[ep]);
        if(is_sample == true){
          had_sample = true;
          tmp_tau = coordinates[sorted_indices[k_tmp]] - lower_coord; 
          if(tmp_tau_remaining < tmp_tau){
            //changed interval end here
            assert(lower_coord <= coordinates[sorted_indices[k_tmp]]);
            tmp_tau       = tmp_tau_remaining;
            delta_tmp_tau = epoch[ep+1] - lower_coord;
            lower_coord  += tmp_tau_remaining;
            assert(lower_coord <= coordinates[sorted_indices[k_tmp]]);
            is_sample     = false;
            has_coal      = true;
            k_tmp--;
          }else{
            //tmp_tau is difference to sample age
            tmp_tau_remaining -= tmp_tau;
            delta_tmp_tau = epoch[ep+1] - lower_coord;
            lower_coord   = coordinates[sorted_indices[k_tmp]];
          }
          assert(tmp_tau >= 0.0);
        }else{

          if(had_sample){
            tmp_tau       = tmp_tau_remaining;
            delta_tmp_tau = epoch[ep+1] - lower_coord;
            has_coal      = true;
            lower_coord  += tmp_tau_remaining;
            assert(lower_coord <= coordinates[sorted_indices[k_tmp]] + delta_tau);
            k_tmp--;
          }else{         
            tmp_tau = coordinates[sorted_indices[k_tmp]] + delta_tau - lower_coord;
            delta_tmp_tau = epoch[ep+1] - lower_coord;
            has_coal      = true;
            lower_coord   = coordinates[sorted_indices[k_tmp]] + delta_tau;
          }
          //lower_coord   = coordinates[sorted_indices[k_tmp]];
        }
      }
      //std::cerr << "debug: " << ep << " " << delta_tmp_tau << " " << tmp_tau << " " << tmp_tau_remaining << std::endl;
      assert(delta_tmp_tau >= 0.0);
      //if epoch[ep+1] begins before the current interval (while num_lineages_tmp remain) ends, enter if statement
      if(delta_tmp_tau <= tmp_tau){

        //add up rate parameters for this interval
        if(coal_rate[ep] > 0.0){
          log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
        }
        tmp_tau                -= delta_tmp_tau;

        ep++;
        delta_tmp_tau           = epoch[ep+1] - epoch[ep];
        while(tmp_tau > delta_tmp_tau && ep < epoch.size() - 1){
          if(coal_rate[ep] > 0.0){
            log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
          }
          tmp_tau              -= delta_tmp_tau;

          ep++;
          delta_tmp_tau         = epoch[ep+1] - epoch[ep];
        }
        assert(tmp_tau >= 0.0);
        if(coal_rate[ep] == 0){
          log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau;
          if(!is_sample) log_likelihood_ratio += log(coal_rate[ep]);
        }

      }else{
        if(coal_rate[ep] == 0){
          log_likelihood_ratio = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau;
          if(!is_sample) log_likelihood_ratio += log(coal_rate[ep]);
        }
      }

    }else{

      if(coal_rate[ep] == 0){
        log_likelihood_ratio = -std::numeric_limits<float>::infinity();
      }else{

        if(has_coal){
          tmp_tau       = coordinates[sorted_indices[k_tmp]] - lower_coord;
          lower_coord   = coordinates[sorted_indices[k_tmp]];
          if(!is_sample){
            tmp_tau     += delta_tau;
            lower_coord += delta_tau;
          }
          //update k_choose_2
        }else{
          assert(coordinates[sorted_indices[k_tmp-1]] >= epoch[ep]);
          if(is_sample == true){
            had_sample = true;
            tmp_tau = coordinates[sorted_indices[k_tmp]] - lower_coord;
            if(tmp_tau_remaining < tmp_tau){
              //changed interval end here
              lower_coord  += tmp_tau_remaining;
              assert(lower_coord <= coordinates[sorted_indices[k_tmp]]);
              tmp_tau       = tmp_tau_remaining;
              is_sample     = false;
              has_coal      = true;
              k_tmp--;
            }else{
              tmp_tau_remaining -= tmp_tau;
              lower_coord   = coordinates[sorted_indices[k_tmp]];
            }
            assert(tmp_tau >= 0.0);
          }else{
            if(had_sample){
              tmp_tau       = tmp_tau_remaining;
              has_coal      = true;
              lower_coord  += tmp_tau_remaining;
              assert(lower_coord <= coordinates[sorted_indices[k_tmp]] + delta_tau);
              k_tmp--;
            }else{         
              tmp_tau = coordinates[sorted_indices[k_tmp]] + delta_tau - lower_coord;
              has_coal      = true;
              lower_coord   = coordinates[sorted_indices[k_tmp]] + delta_tau;
            }
          }
        }

        log_likelihood_ratio -= k_choose_2_tmp * coal_rate[ep] * tmp_tau;
        if(!is_sample) log_likelihood_ratio += log(coal_rate[ep]);

      }

    }

  }

  if(log_likelihood_ratio != -std::numeric_limits<float>::infinity()){

    ep            = ep_begin;
    tmp_tau       = tau_old;
    lower_coord = coordinates[node];

    int k_max = k_tmp;
    k_tmp = k;
    num_lineages_tmp = num_lineages[sorted_indices[k_tmp]];
    assert(coordinates[sorted_indices[k_tmp+1]] - lower_coord == tau_old);
    is_sample = false;
    has_coal = false;
    end_ep = -1;
    // For k_tmp == k, I am actually changing the time while num_lineages ancestors remain.
    // For k_tmp > k, I need to check whether these were pushed to older epochs by the update, and adjust coalescent times.
    while(k_tmp < k_max){

      k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
      assert(num_lineages_tmp >= 1);

      k_tmp++;
      is_sample = false;
      if(sorted_indices[k_tmp] < N){
        age = sample_age[sorted_indices[k_tmp]];
        if(sorted_indices[k_tmp-1] < N) assert(age > coordinates[sorted_indices[k_tmp-1]]);
        while(sorted_indices[k_tmp] < N){ 
          k_tmp++;
          if(sorted_indices[k_tmp] < N){
            if(sample_age[sorted_indices[k_tmp]] != age) break;
          }
        }
        k_tmp--;
        if(sorted_indices[k_tmp] < N) is_sample = true;
      }
      num_lineages_tmp = num_lineages[sorted_indices[k_tmp]];

      //std::cerr << "old: " << is_sample << " " << has_coal << " " << ep << " " << epoch[ep] << " " << epoch[ep+1] << " " << coordinates[sorted_indices[k_tmp-1]] << " " << lower_coord << std::endl;
      if(ep < epoch.size() - 1){

        if(has_coal){
          tmp_tau       = coordinates[sorted_indices[k_tmp]] - lower_coord;
          delta_tmp_tau = epoch[ep+1] - lower_coord;
          if(delta_tmp_tau < 0.0){
            std::cerr << k << " " << k_tmp << " " << epoch[ep+1] << " " << lower_coord << " " << coordinates[sorted_indices[k_tmp-2]] << " " << coordinates[sorted_indices[k_tmp-1]] << " " << coordinates[sorted_indices[k]] << " " << node << " " << sorted_indices[k_tmp-2] << " " << sorted_indices[k_tmp-1] << " " << tau_old << std::endl;

            for(int i = 0; i < N_total; i++){
              std::cerr << i << " " << coordinates[sorted_indices[i]] << std::endl;
            }
            std::cerr << std::endl << std::endl;


          }
          lower_coord   = coordinates[sorted_indices[k_tmp]];
          //update k_choose_2
        }else{
          assert(coordinates[sorted_indices[k_tmp-1]] >= epoch[ep]);
          if(is_sample == true){
            tmp_tau = coordinates[sorted_indices[k_tmp]] - lower_coord;
            //tmp_tau is difference to sample age
            delta_tmp_tau = epoch[ep+1] - lower_coord;
            lower_coord   = coordinates[sorted_indices[k_tmp]];
            assert(tmp_tau >= 0.0);
          }else{
            tmp_tau       = coordinates[sorted_indices[k_tmp]] - lower_coord;
            delta_tmp_tau = epoch[ep+1] - lower_coord;
            lower_coord   = coordinates[sorted_indices[k_tmp]];
            has_coal      = true;
          }
        }

        assert(delta_tmp_tau >= 0.0);
        //if epoch[ep+1] begins before the current interval (while num_lineages_tmp remain) ends, enter if statement
        if(delta_tmp_tau <= tmp_tau){
          //add up rate parameters for this interval
          if(coal_rate[ep] > 0.0){
            log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
          }
          tmp_tau                -= delta_tmp_tau;

          ep++;
          delta_tmp_tau           = epoch[ep+1] - epoch[ep];
          while(tmp_tau > delta_tmp_tau && ep < epoch.size() - 1){
            if(coal_rate[ep] > 0.0){
              log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
            }
            tmp_tau              -= delta_tmp_tau;

            ep++;
            delta_tmp_tau         = epoch[ep+1] - epoch[ep];
          }
          assert(tmp_tau >= 0.0);
          if(coal_rate[ep] == 0){
            log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          }else{
            log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau;
            if(!is_sample) log_likelihood_ratio -= log(coal_rate[ep]);
          }

        }else{
          if(coal_rate[ep] == 0){
            log_likelihood_ratio = -std::numeric_limits<float>::infinity();
          }else{
            log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau;
            if(!is_sample) log_likelihood_ratio -= log(coal_rate[ep]);
          }
        }

      }else{

        if(coal_rate[ep] == 0){
          log_likelihood_ratio = -std::numeric_limits<float>::infinity();
        }else{

          if(has_coal){
            tmp_tau       = coordinates[sorted_indices[k_tmp]] - lower_coord;
            lower_coord   = coordinates[sorted_indices[k_tmp]];
            //update k_choose_2
          }else{
            assert(coordinates[sorted_indices[k_tmp-1]] >= epoch[ep]);
            if(is_sample == true){
              tmp_tau = coordinates[sorted_indices[k_tmp]] - lower_coord;
              lower_coord   = coordinates[sorted_indices[k_tmp]];
              assert(tmp_tau >= 0.0);
            }else{
              tmp_tau       = coordinates[sorted_indices[k_tmp]] - lower_coord;
              delta_tmp_tau = epoch[ep+1] - lower_coord;
              lower_coord   = coordinates[sorted_indices[k_tmp]];
              has_coal      = true;
            }
          }

          log_likelihood_ratio += k_choose_2_tmp * coal_rate[ep] * tmp_tau;
          if(!is_sample) log_likelihood_ratio -= log(coal_rate[ep]);

        }

      }

    }

    ////////////////////////////////////////////
    //likelihood P(Data | tb)

    if(log_likelihood_ratio != std::numeric_limits<float>::infinity()){
      //assert(order[node_k] == k);
      int count_number_of_spanning_branches = 0;
      age = coordinates[node];
      it_sorted_indices = std::next(sorted_indices.begin(), k+1);
      for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){

        if(*it_sorted_indices < N){
          n = tree.nodes[*it_sorted_indices];
          tb     = n.branch_length;
          tb_new = tb + delta_tau;

          if(tb == 0.0){
            log_likelihood_ratio  = std::numeric_limits<float>::infinity();
            break;
          }else if(tb_new <= 0.0){
            log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
            break;
          }else{
            log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
            log_likelihood_ratio += n.num_events * log(tb_new/tb);
          }
        }else{
          if(coordinates[(*tree.nodes[*it_sorted_indices].child_left).label] < age){
            count_number_of_spanning_branches++;

            n = *tree.nodes[*it_sorted_indices].child_left;
            assert(order[n.label] < k+1);
            tb     = n.branch_length;
            tb_new = tb + delta_tau;

            if(tb == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
              break;
            }else if(tb_new <= 0.0){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
              break;
            }else{
              log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
              log_likelihood_ratio += n.num_events * log(tb_new/tb);
            }
          }

          if(coordinates[(*tree.nodes[*it_sorted_indices].child_right).label] < age){
            count_number_of_spanning_branches++;

            n = *tree.nodes[*it_sorted_indices].child_right;
            assert(order[n.label] < k+1);
            tb     = n.branch_length;
            tb_new = tb + delta_tau;
            //assert(tb_new >= 0.0);

            if(tb == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
              break;
            }else if(tb_new <= 0.0){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
              break;
            }else{
              log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
              log_likelihood_ratio += n.num_events * log(tb_new/tb);
            }
          }

        }
        if(count_number_of_spanning_branches == num_lineages[node]) break;
      }
      //assert(count_number_of_spanning_branches == num_lineages);
    }
  }

  //decide whether to accept proposal or not.
  accept = true;
  if(log_likelihood_ratio < 0.0){
    //accept with probability exp(log_p)
    if(dist_unif(rng) > exp(log_likelihood_ratio)){
      accept = false;
    }
  }

  //update coordinates and sorted_indices
  if(accept){
    //calculate new branch lengths
    int k_init = k;
    k++;  
    update_node1 = k;
    //sorted_indices, order, num_lineages will change here
    for(it_sorted_indices = std::next(sorted_indices.begin(), k); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        coordinates[*it_sorted_indices] += delta_tau;
        child_left_label                                 = (*tree.nodes[*it_sorted_indices].child_left).label;
        tree.nodes[child_left_label].branch_length       = coordinates[*it_sorted_indices] - coordinates[child_left_label];
        child_right_label                                = (*tree.nodes[*it_sorted_indices].child_right).label;
        tree.nodes[child_right_label].branch_length      = coordinates[*it_sorted_indices] - coordinates[child_right_label];
        assert(tree.nodes[child_left_label].branch_length >= 0.0);
        assert(tree.nodes[child_right_label].branch_length >= 0.0);
      }
    } 

    std::size_t m1(0);
    std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
    //std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) { return coordinates[i1] < coordinates[i2]; } );
    std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
        return std::tie(coordinates[i1],i1) < std::tie(coordinates[i2],i2); } );

    //obtain order of coalescent events
    std::fill(order.begin(), order.end(), 0);
    std::size_t m2(0);
    std::generate(std::begin(order), std::end(order), [&]{ return m2++; });
    std::sort(std::begin(order), std::end(order), [&](int i1, int i2) { return sorted_indices[i1] < sorted_indices[i2]; } );

    int num_lins = 0;
    double ages = sample_age[*sorted_indices.begin()];
    std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
    for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages[*it_sorted_indices_start] = num_lins; 
        }
        num_lins--;
        num_lineages[*it_sorted_indices] = num_lins;
        it_sorted_indices_start++;
      }else if(ages < sample_age[*it_sorted_indices]){      
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages[*it_sorted_indices_start] = num_lins; 
        }
        ages = sample_age[*it_sorted_indices];
        num_lins++;
      }else{
        num_lins++;
      }
    }


  }

}


//Calculate coalescent prior
double
EstimateBranchLengthsWithSampleAge::CalculatePrior(std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages){

  double log_likelihood = 0.0;
  int k_tmp  = 0;
  int node   = p_sorted_indices[k_tmp];
  double age = p_coordinates[node];
  if(node < N){
    assert(age == sample_age[node]);
    while(p_sorted_indices[k_tmp] < N){ 
      k_tmp++;
      if(p_sorted_indices[k_tmp] < N){
        if(sample_age[p_sorted_indices[k_tmp]] != age) break;
      }
    }
    k_tmp--;
  }

  double lower_coord = p_coordinates[p_sorted_indices[k_tmp]], tmp_tau, delta_tmp_tau;
  int num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]], k_choose_2_tmp;
  bool is_sample = false;

  while(k_tmp < 2*N-2){

    k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
    assert(num_lineages_tmp >= 1);

    k_tmp++;
    is_sample = false;
    if(p_sorted_indices[k_tmp] < N){
      age = sample_age[p_sorted_indices[k_tmp]];
      if(p_sorted_indices[k_tmp-1] < N) assert(age > p_coordinates[p_sorted_indices[k_tmp-1]]);
      while(p_sorted_indices[k_tmp] < N){ 
        k_tmp++;
        if(p_sorted_indices[k_tmp] < N){
          if(sample_age[p_sorted_indices[k_tmp]] != age) break;
        }
      }
      k_tmp--;
      if(p_sorted_indices[k_tmp] < N) is_sample = true;
    }
    num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]];

    if(is_sample == true){
      tmp_tau = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
      lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
      if(!(tmp_tau >= 0.0)) std::cerr << tmp_tau << std::endl;
      assert(tmp_tau >= 0.0);
    }else{
      tmp_tau       = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
      lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
    }

    log_likelihood -= k_choose_2_tmp * tmp_tau;

  }

  return(log_likelihood);

}

double
EstimateBranchLengthsWithSampleAge::CalculatePrior(int k_start, int k_end, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages){

  double log_likelihood = 0.0;
  int k_tmp  = k_start;
  int node   = p_sorted_indices[k_tmp];
  double age = p_coordinates[node];
  if(node < N){
    assert(age == sample_age[node]);
    while(p_sorted_indices[k_tmp] < N){ 
      k_tmp++;
      if(p_sorted_indices[k_tmp] < N){
        if(sample_age[p_sorted_indices[k_tmp]] != age) break;
      }
    }
    k_tmp--;
  }

  double lower_coord = p_coordinates[p_sorted_indices[k_tmp]], tmp_tau, delta_tmp_tau;
  int num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]], k_choose_2_tmp;
  bool is_sample = false;

  while(k_tmp < k_end){

    k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
    assert(num_lineages_tmp >= 1);

    k_tmp++;
    is_sample = false;
    if(p_sorted_indices[k_tmp] < N){
      age = sample_age[p_sorted_indices[k_tmp]];
      if(p_sorted_indices[k_tmp-1] < N) assert(age > p_coordinates[p_sorted_indices[k_tmp-1]]);
      while(p_sorted_indices[k_tmp] < N){ 
        k_tmp++;
        if(k_tmp == k_end) break;
        if(p_sorted_indices[k_tmp] < N){
          if(sample_age[p_sorted_indices[k_tmp]] != age) break;
        }
      }
      k_tmp--;
      if(p_sorted_indices[k_tmp] < N) is_sample = true;
    }
    num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]];

    if(is_sample == true){
      tmp_tau = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
      lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
      assert(tmp_tau >= 0.0);
    }else{
      tmp_tau       = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
      lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
    }

    log_likelihood -= k_choose_2_tmp * tmp_tau;

  }

  return(log_likelihood);

}

double
EstimateBranchLengthsWithSampleAge::CalculatePrior(const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages){

  double log_likelihood = 0.0;
  int ep     = 0; 
  int k_tmp  = 0;
  int node   = p_sorted_indices[k_tmp];
  double age = p_coordinates[node];
  if(node < N){
    assert(age == sample_age[node]);
    while(p_sorted_indices[k_tmp] < N){ 
      k_tmp++;
      if(p_sorted_indices[k_tmp] < N){
        if(sample_age[p_sorted_indices[k_tmp]] != age) break;
      }
    }
    k_tmp--;
  }

  double lower_coord = p_coordinates[p_sorted_indices[k_tmp]], tmp_tau, delta_tmp_tau;
  int num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]], k_choose_2_tmp;
  bool is_sample = false;

  while(k_tmp < 2*N-2){

    k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
    assert(num_lineages_tmp >= 1);

    k_tmp++;
    is_sample = false;
    if(p_sorted_indices[k_tmp] < N){
      age = sample_age[p_sorted_indices[k_tmp]];
      if(p_sorted_indices[k_tmp-1] < N) assert(age > p_coordinates[p_sorted_indices[k_tmp-1]]);
      while(p_sorted_indices[k_tmp] < N){ 
        k_tmp++;
        if(p_sorted_indices[k_tmp] < N){
          if(sample_age[p_sorted_indices[k_tmp]] != age) break;
        }
      }
      k_tmp--;
      is_sample = true;
    }
    num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]];

    if(ep < epoch.size() - 1){

      assert(p_coordinates[p_sorted_indices[k_tmp-1]] >= epoch[ep]);
      if(is_sample == true){
        tmp_tau = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
        //tmp_tau is difference to sample age
        delta_tmp_tau = epoch[ep+1] - lower_coord;
        lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
        assert(tmp_tau >= 0.0);
      }else{
        tmp_tau       = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
        delta_tmp_tau = epoch[ep+1] - lower_coord;
        lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
      }
      assert(delta_tmp_tau >= 0.0);
      //if epoch[ep+1] begins before the current interval (while num_lineages_tmp remain) ends, enter if statement
      if(delta_tmp_tau <= tmp_tau){
        //add up rate parameters for this interval
        if(coal_rate[ep] > 0.0){
          log_likelihood -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
        }
        tmp_tau                -= delta_tmp_tau;

        ep++;
        delta_tmp_tau           = epoch[ep+1] - epoch[ep];
        while(tmp_tau > delta_tmp_tau && ep < epoch.size() - 1){
          if(coal_rate[ep] > 0.0){
            log_likelihood -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
          }
          tmp_tau              -= delta_tmp_tau;

          ep++;
          delta_tmp_tau         = epoch[ep+1] - epoch[ep];
        }
        assert(tmp_tau >= 0.0);
        if(coal_rate[ep] == 0){
          log_likelihood  = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood -= k_choose_2_tmp*coal_rate[ep] * tmp_tau;
          if(!is_sample) log_likelihood += log(coal_rate[ep]);
        }

      }else{
        if(coal_rate[ep] == 0){
          log_likelihood = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood -= k_choose_2_tmp*coal_rate[ep] * tmp_tau;
          if(!is_sample) log_likelihood += log(coal_rate[ep]);
        }
      }

    }else{

      if(coal_rate[ep] == 0){
        log_likelihood = -std::numeric_limits<float>::infinity();
      }else{

        assert(p_coordinates[p_sorted_indices[k_tmp-1]] >= epoch[ep]);
        if(is_sample == true){
          tmp_tau = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
          lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
          assert(tmp_tau >= 0.0);
        }else{
          tmp_tau       = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
          delta_tmp_tau = epoch[ep+1] - lower_coord;
          lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
        }

        log_likelihood -= k_choose_2_tmp * coal_rate[ep] * tmp_tau;
        if(!is_sample) log_likelihood += log(coal_rate[ep]);

      }
    }
  }

  return(log_likelihood);

}

double
EstimateBranchLengthsWithSampleAge::CalculatePrior(int k_start, int k_end, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages){

  double log_likelihood = 0.0;
  int k_tmp  = k_start;
  int node   = p_sorted_indices[k_tmp];
  double age = p_coordinates[node];
  if(node < N){
    assert(age == sample_age[node]);
    while(p_sorted_indices[k_tmp] < N){ 
      k_tmp++;
      if(p_sorted_indices[k_tmp] < N){
        if(sample_age[p_sorted_indices[k_tmp]] != age) break;
      }
    }
    k_tmp--;
  }

  //coalescent prior  
  int ep = 0;
  while(p_coordinates[p_sorted_indices[k_tmp]] >= epoch[ep]){
    ep++;
    if(ep == (int)epoch.size()) break;
  }
  ep--;
  assert(ep > -1);
  assert(p_coordinates[p_sorted_indices[k_tmp]] >= epoch[ep]);
  //std::cerr << ep << " " << epoch[ep] << " " << k_tmp << " " << p_sorted_indices[k_tmp] << " " << p_coordinates[p_sorted_indices[k_tmp]] << " " << p_coordinates[p_sorted_indices[k_tmp + 1]] << std::endl;

  double lower_coord = p_coordinates[p_sorted_indices[k_tmp]], tmp_tau, delta_tmp_tau;
  int num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]], k_choose_2_tmp;
  bool is_sample = false;

  while(k_tmp < k_end){

    k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
    if(num_lineages_tmp < 1) std::cerr << num_lineages_tmp << std::endl;
    assert(num_lineages_tmp >= 1);

    k_tmp++;
    is_sample = false;
    if(p_sorted_indices[k_tmp] < N){
      age = sample_age[p_sorted_indices[k_tmp]];
      if(p_sorted_indices[k_tmp-1] < N) assert(age > p_coordinates[p_sorted_indices[k_tmp-1]]);
      while(p_sorted_indices[k_tmp] < N){ 
        k_tmp++;
        if(p_sorted_indices[k_tmp] < N){
          if(sample_age[p_sorted_indices[k_tmp]] != age) break;
        }
      }
      k_tmp--;
      if(p_sorted_indices[k_tmp] < N) is_sample = true;
    }
    num_lineages_tmp = p_num_lineages[p_sorted_indices[k_tmp]];

    if(ep < epoch.size() - 1){

      assert(p_coordinates[p_sorted_indices[k_tmp-1]] >= epoch[ep]);
      if(is_sample == true){
        tmp_tau = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
        //tmp_tau is difference to sample age
        delta_tmp_tau = epoch[ep+1] - lower_coord;
        lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
        assert(tmp_tau >= 0.0);
      }else{
        tmp_tau       = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
        delta_tmp_tau = epoch[ep+1] - lower_coord;
        lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
      }
      assert(delta_tmp_tau >= 0.0);
      //if epoch[ep+1] begins before the current interval (while num_lineages_tmp remain) ends, enter if statement
      if(delta_tmp_tau <= tmp_tau){
        //add up rate parameters for this interval
        if(coal_rate[ep] > 0.0){
          log_likelihood   -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
        }
        tmp_tau                -= delta_tmp_tau;

        ep++;
        delta_tmp_tau           = epoch[ep+1] - epoch[ep];
        while(tmp_tau > delta_tmp_tau && ep < epoch.size() - 1){
          if(coal_rate[ep] > 0.0){
            log_likelihood -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
          }
          tmp_tau              -= delta_tmp_tau;

          ep++;
          delta_tmp_tau         = epoch[ep+1] - epoch[ep];
        }
        assert(tmp_tau >= 0.0);
        if(coal_rate[ep] == 0){
          log_likelihood  = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood  -= k_choose_2_tmp*coal_rate[ep] * tmp_tau;
          if(!is_sample) log_likelihood += log(coal_rate[ep]);
        }

      }else{
        if(coal_rate[ep] == 0){
          log_likelihood = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood  -= k_choose_2_tmp*coal_rate[ep] * tmp_tau;
          if(!is_sample) log_likelihood += log(coal_rate[ep]);
        }
      }

    }else{

      if(coal_rate[ep] == 0){
        log_likelihood = -std::numeric_limits<float>::infinity();
      }else{

        assert(p_coordinates[p_sorted_indices[k_tmp-1]] >= epoch[ep]);
        if(is_sample == true){
          tmp_tau = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
          lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
          assert(tmp_tau >= 0.0);
        }else{
          tmp_tau       = p_coordinates[p_sorted_indices[k_tmp]] - lower_coord;
          delta_tmp_tau = epoch[ep+1] - lower_coord;
          lower_coord   = p_coordinates[p_sorted_indices[k_tmp]];
        }

        log_likelihood -= k_choose_2_tmp * coal_rate[ep] * tmp_tau;
        if(!is_sample) log_likelihood += log(coal_rate[ep]);

      }
    }
  }

  return(log_likelihood);

}


////////////////////// constant NE /////////////////////
//propose a new time and change events only up to t+bound
void
EstimateBranchLengthsWithSampleAge::ChangeTimeWhilekAncestors_new(Tree& tree, int node, std::uniform_real_distribution<double>& dist_unif){

  //This step is O(N)
  int k = order[node];

  assert(node == sorted_indices[k]);
  //coorindates[sorted_indices[k-1]] determines the epoch at the lower end.
  //from there, I have to propose a new time by drawing a time for the first epoch, and if it is exceeding the epoch, for the next epoch etc.

  double age = coordinates[node];
  if(sorted_indices[k] < N){
    //std::cerr << node << " " << age << " " << sample_age[node] << std::endl;
    assert(age == sample_age[node]);
    while(sorted_indices[k] < N){ 
      k++;
      if(sorted_indices[k] < N){
        if(sample_age[sorted_indices[k]] != age) break;
      }
    }
    k--;
  }
  assert(age == coordinates[sorted_indices[k]]);
  node      = sorted_indices[k];
  tau_old   = coordinates[sorted_indices[k+1]] - age;
  k_choose_2 = num_lineages[node] * (num_lineages[node]-1.0)/2.0;
  log_likelihood_ratio = 0.0;

  double min_tip = std::numeric_limits<float>::infinity();
  for(it_order = order.begin(); it_order != std::next(order.begin(), N); it_order++){
    if(*it_order > k){
      double bl = tree.nodes[sorted_indices[*it_order]].branch_length;
      if(bl < min_tip) min_tip = bl;
    }
  }

  if(min_tip > 0 && tau_old > 0){

    if(tau_old > 0.0){

      if(tau_old <= min_tip){
        tau_new   = -log(dist_unif(rng)) * tau_old;
        delta_tau = tau_new - tau_old;
        //calculate ratio of proposals
        log_likelihood_ratio = log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);   
      }else{
        tau_new   = -log(dist_unif(rng)) * (min_tip) + tau_old - min_tip;
        delta_tau = tau_new - tau_old;
        assert(min_tip+delta_tau >= 0.0);
        //calculate ratio of proposals
        log_likelihood_ratio = log(min_tip/(min_tip+delta_tau)) + ((min_tip + delta_tau)/min_tip) - (min_tip/(min_tip+delta_tau));
      }
      assert(tau_new > 0.0);

    }else{
      assert(0);
    }

    ////////////////////////////////////////////////
    //coalescent prior  

    //calculate new coordinates, order, sorted_indices, num_lineages

    double log_likelihood;

    std::vector<int> sorted_indices_new = sorted_indices, order_new = order, num_lineages_new = num_lineages;
    std::vector<double> coordinates_new = coordinates;

    for(it_sorted_indices = std::next(sorted_indices_new.begin(), k+1); it_sorted_indices != sorted_indices_new.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        coordinates_new[*it_sorted_indices] += delta_tau;
      }
    }

    std::size_t m1(0);
    std::generate(std::begin(sorted_indices_new), std::end(sorted_indices_new), [&]{ return m1++; });
    std::sort(std::begin(sorted_indices_new), std::end(sorted_indices_new), [&](int i1, int i2) {
        return std::tie(coordinates_new[i1],i1) < std::tie(coordinates_new[i2],i2); } );

    //obtain order of coalescent events
    std::fill(order_new.begin(), order_new.end(), 0);
    std::size_t m2(0);
    std::generate(std::begin(order_new), std::end(order_new), [&]{ return m2++; });
    std::sort(std::begin(order_new), std::end(order_new), [&](int i1, int i2) { return sorted_indices_new[i1] < sorted_indices_new[i2]; } );

    int num_lins = 0;
    double ages = sample_age[*sorted_indices_new.begin()];
    std::vector<int>::iterator it_sorted_indices_start = sorted_indices_new.begin();
    for(it_sorted_indices = sorted_indices_new.begin(); it_sorted_indices != sorted_indices_new.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages_new[*it_sorted_indices_start] = num_lins; 
        }
        num_lins--;
        num_lineages_new[*it_sorted_indices] = num_lins;
        it_sorted_indices_start++;
      }else if(ages < sample_age[*it_sorted_indices]){      
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages_new[*it_sorted_indices_start] = num_lins; 
        }
        ages = sample_age[*it_sorted_indices];
        num_lins++;
      }else{
        num_lins++;
      }
    }



    log_likelihood = CalculatePrior(coordinates_new, sorted_indices_new, num_lineages_new);
    if(log_likelihood != -std::numeric_limits<float>::infinity()){
      log_likelihood_ratio += log_likelihood;
      log_likelihood = CalculatePrior(coordinates, sorted_indices, num_lineages);
      if(log_likelihood != -std::numeric_limits<float>::infinity()) log_likelihood_ratio -= log_likelihood;
    }


    if(log_likelihood_ratio != -std::numeric_limits<float>::infinity()){

      ////////////////////////////////////////////
      //likelihood P(Data | tb)

      if(log_likelihood_ratio != std::numeric_limits<float>::infinity()){
        //assert(order[node_k] == k);
        int count_number_of_spanning_branches = 0;
        age = coordinates[node];
        it_sorted_indices = std::next(sorted_indices.begin(), k+1);
        for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){

          if(*it_sorted_indices < N){
            n = tree.nodes[*it_sorted_indices];
            tb     = n.branch_length;
            tb_new = tb + delta_tau;

            if(tb == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
              break;
            }else if(tb_new <= 0.0){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
              break;
            }else{
              log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
              log_likelihood_ratio += n.num_events * log(tb_new/tb);
            }
          }else{
            if(coordinates[(*tree.nodes[*it_sorted_indices].child_left).label] < age){
              count_number_of_spanning_branches++;

              n = *tree.nodes[*it_sorted_indices].child_left;
              assert(order[n.label] < k+1);
              tb     = n.branch_length;
              tb_new = tb + delta_tau;

              if(tb == 0.0){
                log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                break;
              }else if(tb_new <= 0.0){
                log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
                break;
              }else{
                log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
                log_likelihood_ratio += n.num_events * log(tb_new/tb);
              }
            }

            if(coordinates[(*tree.nodes[*it_sorted_indices].child_right).label] < age){
              count_number_of_spanning_branches++;

              n = *tree.nodes[*it_sorted_indices].child_right;
              assert(order[n.label] < k+1);
              tb     = n.branch_length;
              tb_new = tb + delta_tau;
              //assert(tb_new >= 0.0);

              if(tb == 0.0){
                log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                break;
              }else if(tb_new <= 0.0){
                log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
                break;
              }else{
                log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
                log_likelihood_ratio += n.num_events * log(tb_new/tb);
              }
            }

          }
          if(count_number_of_spanning_branches == num_lineages[node]) break;
        }
        //assert(count_number_of_spanning_branches == num_lineages);
      }
    }

    //decide whether to accept proposal or not.
    accept = true;
    if(log_likelihood_ratio < 0.0){
      //accept with probability exp(log_p)
      if(dist_unif(rng) > exp(log_likelihood_ratio)){
        accept = false;
      }
    }

    //update coordinates and sorted_indices
    if(accept){
      //calculate new branch lengths
      int k_init = k;
      k++;  
      update_node1 = k;
      coordinates = coordinates_new;

      //sorted_indices, order, num_lineages will change here
      for(it_sorted_indices = std::next(sorted_indices.begin(), k); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
        if(*it_sorted_indices >= N){
          child_left_label                                 = (*tree.nodes[*it_sorted_indices].child_left).label;
          tree.nodes[child_left_label].branch_length       = coordinates[*it_sorted_indices] - coordinates[child_left_label];
          child_right_label                                = (*tree.nodes[*it_sorted_indices].child_right).label;
          tree.nodes[child_right_label].branch_length      = coordinates[*it_sorted_indices] - coordinates[child_right_label];
          assert(tree.nodes[child_left_label].branch_length >= 0.0);
          assert(tree.nodes[child_right_label].branch_length >= 0.0);
        }
      }

      sorted_indices = sorted_indices_new;
      order = order_new;
      num_lineages = num_lineages_new;
    }

  }else{
    update_node1 = k;
  }

}

//This changes the time of one event, with a beta proposal within the time while k ancestors remain
void
EstimateBranchLengthsWithSampleAge::UpdateOneEvent(Tree& tree, int node_k, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif){

  accept = true;
  log_likelihood_ratio = 0.0;

  //propose to change time of node_k between time of older daughter and parent

  if(tree.nodes[node_k].parent == NULL){
    //propose with exponential distribution  

    tau_old   = coordinates[node_k] - coordinates[(*tree.nodes[node_k].child_left).label];
    if(tau_old > coordinates[node_k] - coordinates[(*tree.nodes[node_k].child_right).label]){
      tau_old = coordinates[node_k] - coordinates[(*tree.nodes[node_k].child_right).label];
    }

    if(tau_old > 0.0){
      tau_new   = -fast_log(dist_unif(rng)) * tau_old;
      delta_tau = tau_new - tau_old;
      assert(tau_new >= 0.0);

      //now decide whether to accept delta_tau:
      //calculate likelihood ratio

      //proposal likelihood ratio
      log_likelihood_ratio = fast_log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);
    }else{
      k_choose_2 = num_lineages[node_k] * (num_lineages[node_k]+1.0)/2.0;
      tau_new   = -log(dist_unif(rng)) * 1.0/k_choose_2;
      tau_old   = 0.0;
      delta_tau = tau_new;    
      //calculate ratio of proposals
      log_likelihood_ratio = fast_log(1.0/(tau_new*k_choose_2)) + tau_new*k_choose_2;
    }
    //coalescent prior

    //this case is unaffected by ancient DNA
    log_likelihood_ratio -= delta_tau;

    //likelihood function, likelihood ratio  
    child_left_label  = (*tree.nodes[node_k].child_left).label;
    child_right_label = (*tree.nodes[node_k].child_right).label;

    child_left_num_events  = tree.nodes[child_left_label].num_events;
    child_right_num_events = tree.nodes[child_right_label].num_events;

    tb_child_left      = tree.nodes[child_left_label].branch_length;
    tb_child_right     = tree.nodes[child_right_label].branch_length;

    //mutation and recombination part

    if(tb_child_left == 0.0){
      log_likelihood_ratio  = std::numeric_limits<float>::infinity();
    }else if(tb_child_left <= -delta_tau){
      log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
    }else{

      if(tb_child_right == 0.0){
        log_likelihood_ratio  = std::numeric_limits<float>::infinity();
      }else if(tb_child_right <= -delta_tau){
        log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
      }else{

        log_likelihood_ratio += (- mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
        if(child_right_num_events >= 1.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
        if(child_left_num_events >= 1.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

      }
    }

    //decide whether to accept proposal or not.
    accept = true;
    if(log_likelihood_ratio < 0.0){
      //accept with probability exp(log_p)
      if(dist_unif(rng) > exp(log_likelihood_ratio)){
        accept = false;
      }
    }

    //update coordinates and sorted_indices
    if(accept){
      //calculate new branch lengths
      update_node3         = node_k;
      update_node1         = node_k;
      coordinates[node_k] += delta_tau;
      tree.nodes[child_left_label].branch_length  = coordinates[node_k] - coordinates[child_left_label];
      tree.nodes[child_right_label].branch_length = coordinates[node_k] - coordinates[child_right_label];
      //no change in num_lineages, order, sorted_indices
    }

  }else{ 

    //tau_below    = coordinates[node_k] - coordinates[sorted_indices[k-1]];
    //tau_above    = coordinates[sorted_indices[k+1]] - coordinates[node_k];

    child_left_label  = (*tree.nodes[node_k].child_left).label;
    child_right_label = (*tree.nodes[node_k].child_right).label;
    parent_label      = (*tree.nodes[node_k].parent).label;

    tb_child_left     = tree.nodes[child_left_label].branch_length;
    tb_child_right    = tree.nodes[child_right_label].branch_length;
    tb                = tree.nodes[node_k].branch_length;

    tau_below    = tb_child_left;
    if(tb_child_right < tau_below) tau_below = tb_child_right;
    tau_above    = tb;
    T            = tau_below + tau_above;

    int k_start, k_end, k = order[node_k];
    if(tau_below >= 0 && tau_above >= 0){

      //sample a new tau_below
      tau_new_below        = dist_unif(rng);
      tau_new_below       *= T; 

      delta_tau            = tau_new_below - tau_below;
      tau_new_above        = T - tau_new_below;

      if(tau_new_above >= 0.0 && tau_new_below >= 0.0){

        //now decide whether to accept delta_tau:

        //calculate likelihood ratio
        //proposal likelihood ratio
        log_likelihood_ratio = 0.0;

        //calculate new coordinates, order, sorted_indices, num_lineages
        double log_likelihood;
        k_end = order[parent_label];

        double coords = coordinates[node_k];
        double coords_new = coords + delta_tau;
        if(coords_new > coordinates[parent_label]) coords_new = coordinates[parent_label];
        if(coords_new < coordinates[child_left_label]) coords_new = coordinates[child_left_label];
        if(coords_new < coordinates[child_right_label]) coords_new = coordinates[child_right_label];

        //anything outside k_start and k_end is identical, so only need to update this bit
        if(delta_tau > 0){
          k_start = k;

          sorted_indices_new[k_start-1]               = sorted_indices[k_start-1];
          num_lineages_new[sorted_indices[k_start-1]] = num_lineages[sorted_indices[k_start-1]];
          sorted_indices_new[k_start]               = sorted_indices[k_start];
          num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];

          int node_tmp;
          double age = coords_new;
          for(int k_tmp = k_start; k_tmp < k_end; k_tmp++){
            node_tmp = sorted_indices[k_tmp + 1];
            if(age > coordinates[node_tmp]){
              sorted_indices_new[k_tmp]  = node_tmp;
              order_new[node_tmp]        = k_tmp;
              num_lineages_new[node_tmp] = num_lineages[node_tmp] + 1;
            }else{
              sorted_indices_new[k_tmp] = node_k;
              order_new[node_k]         = k_tmp;
              num_lineages_new[node_k]  = num_lineages_new[sorted_indices_new[k_tmp-1]]-1; 
              k_start = k - 1;
              k_end   = k_tmp + 1;
              num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];
              sorted_indices_new[k_start]               = sorted_indices[k_start];
              order_new[sorted_indices[k_start]]        = k_start;
              num_lineages_new[sorted_indices[k_end]] = num_lineages[sorted_indices[k_end]];
              sorted_indices_new[k_end]               = sorted_indices[k_end];
              order_new[sorted_indices[k_end]]        = k_end;
              break;
            }
          } 

        }else{      
          k_end = k;
          k_start = order[child_left_label];
          if(k_start < order[child_right_label]) k_start = order[child_right_label];

          sorted_indices_new[k_start-1]               = sorted_indices[k_start-1];
          num_lineages_new[sorted_indices[k_start-1]] = num_lineages[sorted_indices[k_start-1]];
          sorted_indices_new[k_start]               = sorted_indices[k_start];
          num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];

          int node_tmp;
          double age = coords_new;
          for(int k_tmp = k_end; k_tmp > k_start; k_tmp--){
            node_tmp = sorted_indices[k_tmp - 1];
            if(age < coordinates[node_tmp]){
              sorted_indices_new[k_tmp]  = node_tmp;
              order_new[node_tmp]        = k_tmp;
              num_lineages_new[node_tmp] = num_lineages[node_tmp] - 1;
            }else{
              sorted_indices_new[k_tmp] = node_k;
              order_new[node_k]         = k_tmp;
              num_lineages_new[node_k]  = num_lineages[sorted_indices[k_tmp-1]]-1;
              k_start = k_tmp - 1;
              k_end   = k + 1;
              num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];
              sorted_indices_new[k_start]               = sorted_indices[k_start];
              order_new[sorted_indices[k_start]]        = k_start;
              num_lineages_new[sorted_indices[k_end]] = num_lineages[sorted_indices[k_end]];
              sorted_indices_new[k_end]               = sorted_indices[k_end];
              order_new[sorted_indices[k_end]]        = k_end;
              break;
            }
          }

          if(0){
            std::cerr << k_start << " " << k_end << " " << delta_tau << " " << age << std::endl;
            for(int kk = k_start; kk < k_end; kk++){
              std::cerr << num_lineages[sorted_indices[kk]] << " ";
            }
            std::cerr << std::endl;

            for(int kk = k_start; kk < k_end; kk++){
              std::cerr << num_lineages_new[sorted_indices[kk]] << " ";
            }
            std::cerr << std::endl;

            for(int kk = k_start; kk < k_end; kk++){
              std::cerr << sorted_indices[kk] << " ";
            }
            std::cerr << std::endl;

            for(int kk = k_start; kk < k_end; kk++){
              std::cerr << coordinates[sorted_indices[kk]] << " ";
            }
            std::cerr << std::endl;
          }

        }

        //for debugging
        if(0){

          std::vector<double> coordinates_new = coordinates;
          coordinates_new[node_k] += delta_tau; 
          std::vector<int> sorted_indices_new2 = sorted_indices, order_new2 = order, num_lineages_new2 = num_lineages;

          std::size_t m1(0);
          std::generate(std::begin(sorted_indices_new2), std::end(sorted_indices_new2), [&]{ return m1++; });
          std::sort(std::begin(sorted_indices_new2), std::end(sorted_indices_new2), [&](int i1, int i2) {
              return std::tie(coordinates_new[i1],i1) < std::tie(coordinates_new[i2],i2); } );

          //obtain order of coalescent events
          std::fill(order_new2.begin(), order_new2.end(), 0);
          std::size_t m2(0);
          std::generate(std::begin(order_new2), std::end(order_new2), [&]{ return m2++; });
          std::sort(std::begin(order_new2), std::end(order_new2), [&](int i1, int i2) { return sorted_indices_new2[i1] < sorted_indices_new2[i2]; } );

          int num_lins = 0;
          double ages = sample_age[*sorted_indices_new2.begin()];
          std::vector<int>::iterator it_sorted_indices_start = sorted_indices_new2.begin();
          for(it_sorted_indices = sorted_indices_new2.begin(); it_sorted_indices != sorted_indices_new2.end(); it_sorted_indices++){
            if(*it_sorted_indices >= N){
              for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
                num_lineages_new2[*it_sorted_indices_start] = num_lins; 
              }
              num_lins--;
              num_lineages_new2[*it_sorted_indices] = num_lins;
              it_sorted_indices_start++;
            }else if(ages < sample_age[*it_sorted_indices]){      
              for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
                num_lineages_new2[*it_sorted_indices_start] = num_lins; 
              }
              ages = sample_age[*it_sorted_indices];
              num_lins++;
            }else{
              num_lins++;
            }
          }



          for(int i = k_start; i <= k_end; i++){
            assert(sorted_indices_new[i] == sorted_indices_new2[i]);
            int node_tmp = sorted_indices_new[i];
            assert(order_new[node_tmp] == order_new2[node_tmp]);
          } 

        }

        coordinates[node_k] = coords_new;
        log_likelihood = CalculatePrior(k_start, k_end, coordinates, sorted_indices_new, num_lineages_new);
        coordinates[node_k] = coords;
        if(log_likelihood != -std::numeric_limits<float>::infinity()){
          //log_likelihood_ratio += log_likelihood;
          log_likelihood -= CalculatePrior(k_start, k_end, coordinates, sorted_indices, num_lineages);
          if(log_likelihood != -std::numeric_limits<float>::infinity()) log_likelihood_ratio += log_likelihood;
        }

        ///////////////////////////

        //likelihood function, likelihood ratio  
        n_num_events           = tree.nodes[node_k].num_events;
        child_left_num_events  = tree.nodes[child_left_label].num_events;
        child_right_num_events = tree.nodes[child_right_label].num_events;

        tb_child_left      = tree.nodes[child_left_label].branch_length;
        tb_child_right     = tree.nodes[child_right_label].branch_length;

        //mutation and recombination part
        if(tb == 0.0){
          log_likelihood_ratio  = std::numeric_limits<float>::infinity();
        }else if(tb <= delta_tau){
          log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
        }else{

          if(tb_child_left == 0.0){
            log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          }else if(tb_child_left <= -delta_tau){
            log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
          }else{

            if(tb_child_right == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
            }else if(tb_child_right <= -delta_tau){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
            }else{

              log_likelihood_ratio += (mut_rate[node_k] - mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
              if(n_num_events >= 1.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
              if(child_right_num_events >= 1.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
              if(child_left_num_events >= 1.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

            }
          }
        }
      }else{
        log_likelihood_ratio = 0.0;
        k_start = k;
        k_end = k_start;
        delta_tau = 0.0;
      }

      //decide whether to accept proposal or not.
      accept = true;
      if(log_likelihood_ratio < 0.0){
        //accept with probability exp(log_p)
        if(dist_unif(rng) > exp(log_likelihood_ratio)){
          accept = false;
        }
      }

      //update coordinates and sorted_indices
      if(accept){
        //calculate new branch lengths
        update_node3         = node_k;
        update_node1         = node_k;
        coordinates[node_k] += delta_tau;
        tree.nodes[child_left_label].branch_length  = coordinates[node_k] - coordinates[child_left_label];
        tree.nodes[child_right_label].branch_length = coordinates[node_k] - coordinates[child_right_label];
        tree.nodes[node_k].branch_length            = coordinates[parent_label] - coordinates[node_k];  

        for(int k_tmp = k_start + 1; k_tmp < k_end; k_tmp++){
          sorted_indices[k_tmp]               = sorted_indices_new[k_tmp];
          order[sorted_indices[k_tmp]]        = order_new[sorted_indices[k_tmp]];
          num_lineages[sorted_indices[k_tmp]] = num_lineages_new[sorted_indices[k_tmp]];
        }
      }

    }

  }

}


///////////////////// variable NE //////////////////////
//propose a new time and change events only up to t+bound
void
EstimateBranchLengthsWithSampleAge::ChangeTimeWhilekAncestorsVP_new(Tree& tree, int node, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif){

  //This step is O(N)
  int k = order[node];

  assert(node == sorted_indices[k]);
  //coorindates[sorted_indices[k-1]] determines the epoch at the lower end.
  //from there, I have to propose a new time by drawing a time for the first epoch, and if it is exceeding the epoch, for the next epoch etc.

  double age = coordinates[node];
  if(sorted_indices[k] < N){
    //std::cerr << node << " " << age << " " << sample_age[node] << std::endl;
    assert(age == sample_age[node]);
    while(sorted_indices[k] < N){ 
      k++;
      if(sorted_indices[k] < N){
        if(sample_age[sorted_indices[k]] != age) break;
      }
    }
    k--;
    assert(sorted_indices[k] < N);
  }
  assert(age == coordinates[sorted_indices[k]]);
  node      = sorted_indices[k];
  tau_old   = coordinates[sorted_indices[k+1]] - age;
  k_choose_2 = num_lineages[node] * (num_lineages[node]-1.0)/2.0;
  log_likelihood_ratio = 0.0;

  double min_tip = std::numeric_limits<float>::infinity();
  for(it_order = order.begin(); it_order != std::next(order.begin(), N); it_order++){
    if(*it_order > k){
      double bl = tree.nodes[sorted_indices[*it_order]].branch_length;
      if(bl < min_tip) min_tip = bl;
    }
  }

  if(min_tip > 0 && tau_old > 0.0){

    if(tau_old > 0.0){

      if(tau_old <= min_tip){
        tau_new   = -log(dist_unif(rng)) * tau_old;
        delta_tau = tau_new - tau_old;
        //calculate ratio of proposals
        log_likelihood_ratio = log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);   
      }else{
        tau_new   = -log(dist_unif(rng)) * (min_tip) + tau_old - min_tip;
        delta_tau = tau_new - tau_old;
        assert(min_tip+delta_tau >= 0.0);
        //calculate ratio of proposals
        log_likelihood_ratio = log(min_tip/(min_tip+delta_tau)) + ((min_tip + delta_tau)/min_tip) - (min_tip/(min_tip+delta_tau));
      }
      assert(tau_new > 0.0);

    }

    ////////////////////////////////////////////////
    //coalescent prior  

    //calculate new coordinates, order, sorted_indices, num_lineages

    double log_likelihood;

    std::vector<int> sorted_indices_new = sorted_indices, order_new = order, num_lineages_new = num_lineages;
    std::vector<double> coordinates_new = coordinates;

    for(it_sorted_indices = std::next(sorted_indices_new.begin(), k+1); it_sorted_indices != sorted_indices_new.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        coordinates_new[*it_sorted_indices] += delta_tau;
      }
    }

    std::size_t m1(0);
    std::generate(std::begin(sorted_indices_new), std::end(sorted_indices_new), [&]{ return m1++; });
    std::sort(std::begin(sorted_indices_new), std::end(sorted_indices_new), [&](int i1, int i2) {
        return std::tie(coordinates_new[i1],i1) < std::tie(coordinates_new[i2],i2); } );

    //obtain order of coalescent events
    std::fill(order_new.begin(), order_new.end(), 0);
    std::size_t m2(0);
    std::generate(std::begin(order_new), std::end(order_new), [&]{ return m2++; });
    std::sort(std::begin(order_new), std::end(order_new), [&](int i1, int i2) { return sorted_indices_new[i1] < sorted_indices_new[i2]; } );

    int num_lins = 0;
    double ages = sample_age[*sorted_indices_new.begin()];
    std::vector<int>::iterator it_sorted_indices_start = sorted_indices_new.begin();
    for(it_sorted_indices = sorted_indices_new.begin(); it_sorted_indices != sorted_indices_new.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages_new[*it_sorted_indices_start] = num_lins; 
        }
        num_lins--;
        num_lineages_new[*it_sorted_indices] = num_lins;
        it_sorted_indices_start++;
      }else if(ages < sample_age[*it_sorted_indices]){      
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages_new[*it_sorted_indices_start] = num_lins; 
        }
        ages = sample_age[*it_sorted_indices];
        num_lins++;
      }else{
        num_lins++;
      }
    }


    log_likelihood = CalculatePrior(epoch, coal_rate, coordinates_new, sorted_indices_new, num_lineages_new);
    if(log_likelihood != -std::numeric_limits<float>::infinity()){
      log_likelihood_ratio += log_likelihood;
      log_likelihood = CalculatePrior(epoch, coal_rate, coordinates, sorted_indices, num_lineages);
      if(log_likelihood != -std::numeric_limits<float>::infinity()) log_likelihood_ratio -= log_likelihood;
    }


    if(log_likelihood_ratio != -std::numeric_limits<float>::infinity()){

      ////////////////////////////////////////////
      //likelihood P(Data | tb)

      if(log_likelihood_ratio != std::numeric_limits<float>::infinity()){
        //assert(order[node_k] == k);
        int count_number_of_spanning_branches = 0;
        age = coordinates[node];
        it_sorted_indices = std::next(sorted_indices.begin(), k+1);
        for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){

          if(*it_sorted_indices < N){
            n = tree.nodes[*it_sorted_indices];
            tb     = n.branch_length;
            tb_new = tb + delta_tau;

            if(tb == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
              break;
            }else if(tb_new <= 0.0){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
              break;
            }else{
              log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
              log_likelihood_ratio += n.num_events * log(tb_new/tb);
            }
          }else{
            if(coordinates[(*tree.nodes[*it_sorted_indices].child_left).label] < age){
              count_number_of_spanning_branches++;

              n = *tree.nodes[*it_sorted_indices].child_left;
              assert(order[n.label] < k+1);
              tb     = n.branch_length;
              tb_new = tb + delta_tau;

              if(tb == 0.0){
                log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                break;
              }else if(tb_new <= 0.0){
                log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
                break;
              }else{
                log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
                log_likelihood_ratio += n.num_events * log(tb_new/tb);
              }
            }

            if(coordinates[(*tree.nodes[*it_sorted_indices].child_right).label] < age){
              count_number_of_spanning_branches++;

              n = *tree.nodes[*it_sorted_indices].child_right;
              assert(order[n.label] < k+1);
              tb     = n.branch_length;
              tb_new = tb + delta_tau;
              //assert(tb_new >= 0.0);

              if(tb == 0.0){
                log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                break;
              }else if(tb_new <= 0.0){
                log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
                break;
              }else{
                log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
                log_likelihood_ratio += n.num_events * log(tb_new/tb);
              }
            }

          }
          if(count_number_of_spanning_branches == num_lineages[node]) break;
        }
        //assert(count_number_of_spanning_branches == num_lineages);
      }
    }

    //decide whether to accept proposal or not.
    accept = true;
    if(log_likelihood_ratio < 0.0){
      //accept with probability exp(log_p)
      if(dist_unif(rng) > exp(log_likelihood_ratio)){
        accept = false;
      }
    }

    //update coordinates and sorted_indices
    if(accept){
      //calculate new branch lengths
      int k_init = k;
      k++;  
      update_node1 = k;
      coordinates = coordinates_new;

      //sorted_indices, order, num_lineages will change here
      for(it_sorted_indices = std::next(sorted_indices.begin(), k); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
        if(*it_sorted_indices >= N){
          child_left_label                                 = (*tree.nodes[*it_sorted_indices].child_left).label;
          tree.nodes[child_left_label].branch_length       = std::max(0.0, coordinates[*it_sorted_indices] - coordinates[child_left_label]);
          child_right_label                                = (*tree.nodes[*it_sorted_indices].child_right).label;
          tree.nodes[child_right_label].branch_length      = std::max(0.0, coordinates[*it_sorted_indices] - coordinates[child_right_label]);
          assert(tree.nodes[child_left_label].branch_length >= 0.0);
          assert(tree.nodes[child_right_label].branch_length >= 0.0);
        }
      }

      sorted_indices = sorted_indices_new;
      order = order_new;
      num_lineages = num_lineages_new;
    }

  }else{
    update_node1 = k;
  }

}

//This changes the time of one event, with a beta proposal within the time while k ancestors remain
void
EstimateBranchLengthsWithSampleAge::UpdateOneEventVP(Tree& tree, int node_k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif){

  accept = true;
  log_likelihood_ratio = 0.0;

  //propose to change time of node_k between time of older daughter and parent
  if(tree.nodes[node_k].parent == NULL){

    //propose with exponential distribution  
    tau_old   = coordinates[node_k] - coordinates[(*tree.nodes[node_k].child_left).label];
    if(tau_old > coordinates[node_k] - coordinates[(*tree.nodes[node_k].child_right).label]){
      tau_old = coordinates[node_k] - coordinates[(*tree.nodes[node_k].child_right).label];
    }
    if(tau_old > 0.0){
      tau_new   = -fast_log(dist_unif(rng)) * tau_old;
      delta_tau = tau_new - tau_old;
      assert(tau_new >= 0.0);

      //now decide whether to accept delta_tau:
      //calculate likelihood ratio

      //proposal likelihood ratio
      log_likelihood_ratio = fast_log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);
    }else{
      k_choose_2 = num_lineages[node_k] * (num_lineages[node_k]+1.0)/2.0;
      tau_new   = -log(dist_unif(rng)) * 1.0/k_choose_2;
      tau_old   = 0.0;
      delta_tau = tau_new;   
      //calculate ratio of proposals
      log_likelihood_ratio = fast_log(1.0/(tau_new*k_choose_2)) + tau_new*k_choose_2;
    }

    //coalescent prior

    //this case is unaffected by ancient DNA
    coordinates[node_k] += delta_tau;
    assert(order[node_k] == order.size() - 1);
    int k_end = order.size() - 1;
    int k_start = order.size() - 2;
    double log_likelihood = CalculatePrior(k_start, k_end, epoch, coal_rate, coordinates, sorted_indices, num_lineages);
    coordinates[node_k] -= delta_tau;
    if(log_likelihood != -std::numeric_limits<float>::infinity()){
      //log_likelihood_ratio += log_likelihood;
      log_likelihood -= CalculatePrior(k_start, k_end, epoch, coal_rate, coordinates, sorted_indices, num_lineages);
      if(log_likelihood != -std::numeric_limits<float>::infinity()) log_likelihood_ratio += log_likelihood;
    }

    //likelihood function, likelihood ratio  
    child_left_label  = (*tree.nodes[node_k].child_left).label;
    child_right_label = (*tree.nodes[node_k].child_right).label;

    child_left_num_events  = tree.nodes[child_left_label].num_events;
    child_right_num_events = tree.nodes[child_right_label].num_events;

    tb_child_left      = tree.nodes[child_left_label].branch_length;
    tb_child_right     = tree.nodes[child_right_label].branch_length;

    //mutation and recombination part

    if(tb_child_left == 0.0){
      log_likelihood_ratio  = std::numeric_limits<float>::infinity();
    }else if(tb_child_left <= -delta_tau){
      log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
    }else{

      if(tb_child_right == 0.0){
        log_likelihood_ratio  = std::numeric_limits<float>::infinity();
      }else if(tb_child_right <= -delta_tau){
        log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
      }else{

        log_likelihood_ratio += (- mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
        if(child_right_num_events >= 1.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
        if(child_left_num_events >= 1.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

      }
    }

    //decide whether to accept proposal or not.
    accept = true;
    if(log_likelihood_ratio < 0.0){
      //accept with probability exp(log_p)
      if(dist_unif(rng) > exp(log_likelihood_ratio)){
        accept = false;
      }
    }

    //update coordinates and sorted_indices
    if(accept){
      //calculate new branch lengths
      update_node3         = node_k;
      update_node1         = node_k;
      coordinates[node_k] += delta_tau;
      tree.nodes[child_left_label].branch_length  = coordinates[node_k] - coordinates[child_left_label];
      tree.nodes[child_right_label].branch_length = coordinates[node_k] - coordinates[child_right_label];
      //no change in num_lineages, order, sorted_indices
    }

  }else{

    //tau_below    = coordinates[node_k] - coordinates[sorted_indices[k-1]];
    //tau_above    = coordinates[sorted_indices[k+1]] - coordinates[node_k];

    child_left_label  = (*tree.nodes[node_k].child_left).label;
    child_right_label = (*tree.nodes[node_k].child_right).label;
    parent_label      = (*tree.nodes[node_k].parent).label;

    tb_child_left     = tree.nodes[child_left_label].branch_length;
    tb_child_right    = tree.nodes[child_right_label].branch_length;
    tb                = tree.nodes[node_k].branch_length;

    tau_below    = tb_child_left;
    if(tb_child_right < tau_below) tau_below = tb_child_right;
    tau_above    = tb;
    T            = tau_below + tau_above;

    int k_start, k_end, k = order[node_k];
    if(tau_above >= 0.0 || tau_below >= 0.0){

      //sample a new tau_below 
      tau_new_below        = dist_unif(rng) * T;
      delta_tau            = tau_new_below - tau_below;
      tau_new_above        = T - tau_new_below;

      if(tau_new_above >= 0.0 || tau_new_below >= 0.0){

        //now decide whether to accept delta_tau:

        //calculate likelihood ratio
        //proposal likelihood ratio 
        //P(old | new bl)/P(new | old bl)
        log_likelihood_ratio = 0.0;

        //calculate new coordinates, order, sorted_indices, num_lineages
        double log_likelihood;
        k_end = order[parent_label];

        double coords = coordinates[node_k];
        double coords_new = coords + delta_tau;
        if(coords_new > coordinates[parent_label]) coords_new = coordinates[parent_label];
        if(coords_new < coordinates[child_left_label]) coords_new = coordinates[child_left_label];
        if(coords_new < coordinates[child_right_label]) coords_new = coordinates[child_right_label];

        bool flag = false;

        if(1){

          //anything outside k_start and k_end is identical, so only need to update this bit
          if(delta_tau > 0){
            k_start = k;

            int node_tmp;
            double age = coords_new;

            sorted_indices_new[k_start-1]               = sorted_indices[k_start-1];
            num_lineages_new[sorted_indices[k_start-1]] = num_lineages[sorted_indices[k_start-1]];
            sorted_indices_new[k_start]                 = sorted_indices[k_start];
            num_lineages_new[sorted_indices[k_start]]   = num_lineages[sorted_indices[k_start]];

            for(int k_tmp = k_start; k_tmp < k_end; k_tmp++){
              node_tmp = sorted_indices[k_tmp + 1];
              if(age > coordinates[node_tmp]){
                sorted_indices_new[k_tmp]  = node_tmp;
                order_new[node_tmp]        = k_tmp;
                num_lineages_new[node_tmp] = num_lineages[node_tmp] + 1;
              }else{
                sorted_indices_new[k_tmp] = node_k;
                order_new[node_k]         = k_tmp;
                num_lineages_new[node_k]  = num_lineages_new[sorted_indices_new[k_tmp-1]]-1; 

                k_start = k - 1;
                k_end   = k_tmp + 1;
                num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];
                sorted_indices_new[k_start]               = sorted_indices[k_start];
                order_new[sorted_indices[k_start]]        = k_start;
                num_lineages_new[sorted_indices[k_end]] = num_lineages[sorted_indices[k_end]];
                sorted_indices_new[k_end]               = sorted_indices[k_end];
                order_new[sorted_indices[k_end]]        = k_end;
                break;
              }
            }
          }else{      
            k_end = k;
            k_start = order[child_left_label];
            if(k_start < order[child_right_label]) k_start = order[child_right_label]; 

            sorted_indices_new[k_start-1]               = sorted_indices[k_start-1];
            num_lineages_new[sorted_indices[k_start-1]] = num_lineages[sorted_indices[k_start-1]];
            sorted_indices_new[k_start]               = sorted_indices[k_start];
            num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];

            int node_tmp;
            double age = coords_new;
            for(int k_tmp = k_end; k_tmp > k_start; k_tmp--){
              node_tmp = sorted_indices[k_tmp - 1];
              if(age < coordinates[node_tmp]){
                sorted_indices_new[k_tmp]  = node_tmp;
                order_new[node_tmp]        = k_tmp;
                num_lineages_new[node_tmp] = num_lineages[node_tmp] - 1;
              }else{
                sorted_indices_new[k_tmp] = node_k;
                order_new[node_k]         = k_tmp;
                num_lineages_new[node_k]  = num_lineages[sorted_indices[k_tmp-1]]-1;
                k_start = k_tmp - 1;
                k_end   = k + 1;
                num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];
                sorted_indices_new[k_start]               = sorted_indices[k_start];
                order_new[sorted_indices[k_start]]        = k_start;
                num_lineages_new[sorted_indices[k_end]] = num_lineages[sorted_indices[k_end]];
                sorted_indices_new[k_end]               = sorted_indices[k_end];
                order_new[sorted_indices[k_end]]        = k_end;
                break;
              }
            }

          }

          coordinates[node_k] = coords_new;
          log_likelihood = CalculatePrior(k_start, k_end, epoch, coal_rate, coordinates, sorted_indices_new, num_lineages_new);
          coordinates[node_k] = coords;
          if(log_likelihood != -std::numeric_limits<float>::infinity()){
            //log_likelihood_ratio += log_likelihood;
            log_likelihood -= CalculatePrior(k_start, k_end, epoch, coal_rate, coordinates, sorted_indices, num_lineages);
            if(log_likelihood != -std::numeric_limits<float>::infinity()) log_likelihood_ratio += log_likelihood;
          }
        }else{
          //need sorted_indices_new, num_lineages_new, order_new

          coordinates[node_k] = coords_new;
          sorted_indices_new = sorted_indices;
          num_lineages_new   = num_lineages;
          order_new          = order;

          std::size_t m1(0);
          std::generate(std::begin(sorted_indices_new), std::end(sorted_indices_new), [&]{ return m1++; });
          std::sort(std::begin(sorted_indices_new), std::end(sorted_indices_new), [&](int i1, int i2) {
              return std::tie(coordinates[i1],i1) < std::tie(coordinates[i2],i2); } );

          //obtain order of coalescent events
          std::fill(order_new.begin(), order_new.end(), 0);
          std::size_t m2(0);
          std::generate(std::begin(order_new), std::end(order_new), [&]{ return m2++; });
          std::sort(std::begin(order_new), std::end(order_new), [&](int i1, int i2) { return sorted_indices_new[i1] < sorted_indices_new[i2]; } );

          ////////////////////////////////
          int num_lins = 0;
          double ages = sample_age[*sorted_indices_new.begin()];
          std::vector<int>::iterator it_sorted_indices_start = sorted_indices_new.begin();
          for(it_sorted_indices = sorted_indices_new.begin(); it_sorted_indices != sorted_indices_new.end(); it_sorted_indices++){
            if(*it_sorted_indices >= N){
              for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
                num_lineages_new[*it_sorted_indices_start] = num_lins; 
              }
              num_lins--;
              num_lineages_new[*it_sorted_indices] = num_lins;
              it_sorted_indices_start++;
            }else if(ages < sample_age[*it_sorted_indices]){      
              for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
                num_lineages_new[*it_sorted_indices_start] = num_lins; 
              }
              ages = sample_age[*it_sorted_indices];
              num_lins++;
            }else{
              num_lins++;
            }
          }

          log_likelihood = CalculatePrior(epoch, coal_rate, coordinates, sorted_indices_new, num_lineages_new);
          coordinates[node_k] = coords;
          if(log_likelihood != -std::numeric_limits<float>::infinity()){
            //log_likelihood_ratio += log_likelihood;
            log_likelihood -= CalculatePrior(epoch, coal_rate, coordinates, sorted_indices, num_lineages);
            if(log_likelihood != -std::numeric_limits<float>::infinity()) log_likelihood_ratio += log_likelihood;
          } 

        }

        ///////////////////////////

        if(1){

          //likelihood function, likelihood ratio  
          n_num_events           = tree.nodes[node_k].num_events;
          child_left_num_events  = tree.nodes[child_left_label].num_events;
          child_right_num_events = tree.nodes[child_right_label].num_events;

          tb_child_left      = tree.nodes[child_left_label].branch_length;
          tb_child_right     = tree.nodes[child_right_label].branch_length;

          //mutation and recombination part
          if(tb == 0.0){
            log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          }else if(tb <= delta_tau){
            log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          }else{

            if(tb_child_left == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
            }else if(tb_child_left <= -delta_tau){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
            }else{

              if(tb_child_right == 0.0){
                log_likelihood_ratio  = std::numeric_limits<float>::infinity();
              }else if(tb_child_right <= -delta_tau){
                log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
              }else{

                log_likelihood_ratio += (mut_rate[node_k] - mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
                if(n_num_events >= 1.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
                if(child_right_num_events >= 1.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
                if(child_left_num_events >= 1.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

              }
            }
          }

        }

      }else{
        log_likelihood_ratio = 0.0;
        k_start = k;
        k_end = k_start;
        delta_tau = 0.0;
      }

      //decide whether to accept proposal or not.
      accept = true;
      if(log_likelihood_ratio < 0.0){
        //accept with probability exp(log_p)
        if(dist_unif(rng) > exp(log_likelihood_ratio)){
          accept = false;
        }
      }

      //update coordinates and sorted_indices
      if(accept){
        //calculate new branch lengths
        update_node3         = node_k;
        update_node1         = node_k;
        coordinates[node_k] += delta_tau;
        tree.nodes[child_left_label].branch_length  = coordinates[node_k] - coordinates[child_left_label];
        tree.nodes[child_right_label].branch_length = coordinates[node_k] - coordinates[child_right_label];
        tree.nodes[node_k].branch_length            = coordinates[parent_label] - coordinates[node_k];  

        for(int k_tmp = k_start; k_tmp < k_end; k_tmp++){
          sorted_indices[k_tmp]               = sorted_indices_new[k_tmp];
          order[sorted_indices[k_tmp]]        = order_new[sorted_indices[k_tmp]];
          num_lineages[sorted_indices[k_tmp]] = num_lineages_new[sorted_indices[k_tmp]];
        }
      }

    }

  }

}


///////////////////////////////////////////

void
EstimateBranchLengthsWithSampleAge::GetCoordinates(Node& n, std::vector<double>& coords){

  if(n.child_left != NULL){
    GetCoordinates(*n.child_left, coords);
    GetCoordinates(*n.child_right, coords);

    assert((*n.child_left).branch_length >= 0.0);

    coords[n.label] = std::max(coords[(*n.child_right).label] + (*n.child_right).branch_length, coords[(*n.child_left).label] + (*n.child_left).branch_length);
  }else{
    assert(n.label < (coords.size() + 1)/2.0);
    coords[n.label] = sample_age[n.label];
  }

}

void
EstimateBranchLengthsWithSampleAge::MCMC(const Data& data, Tree& tree, const int seed){

  //count_accept = 0;
  //count_proposal = 0;

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 0.6;

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  InitializeMCMC(data, tree); //Initialize using coalescent prior 
  InitializeOrder(tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) 10*data.N * data.N; j++){
    RandomSwitchOrder(tree, dist_n(rng), dist_unif);
  }
  //Initialise branch lengths
  InitializeBranchLengths(tree);

  for(int i = 0; i < N_total-1; i++){
    assert(tree.nodes[i].branch_length >= 0.0);
    assert(order[sorted_indices[i]] == i);
    assert(order[i] < order[(*tree.nodes[i].parent).label]);
    assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
  }

  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 1000*delta; count++){

    for(int i = 0; i < N_total-1; i++){
      assert(tree.nodes[i].branch_length >= 0.0);
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
      assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
    }

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p1/N){
      //std::cerr << "v1" << std::endl;
      ChangeTimeWhilekAncestors_new(tree, dist_tip(rng), dist_unif);
    }else if(uniform_rng <= p1){
      //std::cerr << "v2" << std::endl;
      ChangeTimeWhilekAncestors_new(tree, dist_n(rng), dist_unif);
    }else if(uniform_rng <= p2){
      //std::cerr << "v3" << std::endl;
      UpdateOneEvent(tree, dist_oneevent(rng), dist_gamma, dist_unif);
    }else{ 
      //std::cerr << "v4" << std::endl;
      SwitchOrder(tree, dist_n(rng), dist_unif);
    }      

    for(int k = 0; k < N_total; k++){
      assert(tree.nodes[k].branch_length >= 0.0);
    }

    for(int k = N; k < N_total; k++){
      assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_left).label]));
      assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_right).label]));
      assert(sorted_indices[order[k]] == k);
    }
  }

  avg              = coordinates;
  last_coordinates = coordinates;
  last_update.resize(N_total);
  std::fill(last_update.begin(), last_update.end(), 1);
  count = 1;

  int num_iterations = 0, iterations_threshold = 500*log(data.N);
  bool is_count_threshold = false;
  std::vector<int> count_proposals(N_total-N, 0);
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < p1/N){
        //std::cerr << "v1" << std::endl;
        //update_time_while_k_anc++;
        int k_candidate = dist_tip(rng);
        ChangeTimeWhilekAncestors_new(tree, k_candidate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng < p1){
        //std::cerr << "v2" << std::endl;
        //update_time_while_k_anc++;
        int k_candidate = dist_n(rng);
        ChangeTimeWhilekAncestors_new(tree, k_candidate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        int k_candidate = dist_oneevent(rng);
        count_proposals[k_candidate-N]++;
        UpdateOneEvent(tree, k_candidate, dist_gamma, dist_unif);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
        UpdateAvg(tree);
      }

      for(int k = 0; k < N_total; k++){
        assert(tree.nodes[k].branch_length >= 0.0);
      }

      for(int k = N; k < N_total; k++){
        assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_left).label]));
        assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_right).label]));
        assert(sorted_indices[order[k]] == k);
      }

    }while(count % delta != 0 );

    //std::cerr << num_iterations << std::endl;
    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    if(1){
      is_avg_increasing = true;
      if(!is_count_threshold){
        //assert(is_avg_increasing);
        for(std::vector<int>::iterator it_count = count_proposals.begin(); it_count != count_proposals.end(); it_count++){
          if(*it_count < 100){
            is_avg_increasing = false;
            break;
          }
        }
        if(is_avg_increasing) is_count_threshold = true;
      }
    }else{
      is_avg_increasing = true; 
      if(num_iterations < iterations_threshold) is_avg_increasing = false;     
    }

    //std::cerr << num_iterations << " " << is_avg_increasing << std::endl;

    if(is_avg_increasing){
      //update all nodes in avg 
      it_avg = std::next(avg.begin(), N);
      it_coords = std::next(coordinates.begin(), N);
      it_last_update = std::next(last_update.begin(), N);
      it_last_coords = std::next(last_coordinates.begin(), N);
      //update avg
      for(int ell = N; ell < N_total; ell++){
        *it_avg  += ((count - *it_last_update) * (*it_last_coords - *it_avg))/count;
        *it_last_update = count;
        *it_last_coords = *it_coords;
        it_avg++;
        it_coords++;
        it_last_update++;
        it_last_coords++;
      }
      //check if coalescent ages are non-decreasing in avg
      int ell = N;
      for(it_avg = std::next(avg.begin(),N); it_avg != avg.end(); it_avg++){
        if(ell < root){
          if(*it_avg > avg[(*tree.nodes[ell].parent).label]){
            is_avg_increasing = false;
            break;
          }
        }
        ell++;
      }
    }

  }

  for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
    (*it_n).branch_length = ((double) Ne) * (avg[(*(*it_n).parent).label] - avg[(*it_n).label]);
  }

  /*
     for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
     (*it_n).branch_length = ((double) Ne) * (*it_n).branch_length;
     }
     */

}  

void
EstimateBranchLengthsWithSampleAge::MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 0.6;

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  InitializeMCMC(data, tree); //Initialize using coalescent prior 

  if(0){

    InitializeOrder(tree); 

    //Randomly switch around order of coalescences
    for(int j = 0; j < (int) 10*data.N * data.N; j++){
      RandomSwitchOrder(tree, dist_n(rng), dist_unif);
    }
    //Initialise branch lengths
    InitializeBranchLengths(tree);

    for(int i = 0; i < N_total-1; i++){
      assert(tree.nodes[i].branch_length >= 0.0);
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
      assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
      assert(num_lineages[sorted_indices[i]] >= 1);
    }

  }else{

    for(std::vector<Node>::iterator it_node = tree.nodes.begin(); it_node != tree.nodes.end(); it_node++){
      (*it_node).branch_length /= data.Ne;
    }

    coordinates.resize(N_total);
    GetCoordinates(tree.nodes[root], coordinates);
    //avg   = coordinates;

    //as a sanity check that branch lengths are consistent with sample_ages, calculate branch lengths to root + sample_age
    double dist_to_root, dist_to_root_0 = 0.0;
    dist_to_root = sample_age[0];
    int n = 0;
    while(tree.nodes[n].parent != NULL){
      dist_to_root_0 += tree.nodes[n].branch_length;
      n = (*tree.nodes[n].parent).label;
    }
    for(int i = 1; i < N; i++){
      dist_to_root = sample_age[i];
      n = i;
      while(tree.nodes[n].parent != NULL){
        dist_to_root += tree.nodes[n].branch_length;
        n = (*tree.nodes[n].parent).label;
      }
      assert(dist_to_root - dist_to_root_0 < 1e-1);
    }


    std::size_t m1(0);
    std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
    std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
        return std::tie(coordinates[i1],i1) < std::tie(coordinates[i2],i2); } );

    //obtain order of coalescent events
    std::fill(order.begin(), order.end(), 0);
    std::size_t m2(0);
    std::generate(std::begin(order), std::end(order), [&]{ return m2++; });
    std::sort(std::begin(order), std::end(order), [&](int i1, int i2) { return sorted_indices[i1] < sorted_indices[i2]; } );

    ////////////////////////////////

    int num_lins = 0;
    double ages = sample_age[*sorted_indices.begin()];
    std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
    for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages[*it_sorted_indices_start] = num_lins; 
        }
        num_lins--;
        num_lineages[*it_sorted_indices] = num_lins;
        it_sorted_indices_start++;
      }else if(ages < sample_age[*it_sorted_indices]){      
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages[*it_sorted_indices_start] = num_lins; 
        }
        ages = sample_age[*it_sorted_indices];
        num_lins++;
      }else{
        num_lins++;
      }
    }

    sorted_indices_new = sorted_indices;
    order_new          = order;
    num_lineages_new   = num_lineages;

    //debug
    for(int i = 0; i < N_total-1; i++){
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
      assert(num_lineages[sorted_indices[i]] >= 1);
    }

  }



  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 1000*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p1/N){
      //std::cerr << "v1" << std::endl;
      ChangeTimeWhilekAncestorsVP_new(tree, dist_tip(rng), epoch, coal_rate, dist_unif);
      //ChangeTimeWhilekAncestors_new(tree, dist_tip(rng), dist_unif);
    }else if(uniform_rng <= p1){
      //std::cerr << "v2" << std::endl;
      ChangeTimeWhilekAncestorsVP_new(tree, dist_n(rng), epoch, coal_rate, dist_unif);
      //ChangeTimeWhilekAncestors_new(tree, dist_n(rng), dist_unif);
    }else if(uniform_rng <= p2){
      //std::cerr << "v3" << std::endl;
      UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
      //UpdateOneEvent(tree, dist_oneevent(rng), dist_gamma, dist_unif);
    }else{ 
      //std::cerr << "v4" << std::endl;
      SwitchOrder(tree, dist_n(rng), dist_unif);
    }    

    for(int k = 0; k < N_total; k++){
      assert(tree.nodes[k].branch_length >= 0.0);
    }

    for(int k = N; k < N_total; k++){
      assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_left).label]));
      assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_right).label]));
      //std::cerr << k << " " << order[k] << " " << sorted_indices[order[k]] << " " << coordinates[k] << " " << coordinates[sorted_indices[order[k]]] << std::endl;
      assert(sorted_indices[order[k]] == k);
    }
  }

  avg              = coordinates;
  last_coordinates = coordinates;
  last_update.resize(N_total);
  std::fill(last_update.begin(), last_update.end(), 1);
  count = 1;

  int num_iterations = 0, iterations_threshold = 500*log(data.N);
  bool is_count_threshold = false;
  std::vector<int> count_proposals(N_total-N, 0);
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < p1/N){
        //std::cerr << "v1" << std::endl;
        //update_time_while_k_anc++;
        int k_candidate = dist_tip(rng);
        ChangeTimeWhilekAncestorsVP_new(tree, k_candidate, epoch, coal_rate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng < p1){
        //std::cerr << "v2" << std::endl;
        //update_time_while_k_anc++;
        int k_candidate = dist_n(rng);
        ChangeTimeWhilekAncestorsVP_new(tree, k_candidate, epoch, coal_rate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        int k_candidate = dist_oneevent(rng);
        count_proposals[k_candidate-N]++;
        UpdateOneEventVP(tree, k_candidate, epoch, coal_rate, dist_gamma, dist_unif);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
        UpdateAvg(tree);
      }

      for(int k = 0; k < N_total; k++){
        assert(tree.nodes[k].branch_length >= 0.0);
      }

      for(int k = N; k < N_total; k++){
        assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_left).label]));
        assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_right).label]));
        //std::cerr << k << " " << order[k] << " " << sorted_indices[order[k]] << " " << coordinates[k] << " " << coordinates[sorted_indices[order[k]]] << std::endl;
        assert(sorted_indices[order[k]] == k);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    if(1){
      is_avg_increasing = true;
      if(!is_count_threshold){
        //assert(is_avg_increasing);
        for(std::vector<int>::iterator it_count = count_proposals.begin(); it_count != count_proposals.end(); it_count++){
          if(*it_count < 100){
            is_avg_increasing = false;
            break;
          }
        }
        if(is_avg_increasing) is_count_threshold = true;
      }
    }else{
      is_avg_increasing = true; 
      if(num_iterations < iterations_threshold) is_avg_increasing = false;     
    }

    if(is_avg_increasing){
      //update all nodes in avg 
      it_avg = std::next(avg.begin(), N);
      it_coords = std::next(coordinates.begin(), N);
      it_last_update = std::next(last_update.begin(), N);
      it_last_coords = std::next(last_coordinates.begin(), N);
      //update avg
      for(int ell = N; ell < N_total; ell++){
        *it_avg  += ((count - *it_last_update) * (*it_last_coords - *it_avg))/count;
        *it_last_update = count;
        *it_last_coords = *it_coords;
        it_avg++;
        it_coords++;
        it_last_update++;
        it_last_coords++;
      }
      //check if coalescent ages are non-decreasing in avg
      int ell = N;
      for(it_avg = std::next(avg.begin(),N); it_avg != avg.end(); it_avg++){
        if(ell < root){
          if(*it_avg > avg[(*tree.nodes[ell].parent).label]){
            is_avg_increasing = false;
            break;
          }
        }
        ell++;
      }
    }

  }

  for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
    (*it_n).branch_length = ((double) Ne) * (avg[(*(*it_n).parent).label] - avg[(*it_n).label]);
  }

  /*
     for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
     (*it_n).branch_length = ((double) Ne) * (*it_n).branch_length;
     }
     */

}  

void
EstimateBranchLengthsWithSampleAge::MCMCVariablePopulationSizeForRelate(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 0.6;

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  InitializeMCMC(data, tree); //Initialize using coalescent prior 
  InitializeOrder(tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) 10*data.N * data.N; j++){
    RandomSwitchOrder(tree, dist_n(rng), dist_unif);
  }
  //Initialise branch lengths
  InitializeBranchLengths(tree);

  for(int i = 0; i < N_total-1; i++){
    assert(tree.nodes[i].branch_length >= 0.0);
    assert(order[sorted_indices[i]] == i);
    assert(order[i] < order[(*tree.nodes[i].parent).label]);
    assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
  }

  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 1000*delta; count++){

    for(int i = 0; i < N_total-1; i++){
      assert(tree.nodes[i].branch_length >= 0.0);
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
      assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
    }

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p1/N){
      //std::cerr << "v1" << std::endl;
      ChangeTimeWhilekAncestorsVP_new(tree, dist_tip(rng), epoch, coal_rate, dist_unif);
    }else if(uniform_rng <= p1){
      //std::cerr << "v2" << std::endl;
      ChangeTimeWhilekAncestorsVP_new(tree, dist_n(rng), epoch, coal_rate, dist_unif);
    }else if(uniform_rng <= p2){
      //std::cerr << "v3" << std::endl;
      UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
    }else{ 
      //std::cerr << "v4" << std::endl;
      SwitchOrder(tree, dist_n(rng), dist_unif);
    }    

    for(int k = 0; k < N_total; k++){
      assert(tree.nodes[k].branch_length >= 0.0);
    }

    for(int k = N; k < N_total; k++){
      //std::cerr << k << " " << order[k] << " " << sorted_indices[order[k]] << " " << coordinates[k] << " " << coordinates[sorted_indices[order[k]]] << " " <<  << std::endl;
      assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_left).label]));
      assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_right).label]));
      assert(sorted_indices[order[k]] == k);
    }
  }

  avg              = coordinates;
  last_coordinates = coordinates;
  last_update.resize(N_total);
  std::fill(last_update.begin(), last_update.end(), 1);
  count = 1;

  int num_iterations = 0, iterations_threshold = 500*log(data.N);
  bool is_count_threshold = false;
  std::vector<int> count_proposals(N_total-N, 0);
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < p1/N){
        //std::cerr << "v1" << std::endl;
        //update_time_while_k_anc++;
        int k_candidate = dist_tip(rng);
        ChangeTimeWhilekAncestorsVP_new(tree, k_candidate, epoch, coal_rate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng < p1){
        //std::cerr << "v2" << std::endl;
        //update_time_while_k_anc++;
        int k_candidate = dist_n(rng);
        ChangeTimeWhilekAncestorsVP_new(tree, k_candidate, epoch, coal_rate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        int k_candidate = dist_oneevent(rng);
        count_proposals[k_candidate-N]++;
        UpdateOneEventVP(tree, k_candidate, epoch, coal_rate, dist_gamma, dist_unif);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
        UpdateAvg(tree);
      }

      for(int k = 0; k < N_total; k++){
        assert(tree.nodes[k].branch_length >= 0.0);
      }

      for(int k = N; k < N_total; k++){
        assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_left).label]));
        assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_right).label]));
        //std::cerr << k << " " << order[k] << " " << sorted_indices[order[k]] << " " << coordinates[k] << " " << coordinates[sorted_indices[order[k]]] << std::endl;
        assert(sorted_indices[order[k]] == k);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    if(1){
      is_avg_increasing = true;
      if(!is_count_threshold){
        //assert(is_avg_increasing);
        for(std::vector<int>::iterator it_count = count_proposals.begin(); it_count != count_proposals.end(); it_count++){
          if(*it_count < 100){
            is_avg_increasing = false;
            break;
          }
        }
        if(is_avg_increasing) is_count_threshold = true;
      }
    }else{
      is_avg_increasing = true; 
      if(num_iterations < iterations_threshold) is_avg_increasing = false;     
    }


    if(is_avg_increasing){
      //update all nodes in avg 
      it_avg = std::next(avg.begin(), N);
      it_coords = std::next(coordinates.begin(), N);
      it_last_update = std::next(last_update.begin(), N);
      it_last_coords = std::next(last_coordinates.begin(), N);
      //update avg
      for(int ell = N; ell < N_total; ell++){
        *it_avg  += ((count - *it_last_update) * (*it_last_coords - *it_avg))/count;
        *it_last_update = count;
        *it_last_coords = *it_coords;
        it_avg++;
        it_coords++;
        it_last_update++;
        it_last_coords++;
      }
      //check if coalescent ages are non-decreasing in avg
      int ell = N;
      for(it_avg = std::next(avg.begin(),N); it_avg != avg.end(); it_avg++){
        if(ell < root){
          if(*it_avg > avg[(*tree.nodes[ell].parent).label]){
            is_avg_increasing = false;
            break;
          }
        }
        ell++;
      }
    }

  }

  for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
    (*it_n).branch_length = ((double) Ne) * (avg[(*(*it_n).parent).label] - avg[(*it_n).label]);
  }

  /*
     for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
     (*it_n).branch_length = ((double) Ne) * (*it_n).branch_length;
     }
     */

}  

void
EstimateBranchLengthsWithSampleAge::MCMCVariablePopulationSizeSample(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, int num_proposals, const bool init, const int seed){

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 0.6;

  if(init == 1){
    rng.seed(seed);
    root = N_total - 1;

    InitializeMCMC(data, tree); 
    coordinates.resize(N_total);

    //as a sanity check that branch lengths are consistent with sample_ages, calculate branch lengths to root + sample_age
    double dist_to_root, dist_to_root_0 = 0.0;
    dist_to_root = sample_age[0];
    int n = 0;
    while(tree.nodes[n].parent != NULL){
      dist_to_root_0 += tree.nodes[n].branch_length;
      n = (*tree.nodes[n].parent).label;
    }
    for(int i = 1; i < N; i++){
      dist_to_root = sample_age[i];
      n = i;
      while(tree.nodes[n].parent != NULL){
        dist_to_root += tree.nodes[n].branch_length;
        n = (*tree.nodes[n].parent).label;
      }
      assert(dist_to_root - dist_to_root_0 < 1e-1);
    }

    GetCoordinates(tree.nodes[root], coordinates);

    std::size_t m1(0);
    std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
    std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
        return std::tie(coordinates[i1],i1) < std::tie(coordinates[i2],i2); } );

    //obtain order of coalescent events
    std::fill(order.begin(), order.end(), 0);
    std::size_t m2(0);
    std::generate(std::begin(order), std::end(order), [&]{ return m2++; });
    std::sort(std::begin(order), std::end(order), [&](int i1, int i2) { return sorted_indices[i1] < sorted_indices[i2]; } );

    ////////////////////////////////
    int num_lins = 0;
    double ages = sample_age[*sorted_indices.begin()];
    std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
    for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      if(*it_sorted_indices >= N){
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages[*it_sorted_indices_start] = num_lins; 
        }
        num_lins--;
        num_lineages[*it_sorted_indices] = num_lins;
        it_sorted_indices_start++;
      }else if(ages < sample_age[*it_sorted_indices]){      
        for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
          num_lineages[*it_sorted_indices_start] = num_lins; 
        }
        ages = sample_age[*it_sorted_indices];
        num_lins++;
      }else{
        num_lins++;
      }
    }


    sorted_indices_new = sorted_indices;
    order_new          = order;
    num_lineages_new   = num_lineages;

    //debug
    for(int i = 0; i < N_total-1; i++){
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
    }
  }

  if(0){
    ///////////////
    int delta = std::max(data.N/10.0, 10.0);
    root = N_total - 1;

    InitializeMCMC(data, tree); 
    InitializeOrder(tree); 

    if(1){
      //Randomly switch around order of coalescences
      for(int j = 0; j < (int) 10*data.N * data.N; j++){
        RandomSwitchOrder(tree, dist_n(rng), dist_unif);
      }
    }
    //Initialise branch lengths
    InitializeBranchLengths(tree);
  }

  ////////////////// Sample branch lengths /////////////////

  count = 0;
  for(; count < num_proposals; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= 0.5*p1){
      ChangeTimeWhilekAncestorsVP_new(tree, dist_tip(rng), epoch, coal_rate, dist_unif);
    }else if(uniform_rng <= p1){
      ChangeTimeWhilekAncestorsVP_new(tree, dist_n(rng), epoch, coal_rate, dist_unif);
    }else if(uniform_rng <= p2){
      UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
    }else{ 
      SwitchOrder(tree, dist_n(rng), dist_unif);
    }

    if(0){
      std::vector<int> sorted_indices_foo = sorted_indices;
      std::vector<int> order_foo = order;
      std::vector<int> num_lineages_foo = num_lineages;

      std::size_t m1(0);
      std::generate(std::begin(sorted_indices_foo), std::end(sorted_indices_foo), [&]{ return m1++; });
      std::sort(std::begin(sorted_indices_foo), std::end(sorted_indices_foo), [&](int i1, int i2) {
          return std::tie(coordinates[i1],i1) < std::tie(coordinates[i2],i2); } );

      //obtain order of coalescent events
      std::fill(order_foo.begin(), order_foo.end(), 0);
      std::size_t m2(0);
      std::generate(std::begin(order_foo), std::end(order_foo), [&]{ return m2++; });
      std::sort(std::begin(order_foo), std::end(order_foo), [&](int i1, int i2) { return sorted_indices_foo[i1] < sorted_indices_foo[i2]; } );

      ////////////////////////////////
      int num_lins = 0;
      double ages = sample_age[*sorted_indices_foo.begin()];
      std::vector<int>::iterator it_sorted_indices_start = sorted_indices_foo.begin();
      for(it_sorted_indices = sorted_indices_foo.begin(); it_sorted_indices != sorted_indices_foo.end(); it_sorted_indices++){
        if(*it_sorted_indices >= N){
          for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
            num_lineages_foo[*it_sorted_indices_start] = num_lins; 
          }
          num_lins--;
          num_lineages_foo[*it_sorted_indices] = num_lins;
          it_sorted_indices_start++;
        }else if(ages < sample_age[*it_sorted_indices]){      
          for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
            num_lineages_foo[*it_sorted_indices_start] = num_lins; 
          }
          ages = sample_age[*it_sorted_indices];
          num_lins++;
        }else{
          num_lins++;
        }
      }

      for(int i = 0; i < N_total; i++){
        if(!(order[i] == order_foo[i])) std::cerr << order[i] << " " << order_foo[i] << std::endl;
        assert(sorted_indices[i] == sorted_indices_foo[i]);
        assert(order[i] == order_foo[i]);
        assert(num_lineages[i] == num_lineages_foo[i]);
      }

      sorted_indices = sorted_indices_foo;
      order = order_foo;
      num_lineages = num_lineages_foo;
    }

  }

}  


