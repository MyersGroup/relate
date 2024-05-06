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

EstimateBranchLengthsWithSampleAge::EstimateBranchLengthsWithSampleAge(const Data& data){
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

  sample_age.resize(N);
  std::fill(sample_age.begin(), sample_age.end(), 0);

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
        coordinates[node_i] = (tmp - coordinates[sorted_indices[i-1]])/10.0 + coordinates[sorted_indices[i-1]];
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
    //std::cerr << i << "," << coordinates[sorted_indices[i+1]] << " ";
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
  epsilon /= 10.0;
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

    assert(snp_end < data.dist.size());
    mut_rate[i]            = 0.0;
    for(int snp = snp_begin; snp < snp_end; snp++){
      mut_rate[i]         += data.dist[snp];
    }

    if(snp_begin > 0){
      snp_begin--;
      mut_rate[i]         += 0.5 * data.dist[snp_begin];
    }
    if(snp_end < data.L-1){
      mut_rate[i]         += 0.5 * data.dist[snp_end];
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


void
EstimateBranchLengthsWithSampleAge::SwitchTopo(Tree& tree, std::vector<Leaves>& desc, const std::vector<double>& epoch, std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, int node_k, std::uniform_real_distribution<double>& dist_unif){

  accept = true;
  int node_swap_k;
  log_likelihood_ratio = 0.0;

  int k = order[node_k];

  int parent = (*tree.nodes[node_k].parent).label;
  int sib  = (*tree.nodes[parent].child_left).label;
  if(sib == node_k) sib  = (*tree.nodes[parent].child_right).label;
  int child_left = (*tree.nodes[node_k].child_left).label;
  int child_right = (*tree.nodes[node_k].child_right).label;

  if(node_k >= N && order[sib] < order[node_k] && tree.nodes[node_k].num_events == 0.0){

    float bl_child_left  = tree.nodes[child_left].branch_length;
    float bl_sib         = tree.nodes[sib].branch_length;
    float bl_node_k      = tree.nodes[node_k].branch_length;
    float bl_child_right = tree.nodes[child_right].branch_length;

    //order, sorted_indices are unchanged
    //so just update remaining
    int k_start = order[node_k];
    int k_end   = order[parent];
    assert(k_start < k_end);
    log_likelihood_ratio = -CalculatePrior(k_start, k_end, tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);

    if(dist_unif(rng) < 0.5){

      //merge sib with child_left
      //then child_right
      tree.nodes[child_left].parent = &tree.nodes[node_k];
      tree.nodes[sib].parent = &tree.nodes[node_k];
      tree.nodes[node_k].child_left = &tree.nodes[child_left];
      tree.nodes[node_k].child_right = &tree.nodes[sib];

      tree.nodes[node_k].parent = &tree.nodes[parent];
      tree.nodes[child_right].parent = &tree.nodes[parent];
      tree.nodes[parent].child_left = &tree.nodes[node_k];
      tree.nodes[parent].child_right = &tree.nodes[child_right];

      tree.nodes[child_left].branch_length   = coordinates[node_k] - coordinates[child_left]; 
      tree.nodes[sib].branch_length          = coordinates[node_k] - coordinates[sib]; 
      tree.nodes[node_k].branch_length       = coordinates[parent] - coordinates[node_k]; 
      tree.nodes[child_right].branch_length  = coordinates[parent] - coordinates[child_right];  

      //everyone is unchanged until node_k
      //remove sib and add child_right
      for(int k = k_start; k < k_end; k++){
        remaining_new[sorted_indices[k]] = remaining[sorted_indices[k]];
        for(std::vector<int>::iterator it = remaining_new[sorted_indices[k]].begin(); it != remaining_new[sorted_indices[k]].end(); it++){
          if(*it == sib){
            *it = child_right;
            break;
          }
        } 
      }

      //update coal_rate_pair;
      //coal_rate_pair[e][node_k][*]
      //and coal_rate_pair[e][*][node_k]
      for(int e = 0; e < epoch.size(); e++){
        for(int i = 0; i < tree.nodes.size(); i++){
          if(node_k != i){
            coal_rate_pair[e][node_k][i] = (desc[child_left].num_leaves * coal_rate_pair[e][child_left][i] + desc[sib].num_leaves * coal_rate_pair[e][sib][i])/(desc[child_left].num_leaves + desc[sib].num_leaves);
            coal_rate_pair[e][node_k][i] = coal_rate_pair[e][i][node_k];
          }
        }
      }

    }else{

      //merge sib with child_right
      //then child_left
      tree.nodes[child_right].parent = &tree.nodes[node_k];
      tree.nodes[sib].parent = &tree.nodes[node_k];
      tree.nodes[node_k].child_right = &tree.nodes[child_right];
      tree.nodes[node_k].child_left = &tree.nodes[sib];

      tree.nodes[node_k].parent = &tree.nodes[parent];
      tree.nodes[child_left].parent = &tree.nodes[parent];
      tree.nodes[parent].child_right = &tree.nodes[node_k];
      tree.nodes[parent].child_left = &tree.nodes[child_left];

      tree.nodes[child_right].branch_length   = coordinates[node_k] - coordinates[child_right]; 
      tree.nodes[sib].branch_length          = coordinates[node_k] - coordinates[sib]; 
      tree.nodes[node_k].branch_length       = coordinates[parent] - coordinates[node_k]; 
      tree.nodes[child_left].branch_length  = coordinates[parent] - coordinates[child_left]; 

      //order, sorted_indices are unchanged
      //so just update remaining

      //everyone is unchanged until node_k
      //remove sib and add child_left
      for(int k = k_start; k < k_end; k++){
        remaining_new[sorted_indices[k]] = remaining[sorted_indices[k]];
        for(std::vector<int>::iterator it = remaining_new[sorted_indices[k]].begin(); it != remaining_new[sorted_indices[k]].end(); it++){
          if(*it == sib){
            *it = child_left;
            break;
          }
        }
      }

      //update coal_rate_pair;
      //coal_rate_pair[e][node_k][*]
      //and coal_rate_pair[e][*][node_k]
      for(int e = 0; e < epoch.size(); e++){
        for(int i = 0; i < tree.nodes.size(); i++){
          if(node_k != i){
            coal_rate_pair[e][node_k][i] = (desc[child_right].num_leaves * coal_rate_pair[e][child_right][i] + desc[sib].num_leaves * coal_rate_pair[e][sib][i])/(desc[child_right].num_leaves + desc[sib].num_leaves);
            coal_rate_pair[e][node_k][i] = coal_rate_pair[e][i][node_k];
          }
        }
      }

    }

    //calc Prior ratio
    log_likelihood_ratio += CalculatePrior(k_start, k_end, tree, epoch, coal_rate_pair, remaining_new, coordinates, sorted_indices, num_lineages);

    //then likelihood ratio
    delta_tau = tree.nodes[child_left].branch_length - bl_child_left;
    log_likelihood_ratio -= delta_tau*mut_rate[child_left];
    if(tree.nodes[child_left].num_events >= 0.0) log_likelihood_ratio += tree.nodes[child_left].num_events * log_deltat(delta_tau/bl_child_left);

    delta_tau = tree.nodes[child_right].branch_length - bl_child_right;
    log_likelihood_ratio -= delta_tau*mut_rate[child_right];
    if(tree.nodes[child_right].num_events >= 0.0) log_likelihood_ratio += tree.nodes[child_right].num_events * log_deltat(delta_tau/bl_child_right);

    delta_tau = tree.nodes[sib].branch_length - bl_sib;
    log_likelihood_ratio -= delta_tau*mut_rate[sib];
    if(tree.nodes[sib].num_events >= 0.0) log_likelihood_ratio += tree.nodes[sib].num_events * log_deltat(delta_tau/bl_sib);

    delta_tau = tree.nodes[node_k].branch_length - bl_node_k;
    log_likelihood_ratio -= delta_tau*mut_rate[node_k];
    if(tree.nodes[node_k].num_events >= 0.0) log_likelihood_ratio += tree.nodes[node_k].num_events * log_deltat(delta_tau/bl_node_k);

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
      for(int k_tmp = k_start; k_tmp < k_end; k_tmp++){
        remaining[sorted_indices[k_tmp]]    = remaining_new[sorted_indices[k_tmp]];
        assert(remaining[sorted_indices[k_tmp]].size() == num_lineages[sorted_indices[k_tmp]]);
      }
    }else{
      //reverse changes

      tree.nodes[child_left].parent = &tree.nodes[node_k];
      tree.nodes[child_right].parent = &tree.nodes[node_k];
      tree.nodes[node_k].child_left = &tree.nodes[child_left];
      tree.nodes[node_k].child_right = &tree.nodes[child_right];

      tree.nodes[sib].parent = &tree.nodes[parent];
      tree.nodes[parent].child_left = &tree.nodes[node_k];
      tree.nodes[parent].child_right = &tree.nodes[sib];

      tree.nodes[child_left].branch_length   = coordinates[node_k] - coordinates[child_left]; 
      tree.nodes[child_right].branch_length  = coordinates[node_k] - coordinates[child_right]; 
      tree.nodes[sib].branch_length          = coordinates[parent] - coordinates[sib]; 
      tree.nodes[node_k].branch_length       = coordinates[parent] - coordinates[node_k]; 

      //update coal_rate_pair;
      //coal_rate_pair[e][node_k][*]
      //and coal_rate_pair[e][*][node_k]
      for(int e = 0; e < epoch.size(); e++){
        for(int i = 0; i < tree.nodes.size(); i++){
          if(node_k != i){
            coal_rate_pair[e][node_k][i] = (desc[child_right].num_leaves * coal_rate_pair[e][child_right][i] + desc[child_left].num_leaves * coal_rate_pair[e][child_left][i])/(desc[child_right].num_leaves + desc[child_left].num_leaves);
            coal_rate_pair[e][node_k][i] = coal_rate_pair[e][i][node_k];
          }
        }
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

			//std::cerr << ep << " " << epoch[ep] << " " << p_coordinates[p_sorted_indices[k_tmp-1]] << std::endl;
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

double
EstimateBranchLengthsWithSampleAge::CalculatePrior(const Tree& tree, const std::vector<double>& epoch, const std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages){

  double log_likelihood = 0.0;
  int ep     = 0; 
  int k_tmp  = 0, k_prev = 0;
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
  double coal = 0.0;

  while(k_tmp < 2*N-2){

    k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
    assert(num_lineages_tmp >= 1);
    k_prev = k_tmp;

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

      //std::cerr << k_prev << " " << ep << " " << p_coordinates[p_sorted_indices[k_prev]] << " " << epoch[ep] << std::endl;
      assert(p_coordinates[p_sorted_indices[k_prev]] >= epoch[ep]);
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
        for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
          for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
            if(*itr1 != *itr2){
              coal += coal_rate_pair[ep][*itr1][*itr2];
            }
          }
        }
        if(coal > 0){
          log_likelihood -= coal * delta_tmp_tau; 
        }
        tmp_tau                -= delta_tmp_tau;

        ep++;
        delta_tmp_tau           = epoch[ep+1] - epoch[ep];
        while(tmp_tau > delta_tmp_tau && ep < epoch.size() - 1){
          coal = 0.0;
          for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
            for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
              if(*itr1 != *itr2){
                coal += coal_rate_pair[ep][*itr1][*itr2];
              }
            }
          }
          if(coal > 0){;
            log_likelihood -= coal * delta_tmp_tau; 
          }
          tmp_tau              -= delta_tmp_tau;

          ep++;
          delta_tmp_tau         = epoch[ep+1] - epoch[ep];
        }
        assert(tmp_tau >= 0.0);

        coal = 0.0;
        for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
          for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
            if(*itr1 != *itr2){
              coal += coal_rate_pair[ep][*itr1][*itr2];
            }
          }
        }

        if(coal == 0){
          log_likelihood  = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood -= coal * tmp_tau;
          if(!is_sample){
            log_likelihood += log(coal_rate_pair[ep][(*tree.nodes[p_sorted_indices[k_tmp]].child_left).label][(*tree.nodes[p_sorted_indices[k_tmp]].child_right).label]);
          }
        }

      }else{

        coal = 0.0;
        for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
          for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
            if(*itr1 != *itr2){
              coal += coal_rate_pair[ep][*itr1][*itr2];
            }
          }
        }

        if(coal == 0){
          log_likelihood  = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood -= coal * tmp_tau;
          if(!is_sample){
            log_likelihood += log(coal_rate_pair[ep][(*tree.nodes[p_sorted_indices[k_tmp]].child_left).label][(*tree.nodes[p_sorted_indices[k_tmp]].child_right).label]);
          }
        }

      }

    }else{

      coal = 0.0;
      for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
        for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
          if(*itr1 != *itr2){
            coal += coal_rate_pair[ep][*itr1][*itr2];
          }
        }
      }

      if(coal == 0){
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

        log_likelihood -= coal * tmp_tau;
        if(!is_sample){
          log_likelihood += log(coal_rate_pair[ep][(*tree.nodes[p_sorted_indices[k_tmp]].child_left).label][(*tree.nodes[p_sorted_indices[k_tmp]].child_right).label]);
        }

      }
    }
  }

  return(log_likelihood);

}

double
EstimateBranchLengthsWithSampleAge::CalculatePrior(int k_start, int k_end, const Tree& tree, const std::vector<double>& epoch, const std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, std::vector<double>& p_coordinates, std::vector<int>& p_sorted_indices, std::vector<int>& p_num_lineages){

  double log_likelihood = 0.0;
  int k_tmp  = k_start, k_prev;
  int node   = p_sorted_indices[k_tmp];
  double age = p_coordinates[node];
  double coal = 0.0;

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

  if(k_start == order.size()-2){
    if(remaining[p_sorted_indices[k_start]].size() != 2){
      std::cerr << remaining[p_sorted_indices[k_start]].size() << " " << p_sorted_indices[k_start] << std::endl;
    }
    assert(num_lineages_tmp == 2);
    assert(remaining[p_sorted_indices[k_start]].size() == 2);
  }

  while(k_tmp < k_end){

    //if(k_start == order.size()-2) std::cerr << num_lineages_tmp << std::endl;

    k_choose_2_tmp   = (num_lineages_tmp * (num_lineages_tmp - 1.0))/2.0;
    if(num_lineages_tmp < 1) std::cerr << num_lineages_tmp << std::endl;
    k_prev = k_tmp;
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
        coal = 0.0;
        for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
          for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
            if(*itr1 != *itr2){
              coal += coal_rate_pair[ep][*itr1][*itr2];
            }
          }
        } 

        if(coal > 0){
          log_likelihood -= coal * delta_tmp_tau; 
        }

        tmp_tau                -= delta_tmp_tau;

        ep++;
        delta_tmp_tau           = epoch[ep+1] - epoch[ep];
        while(tmp_tau > delta_tmp_tau && ep < epoch.size() - 1){

          coal = 0.0;
          for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
            for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
              if(*itr1 != *itr2){
                coal += coal_rate_pair[ep][*itr1][*itr2];
              }
            }
          }
          if(coal > 0){
            log_likelihood -= coal * delta_tmp_tau; 
          }
          tmp_tau              -= delta_tmp_tau;

          ep++;
          delta_tmp_tau         = epoch[ep+1] - epoch[ep];
        }
        assert(tmp_tau >= 0.0);

        coal = 0.0;
        for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
          for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
            if(*itr1 != *itr2){
              coal += coal_rate_pair[ep][*itr1][*itr2];
            }
          }
        }

        if(coal == 0){
          log_likelihood  = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood -= coal * tmp_tau; 
          if(!is_sample) log_likelihood += log(coal_rate_pair[ep][(*tree.nodes[p_sorted_indices[k_tmp]].child_left).label][(*tree.nodes[p_sorted_indices[k_tmp]].child_right).label]);
        }

      }else{

        coal = 0.0;
        for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
          for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
            if(*itr1 != *itr2){
              coal += coal_rate_pair[ep][*itr1][*itr2];
            }
          }
        }

        if(coal == 0){
          log_likelihood  = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood -= coal * tmp_tau; 
          if(!is_sample) log_likelihood += log(coal_rate_pair[ep][(*tree.nodes[p_sorted_indices[k_tmp]].child_left).label][(*tree.nodes[p_sorted_indices[k_tmp]].child_right).label]);
        }

      }

    }else{

      coal = 0.0;
      for(std::vector<int>::iterator itr1 = remaining[p_sorted_indices[k_prev]].begin(); itr1 != remaining[p_sorted_indices[k_prev]].end(); itr1++){
        for(std::vector<int>::iterator itr2 = remaining[p_sorted_indices[k_prev]].begin(); itr2 != itr1; itr2++){
          if(*itr1 != *itr2){
            coal += coal_rate_pair[ep][*itr1][*itr2];
          }
        }
      }

      if(coal == 0){
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

        log_likelihood -= coal * tmp_tau;
        if(!is_sample) log_likelihood += log(coal_rate_pair[ep][(*tree.nodes[p_sorted_indices[k_tmp]].child_left).label][(*tree.nodes[p_sorted_indices[k_tmp]].child_right).label]);

      }
    }
  }

  return(log_likelihood);

}




////////////////////// constant NE /////////////////////
//propose a new time and change events only up to t+bound

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

    if(1){
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

void
EstimateBranchLengthsWithSampleAge::UpdateOneEventVP(Tree& tree, int node_k, const std::vector<double>& epoch, const std::vector<std::vector<std::vector<float>>>& coal_rate_pair, std::vector<std::vector<int>>& remaining, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif){

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

    double log_likelihood = CalculatePrior(k_start, k_end, tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);
    //double log_likelihood = CalculatePrior(tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);
    coordinates[node_k] -= delta_tau;
    if(log_likelihood != -std::numeric_limits<float>::infinity()){
      //log_likelihood_ratio += log_likelihood;
      log_likelihood -= CalculatePrior(k_start, k_end, tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);
      //log_likelihood -= CalculatePrior(tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);
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

    if(1){
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
    }

    //decide whether to accept proposal or not.
    accept = true;
    if(log_likelihood_ratio < 0.0){
      //accept with probability exp(log_p)
      if(dist_unif(rng) > exp(log_likelihood_ratio)){
        accept = false;
      }
    }

    //accept = false;

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
            remaining_new[sorted_indices[k_start-1]]    = remaining[sorted_indices[k_start-1]];
            sorted_indices_new[k_start]                 = sorted_indices[k_start];
            num_lineages_new[sorted_indices[k_start]]   = num_lineages[sorted_indices[k_start]];
            remaining_new[sorted_indices[k_start]]      = remaining[sorted_indices[k_start]];

            for(int k_tmp = k_start; k_tmp < k_end; k_tmp++){
              node_tmp = sorted_indices[k_tmp + 1];
              if(age > coordinates[node_tmp]){
                sorted_indices_new[k_tmp]  = node_tmp;
                order_new[node_tmp]        = k_tmp;
                num_lineages_new[node_tmp] = num_lineages[node_tmp] + 1;

                remaining_new[node_tmp]    = remaining[node_tmp];
                for(std::vector<int>::iterator it = remaining_new[node_tmp].begin(); it != remaining_new[node_tmp].end(); it++){
                  if(*it == node_k) *it = child_left_label;
                }
                remaining_new[node_tmp].push_back(child_right_label);
              }else{
                sorted_indices_new[k_tmp] = node_k;
                order_new[node_k]         = k_tmp;
                num_lineages_new[node_k]  = num_lineages_new[sorted_indices_new[k_tmp-1]]-1; 
                remaining_new[node_k]     = remaining_new[sorted_indices_new[k_tmp-1]];
                for(std::vector<int>::iterator it = remaining_new[node_k].begin(); it != remaining_new[node_k].end(); it++){
                  if(*it == child_left_label){
                    *it = node_k;
                    break;
                  }
                }
                for(std::vector<int>::iterator it = remaining_new[node_k].begin(); it != remaining_new[node_k].end(); it++){
                  if(*it == child_right_label){
                    *it = remaining_new[node_k][remaining_new[node_k].size()-1];
                    break;
                  }
                }
                remaining_new[node_k].pop_back();

                k_start = k - 1;
                k_end   = k_tmp + 1;
                num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];
                sorted_indices_new[k_start]               = sorted_indices[k_start];
                order_new[sorted_indices[k_start]]        = k_start;
                remaining_new[sorted_indices[k_start]]    = remaining[sorted_indices[k_start]];
                num_lineages_new[sorted_indices[k_end]] = num_lineages[sorted_indices[k_end]];
                sorted_indices_new[k_end]               = sorted_indices[k_end];
                order_new[sorted_indices[k_end]]        = k_end;
                remaining_new[sorted_indices[k_end]]    = remaining[sorted_indices[k_end]];
                break;
              }
            }

          }else{      
            k_end = k;
            k_start = order[child_left_label];
            if(k_start < order[child_right_label]) k_start = order[child_right_label]; 

            sorted_indices_new[k_start-1]               = sorted_indices[k_start-1];
            num_lineages_new[sorted_indices[k_start-1]] = num_lineages[sorted_indices[k_start-1]];
            remaining_new[sorted_indices[k_start-1]]    = remaining[sorted_indices[k_start-1]];
            sorted_indices_new[k_start]               = sorted_indices[k_start];
            num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];
            remaining_new[sorted_indices[k_start]]    = remaining[sorted_indices[k_start]];

            int node_tmp;
            double age = coords_new;
            for(int k_tmp = k_end; k_tmp > k_start; k_tmp--){
              node_tmp = sorted_indices[k_tmp - 1];
              if(age < coordinates[node_tmp]){
                sorted_indices_new[k_tmp]  = node_tmp;
                order_new[node_tmp]        = k_tmp;
                num_lineages_new[node_tmp] = num_lineages[node_tmp] - 1;

                remaining_new[node_tmp]    = remaining[node_tmp];
                for(std::vector<int>::iterator it = remaining_new[node_tmp].begin(); it != remaining_new[node_tmp].end(); it++){
                  if(*it == child_left_label){
                    *it = node_k;
                    break;
                  }
                }
                for(std::vector<int>::iterator it = remaining_new[node_tmp].begin(); it != remaining_new[node_tmp].end(); it++){
                  if(*it == child_right_label){
                    *it = remaining_new[node_tmp][remaining_new[node_tmp].size()-1];
                    break;
                  }
                }
                remaining_new[node_tmp].pop_back();
              }else{
                sorted_indices_new[k_tmp] = node_k;
                order_new[node_k]         = k_tmp;
                num_lineages_new[node_k]  = num_lineages[sorted_indices[k_tmp-1]]-1;
                remaining_new[node_k]     = remaining[sorted_indices[k_tmp-1]];
                for(std::vector<int>::iterator it = remaining_new[node_k].begin(); it != remaining_new[node_k].end(); it++){
                  if(*it == child_left_label){
                    *it = node_k;
                    break;
                  }
                }
                for(std::vector<int>::iterator it = remaining_new[node_k].begin(); it != remaining_new[node_k].end(); it++){
                  if(*it == child_right_label){
                    *it = remaining_new[node_k][remaining_new[node_k].size()-1];
                    break;
                  }
                }
                remaining_new[node_k].pop_back();

                k_start = k_tmp - 1;
                k_end   = k + 1;
                sorted_indices_new[k_start]               = sorted_indices[k_start];
                order_new[sorted_indices[k_start]]        = k_start;
                num_lineages_new[sorted_indices[k_start]] = num_lineages[sorted_indices[k_start]];
                remaining_new[sorted_indices[k_start]]    = remaining[sorted_indices[k_start]];

                sorted_indices_new[k_end]               = sorted_indices[k_end];
                order_new[sorted_indices[k_end]]        = k_end;
                num_lineages_new[sorted_indices[k_end]] = num_lineages[sorted_indices[k_end]];
                remaining_new[sorted_indices[k_end]]    = remaining[sorted_indices[k_end]];

                break;
              }
            }

          }

          for(int k_tmp = k_start; k_tmp < k_end; k_tmp++){
            //std::cerr << k_tmp << " " << k_start << " " << k_end << " " << remaining.size() << " " << sorted_indices[k_tmp] << " " << remaining[sorted_indices[k_tmp]].size() << " " << remaining_new[sorted_indices_new[k_tmp]].size() << " " << num_lineages[sorted_indices[k_tmp]] << " " << num_lineages_new[sorted_indices_new[k_tmp]] << std::endl;
            assert(remaining[sorted_indices[k_tmp]].size() ==  num_lineages[sorted_indices[k_tmp]]);
            assert(remaining_new[sorted_indices_new[k_tmp]].size() ==  num_lineages_new[sorted_indices_new[k_tmp]]);
          }

          coordinates[node_k] = coords_new;

          log_likelihood = CalculatePrior(k_start, k_end, tree, epoch, coal_rate_pair, remaining_new, coordinates, sorted_indices_new, num_lineages_new);
          //log_likelihood = CalculatePrior(tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices_new, num_lineages_new);
          coordinates[node_k] = coords;
          if(log_likelihood != -std::numeric_limits<float>::infinity()){
            //log_likelihood_ratio += log_likelihood;
            log_likelihood -= CalculatePrior(k_start, k_end, tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);
            //log_likelihood = CalculatePrior(tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);
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

          //Here doing from scratch, more efficient if I set active = remaining[sorted_indices[k_start]] and update from there
          std::vector<int> active;
          ages = sample_age[*sorted_indices.begin()];
          it_sorted_indices_start = sorted_indices.begin();
          for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
            if(*it_sorted_indices >= N){
              for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
                remaining_new[*it_sorted_indices_start] = active;
              }
              //remove two children and replace by *it_sorted_indices
              int ind1, ind2;
              int c = 0;
              for(std::vector<int>::iterator it_act = active.begin(); it_act != active.end(); it_act++){
                if(*it_act == (*tree.nodes[*it_sorted_indices].child_left).label) ind1 = c;
                if(*it_act == (*tree.nodes[*it_sorted_indices].child_right).label) ind1 = c;
                c++;
              }
              active[ind1] = *it_sorted_indices;
              active[ind2] = active[active.size()-1];
              active.pop_back();
              remaining_new[*it_sorted_indices] = active;
              it_sorted_indices_start++;
            }else if(ages < sample_age[*it_sorted_indices]){      
              for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
                remaining_new[*it_sorted_indices_start] = active;
              }
              ages = sample_age[*it_sorted_indices];
              active.push_back(*it_sorted_indices);
            }else{
              active.push_back(*it_sorted_indices);
            }
          }

          log_likelihood = CalculatePrior(tree, epoch, coal_rate_pair, remaining_new, coordinates, sorted_indices_new, num_lineages_new);
          coordinates[node_k] = coords;
          if(log_likelihood != -std::numeric_limits<float>::infinity()){
            //log_likelihood_ratio += log_likelihood;
            log_likelihood -= CalculatePrior(tree, epoch, coal_rate_pair, remaining, coordinates, sorted_indices, num_lineages);
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

        if(1){
          for(int k_tmp = k_start; k_tmp < k_end; k_tmp++){
            sorted_indices[k_tmp]               = sorted_indices_new[k_tmp];
            order[sorted_indices[k_tmp]]        = order_new[sorted_indices[k_tmp]];
            num_lineages[sorted_indices[k_tmp]] = num_lineages_new[sorted_indices[k_tmp]];
            remaining[sorted_indices[k_tmp]]    = remaining_new[sorted_indices[k_tmp]];
            assert(remaining[sorted_indices[k_tmp]].size() == num_lineages[sorted_indices[k_tmp]]);
          }
        }else{
          sorted_indices = sorted_indices_new;
          order = order_new;
          num_lineages = num_lineages_new;
          remaining = remaining_new;
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

  bool debug = false;

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
  p2 = 0.7;

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  InitializeMCMC(data, tree); //Initialize using coalescent prior

  std::vector<double> sample_age_tmp;

  bool is_ancient = false;
  for(int i = 0; i < sample_age.size(); i++){
    if(sample_age[i] > 0){
      is_ancient = true;
      break;
    }
  }

  if(is_ancient){
    sample_age_tmp = sample_age;
    std::fill(sample_age.begin(), sample_age.end(), 0.0);
  }

  InitializeOrder(tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) 2*data.N * data.N; j++){
    RandomSwitchOrder(tree, dist_n(rng), dist_unif);
  }
  //Initialise branch lengths
  InitializeBranchLengths(tree);

  if(debug){
    for(int i = 0; i < N_total-1; i++){
      assert(tree.nodes[i].branch_length >= 0.0);
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
      assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
    }
  }

  if(is_ancient){
    count = 0;
    for(; count < 50*delta; count++){
      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        UpdateOneEvent(tree, dist_oneevent(rng), dist_gamma, dist_unif);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
      }      

      for(int k = 0; k < N_total; k++){
        assert(tree.nodes[k].branch_length >= 0.0);
      }
    }

    GetCoordinates(tree.nodes[root], coordinates);

    sample_age = sample_age_tmp;
    double min_sample_age = sample_age[0];
    for(int i = 0; i < N; i++){
      if(min_sample_age > sample_age[i]){
        min_sample_age = sample_age[i];
      }
    }
    if(min_sample_age > 0){
      for(std::vector<double>::iterator it_coords = coordinates.begin(); it_coords != coordinates.end(); it_coords++){
        *it_coords += min_sample_age;
      }
    }

    for(int i = 0; i < N; i++){
      if(sample_age[i] > 0){

        int n = (*tree.nodes[i].parent).label;
        if(coordinates[n] > sample_age[i]){
          //tree.nodes[n].branch_length -= sample_age[i];
          coordinates[i] = sample_age[i];
        }else{
          coordinates[i] = sample_age[i];
          float prev_coords = coordinates[i];
          coordinates[n] += sample_age[i];
          prev_coords = coordinates[n];
          while(tree.nodes[n].parent != NULL){
            n = (*tree.nodes[n].parent).label;
            if(coordinates[n] <= prev_coords){
              coordinates[n] += sample_age[i];
              prev_coords = coordinates[n];
            }else{
              break;
            }
          }
        }

      }
    }

    for(int i = 0; i < N_total-1; i++){
      tree.nodes[i].branch_length = coordinates[(*tree.nodes[i].parent).label] - coordinates[i];
    }
    order.clear();
    sorted_indices.clear();
    order.resize(N_total);
    sorted_indices.resize(N_total);

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

  }

  sorted_indices_new = sorted_indices;
  order_new          = order;
  num_lineages_new   = num_lineages;

  for(int i = 0; i < N_total-1; i++){
    assert(tree.nodes[i].branch_length >= 0.0);
    assert(order[sorted_indices[i]] == i);
    assert(order[i] < order[(*tree.nodes[i].parent).label]);
    assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
  }


  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 50*delta; count++){

    for(int i = 0; i < N_total-1; i++){
      assert(tree.nodes[i].branch_length >= 0.0);
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
      assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
    }

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p2){
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
      if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        int k_candidate = dist_oneevent(rng);
        count_proposals[k_candidate-N]++;
        UpdateOneEvent(tree, k_candidate, dist_gamma, dist_unif);
        UpdateAvg(tree);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
        UpdateAvg(tree);
      }

      if(debug){
        for(int k = 0; k < N_total; k++){
          assert(tree.nodes[k].branch_length >= 0.0);
        }
        for(int k = N; k < N_total; k++){
          assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_left).label]));
          assert((coordinates[k] >= coordinates[(*tree.nodes[k].child_right).label]));
          assert(sorted_indices[order[k]] == k);
        }
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
          if(*it_count < 50){
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

//TODO: undo some changes
void
EstimateBranchLengthsWithSampleAge::MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  bool debug = false;

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 0.7;

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  InitializeMCMC(data, tree); //Initialize using coalescent prior 

  double total_bl = 0.0;
  for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != tree.nodes.end(); it_n++){
    total_bl += (*it_n).branch_length;
  }

  if(total_bl != 0){
    //as a sanity check that branch lengths are consistent with sample_ages, calculate branch lengths to root + sample_age
    double dist_to_root, dist_to_root_0 = 0.0;
    dist_to_root = sample_age[0];
    int n = 0;
    while(tree.nodes[n].parent != NULL){
      dist_to_root_0 += tree.nodes[n].branch_length/Ne;
      n = (*tree.nodes[n].parent).label;
    }
    for(int i = 1; i < N; i++){
      dist_to_root = sample_age[i];
      n = i;
      while(tree.nodes[n].parent != NULL){
        dist_to_root += tree.nodes[n].branch_length/Ne;
        n = (*tree.nodes[n].parent).label;
      }
      if(std::fabs(dist_to_root - dist_to_root_0) > 1e-1){
        total_bl = 0;
        break;
      }
    }
  }

  if(total_bl == 0){

    std::vector<double> sample_age_tmp;

    bool is_ancient = false;
    for(int i = 0; i < sample_age.size(); i++){
      if(sample_age[i] > 0){
        is_ancient = true;
        break;
      }
    }

    if(is_ancient){
      sample_age_tmp = sample_age;
      std::fill(sample_age.begin(), sample_age.end(), 0.0);
    }

    InitializeOrder(tree); 

    //Randomly switch around order of coalescences
    for(int j = 0; j < (int) 2*data.N * data.N; j++){
      RandomSwitchOrder(tree, dist_n(rng), dist_unif);
    }
    //Initialise branch lengths
    InitializeBranchLengths(tree);

    if(debug){
      for(int i = 0; i < N_total-1; i++){
        assert(tree.nodes[i].branch_length >= 0.0);
        assert(order[sorted_indices[i]] == i);
        assert(order[i] < order[(*tree.nodes[i].parent).label]);
        assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
      }
    }

    if(is_ancient){
      count = 0;
      for(; count < 50*delta; count++){
        //Either switch order of coalescent event or extent time while k ancestors 
        uniform_rng = dist_unif(rng);
        if(uniform_rng <= p2){
          //std::cerr << "v3" << std::endl;
          UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
        }else{ 
          //std::cerr << "v4" << std::endl;
          SwitchOrder(tree, dist_n(rng), dist_unif);
        }    

        if(debug){
          for(int k = 0; k < N_total; k++){
            assert(tree.nodes[k].branch_length >= 0.0);
          }
        }
      }

      GetCoordinates(tree.nodes[root], coordinates);

      sample_age = sample_age_tmp;
      double min_sample_age = sample_age[0];
      for(int i = 0; i < N; i++){
        if(min_sample_age > sample_age[i]){
          min_sample_age = sample_age[i];
        }
      }
      if(min_sample_age > 0){
        for(std::vector<double>::iterator it_coords = coordinates.begin(); it_coords != coordinates.end(); it_coords++){
          *it_coords += min_sample_age;
        }
      }

      for(int i = 0; i < N; i++){
        if(sample_age[i] > 0){

          int n = (*tree.nodes[i].parent).label;
          if(coordinates[n] > sample_age[i]){
            //tree.nodes[n].branch_length -= sample_age[i];
            coordinates[i] = sample_age[i];
          }else{
            coordinates[i] = sample_age[i];
            float prev_coords = coordinates[i];
            coordinates[n] += sample_age[i];
            prev_coords = coordinates[n];
            while(tree.nodes[n].parent != NULL){
              n = (*tree.nodes[n].parent).label;
              if(coordinates[n] <= prev_coords){
                coordinates[n] += sample_age[i];
                prev_coords = coordinates[n];
              }else{
                break;
              }
            }
          }

        }
      }

      for(int i = 0; i < N_total-1; i++){
        tree.nodes[i].branch_length = coordinates[(*tree.nodes[i].parent).label] - coordinates[i];
      }
      order.clear();
      sorted_indices.clear();
      order.resize(N_total);
      sorted_indices.resize(N_total);

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

      if(debug){
        for(int i = 0; i < N_total-1; i++){
          assert(tree.nodes[i].branch_length >= 0.0);
          assert(order[sorted_indices[i]] == i);
          assert(order[i] < order[(*tree.nodes[i].parent).label]);
          assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
          assert(num_lineages[sorted_indices[i]] >= 1);
        }
      }

    }

  }else{

    for(std::vector<Node>::iterator it_node = tree.nodes.begin(); it_node != tree.nodes.end(); it_node++){
      (*it_node).branch_length /= data.Ne;
    }

    coordinates.resize(N_total);
    GetCoordinates(tree.nodes[root], coordinates);
    //avg   = coordinates;

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
    if(debug){
      for(int i = 0; i < N_total-1; i++){
        assert(order[sorted_indices[i]] == i);
        assert(order[i] < order[(*tree.nodes[i].parent).label]);
        assert(num_lineages[sorted_indices[i]] >= 1);
      }
    }

  }

  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 50*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p2){
      //std::cerr << "v3" << std::endl;
      UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
      //UpdateOneEvent(tree, dist_oneevent(rng), dist_gamma, dist_unif);
    }else{ 
      //std::cerr << "v4" << std::endl;
      SwitchOrder(tree, dist_n(rng), dist_unif);
    }    

    if(debug){
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
      if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        int k_candidate = dist_oneevent(rng);
        count_proposals[k_candidate-N]++;
        UpdateOneEventVP(tree, k_candidate, epoch, coal_rate, dist_gamma, dist_unif);
        UpdateAvg(tree);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
        UpdateAvg(tree);
      }

      if(debug){
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
          if(*it_count < 50){
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

  if(1){
    for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
      (*it_n).branch_length = ((double) Ne) * (avg[(*(*it_n).parent).label] - avg[(*it_n).label]);
    }
  }else{
    for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
      (*it_n).branch_length *= ((double) Ne);
    }
  }

}  

void
EstimateBranchLengthsWithSampleAge::MCMCVariablePopulationSizeForRelate(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  bool debug = false;

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 0.7;

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  ///////////////////////////////////////////////////////////////////////////////////

  InitializeMCMC(data, tree); //Initialize using coalescent prior 

  std::vector<double> sample_age_tmp;

  bool is_ancient = false;
  for(int i = 0; i < sample_age.size(); i++){
    if(sample_age[i] > 0){
      is_ancient = true;
      break;
    }
  }

  if(is_ancient){
    sample_age_tmp = sample_age;
    std::fill(sample_age.begin(), sample_age.end(), 0.0);
  }
  
  InitializeOrder(tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) 2*data.N * data.N; j++){
    RandomSwitchOrder(tree, dist_n(rng), dist_unif);
  }
  //Initialise branch lengths
  InitializeBranchLengths(tree);

  if(debug){
    for(int i = 0; i < N_total-1; i++){
      assert(tree.nodes[i].branch_length >= 0.0);
      assert(order[sorted_indices[i]] == i);
      assert(order[i] < order[(*tree.nodes[i].parent).label]);
      assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
    }
  }

  if(is_ancient){

    count = 0;
    for(; count < 50*delta; count++){
      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
      }    

      if(debug){
        for(int k = 0; k < N_total; k++){
          assert(tree.nodes[k].branch_length >= 0.0);
        }
      }
    }

    if(1){
      GetCoordinates(tree.nodes[root], coordinates);

      sample_age = sample_age_tmp;
      double min_sample_age = sample_age[0];
      for(int i = 0; i < N; i++){
        if(min_sample_age > sample_age[i]){
          min_sample_age = sample_age[i];
        }
      }
      if(min_sample_age > 0){
        for(std::vector<double>::iterator it_coords = coordinates.begin(); it_coords != coordinates.end(); it_coords++){
          *it_coords += min_sample_age;
        }
      }

      for(int i = 0; i < N; i++){
        if(sample_age[i] > 0){

          int n = (*tree.nodes[i].parent).label;
          if(coordinates[n] > sample_age[i]){
            //tree.nodes[n].branch_length -= sample_age[i];
            coordinates[i] = sample_age[i];
          }else{
            coordinates[i] = sample_age[i];
            float prev_coords = coordinates[i];
            coordinates[n] += sample_age[i];
            prev_coords = coordinates[n];
            while(tree.nodes[n].parent != NULL){
              n = (*tree.nodes[n].parent).label;
              if(coordinates[n] <= prev_coords){
                coordinates[n] += sample_age[i];
                prev_coords = coordinates[n];
              }else{
                break;
              }
            }
          }

        }
      }
    }

    for(int i = 0; i < N_total-1; i++){
      tree.nodes[i].branch_length = coordinates[(*tree.nodes[i].parent).label] - coordinates[i];
    }
    order.clear();
    sorted_indices.clear();
    order.resize(N_total);
    sorted_indices.resize(N_total);

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

  }
  /////////////////////////////////////////////////////////

  if(1){

    ////////////////// Transient /////////////////

    count = 0;
    for(; count < 50*delta; count++){

      if(debug){
        for(int i = 0; i < N_total-1; i++){
          assert(tree.nodes[i].branch_length >= 0.0);
          assert(order[sorted_indices[i]] == i);
          assert(order[i] < order[(*tree.nodes[i].parent).label]);
          assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
        }
      }

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
      }    

      if(debug){
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
        if(uniform_rng <= p2){
          //std::cerr << "v3" << std::endl;
          int k_candidate = dist_oneevent(rng);
          count_proposals[k_candidate-N]++;
          UpdateOneEventVP(tree, k_candidate, epoch, coal_rate, dist_gamma, dist_unif);
          UpdateAvg(tree);
        }else{ 
          //std::cerr << "v4" << std::endl;
          SwitchOrder(tree, dist_n(rng), dist_unif);
          UpdateAvg(tree);
        }

        if(debug){
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
            if(*it_count < 50){
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

  }else{
    GetCoordinates(tree.nodes[root], coordinates);
    avg              = coordinates;
  }

  //UNDO
  if(1){
    for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
      (*it_n).branch_length = ((double) Ne) * (avg[(*(*it_n).parent).label] - avg[(*it_n).label]);
    }
  }else{
    for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
      (*it_n).branch_length *= ((double) Ne);
    }
  }

}  

void
EstimateBranchLengthsWithSampleAge::MCMCCoalRatesForRelate(const Data& data, Tree& tree, const std::vector<int>& membership, const std::vector<double>& epoch, std::vector<std::vector<std::vector<double>>>& coal_rate, const int seed){

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

  ///////////////////////////////////////////////////////////////////////////////////

  InitializeMCMC(data, tree); //Initialize using coalescent prior 

  std::vector<double> sample_age_tmp;

  bool is_ancient = false;
  for(int i = 0; i < sample_age.size(); i++){
    if(sample_age[i] > 0){
      is_ancient = true;
      break;
    }
  }

  if(is_ancient){
    sample_age_tmp = sample_age;
    std::fill(sample_age.begin(), sample_age.end(), 0.0);
  }

  InitializeOrder(tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) 2*data.N * data.N; j++){
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

  if(is_ancient){

    count = 0;
    for(; count < 10*delta; count++){
      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng <= p2){
        //std::cerr << "v3" << std::endl;
        UpdateOneEvent(tree, dist_oneevent(rng), dist_gamma, dist_unif);
      }else{ 
        //std::cerr << "v4" << std::endl;
        SwitchOrder(tree, dist_n(rng), dist_unif);
      }    

      for(int k = 0; k < N_total; k++){
        assert(tree.nodes[k].branch_length >= 0.0);
      }
    }

    GetCoordinates(tree.nodes[root], coordinates);

    sample_age = sample_age_tmp;
    double min_sample_age = sample_age[0];
    for(int i = 0; i < N; i++){
      if(min_sample_age > sample_age[i]){
        min_sample_age = sample_age[i];
      }
    }
    if(min_sample_age > 0){
      for(std::vector<double>::iterator it_coords = coordinates.begin(); it_coords != coordinates.end(); it_coords++){
        *it_coords += min_sample_age;
      }
    }


    for(int i = 0; i < N; i++){
      if(sample_age[i] > 0){

        int n = (*tree.nodes[i].parent).label;
        if(coordinates[n] > sample_age[i]){
          //tree.nodes[n].branch_length -= sample_age[i];
          coordinates[i] = sample_age[i];
        }else{
          coordinates[i] = sample_age[i];
          float prev_coords = coordinates[i];
          coordinates[n] += sample_age[i];
          prev_coords = coordinates[n];
          while(tree.nodes[n].parent != NULL){
            n = (*tree.nodes[n].parent).label;
            if(coordinates[n] <= prev_coords){
              coordinates[n] += sample_age[i];
              prev_coords = coordinates[n];
            }else{
              break;
            }
          }
        }

      }
    }


    for(int i = 0; i < N_total-1; i++){
      tree.nodes[i].branch_length = coordinates[(*tree.nodes[i].parent).label] - coordinates[i];
    }
    order.clear();
    sorted_indices.clear();
    order.resize(N_total);
    sorted_indices.resize(N_total);

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

  }
  /////////////////////////////////////////////////////////

  std::vector<Leaves> desc;
  tree.FindAllLeaves(desc);

  //I need which lineages are remaining
  //std::vector<std::vector<float>> num_members(tree.nodes.size());
  //std::vector<float> active_num(sample.groups.size(), 0);
  //for(int i = 0; i < tree.nodes.size(); i++){
  //  num_members.resize(sample.groups.size());
  //}

  remaining.resize(tree.nodes.size());
  remaining_new.resize(tree.nodes.size());
  std::vector<int> active;
  double ages = sample_age[*sorted_indices.begin()];
  std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
  int k = 0;
  for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    if(*it_sorted_indices >= N){
      for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        remaining[*it_sorted_indices_start] = active;
      }
      //remove two children and replace by *it_sorted_indices
      int ind1, ind2;
      int c = 0;
      for(std::vector<int>::iterator it_act = active.begin(); it_act != active.end(); it_act++){
        if(*it_act == (*tree.nodes[*it_sorted_indices].child_left).label) ind1 = c;
        if(*it_act == (*tree.nodes[*it_sorted_indices].child_right).label) ind2 = c;
        c++;
      }
      active[ind1] = *it_sorted_indices;
      active[ind2] = active[active.size()-1];
      active.pop_back();

      int c1 = (*tree.nodes[*it_sorted_indices].child_left).label;
      int c2 = (*tree.nodes[*it_sorted_indices].child_right).label;
      float frac = 1.0/desc[c1].num_leaves;
      //for(std::vector<int>::iterator it_mem = desc[c1].member.begin(); it_mem != desc[c1].member.end(); it_mem++){
      //  active_num[membership[*it_mem]] -= frac;
      //}
      frac = 1.0/desc[c2].num_leaves;
      //for(std::vector<int>::iterator it_mem = desc[c2].member.begin(); it_mem != desc[c2].member.end(); it_mem++){
      //  active_num[membership[*it_mem]] -= frac;
      //}
      frac = 1.0/desc[*it_sorted_indices].num_leaves;
      //for(std::vector<int>::iterator it_mem = desc[*it_sorted_indices].member.begin(); it_mem != desc[*it_sorted_indices].member.end(); it_mem++){
      //  active_num[membership[*it_mem]] += frac;
      //}
      //num_member[*it_sorted_indices] = active_num;

      remaining[*it_sorted_indices] = active;
      it_sorted_indices_start++;
    }else if(ages < sample_age[*it_sorted_indices]){      
      for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        remaining[*it_sorted_indices_start] = active;
        //num_member[*it_sorted_indices_start] = active_num;
      }
      ages = sample_age[*it_sorted_indices];
      active.push_back(*it_sorted_indices);
      //active_num[membership[*it_sorted_indices]] += 1.0;
    }else{
      active.push_back(*it_sorted_indices);
      //active_num[membership[*it_sorted_indices]] += 1.0;
    }
  }

  //I need to calculate E x #nodes for coalescence rate
  std::vector<std::vector<std::vector<float>>> coal_rate_pair(epoch.size());
  for(int e = 0; e < epoch.size(); e++){
    coal_rate_pair[e].resize(tree.nodes.size());
    for(int i = 0; i < tree.nodes.size(); i++){
      coal_rate_pair[e][i].resize(tree.nodes.size());
      coal_rate_pair[e][i][i] = 0.0;
      for(int j = 0; j < i; j++){
        //calculate coal rate between i and j in epoch e
        coal_rate_pair[e][i][j] = 0.0;
        for(std::vector<int>::iterator it1 = desc[i].member.begin(); it1 != desc[i].member.end(); it1++){
          for(std::vector<int>::iterator it2 = desc[j].member.begin(); it2 != desc[j].member.end(); it2++){
            coal_rate_pair[e][i][j] += coal_rate[e][membership[*it1]][membership[*it2]];
          }
        }        
        coal_rate_pair[e][i][j] /= (desc[i].num_leaves * desc[j].num_leaves);
        coal_rate_pair[e][j][i] = coal_rate_pair[e][i][j];
      }
    }
  }


  if(1){

    p2 = 1.0;
    ////////////////// Transient /////////////////

    count = 0;
    for(; count < 10*delta; count++){

      for(int i = 0; i < N_total-1; i++){
        assert(tree.nodes[i].branch_length >= 0.0);
        assert(order[sorted_indices[i]] == i);
        assert(order[i] < order[(*tree.nodes[i].parent).label]);
        assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
      }

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng <= p2){
        //std::vector<double> coal_rate_one(coal_rate_pair.size(), 1.5);
        //UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate_one, dist_gamma, dist_unif);
        UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate_pair, remaining, dist_gamma, dist_unif);
      }else{ 
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
        if(uniform_rng <= p2){
          int k_candidate = dist_oneevent(rng);
          count_proposals[k_candidate-N]++;
          //std::vector<double> coal_rate_one(coal_rate_pair.size(), 1.5);
          //UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate_one, dist_gamma, dist_unif);
          UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate_pair, remaining, dist_gamma, dist_unif);
          UpdateAvg(tree);
        }else{ 
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
            if(*it_count < 10){
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

    //UNDO
    if(1){
      for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
        (*it_n).branch_length = ((double) Ne) * (avg[(*(*it_n).parent).label] - avg[(*it_n).label]);
      }
    }else{
      for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
        (*it_n).branch_length *= ((double) Ne);
      }
    }

  }

}  


//Make a version that uses pairwise coal rates
void
EstimateBranchLengthsWithSampleAge::MCMCVariablePopulationSizeSample(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, int num_proposals, const bool init, const int seed){

  bool debug = false;

  float uniform_rng;
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 0.7;

  double total_bl = 0.0;
  for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != tree.nodes.end(); it_n++){
    total_bl += (*it_n).branch_length;
  }

  if(init == 1 && total_bl > 0){
    rng.seed(seed);
    root = N_total - 1;

    InitializeMCMC(data, tree); 
    coordinates.resize(N_total);

    GetCoordinates(tree.nodes[root], coordinates);

    //as a sanity check that branch lengths are consistent with sample_ages, calculate branch lengths to root + sample_age
    double dist_to_root, dist_to_root_0 = sample_age[0];
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
    }
  }else if(total_bl == 0){
    ///////////////
    int delta = std::max(data.N/10.0, 10.0);
    root = N_total - 1;

    InitializeMCMC(data, tree);

    std::vector<double> sample_age_tmp;

    bool is_ancient = false;
    for(int i = 0; i < sample_age.size(); i++){
      if(sample_age[i] > 0){
        is_ancient = true;
        break;
      }
    }

    if(is_ancient){
      sample_age_tmp = sample_age;
      std::fill(sample_age.begin(), sample_age.end(), 0.0);
    }

    InitializeOrder(tree); 

    //Randomly switch around order of coalescences
    for(int j = 0; j < (int) 2*data.N * data.N; j++){
      RandomSwitchOrder(tree, dist_n(rng), dist_unif);
    }
    //Initialise branch lengths
    InitializeBranchLengths(tree);

    if(debug){
      for(int i = 0; i < N_total-1; i++){
        assert(tree.nodes[i].branch_length >= 0.0);
        assert(order[sorted_indices[i]] == i);
        assert(order[i] < order[(*tree.nodes[i].parent).label]);
        assert(coordinates[sorted_indices[i]] <= coordinates[sorted_indices[i+1]]);
      }
    }

    if(is_ancient){

      count = 0;
      for(; count < 50*delta; count++){
        //Either switch order of coalescent event or extent time while k ancestors 
        uniform_rng = dist_unif(rng);
        if(uniform_rng <= p2){
          //std::cerr << "v3" << std::endl;
          UpdateOneEventVP(tree, dist_oneevent(rng), epoch, coal_rate, dist_gamma, dist_unif);
        }else{ 
          //std::cerr << "v4" << std::endl;
          SwitchOrder(tree, dist_n(rng), dist_unif);
        }    

        if(debug){
          for(int k = 0; k < N_total; k++){
            assert(tree.nodes[k].branch_length >= 0.0);
          }
        }
      }

      GetCoordinates(tree.nodes[root], coordinates);

      sample_age = sample_age_tmp;
      double min_sample_age = sample_age[0];
      for(int i = 0; i < N; i++){
        if(min_sample_age > sample_age[i]){
          min_sample_age = sample_age[i];
        }
      }
      if(min_sample_age > 0){
        for(std::vector<double>::iterator it_coords = coordinates.begin(); it_coords != coordinates.end(); it_coords++){
          *it_coords += min_sample_age;
        }
      }

      for(int i = 0; i < N; i++){
        if(sample_age[i] > 0){

          int n = (*tree.nodes[i].parent).label;
          if(coordinates[n] > sample_age[i]){
            //tree.nodes[n].branch_length -= sample_age[i];
            coordinates[i] = sample_age[i];
          }else{
            coordinates[i] = sample_age[i];
            float prev_coords = coordinates[i];
            coordinates[n] += sample_age[i];
            prev_coords = coordinates[n];
            while(tree.nodes[n].parent != NULL){
              n = (*tree.nodes[n].parent).label;
              if(coordinates[n] <= prev_coords){
                coordinates[n] += sample_age[i];
                prev_coords = coordinates[n];
              }else{
                break;
              }
            }
          }

        }
      }

      for(int i = 0; i < N_total-1; i++){
        tree.nodes[i].branch_length = coordinates[(*tree.nodes[i].parent).label] - coordinates[i];
      }
      order.clear();
      sorted_indices.clear();
      order.resize(N_total);
      sorted_indices.resize(N_total);

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

    }
    /////////////////////////////////////////////////////////

  }

  ////////////////// Sample branch lengths /////////////////

  count = 0;
  for(; count < num_proposals; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p2){
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

      if(debug){
        for(int i = 0; i < N_total; i++){
          if(!(order[i] == order_foo[i])) std::cerr << order[i] << " " << order_foo[i] << std::endl;
          assert(sorted_indices[i] == sorted_indices_foo[i]);
          assert(order[i] == order_foo[i]);
          assert(num_lineages[i] == num_lineages_foo[i]);
        }
      }

      sorted_indices = sorted_indices_foo;
      order = order_foo;
      num_lineages = num_lineages_foo;
    }

  }

}  

void
EstimateBranchLengthsWithSampleAge::MCMCCoalRatesSample(const Data& data, Tree& tree, const std::vector<int>& membership, const std::vector<double>& epoch, std::vector<std::vector<std::vector<double>>>& coal_rate, int num_proposals, const bool init, const int seed){

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_tip(0,N-1);
  std::uniform_int_distribution<int> dist_n(N,N_total-2);
  std::uniform_int_distribution<int> dist_oneevent(N,N_total-1);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1), p2;
  p1 = 0.0;
  p2 = 1.0;

  double total_bl = 0.0;
  for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != tree.nodes.end(); it_n++){
    total_bl += (*it_n).branch_length;
  }

  if(init == 1 && total_bl > 0){
    rng.seed(seed);
    root = N_total - 1;

    InitializeMCMC(data, tree); 
    coordinates.resize(N_total);

    GetCoordinates(tree.nodes[root], coordinates);

    //as a sanity check that branch lengths are consistent with sample_ages, calculate branch lengths to root + sample_age
    double dist_to_root, dist_to_root_0 = sample_age[0];
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
    }
  }else if(total_bl == 0){
    ///////////////
    int delta = std::max(data.N/10.0, 10.0);
    root = N_total - 1;

    InitializeMCMC(data, tree);

    std::vector<double> sample_age_tmp;

    bool is_ancient = false;
    for(int i = 0; i < sample_age.size(); i++){
      if(sample_age[i] > 0){
        is_ancient = true;
        break;
      }
    }

    if(is_ancient){
      sample_age_tmp = sample_age;
      std::fill(sample_age.begin(), sample_age.end(), 0.0);
    }

    InitializeOrder(tree); 

    //Randomly switch around order of coalescences
    for(int j = 0; j < (int) 2*data.N * data.N; j++){
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

    if(is_ancient){

      count = 0;
      for(; count < 50*delta; count++){
        //Either switch order of coalescent event or extent time while k ancestors 
        uniform_rng = dist_unif(rng);
        if(uniform_rng <= p2){
          //std::cerr << "v3" << std::endl;
          UpdateOneEvent(tree, dist_oneevent(rng), dist_gamma, dist_unif);
        }else{ 
          //std::cerr << "v4" << std::endl;
          SwitchOrder(tree, dist_n(rng), dist_unif);
        }    

        for(int k = 0; k < N_total; k++){
          assert(tree.nodes[k].branch_length >= 0.0);
        }
      }

      GetCoordinates(tree.nodes[root], coordinates);

      sample_age = sample_age_tmp;
      double min_sample_age = sample_age[0];
      for(int i = 0; i < N; i++){
        if(min_sample_age > sample_age[i]){
          min_sample_age = sample_age[i];
        }
      }
      if(min_sample_age > 0){
        for(std::vector<double>::iterator it_coords = coordinates.begin(); it_coords != coordinates.end(); it_coords++){
          *it_coords += min_sample_age;
        }
      }

      for(int i = 0; i < N; i++){
        if(sample_age[i] > 0){

          int n = (*tree.nodes[i].parent).label;
          if(coordinates[n] > sample_age[i]){
            //tree.nodes[n].branch_length -= sample_age[i];
            coordinates[i] = sample_age[i];
          }else{
            coordinates[i] = sample_age[i];
            float prev_coords = coordinates[i];
            coordinates[n] += sample_age[i];
            prev_coords = coordinates[n];
            while(tree.nodes[n].parent != NULL){
              n = (*tree.nodes[n].parent).label;
              if(coordinates[n] <= prev_coords){
                coordinates[n] += sample_age[i];
                prev_coords = coordinates[n];
              }else{
                break;
              }
            }
          }

        }
      }

      for(int i = 0; i < N_total-1; i++){
        tree.nodes[i].branch_length = coordinates[(*tree.nodes[i].parent).label] - coordinates[i];
      }
      order.clear();
      sorted_indices.clear();
      order.resize(N_total);
      sorted_indices.resize(N_total);

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

    }
    /////////////////////////////////////////////////////////

  }

  /////////////////////////////////////////////////////////

  std::vector<Leaves> desc;
  tree.FindAllLeaves(desc);

  //I need which lineages are remaining
  //std::vector<std::vector<float>> num_members(tree.nodes.size());
  //std::vector<float> active_num(sample.groups.size(), 0);
  //for(int i = 0; i < tree.nodes.size(); i++){
  //  num_members.resize(sample.groups.size());
  //}

  remaining.resize(tree.nodes.size());
  remaining_new.resize(tree.nodes.size());
  std::vector<int> active;
  double ages = sample_age[*sorted_indices.begin()];
  std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
  int k = 0;
  for(it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    if(*it_sorted_indices >= N){
      for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        remaining[*it_sorted_indices_start] = active;
      }
      //remove two children and replace by *it_sorted_indices
      int ind1, ind2;
      int c = 0;
      for(std::vector<int>::iterator it_act = active.begin(); it_act != active.end(); it_act++){
        if(*it_act == (*tree.nodes[*it_sorted_indices].child_left).label) ind1 = c;
        if(*it_act == (*tree.nodes[*it_sorted_indices].child_right).label) ind2 = c;
        c++;
      }
      active[ind1] = *it_sorted_indices;
      active[ind2] = active[active.size()-1];
      active.pop_back();

      int c1 = (*tree.nodes[*it_sorted_indices].child_left).label;
      int c2 = (*tree.nodes[*it_sorted_indices].child_right).label;
      float frac = 1.0/desc[c1].num_leaves;
      //for(std::vector<int>::iterator it_mem = desc[c1].member.begin(); it_mem != desc[c1].member.end(); it_mem++){
      //  active_num[membership[*it_mem]] -= frac;
      //}
      frac = 1.0/desc[c2].num_leaves;
      //for(std::vector<int>::iterator it_mem = desc[c2].member.begin(); it_mem != desc[c2].member.end(); it_mem++){
      //  active_num[membership[*it_mem]] -= frac;
      //}
      frac = 1.0/desc[*it_sorted_indices].num_leaves;
      //for(std::vector<int>::iterator it_mem = desc[*it_sorted_indices].member.begin(); it_mem != desc[*it_sorted_indices].member.end(); it_mem++){
      //  active_num[membership[*it_mem]] += frac;
      //}
      //num_member[*it_sorted_indices] = active_num;

      remaining[*it_sorted_indices] = active;
      it_sorted_indices_start++;
    }else if(ages < sample_age[*it_sorted_indices]){      
      for(;it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        remaining[*it_sorted_indices_start] = active;
        //num_member[*it_sorted_indices_start] = active_num;
      }
      ages = sample_age[*it_sorted_indices];
      active.push_back(*it_sorted_indices);
      //active_num[membership[*it_sorted_indices]] += 1.0;
    }else{
      active.push_back(*it_sorted_indices);
      //active_num[membership[*it_sorted_indices]] += 1.0;
    }
  }

  //I need to calculate E x #nodes for coalescence rate
  std::vector<std::vector<std::vector<float>>> coal_rate_pair(epoch.size());
  for(int e = 0; e < epoch.size(); e++){
    coal_rate_pair[e].resize(tree.nodes.size());
    for(int i = 0; i < tree.nodes.size(); i++){
      coal_rate_pair[e][i].resize(tree.nodes.size());
      coal_rate_pair[e][i][i] = 0.0;
      for(int j = 0; j < i; j++){
        //calculate coal rate between i and j in epoch e
        //TODO: use children to speed this up
        coal_rate_pair[e][i][j] = 0.0;
        for(std::vector<int>::iterator it1 = desc[i].member.begin(); it1 != desc[i].member.end(); it1++){
          for(std::vector<int>::iterator it2 = desc[j].member.begin(); it2 != desc[j].member.end(); it2++){
            coal_rate_pair[e][i][j] += coal_rate[e][membership[*it1]][membership[*it2]];
          }
        }        
        coal_rate_pair[e][i][j] /= (desc[i].num_leaves * desc[j].num_leaves);
        coal_rate_pair[e][j][i] = coal_rate_pair[e][i][j];
      }
    }
  }


  float frac = 0.0;
  float count = 0.0;
  std::vector<int> swap_nodes;
  for(int i = data.N; i < tree.nodes.size()-1; i++){
    if(tree.nodes[i].num_events == 0.0){
      frac  += (mut_rate[i] < 2);
      count += 1.0;
      if(mut_rate[i] < 2) swap_nodes.push_back(i);
    }
  }
  std::uniform_int_distribution<int> dist_swap(0,swap_nodes.size()-1);
  //std::cerr << count << " " << frac/count << std::endl;

  ////////////////// Sample branch lengths /////////////////

  if(init == 1 && swap_nodes.size() > 0){
    count = 0;
    for(; count < num_proposals/10.0; count++){
      int n = swap_nodes[dist_swap(rng)];
      SwitchTopo(tree, desc, epoch, coal_rate_pair, remaining, n, dist_unif); //implement switch topology
    }
  }

  if(swap_nodes.size() == 0){
    p2 = 1.0; 
  }else{
    p2 = 0.5;
  }

  count = 0;
  for(; count < num_proposals; count++){

    //std::cerr << count << std::endl;
    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p2){
      int n = dist_oneevent(rng);
      UpdateOneEventVP(tree, n, epoch, coal_rate_pair, remaining, dist_gamma, dist_unif);
    }else{ 
      int n = swap_nodes[dist_swap(rng)];
      SwitchTopo(tree, desc, epoch, coal_rate_pair, remaining, n, dist_unif); //implement switch topology
    }

    for(int i = 0; i < tree.nodes.size()-1; i++){
      //if( order[i] > order[(*tree.nodes[i].parent).label] ){
      //  std::cerr << order[i] << " " << order[(*tree.nodes[i].parent).label] << " " << i << " " << (*tree.nodes[i].parent).label << std::endl;
      //}
      assert( order[i] < order[(*tree.nodes[i].parent).label] );
    }

  }

}  


