#include "branch_length_estimator.hpp"

//////////////////////////////////////////


EstimateBranchLengths::EstimateBranchLengths(const Data& data){
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

  coordinates.resize(N_total);
  sorted_indices.resize(N_total); //node indices in order of coalescent events
  order.resize(N_total); //order of coalescent events
};


//MCMC
void
EstimateBranchLengths::InitializeBranchLengths(Tree& tree){

  int node_i, num_lineages;
  //initialize using coalescent prior
  coordinates.resize(N_total);
  for(int i = 0; i < N; i++){
    coordinates[i] = 0.0;
  }
  for(int i = N; i < N_total; i++){
    num_lineages = 2*N-i;
    node_i = sorted_indices[i];
    coordinates[node_i] = coordinates[sorted_indices[i-1]] + 2.0/( num_lineages * (num_lineages - 1.0) ); // determined by the prior
    (*tree.nodes[node_i].child_left).branch_length  = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_left).label];
    (*tree.nodes[node_i].child_right).branch_length = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_right).label];
  }

}

void
EstimateBranchLengths::InitializeMCMC(const Data& data, Tree& tree){

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
    //mut_rate[i]     = 0.0;
  }

  //initialize
  //1. sort coordinate vector to obtain sorted_indices
  //2. sort sorted_indices to obtain order

  order.resize(N_total);
  sorted_indices.resize(N_total);

  for(int i = 0; i < N_total; i++){
    order[i] = i;
    sorted_indices[i] = i;
  }

  //debug
  //for(int i = 0; i < N_total-1; i++){
  //  assert(order[i] < order[(*tree.nodes[i].parent).label]);
  //}

}

void
EstimateBranchLengths::UpdateAvg(Tree& tree){


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
EstimateBranchLengths::log_deltat(float t){
  if(t >= 0){
    //return(fast_log(1+t));
    if(t < 1){
      return(logt_pos[(int) (t*10000)]);
    }else{
      return(fast_log(1+t));
    }  
  }else{
    //return(fast_log(1+t));
    if(t > -0.1){
      return(logt_neg[(int) (-t*10000)]);
    }else{
      return(fast_log(1+t));
    }  
  }
}


//////////////////////////////////////////

//This function switches the order of coalescent events randomly for initialisation.
//Alternatively, I could just go from bottom to top, and reorder
//TODO: implement version in which I have a list of possible next coalescent events, and choose a random one from there.
void
EstimateBranchLengths::RandomSwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif){

  int node_k, node_swap_k;

  //check order of parent and order of children
  //choose a random number in this range (not including)
  //swap order with that node

  node_k = sorted_indices[k];

  int parent_order    = order[(*tree.nodes[node_k].parent).label];
  int child_order     = order[(*tree.nodes[node_k].child_left).label];
  int child_order_alt = order[(*tree.nodes[node_k].child_right).label];
  if(child_order < child_order_alt) child_order = child_order_alt;
  if(child_order < N) child_order = N-1;
  assert(child_order < parent_order);

  if( parent_order - child_order > 2 ){

    std::uniform_int_distribution<int> d_swap(child_order+1, parent_order-1);
    int new_order = d_swap(rng);

    //calculate odds
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
      }

    }

  }

}

void
EstimateBranchLengths::InitialiseEventOrder(Tree& tree, std::uniform_real_distribution<double>& dist_unif){

  //first I need all events with two descendants,
  //choose one at random and replace by its parent event, but only if the other daughter event is not in the list

  //Binary vector of lengths N-1, specifying is that event is active
  //list with active events

  int num_active_events = 0;
  std::vector<int> events(N-1,0);
  std::vector<int> active_events(N-1,-1);
  std::vector<int> one_side(N-1, 0);

  for(std::vector<Node>::iterator it_node = tree.nodes.begin(); it_node != std::next(tree.nodes.begin(), N); it_node++){
    one_side[(*(*it_node).parent).label - N]++;
  }

  int i = N;
  for(std::vector<int>::iterator it_one_side = one_side.begin(); it_one_side != one_side.end(); it_one_side++){
    if(*it_one_side == 2){
      events[i-N]                      = 1;
      active_events[num_active_events] = i;
      num_active_events++;
    }
    i++;
  }

  std::cerr << num_active_events << std::endl;

  //sorted_indices: kth entry is node label of kth event
  //order: kth entry is order of event k
  int choose_event;
  int i_order = N;
  Node node, parent;
  while(1){

    choose_event                       = (int)(dist_unif(rng) * num_active_events);
    //this is the next event
    sorted_indices[i_order]            = active_events[choose_event];
    order[active_events[choose_event]] = i_order;

    assert(one_side[active_events[choose_event] - N] == 2);
    one_side[active_events[choose_event]-N] = 0;

    //delete active_events[choose_event]
    //either replace by active_events[num_active_events] or by a new event
    node                               = *std::next(tree.nodes.begin(), active_events[choose_event]);
    if(node.parent == NULL){
      break;
    }else{
      parent = *(node.parent);

      if(one_side[parent.label - N] == 0){
        //don't include parent in active events yet
        std::cerr << num_active_events << " " << choose_event << " " << active_events[choose_event] << " " << active_events[num_active_events-1] << std::endl;
        active_events[choose_event]        = active_events[num_active_events-1];
        active_events[num_active_events-1] = -1; 
        events[node.label - N]             = 0; 
        one_side[parent.label - N]++;
        num_active_events--;
      }else if(one_side[parent.label - N] == 1){
        //both daughter events have already been allocated, so include parent
        std::cerr << num_active_events << " " << choose_event << " " << active_events[choose_event] << " " << parent.label << " " << (*(parent).child_left).label << " " << (*(parent).child_right).label << std::endl;
        active_events[choose_event] = parent.label;
        events[node.label - N]      = 0;
        events[parent.label - N]    = 1; 
        one_side[parent.label - N]++;
      } 
    }

    i_order++;

  }

}

void
EstimateBranchLengths::InitialiseEventOrder2(Tree& tree, std::uniform_real_distribution<double>& dist_unif){

  //first I need all events with two descendants,
  //choose one at random and replace by its parent event, but only if the other daughter event is not in the list

  //Binary vector of lengths N-1, specifying is that event is active
  //list with active events

  int num_active_events;
  std::vector<int> active_events(N-1,-1);
  int choose_event;
  int i_order = 2*N-2;
  int child;

  assert(tree.nodes[i_order].parent == NULL);
  order[i_order]        = i_order;
  sorted_indices[i_order] = i_order;

  num_active_events                = 0;
  active_events[0]                 = (*tree.nodes[i_order].child_left).label; 
  if(active_events[0] >= N) num_active_events++;
  active_events[num_active_events] = (*tree.nodes[i_order].child_right).label; 
  if(active_events[num_active_events] >= N) num_active_events++;

  //sorted_indices: kth entry is node label of kth event
  //order: kth entry is order of event k
  i_order--;
  while(num_active_events > 0){

    choose_event                       = (int)(dist_unif(rng) * num_active_events);
    
    //std::cerr << num_active_events << " " << i_order << " " << choose_event << " " << active_events[choose_event] << std::endl;
    //this is the next event
    sorted_indices[i_order]            = active_events[choose_event];
    order[active_events[choose_event]] = i_order;

    //delete active_events[choose_event]
    //either replace by active_events[num_active_events] or by a new event
    n                                  = *std::next(tree.nodes.begin(), active_events[choose_event]);
    child = (*(n.child_left)).label;
    if(child >= N){
      //replace
      active_events[choose_event] = child;
    }else{
      active_events[choose_event]        = active_events[num_active_events-1];
      active_events[num_active_events-1] = -1;
      num_active_events--; 
    }

    child = (*(n.child_right)).label;
    if(child >= N){
      //add
      active_events[num_active_events] = child; 
      num_active_events++;
    }

    i_order--;

  }

}


void
EstimateBranchLengths::SwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif){

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
  int node_k, node_swap_k;
  log_likelihood_ratio = 0.0;

  //check order of parent and order of children
  //choose a random number in this range (not including)
  //swap order with that node

  node_k = sorted_indices[k];

  int parent_order    = order[(*tree.nodes[node_k].parent).label];
  int child_order     = order[(*tree.nodes[node_k].child_left).label];
  int child_order_alt = order[(*tree.nodes[node_k].child_right).label];
  if(child_order < child_order_alt) child_order = child_order_alt;
  if(child_order < N) child_order = N-1;
  if(child_order >= parent_order){
    std::cerr << child_order << " " << k << " " << parent_order << std::endl;
    std::cerr << coordinates[node_k] << " " << coordinates[(*tree.nodes[node_k].child_left).label] << " " << coordinates[(*tree.nodes[node_k].child_right).label] << " " <<  coordinates[(*tree.nodes[node_k].parent).label] << std::endl;
  }
  assert(child_order < parent_order);

  if( parent_order - child_order > 2 ){

    std::uniform_int_distribution<int> d_swap(child_order+1, parent_order-1);
    int new_order = d_swap(rng);

    //calculate odds
    node_swap_k     = sorted_indices[new_order];
    parent_order    = order[(*tree.nodes[node_swap_k].parent).label];
    child_order     = order[(*tree.nodes[node_swap_k].child_left).label];
    child_order_alt = order[(*tree.nodes[node_swap_k].child_right).label];
    if(child_order < child_order_alt) child_order = child_order_alt;
    if(child_order < N) child_order = N-1;

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
            if(n_num_events > 0.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
            if(child_right_num_events > 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
            if(child_left_num_events > 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

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
                  if(n_num_events > 0.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
                  if(child_right_num_events > 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
                  if(child_left_num_events > 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);
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

        double tmp_coords                       = coordinates[node_k];
        coordinates[node_k]                     = coordinates[node_swap_k];
        coordinates[node_swap_k]                = tmp_coords;
        update_node1 = node_k;
        update_node2 = node_swap_k;

        //calculate new branch lengths
        tree.nodes[node_k].branch_length                 = coordinates[(*tree.nodes[node_k].parent).label]  - coordinates[node_k];
        assert(tree.nodes[node_k].branch_length >= 0.0); 
        child_left_label                                 = (*tree.nodes[node_k].child_left).label;
        tree.nodes[child_left_label].branch_length       = coordinates[node_k] - coordinates[child_left_label];
        assert(tree.nodes[child_left_label].branch_length >= 0.0);
        child_right_label                                = (*tree.nodes[node_k].child_right).label;
        tree.nodes[child_right_label].branch_length      = coordinates[node_k] - coordinates[child_right_label];
        assert(tree.nodes[child_right_label].branch_length >= 0.0);

        tree.nodes[node_swap_k].branch_length            = coordinates[(*tree.nodes[node_swap_k].parent).label]  - coordinates[node_swap_k];
        assert(tree.nodes[node_swap_k].branch_length >= 0.0); 
        child_left_label                                 = (*tree.nodes[node_swap_k].child_left).label;
        tree.nodes[child_left_label].branch_length       = coordinates[node_swap_k] - coordinates[child_left_label];
        assert(tree.nodes[child_left_label].branch_length >= 0.0);
        child_right_label                                = (*tree.nodes[node_swap_k].child_right).label;
        tree.nodes[child_right_label].branch_length      = coordinates[node_swap_k] - coordinates[child_right_label];
        assert(tree.nodes[child_right_label].branch_length >= 0.0);

      }
    }

  }

}


//This changes the time of one event, with a beta proposal within the time while k ancestors remain
void
EstimateBranchLengths::UpdateOneEvent(Tree& tree, int k, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif){

  accept = true;
  int node_k;
  log_likelihood_ratio = 0.0;

  node_k = sorted_indices[k];
  //propose to change time of node_k between time of older daughter and parent

  if(tree.nodes[node_k].parent == NULL){

    //propose with exponential distribution  
    num_lineages = 2*N-k;
    assert(num_lineages == 2);

    tau_old   = coordinates[sorted_indices[k]] - coordinates[sorted_indices[k-1]];
    tau_new   = -fast_log(dist_unif(rng)) * tau_old;
    delta_tau = tau_new - tau_old;
    assert(tau_new > 0.0);

    //now decide whether to accept delta_tau:
    //calculate likelihood ratio

    //proposal likelihood ratio
    log_likelihood_ratio = fast_log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);

    //coalescent prior likelihood ratio
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
        if(child_right_num_events > 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
        if(child_left_num_events > 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

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
    }

  }else{

    num_lineages = 2*N-k;
    assert(num_lineages <= N);
    assert(num_lineages > 2);

    tau_below    = coordinates[node_k] - coordinates[sorted_indices[k-1]];
    tau_above    = coordinates[sorted_indices[k+1]] - coordinates[node_k];
    T            = tau_below + tau_above;

    child_left_label  = (*tree.nodes[node_k].child_left).label;
    child_right_label = (*tree.nodes[node_k].child_right).label;
    parent_label      = (*tree.nodes[node_k].parent).label;

    tb_child_left     = tree.nodes[child_left_label].branch_length;
    tb_child_right    = tree.nodes[child_right_label].branch_length;
    tb                = tree.nodes[node_k].branch_length;

    //sample a new tau_below
    x2 = tau_above/tau_below;
    std::gamma_distribution<double> dist_gamma2(2.0 * x2,1.0);

    tau_new_below        = dist_gamma(rng);
    tau_new_below        = tau_new_below/(tau_new_below + dist_gamma2(rng));
    assert(tau_new_below <= 1.0);
    tau_new_below       *= T;
    delta_tau            = tau_new_below - tau_below;
    tau_new_above        = T - tau_new_below;

    if(tau_new_above > 0.0 && tau_new_below > 0.0){

      //now decide whether to accept delta_tau:

      //calculate likelihood ratio

      //proposal likelihood ratio
      x1 = tau_new_above/tau_new_below;
      log_likelihood_ratio = fast_log( x1/x2*(2.0*x1+1.0)/(2.0*x2+1.0)*(tau_below/tau_new_below)) + (2.0*x1 - 1.0) * fast_log(1.0 - 1.0/(1.0+x2)) - (2.0*x2 - 1.0) * fast_log(1.0 - 1.0/(1.0+x1));

      //coalescent prior likelihood ratio
      log_likelihood_ratio -= (num_lineages - 1.0) * delta_tau;

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
            if(n_num_events > 0.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
            if(child_right_num_events > 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
            if(child_left_num_events > 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

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
        tree.nodes[node_k].branch_length            = coordinates[parent_label] - coordinates[node_k];  
      }

    }

  }


}

//This changes the time while k ancestors remain. Want to only make a constant number of updates for each k
void
EstimateBranchLengths::ChangeTimeWhilekAncestors(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif){

  //This step is O(N)
  num_lineages = 2*N-k;
  assert(num_lineages <= N);
  assert(num_lineages >= 2);

  k_choose_2 = num_lineages * (num_lineages-1.0)/2.0;
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

    if(count_number_of_spanning_branches == num_lineages) break;
  }
  //assert(count_number_of_spanning_branches == num_lineages);

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



//This changes the time of one event, with a beta proposal within the time while k ancestors remain
void
EstimateBranchLengths::UpdateOneEventVP(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::gamma_distribution<double>& dist_gamma, std::uniform_real_distribution<double>& dist_unif){

  accept = true;
  int node_k;
  log_likelihood_ratio = 0.0;

  node_k = sorted_indices[k];
  //propose to change time of node_k between time of older daughter and parent

  if(tree.nodes[node_k].parent == NULL){

    //propose with exponential distribution  
    num_lineages = 2*N-k;
    assert(num_lineages == 2);

    tau_old   = coordinates[node_k] - coordinates[sorted_indices[k-1]];
    tau_new   = -fast_log(dist_unif(rng)) * tau_old;
    delta_tau = tau_new - tau_old;
    assert(tau_new > 0.0);

    //now decide whether to accept delta_tau:
    //calculate likelihood ratio

    //proposal likelihood ratio
    log_likelihood_ratio = fast_log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);

    //coalescent prior likelihood ratio
    float new_coordinate = coordinates[node_k] + delta_tau;
    if(delta_tau < 0){

      int ep_begin = 0;
      while(new_coordinate >= epoch[ep_begin]){
        ep_begin++;
        if(ep_begin == (int)epoch.size()) break;
      }
      ep_begin--;
      assert(ep_begin > -1);

      //epoch[ep] <= coordinates[node_k], epoch[ep+1] > coordinates[node_k], unless
      if(new_coordinate >= epoch[ep_begin+1]){
        assert(ep_begin == (int) epoch.size() - 1);
        log_likelihood_ratio -= coal_rate[ep_begin] * delta_tau;
      }else if(coordinates[node_k] <= epoch[ep_begin+1]){
        log_likelihood_ratio -= coal_rate[ep_begin] * delta_tau;
      }else{
        //need to actually calculate things
        float ratio = coal_rate[ep_begin], integral = 0.0;
        while(coordinates[node_k] > epoch[ep_begin+1]){
          integral += coal_rate[ep_begin] * (epoch[ep_begin+1]-epoch[ep_begin]);
          ep_begin++;
        }
        integral += coal_rate[ep_begin] * (coordinates[node_k]-epoch[ep_begin]);
        log_likelihood_ratio += fast_log(ratio/coal_rate[ep_begin]) + integral;
      }

    }else{

      int ep_begin = 0;
      while(coordinates[node_k] >= epoch[ep_begin]){
        ep_begin++;
        if(ep_begin == (int)epoch.size()) break;
      }
      ep_begin--;
      assert(ep_begin > -1);

      //epoch[ep] <= coordinates[node_k], epoch[ep+1] > coordinates[node_k], unless
      if(coordinates[node_k] >= epoch[ep_begin+1]){
        assert(ep_begin == (int) epoch.size() - 1);
        log_likelihood_ratio -= coal_rate[ep_begin] * delta_tau;
      }else if(new_coordinate <= epoch[ep_begin+1]){
        log_likelihood_ratio -= coal_rate[ep_begin] * delta_tau;
      }else{
        //need to actually calculate things
        float ratio = coal_rate[ep_begin], integral = 0.0;
        while(new_coordinate > epoch[ep_begin+1]){
          integral -= coal_rate[ep_begin] * (epoch[ep_begin+1]-epoch[ep_begin]);
          ep_begin++;
        }
        integral -= coal_rate[ep_begin] * (new_coordinate-epoch[ep_begin]);
        log_likelihood_ratio += fast_log(coal_rate[ep_begin]/ratio) + integral;
      }

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
        if(child_right_num_events > 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
        if(child_left_num_events > 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

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
    }

  }else{

    num_lineages = 2*N-k;
    assert(num_lineages <= N);
    assert(num_lineages > 2);

    tau_below    = coordinates[node_k] - coordinates[sorted_indices[k-1]];
    tau_above    = coordinates[sorted_indices[k+1]] - coordinates[node_k];
    T            = tau_below + tau_above;

    child_left_label  = (*tree.nodes[node_k].child_left).label;
    child_right_label = (*tree.nodes[node_k].child_right).label;
    parent_label      = (*tree.nodes[node_k].parent).label;

    tb_child_left     = tree.nodes[child_left_label].branch_length;
    tb_child_right    = tree.nodes[child_right_label].branch_length;
    tb                = tree.nodes[node_k].branch_length;

    //sample a new tau_below
    x2 = tau_above/tau_below;
    std::gamma_distribution<double> dist_gamma2(2.0 * x2,1.0);

    tau_new_below        = dist_gamma(rng);
    tau_new_below        = tau_new_below/(tau_new_below + dist_gamma2(rng));
    assert(tau_new_below <= 1.0);
    tau_new_below       *= T;
    delta_tau            = tau_new_below - tau_below;
    tau_new_above        = T - tau_new_below;

    if(tau_new_above > 0.0 && tau_new_below > 0.0){

      //now decide whether to accept delta_tau:

      //calculate likelihood ratio

      //proposal likelihood ratio
      x1 = tau_new_above/tau_new_below;
      log_likelihood_ratio = fast_log( x1/x2*(2.0*x1+1.0)/(2.0*x2+1.0)*(tau_below/tau_new_below)) + (2.0*x1 - 1.0) * fast_log(1.0 - 1.0/(1.0+x2)) - (2.0*x2 - 1.0) * fast_log(1.0 - 1.0/(1.0+x1));

      //coalescent prior likelihood ratio
      float new_coordinate = coordinates[node_k] + delta_tau;
      if(delta_tau < 0){

        int ep_begin = 0;
        while(new_coordinate >= epoch[ep_begin]){
          ep_begin++;
          if(ep_begin == (int)epoch.size()) break;
        }
        ep_begin--;
        assert(ep_begin > -1);

        //epoch[ep] <= coordinates[node_k], epoch[ep+1] > coordinates[node_k], unless
        if(new_coordinate >= epoch[ep_begin+1]){
          assert(ep_begin == (int) epoch.size() - 1);
          log_likelihood_ratio -= (num_lineages - 1.0) * coal_rate[ep_begin] * delta_tau;
        }else if(coordinates[node_k] <= epoch[ep_begin+1]){
          log_likelihood_ratio -= (num_lineages - 1.0) * coal_rate[ep_begin] * delta_tau;
        }else{
          //need to actually calculate things
          float ratio = coal_rate[ep_begin], integral = 0.0;
          while(coordinates[node_k] > epoch[ep_begin+1]){
            integral += coal_rate[ep_begin] * (epoch[ep_begin+1]-epoch[ep_begin]);
            ep_begin++;
          }
          integral += coal_rate[ep_begin] * (coordinates[node_k]-epoch[ep_begin]);
          log_likelihood_ratio += fast_log(ratio/coal_rate[ep_begin]) + (num_lineages - 1.0) * integral;
        }

      }else{

        int ep_begin = 0;
        while(coordinates[node_k] >= epoch[ep_begin]){
          ep_begin++;
          if(ep_begin == (int)epoch.size()) break;
        }
        ep_begin--;
        assert(ep_begin > -1);

        //epoch[ep] <= coordinates[node_k], epoch[ep+1] > coordinates[node_k], unless
        if(coordinates[node_k] >= epoch[ep_begin+1]){
          assert(ep_begin == (int) epoch.size() - 1);
          log_likelihood_ratio -= (num_lineages - 1.0) * coal_rate[ep_begin] * delta_tau;
        }else if(new_coordinate <= epoch[ep_begin+1]){
          log_likelihood_ratio -= (num_lineages - 1.0) * coal_rate[ep_begin] * delta_tau;
        }else{
          //need to actually calculate things
          float ratio = coal_rate[ep_begin], integral = 0.0;
          while(new_coordinate > epoch[ep_begin+1]){
            integral -= coal_rate[ep_begin] * (epoch[ep_begin+1]-epoch[ep_begin]);
            ep_begin++;
          }
          integral -= coal_rate[ep_begin] * (new_coordinate-epoch[ep_begin]);
          log_likelihood_ratio += fast_log(coal_rate[ep_begin]/ratio) + (num_lineages - 1.0) * integral;
        }

      }

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
            if(n_num_events > 0.0) log_likelihood_ratio += n_num_events * log_deltat(-delta_tau/tb);
            if(child_right_num_events > 0.0) log_likelihood_ratio += child_right_num_events * log_deltat(delta_tau/tb_child_right);
            if(child_left_num_events > 0.0)  log_likelihood_ratio += child_left_num_events * log_deltat(delta_tau/tb_child_left);

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
        tree.nodes[node_k].branch_length            = coordinates[parent_label] - coordinates[node_k];  
      }

    }

  }


}


//propose a new time and change events only up to t+bound
void
EstimateBranchLengths::ChangeTimeWhilekAncestorsVP(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif){

  //This step is O(N)
  num_lineages = 2*N-k;
  assert(num_lineages <= N);
  assert(num_lineages >= 2);

  k_choose_2 = num_lineages * (num_lineages-1.0)/2.0;

  //coorindates[sorted_indices[k-1]] determines the epoch at the lower end.
  //from there, I have to propose a new time by drawing a time for the first epoch, and if it is exceeding the epoch, for the next epoch etc.

  tau_old   = coordinates[sorted_indices[k]] - coordinates[sorted_indices[k-1]];

  log_likelihood_ratio = 0.0;

  if(tau_old > 0.0){
    tau_new   = -log(dist_unif(rng)) * tau_old;
    delta_tau = tau_new - tau_old;
    //calculate ratio of proposals
    log_likelihood_ratio = log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);   
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
  while(coordinates[sorted_indices[k-1]] >= epoch[ep_begin]){
    ep_begin++;
    if(ep_begin == (int)epoch.size()) break;
  }
  ep_begin--;
  assert(ep_begin > -1);
  assert(coordinates[sorted_indices[k-1]] >= epoch[ep_begin]);
  if( coordinates[sorted_indices[k-1]] >= epoch[ep_begin+1]  ){
    assert(ep_begin == (int) epoch.size() - 1);
  }

  int ep;
  double tmp_tau, delta_tmp_tau;
  ep            = ep_begin;
  tmp_tau       = tau_new;

  int k_tmp = k;
  int num_lineages_tmp = num_lineages;
  float k_choose_2_tmp = k_choose_2;
  bool foo = false;
  int end_ep = -1;
  // For k_tmp == k, I am actually changing the time while num_lineages ancestors remain.
  // For k_tmp > k, I need to check whether these were pushed to older epochs by the update, and adjust coalescent times.
  while(k_tmp < 2*N-1){
    if(ep < epoch.size() - 1){

      if(k_tmp > k){
        tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
        delta_tmp_tau = epoch[ep+1] - (coordinates[sorted_indices[k_tmp-1]] + delta_tau);
        //update k_choose_2
        k_choose_2_tmp *= (num_lineages_tmp - 2.0)/num_lineages_tmp;
        num_lineages_tmp--;
        assert(num_lineages_tmp > 1);
      }else{
        assert(coordinates[sorted_indices[k_tmp-1]] >= epoch[ep]);
        delta_tmp_tau = epoch[ep+1] - coordinates[sorted_indices[k_tmp-1]];
      }

      assert(delta_tmp_tau >= 0.0);
      //if epoch[ep+1] begins before the current interval (while num_lineages_tmp remain) ends, enter if statement
      if(delta_tmp_tau <= tmp_tau){

        //add up rate parameters for this interval
        if(coal_rate[ep] > 0.0){
          log_likelihood_ratio   -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
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
          log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
        }

      }else{
        if(coal_rate[ep] == 0){
          log_likelihood_ratio = -std::numeric_limits<float>::infinity();
        }else{
          log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
        }
      }

    }else{

      if(coal_rate[ep] == 0){
        log_likelihood_ratio = -std::numeric_limits<float>::infinity();
      }else{
        if(k_tmp > k){
          tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
        }
        log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
      }
    }
    k_tmp++;

  }

  if(log_likelihood_ratio != -std::numeric_limits<float>::infinity()){

    ep            = ep_begin;
    tmp_tau       = tau_old;

    int k_max = k_tmp;
    //int k_max = 2*N-1;
    k_tmp = k;
    k_choose_2_tmp = k_choose_2;
    num_lineages_tmp = num_lineages;
    while(k_tmp < k_max){

      if(ep < epoch.size() - 1){

        if(k_tmp > k){
          tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
          delta_tmp_tau = epoch[ep+1] - coordinates[sorted_indices[k_tmp-1]];
          //update k_choose_2
          k_choose_2_tmp *= (num_lineages_tmp - 2.0)/num_lineages_tmp;
          num_lineages_tmp--;
        }else{
          delta_tmp_tau = epoch[ep+1] - coordinates[sorted_indices[k_tmp-1]];
        }

        assert(delta_tmp_tau >= 0.0);
        if(delta_tmp_tau <= tmp_tau){
          if(coal_rate[ep] > 0.0){
            log_likelihood_ratio  += k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
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
            log_likelihood_ratio = std::numeric_limits<float>::infinity();
          }else{
            log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
          }

        }else{
          if(coal_rate[ep] == 0){
            log_likelihood_ratio = std::numeric_limits<float>::infinity();
          }else{
            log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
          }
        }

      }else{
        if(coal_rate[ep] == 0){
          log_likelihood_ratio = std::numeric_limits<float>::infinity();
        }else{
          if(k_tmp > k){
            tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
          }
          log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
        }
      }

      k_tmp++;
    }

    ////////////////////////////////////////////
    //likelihood P(Data | tb)

    if(log_likelihood_ratio != std::numeric_limits<float>::infinity()){
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
        if(order[(*tree.nodes[*it_sorted_indices].child_right).label] < k){
          count_number_of_spanning_branches++;

          n = *tree.nodes[*it_sorted_indices].child_right;
          assert(order[n.label] < k);
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
        if(count_number_of_spanning_branches == num_lineages) break;
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
    it_sorted_indices = std::next(sorted_indices.begin(), k);
    update_node1 = k;

    for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      coordinates[*it_sorted_indices]                 += delta_tau;
      assert(coordinates[*it_sorted_indices] >= coordinates[*std::prev(it_sorted_indices,1)]);
      child_left_label                                 = (*tree.nodes[*it_sorted_indices].child_left).label;
      tree.nodes[child_left_label].branch_length       = coordinates[*it_sorted_indices] - coordinates[child_left_label];
      child_right_label                                = (*tree.nodes[*it_sorted_indices].child_right).label;
      tree.nodes[child_right_label].branch_length      = coordinates[*it_sorted_indices] - coordinates[child_right_label];
    }

  }

}


void
EstimateBranchLengths::ChangeTimeWhilekAncestorsVP_new(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif){

  //This step is O(N)
  num_lineages = 2*N-k;
  assert(num_lineages <= N);
  assert(num_lineages >= 2);

  k_choose_2 = num_lineages * (num_lineages-1.0)/2.0;

  //coorindates[sorted_indices[k-1]] determines the epoch at the lower end.
  //from there, I have to propose a new time by drawing a time for the first epoch, and if it is exceeding the epoch, for the next epoch etc.

  int node_k    = sorted_indices[k];
  tau_old   = coordinates[sorted_indices[k]] - coordinates[sorted_indices[k-1]];

  log_likelihood_ratio = 0.0;

  float c = 3.0;
  float prop_c = -1.0/c + c; 

  if(tau_old > 0.0){
    tau_new   = -log(dist_unif(rng)) * tau_old;

    if(tau_new >= c * tau_old){
      tau_new = c * tau_old;
      log_likelihood_ratio = -log(tau_new)+prop_c;   
    }else{
      log_likelihood_ratio = log(tau_old/tau_new) + (tau_new/tau_old - tau_old/tau_new);   
    }
    delta_tau = tau_new - tau_old;
    assert(tau_new > 0.0);
  }else{
    tau_new   = -log(dist_unif(rng)) * 1.0/k_choose_2;
    tau_old   = 0.0;
    if(tau_new >= c/k_choose_2){
      tau_new = c/k_choose_2;
      log_likelihood_ratio = -log(tau_new)-c;
    }else{
      log_likelihood_ratio = log(1.0/(tau_new*k_choose_2)) + tau_new*k_choose_2;
    }
    delta_tau = tau_new;    
    assert(tau_new > 0.0);
  }

  ////////////////////////////////////////////////
  //coalescent prior  
  int num_lineages_tmp = num_lineages;
  float k_choose_2_tmp = k_choose_2;

  //coalescent prior
  //TODO
  //coalescence events up to coordinates[sorted_indices[k-1]] + 2*c*tau_old can change

  float old_coordinate = coordinates[sorted_indices[k]];
  float new_coordinate = old_coordinate + delta_tau;
  float lower_coordinate, upper_coordinate;

  if(delta_tau < 0){

    int ep_begin = 0;
    while(new_coordinate >= epoch[ep_begin]){
      ep_begin++;
      if(ep_begin == (int)epoch.size()) break;
    }
    ep_begin--;
    assert(ep_begin > -1);

    //epoch[ep] <= coordinates[node_k], epoch[ep+1] > coordinates[node_k], unless
    if(new_coordinate >= epoch[ep_begin+1]){
      assert(ep_begin == (int) epoch.size() - 1);
      log_likelihood_ratio -= k_choose_2 * coal_rate[ep_begin] * delta_tau;
    }else{ 
      if(coordinates[node_k] <= epoch[ep_begin+1]){
        log_likelihood_ratio -= k_choose_2 * coal_rate[ep_begin] * delta_tau;
      }else{
        //need to actually calculate things
        float ratio = coal_rate[ep_begin], integral = 0.0;
        while(coordinates[node_k] > epoch[ep_begin+1]){
          integral += coal_rate[ep_begin] * (epoch[ep_begin+1]-epoch[ep_begin]);
          ep_begin++;
        }
        integral += coal_rate[ep_begin] * (coordinates[node_k]-epoch[ep_begin]);
        log_likelihood_ratio += fast_log(ratio/coal_rate[ep_begin]) + k_choose_2 * integral;
      }

      //now calculate the remaining c_k's


    }

  }else{

    int ep_begin = 0;
    int ep, k_tmp;
    double tmp_tau, delta_tmp_tau;

    while(coordinates[node_k] >= epoch[ep_begin]){
      ep_begin++;
      if(ep_begin == (int)epoch.size()) break;
    }
    ep_begin--;
    assert(ep_begin > -1);

    //epoch[ep] <= coordinates[node_k], epoch[ep+1] > coordinates[node_k], unless
    if(coordinates[node_k] >= epoch[ep_begin+1]){
      assert(ep_begin == (int) epoch.size() - 1);
      log_likelihood_ratio -= k_choose_2 * coal_rate[ep_begin] * delta_tau;
    }else{
      if(new_coordinate <= epoch[ep_begin+1]){
        log_likelihood_ratio -= k_choose_2 * coal_rate[ep_begin] * delta_tau;
      }else{
        //need to actually calculate things
        float ratio = coal_rate[ep_begin], integral = 0.0;
        while(new_coordinate > epoch[ep_begin+1]){
          integral -= coal_rate[ep_begin] * (epoch[ep_begin+1]-epoch[ep_begin]);
          ep_begin++;
        }
        integral -= coal_rate[ep_begin] * (new_coordinate-epoch[ep_begin]);
        log_likelihood_ratio += fast_log(coal_rate[ep_begin]/ratio) + k_choose_2 * integral;
      }

      //now calculate the remaining c_k's
      ep = ep_begin;
      k_tmp = k+1;
      num_lineages_tmp = num_lineages;
      k_choose_2_tmp   = k_choose_2;
      while(k_tmp < 2*N-1){

        /*
        //need to sort this out
        if(coordinates[sorted_indices[k_tmp]] > threshold_coordinate) break;
        if(k_tmp < 2*N-2){
          //if next event is above threshold and age to parent is too small
          if(coordinates[sorted_indices[k_tmp+1]] > threshold_coordinate && coordinates[sorted_indices[k_tmp]] + delta_tau >= coordinates[sorted_indices[k_tmp+1]]){
            tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]] - delta_tau;
            break;
          }
        }
        */

        if(ep < epoch.size() - 1){

          tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
          delta_tmp_tau = epoch[ep+1] - (coordinates[sorted_indices[k_tmp-1]] + delta_tau);
          //update k_choose_2
          k_choose_2_tmp *= (num_lineages_tmp - 2.0)/num_lineages_tmp;
          num_lineages_tmp--;
          assert(num_lineages_tmp > 1);

          assert(delta_tmp_tau >= 0.0);
          
          //if epoch[ep+1] begins before the current interval (while num_lineages_tmp remain) ends, enter if statement
          if(delta_tmp_tau <= tmp_tau){

            //add up rate parameters for this interval
            if(coal_rate[ep] > 0.0){
              log_likelihood_ratio   -= k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
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
              log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
            }

          }else{
            if(coal_rate[ep] == 0){
              log_likelihood_ratio = -std::numeric_limits<float>::infinity();
            }else{
              log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
            }
          }

        }else{

          if(coal_rate[ep] == 0){
            log_likelihood_ratio = -std::numeric_limits<float>::infinity();
          }else{
            if(k_tmp > k){
              tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
            }
            log_likelihood_ratio -= k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
          }

        }
        k_tmp++;

      }

      ep            = ep_begin;
      int k_max = k_tmp;
      k_tmp = k+1;
      k_choose_2_tmp = k_choose_2;
      num_lineages_tmp = num_lineages;
      while(k_tmp < k_max){
        if(ep < epoch.size() - 1){

          tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
          delta_tmp_tau = epoch[ep+1] - coordinates[sorted_indices[k_tmp-1]];
          //update k_choose_2
          k_choose_2_tmp *= (num_lineages_tmp - 2.0)/num_lineages_tmp;
          num_lineages_tmp--;

          assert(delta_tmp_tau >= 0.0);
          if(delta_tmp_tau <= tmp_tau){
            if(coal_rate[ep] > 0.0){
              log_likelihood_ratio  += k_choose_2_tmp*coal_rate[ep] * delta_tmp_tau; 
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
              log_likelihood_ratio = std::numeric_limits<float>::infinity();
            }else{
              log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
            }

          }else{
            if(coal_rate[ep] == 0){
              log_likelihood_ratio = std::numeric_limits<float>::infinity();
            }else{
              log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
            }
          }

        }else{
          if(coal_rate[ep] == 0){
            log_likelihood_ratio = std::numeric_limits<float>::infinity();
          }else{
            if(k_tmp > k){
              tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
            }
            log_likelihood_ratio += k_choose_2_tmp*coal_rate[ep] * tmp_tau - log(coal_rate[ep]);
          }
        }

        k_tmp++;
      }

    }

  }

  ////////////////////////////////////////////
  //likelihood P(Data | tb)

  if(log_likelihood_ratio != std::numeric_limits<float>::infinity()){
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
      if(order[(*tree.nodes[*it_sorted_indices].child_right).label] < k){
        count_number_of_spanning_branches++;

        n = *tree.nodes[*it_sorted_indices].child_right;
        assert(order[n.label] < k);
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
      if(count_number_of_spanning_branches == num_lineages) break;
    }
    //assert(count_number_of_spanning_branches == num_lineages);
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
      coordinates[*it_sorted_indices]                 += delta_tau;
      assert(coordinates[*it_sorted_indices] >= coordinates[*std::prev(it_sorted_indices,1)]);
      child_left_label                                 = (*tree.nodes[*it_sorted_indices].child_left).label;
      tree.nodes[child_left_label].branch_length       = coordinates[*it_sorted_indices] - coordinates[child_left_label];
      child_right_label                                = (*tree.nodes[*it_sorted_indices].child_right).label;
      tree.nodes[child_right_label].branch_length      = coordinates[*it_sorted_indices] - coordinates[child_right_label];
    }

  }

}


///////////////////////////////////////////

void
EstimateBranchLengths::GetCoordinates(Node& n, std::vector<double>& coords){

  if(n.child_left != NULL){
    GetCoordinates(*n.child_left, coords);
    GetCoordinates(*n.child_right, coords);

    coords[n.label] = coords[(*n.child_left).label] + (*n.child_left).branch_length;
  }else{
    coords[n.label] = 0.0;
  }

}

void
EstimateBranchLengths::MCMC(const Data& data, Tree& tree, const int seed){

  int delta = std::max(data.N/10.0, 10.0);

  //Random number generators
  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_k(N,N_total-1);
  std::uniform_int_distribution<int> dist_switch(N,N_total-2);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  root = N_total - 1;

  float p1 = std::min(10.0/data.N, 0.1);
  float p2 = p1 + 0.2;  
  //p1       = 0.0;
  //float p2 = 0.8;

  ////////// Initialize MCMC ///////////

  //Initialize MCMC using coalescent prior
  InitializeMCMC(data, tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) data.N * data.N; j++){
    RandomSwitchOrder(tree, dist_switch(rng), dist_unif);
  }  
  
  //InitialiseEventOrder2(tree, dist_unif); 
  //for(int j = 0; j < 10 * (int) data.N; j++){
  //  RandomSwitchOrder(tree, dist_switch(rng), dist_unif);
  //}  

  //Apply EM algorithm to calculate MLE of branch lengths given order of coalescences
  InitializeBranchLengths(tree);

  /////////////// Start MCMC ///////////////

  //transient
  count = 0;
  for(; count < 100*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p1){
      ChangeTimeWhilekAncestors(tree, dist_k(rng), dist_unif);
    }else if(uniform_rng <= p2){
      UpdateOneEvent(tree, dist_k(rng), dist_gamma, dist_unif);
    }else{ 
      SwitchOrder(tree, dist_switch(rng), dist_unif);
    }    

  }

  //Now start estimating branch lenghts. We store the average of coalescent ages in avg calculate branch lengths from here.
  //UpdateAvg is used to update avg.
  //Coalescent ages of the tree are stored in coordinates (this is updated in SwitchOrder and ChangeTimeWhilekAncestors)

  avg              = coordinates;
  last_coordinates = coordinates;
  last_update.resize(N_total);
  std::fill(last_update.begin(), last_update.end(), 1);
  count = 1;

  int num_iterations = 0, iterations_threshold = 500*log(data.N);
  bool is_count_threshold = false;
  is_avg_increasing = false;
  //std::vector<int> counts(data.N-1,0);
  //int update_time_while_k_anc = 0;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < p1){
        //update_time_while_k_anc++;
        int k_candidate = dist_k(rng);
        ChangeTimeWhilekAncestors(tree, k_candidate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng <= p2){
        int k_candidate = dist_k(rng);
        //counts[k_candidate-data.N]++;
        UpdateOneEvent(tree, k_candidate, dist_gamma, dist_unif);
        UpdateAvg(tree);
      }else{ 
        SwitchOrder(tree, dist_switch(rng), dist_unif);
        UpdateAvg(tree);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    is_avg_increasing = true; 
    if(num_iterations < iterations_threshold) is_avg_increasing = false;


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

  /*
     int min_count = counts[0], max_count = counts[0];
     for(std::vector<int>::iterator it_counts = counts.begin(); it_counts != counts.end(); it_counts++){
     if(min_count > *it_counts){
     min_count = *it_counts;
     }
     if(max_count < *it_counts){
     max_count = *it_counts;
     }
     }
     std::cerr << min_count << " " << max_count << " " << update_time_while_k_anc << std::endl;
     */

  //////////// Caluclate branch lengths from avg ////////////

  //don't need to update avg again because I am guaranteed to finish with having updated all nodes
  for(std::vector<Node>::iterator it_n = tree.nodes.begin(); it_n != std::prev(tree.nodes.end(),1); it_n++){
    (*it_n).branch_length = ((double) Ne) * (avg[(*(*it_n).parent).label] - avg[(*it_n).label]);
  }

}  

void
EstimateBranchLengths::MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  //count_accept = 0;
  //count_proposal = 0;

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_k(N,N_total-1);
  std::uniform_int_distribution<int> dist_switch(N,N_total-2);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1);
  float p2 = p1 + 0.2;

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  InitializeMCMC(data, tree); //Initialize using coalescent prior 

  for(std::vector<Node>::iterator it_node = tree.nodes.begin(); it_node != tree.nodes.end(); it_node++){
    (*it_node).branch_length /= data.Ne;
  }

  coordinates.resize(N_total);
  GetCoordinates(tree.nodes[root], coordinates);
  //avg   = coordinates;

  std::size_t m1(0);
  std::generate(std::begin(sorted_indices) + N, std::end(sorted_indices), [&]{ return m1++; });
  std::sort(std::begin(sorted_indices) + N, std::end(sorted_indices), [&](int i1, int i2) { return coordinates[i1 + N] < coordinates[i2 + N]; } );
  for(int i = 0; i < N; i++){
    sorted_indices[i]  = i;
  }
  for(int i = N; i < N_total; i++){
    sorted_indices[i] += N;
  }

  //obtain order of coalescent events
  std::fill(order.begin(), order.end(), 0);
  std::size_t m2(0);
  std::generate(std::begin(order) + N, std::end(order), [&]{ return m2++; });
  std::sort(std::begin(order) + N, std::end(order), [&](int i1, int i2) { return sorted_indices[i1 + N] < sorted_indices[i2 + N]; } );

  for(int i = 0; i < N; i++){
    order[i]  = i;
  }
  for(int i = N; i < N_total; i++){
    order[i] += N;
  }

  //This is for when branch lengths are zero, because then the above does not guarantee parents to be above children
  bool is_topology_violated = true;
  while(is_topology_violated){
    is_topology_violated = false;
    for(int i = N; i < N_total; i++){
      int node_k = sorted_indices[i];
      if( order[(*tree.nodes[node_k].child_left).label] > order[node_k] ){
        int tmp_order = order[node_k];
        order[node_k] = order[(*tree.nodes[node_k].child_left).label];
        order[(*tree.nodes[node_k].child_left).label] = tmp_order;
        sorted_indices[order[node_k]] = node_k;
        sorted_indices[tmp_order] = (*tree.nodes[node_k].child_left).label;
        is_topology_violated = true;
      }
      if( order[(*tree.nodes[node_k].child_right).label] > order[node_k] ){
        int tmp_order = order[node_k];
        order[node_k] = order[(*tree.nodes[node_k].child_right).label];
        order[(*tree.nodes[node_k].child_right).label] = tmp_order;
        sorted_indices[order[node_k]] = node_k;
        sorted_indices[tmp_order] = (*tree.nodes[node_k].child_right).label;
        is_topology_violated = true;
      }
    }
  }

  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 100*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p1){
      ChangeTimeWhilekAncestorsVP(tree, dist_k(rng), epoch, coal_rate, dist_unif);
    }else if(uniform_rng <= p2){
      UpdateOneEventVP(tree, dist_k(rng), epoch, coal_rate, dist_gamma, dist_unif);
    }else{ 
      SwitchOrder(tree, dist_switch(rng), dist_unif);
    }    

  }

  avg              = coordinates;
  last_coordinates = coordinates;
  last_update.resize(N_total);
  std::fill(last_update.begin(), last_update.end(), 1);
  count = 1;

  int num_iterations = 0, iterations_threshold = 500*log(data.N);
  bool is_count_threshold = false;
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < p1){
        //update_time_while_k_anc++;
        int k_candidate = dist_k(rng);
        ChangeTimeWhilekAncestorsVP(tree, k_candidate, epoch, coal_rate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng <= p2){
        int k_candidate = dist_k(rng);
        //counts[k_candidate-data.N]++;
        UpdateOneEventVP(tree, k_candidate, epoch, coal_rate, dist_gamma, dist_unif);
        UpdateAvg(tree);
      }else{ 
        SwitchOrder(tree, dist_switch(rng), dist_unif);
        UpdateAvg(tree);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    is_avg_increasing = true; 
    if(num_iterations < iterations_threshold) is_avg_increasing = false;


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
EstimateBranchLengths::MCMCVariablePopulationSizeForRelate(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  //count_accept = 0;
  //count_proposal = 0;

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_k(N,N_total-1);
  std::uniform_int_distribution<int> dist_switch(N,N_total-2);
  std::gamma_distribution<double> dist_gamma(2.0,1.0);

  float p1 = std::min(10.0/data.N, 0.1);
  float p2 = p1 + 0.2;  

  int delta = std::max(data.N/10.0, 10.0);
  root = N_total - 1;

  ////////// Initialize MCMC ///////////

  //Initialize MCMC using coalescent prior
  InitializeMCMC(data, tree); 

  //Randomly switch around order of coalescences
  //for(int j = 0; j < (int) data.N * data.N; j++){
  //  RandomSwitchOrder(tree, dist_switch(rng), dist_unif);
  //}  
  InitialiseEventOrder2(tree, dist_unif); 
  for(int j = 0; j < (int) data.N; j++){
    RandomSwitchOrder(tree, dist_switch(rng), dist_unif);
  }  

  //Apply EM algorithm to calculate MLE of branch lengths given order of coalescences
  InitializeBranchLengths(tree);

  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 100*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng <= p1){
      ChangeTimeWhilekAncestorsVP(tree, dist_k(rng), epoch, coal_rate, dist_unif);
    }else if(uniform_rng <= p2){
      UpdateOneEventVP(tree, dist_k(rng), epoch, coal_rate, dist_gamma, dist_unif);
    }else{ 
      SwitchOrder(tree, dist_switch(rng), dist_unif);
    }    

  }

  avg              = coordinates;
  last_coordinates = coordinates;
  last_update.resize(N_total);
  std::fill(last_update.begin(), last_update.end(), 1);
  count = 1;

  int num_iterations = 0, iterations_threshold = 500*log(data.N);
  bool is_count_threshold = false;
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < p1){
        //update_time_while_k_anc++;
        int k_candidate = dist_k(rng);
        ChangeTimeWhilekAncestorsVP(tree, k_candidate, epoch, coal_rate, dist_unif);
        UpdateAvg(tree);
      }else if(uniform_rng <= p2){
        int k_candidate = dist_k(rng);
        //counts[k_candidate-data.N]++;
        UpdateOneEventVP(tree, k_candidate, epoch, coal_rate, dist_gamma, dist_unif);
        UpdateAvg(tree);
      }else{ 
        SwitchOrder(tree, dist_switch(rng), dist_unif);
        UpdateAvg(tree);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    is_avg_increasing = true; 
    if(num_iterations < iterations_threshold) is_avg_increasing = false;


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

