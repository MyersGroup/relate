#include "tree_builder.hpp"


bool operator <(const Candidate& a, const Candidate& b){
  return a.dist < b.dist;
}


//////////////////////////////////////////

MinMatch::MinMatch(Data& data){
  N       = data.N;
  N_total = 2*N-1;
  L       = data.L;
  Ne      = data.Ne;
  threshold = -0.2 * log(data.theta/(1.0 - data.theta)); //this is 0.1 of a mutation
  //threshold = 0;

  convert_index.resize(N); //convert_index converts the indices in cluster_index to their actual indices between 0 and N_total-1.
  cluster_size.resize(N);  //stores size of clusters. Accessed using cluster_index.
  min_values.resize(N);
  min_values_sym.resize(N);

  mcandidates.resize(N);
  mcandidates_sym.resize(N);
  updated_cluster.resize(N);
  //min_updated.resize(N);
  //candidates_to_check.resize(N*(N-1)/2); //for candidates where I am not sure if they are valid
  //updated_candidates.resize(N*(N-1)/2); //needed to break ties
};

void
MinMatch::Initialize(CollapsedMatrix<float>& d){

  //I have to go thourgh every pair of clusters and check if they are matching mins.
  //I am also filling in the vector min_values.
  it_min_values_it = min_values.begin(); //iterator for min_values[i]

  std::deque<int>::iterator jt;  //iterator for j
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    mcandidates[*it].dist = std::numeric_limits<float>::infinity();
    d_it = d.rowbegin(*it);
    for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
      if(*it_min_values_it > *d_it && *l != *it){
        *it_min_values_it = *d_it;
      }
      d_it++;
    }
    *it_min_values_it += threshold;
    it_min_values_it++;

  }

  it_min_values_it = min_values.begin();
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    jt                    = std::next(it,1);
    it_d_it_jt            = std::next(d.rowbegin(*it), *jt);
    it_min_values_jt      = std::next(it_min_values_it,1);

    for(; jt != cluster_index.end(); jt++){

      //check if jt is a min value entry, and for each min entry, check if 'it' is also min entry. 
      if( *it_min_values_it >= *it_d_it_jt){        
        if( *it_min_values_jt >= d[*jt][*it]){ //it, jt is a candidate

          sym_dist = *it_d_it_jt + d[*jt][*it];
          if(mcandidates[*it].dist > sym_dist){
            mcandidates[*it].lin1 = *it;
            mcandidates[*it].lin2 = *jt;
            mcandidates[*it].dist = sym_dist;
          }
          if(mcandidates[*jt].dist > sym_dist){
            mcandidates[*jt].lin1 = *it;
            mcandidates[*jt].lin2 = *jt;
            mcandidates[*jt].dist = sym_dist;
          }
          if(best_candidate.dist > mcandidates[*jt].dist){ //only need to check for *jt because if it is the absolute minimum, it would also be *it's candidate
            best_candidate.lin1 = *it;
            best_candidate.lin2 = *jt;
            best_candidate.dist = sym_dist;
          }

        }
      }

      it_d_it_jt++;
      it_min_values_jt++;

    }

    it_min_values_it++;
  }

}

void
MinMatch::InitializeSym(CollapsedMatrix<float>& sym_d, CollapsedMatrix<float>& d){

  std::deque<int>::iterator jt;  //iterator for j
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){
    jt                = std::next(it,1);
    for(; jt != cluster_index.end(); jt++){
      sym_d[*it][*jt] = d[*it][*jt] + d[*jt][*it]; 
      sym_d[*jt][*it] = sym_d[*it][*jt];
    }
  }

  //I have to go thourgh every pair of clusters and check if they are matching mins.
  //I am also filling in the vector min_values.
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    it_min_values_it = std::next(min_values_sym.begin(), *it);
    mcandidates_sym[*it].dist = std::numeric_limits<float>::infinity();
    d_it = sym_d.rowbegin(*it);
    for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
      if(*it_min_values_it > *std::next(d_it,*l) && *l != *it){
        *it_min_values_it = *std::next(d_it,*l);

        if(mcandidates_sym[*it].dist > *it_min_values_it){
          mcandidates_sym[*it].lin1 = *it;
          mcandidates_sym[*it].lin2 = *l;
          mcandidates_sym[*it].dist = *it_min_values_it;
        }
        if(best_sym_candidate.dist > mcandidates_sym[*it].dist){
          best_sym_candidate.lin1 = *it;
          best_sym_candidate.lin2 = *l;
          best_sym_candidate.dist = *it_min_values_it;
        }

      }
    }

  }

}

void
MinMatch::Coalesce(const int i, const int j, CollapsedMatrix<float>& d){

  /////////////
  //now I have to update the distance matrix. I am simultaneously updating min_values[j], where j is the new cluster     
  float added_cluster_size = cluster_size[i] + cluster_size[j];
  float min_value_k, min_value_j = std::numeric_limits<float>::infinity();
  int updated_cluster_size = 0;

  dj_it = d.rowbegin(j);
  di_it = d.rowbegin(i);
  best_candidate.dist = std::numeric_limits<float>::infinity();
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
    if(j != *k && i != *k){
      dk_it = d.rowbegin(*k);
      dkj   = *std::next(dk_it,j);
      dki   = *std::next(dk_it,i);
      dik   = *std::next(di_it,*k);
      djk   = *std::next(dj_it,*k);
      min_value_k = min_values[*k];

      if(dik != djk){ //if dik == djk, the distance of j to k does not change
        *std::next(dj_it,*k) = (cluster_size[i] * dik + cluster_size[j] * djk)/added_cluster_size;
      }
      if(dki != dkj){ //if dki == dkj, the distance of k to j does not change
        *std::next(dk_it,j)  = (cluster_size[i] * dki + cluster_size[j] * dkj)/added_cluster_size;
      }

      bool min_value_changed = false;
      if(dkj != dki){ // if dkj == dki, min_values for k'th row will not change

        if(std::fabs(min_value_k - threshold - dkj) < 1e-4 || std::fabs(min_value_k - threshold - dki) < 1e-4){ //if distance to j or i was minimum distance, we need to update min_values[*k]

          //Note that min_values can only increase. Therefore, if it is the same as before we can break
          min_value_old     = min_value_k - threshold;
          min_value_k       = std::numeric_limits<float>::infinity();
          min_value_changed = true;
          for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
            if(*l != i && *l != *k){  
              if(min_value_k > *std::next(dk_it,*l)){ //update min_values
                min_value_k    = *std::next(dk_it,*l);
                if(min_value_k == min_value_old){
                  break; // same as before, so break 
                }
              }
            }
          }

          min_value_k   += threshold; //add threshold
          min_values[*k] = min_value_k;

        }

      }

      if( dkj != dki || djk != dik ){ //If *std::next(dj_it,*k) and *std::next(dk_it,j) have not changed, the next candidate for *k will be unchanged

        if(min_value_changed || mcandidates[*k].lin1 == j || mcandidates[*k].lin2 == j || mcandidates[*k].lin1 == i || mcandidates[*k].lin2 == i){

          updated_cluster[updated_cluster_size] = *k; //This is to keep track of which candidates have been changed in current iteration. Needed when I decide that candidates for *k don't change (but candidates of *l might change and could be *k) 
          updated_cluster_size++;

          //count2++;
          mcandidates[*k].dist = std::numeric_limits<float>::infinity();

          for(std::deque<int>::iterator l = cluster_index.begin(); l != k; l++){ //only need *l < *k to avoid going through the same pair twice
            if(*std::next(dk_it,*l) <= min_value_k){//this might be a new candidate
              const float min_value_l = min_values[*l];
              if(*l != j && *l != i){

                //add potential candidates
                if(d[*l][*k] <= min_value_l){ 
                  sym_dist = d[*l][*k] + d[*k][*l];
                  if(mcandidates[*k].dist > sym_dist){
                    mcandidates[*k].lin1 = *k;
                    mcandidates[*k].lin2 = *l;
                    mcandidates[*k].dist = sym_dist;
                  }
                  if(mcandidates[*l].dist > sym_dist){
                    mcandidates[*l].lin1 = *k;
                    mcandidates[*l].lin2 = *l;
                    mcandidates[*l].dist = sym_dist;
                  }
                }

              }
            }
          }

        }else{

          //need to check if *k is a candidate for some *l
          for(std::vector<int>::iterator l = updated_cluster.begin(); l != std::next(updated_cluster.begin(),updated_cluster_size); l++){ 
            if(*std::next(dk_it,*l) <= min_value_k){//this might be a new candidate
              const float min_value_l = min_values[*l];

              //add potential candidates
              if(d[*l][*k] <= min_value_l){ 
                sym_dist = d[*l][*k] + d[*k][*l]; 
                if(mcandidates[*l].dist > sym_dist){
                  mcandidates[*l].lin1 = *k;
                  mcandidates[*l].lin2 = *l;
                  mcandidates[*l].dist = sym_dist;
                }
                if(mcandidates[*k].dist > sym_dist){
                  mcandidates[*k].lin1 = *k;
                  mcandidates[*k].lin2 = *l;
                  mcandidates[*k].dist = sym_dist;
                }
              }

            }
          }

        }

      }else{
        //The candidate of *k is unchanged, but it might have been (*k,i), so have to change that to (*k,j)
        if(mcandidates[*k].lin1 == i){
          mcandidates[*k].lin1 = j;
        }
        if(mcandidates[*k].lin2 == i){
          mcandidates[*k].lin2 = j;
        }

        //need to check if *k is a candidate for some *l
        for(std::vector<int>::iterator l = updated_cluster.begin(); l != std::next(updated_cluster.begin(),updated_cluster_size); l++){ 
          if(*std::next(dk_it,*l) <= min_value_k){//this might be a new candidate
            const float min_value_l = min_values[*l];

            //add potential candidates
            if(d[*l][*k] <= min_value_l){ 
              sym_dist = d[*l][*k] + d[*k][*l]; 
              if(mcandidates[*l].dist > sym_dist){
                mcandidates[*l].lin1 = *k;
                mcandidates[*l].lin2 = *l;
                mcandidates[*l].dist = sym_dist;
              }
              if(mcandidates[*k].dist > sym_dist){
                mcandidates[*k].lin1 = *k;
                mcandidates[*k].lin2 = *l;
                mcandidates[*k].dist = sym_dist;
              }
            }

          }
        }

      }

      //I might have updated mcandidates[*l] for *l < *k, but if it was the absolute minimum, it has also updated mcandidates[*k]
      if(best_candidate.dist > mcandidates[*k].dist){ 
        best_candidate   = mcandidates[*k];
      }

      if(*std::next(dj_it,*k) < min_value_j){
        min_value_j   = *std::next(dj_it,*k);
      }
    }
  }
  min_value_j  += threshold;
  min_values[j] = min_value_j;

  //add candidates with new cluster j 
  mcandidates[j].dist = std::numeric_limits<float>::infinity();
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){ 
    if(*std::next(dj_it,*k) <= min_value_j){ 
      if(d[*k][j] <= min_values[*k]){
        if(*k != i && *k != j){

          sym_dist = d[j][*k] + d[*k][j];
          if(mcandidates[*k].dist > sym_dist){
            mcandidates[*k].lin1 = *k;
            mcandidates[*k].lin2 = j;
            mcandidates[*k].dist = sym_dist;
          }
          if(mcandidates[j].dist > sym_dist){
            mcandidates[j].lin1 = *k;
            mcandidates[j].lin2 = j;
            mcandidates[j].dist = sym_dist;
          }

        }
      }
    }
  }

  if(best_candidate.dist > mcandidates[j].dist){
    best_candidate   = mcandidates[j];
  }

}

void
MinMatch::CoalesceSym(const int i, const int j, CollapsedMatrix<float>& sym_d){

  /////////////
  //now I have to update the distance matrix. I am simultaneously updating min_values[j], where j is the new cluster     
  float added_cluster_size = cluster_size[i] + cluster_size[j];
  float min_value_k, min_value_j = std::numeric_limits<float>::infinity();
  int updated_cluster_size = 0;

  dj_it = sym_d.rowbegin(j);
  di_it = sym_d.rowbegin(i);
  best_sym_candidate.dist = std::numeric_limits<float>::infinity();
  mcandidates_sym[j].dist = std::numeric_limits<float>::infinity();
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
    if(j != *k && i != *k){
      dk_it = sym_d.rowbegin(*k);
      dkj   = *std::next(dk_it,j);
      dki   = *std::next(dk_it,i);
      dik   = *std::next(di_it,*k);
      djk   = *std::next(dj_it,*k);
      min_value_k = min_values_sym[*k];

      if(dik != djk){ //if dik == djk, the distance of j to k does not change
        *std::next(dj_it,*k) = (cluster_size[i] * dik + cluster_size[j] * djk)/added_cluster_size;
      }
      if(dki != dkj){ //if dki == dkj, the distance of k to j does not change
        *std::next(dk_it,j)  = (cluster_size[i] * dki + cluster_size[j] * dkj)/added_cluster_size;
      }

      assert(*std::next(dj_it,*k) == *std::next(dk_it,j)); 

      //simplify this

      bool min_value_changed = false;
      if(dkj != dki){ // if dkj == dki, min_values for k'th row will not change

        if(std::fabs(min_value_k - dkj) < 1e-6 || std::fabs(min_value_k - dki) < 1e-6){ //if distance to j or i was minimum distance, we need to update min_values[*k]

          //Note that min_values can only increase. Therefore, if it is the same as before we can break
          min_value_old     = min_value_k;
          min_value_k       = std::numeric_limits<float>::infinity();
          min_value_changed = true;
          mcandidates_sym[*k].dist = std::numeric_limits<float>::infinity();
          for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
            if(*l != i && *l != *k){  
              if(min_value_k > *std::next(dk_it,*l)){ //update min_values
                min_value_k    = *std::next(dk_it,*l);

                if(mcandidates_sym[*k].dist > min_value_k){
                  mcandidates_sym[*k].lin1 = *k;
                  mcandidates_sym[*k].lin2 = *l;
                  mcandidates_sym[*k].dist = min_value_k;
                }

                if(min_value_k == min_value_old){
                  break; // same as before, so break 
                }
              }
            }
          }
          //if(min_value_k < min_value_old){
          //  std::cerr << min_value_k << " " << min_value_old << std::endl;
          //}
          //assert(min_value_k >= min_value_old);
          min_values_sym[*k] = min_value_k;

        }

      }else{
        if(mcandidates_sym[*k].lin1 == i) mcandidates_sym[*k].lin1 = j;
        if(mcandidates_sym[*k].lin2 == i) mcandidates_sym[*k].lin2 = j;
      }

      //I might have updated mcandidates[*l] for *l < *k, but if it was the absolute minimum, it has also updated mcandidates[*k]
      if(best_sym_candidate.dist > mcandidates_sym[*k].dist){ 
        best_sym_candidate   = mcandidates_sym[*k];
      }

      if(*std::next(dj_it,*k) < min_value_j){
        min_value_j   = *std::next(dj_it,*k);
        if(mcandidates_sym[j].dist > *std::next(dj_it,*k)){
          mcandidates_sym[j].lin1 = *k;
          mcandidates_sym[j].lin2 = j;
          mcandidates_sym[j].dist = *std::next(dj_it,*k);
        }
      }

    }
  }
  min_values_sym[j] = min_value_j;

  if(best_sym_candidate.dist > mcandidates_sym[j].dist){
    best_sym_candidate = mcandidates_sym[j];
  }

}

void
MinMatch::QuickBuild(CollapsedMatrix<float>& d, Tree& tree){

  int root = N_total-1;
  tree.nodes.resize(N_total);
  tree.nodes[root].label  = root;

  assert(d.size() > 0);
  assert(d.subVectorSize(0) == d.size());

  //Initialize tree builder

  cluster_index.resize(N); //the cluster index will always stay between 0 - N-1 so that I can access the distance matrix.
  std::deque<int>::iterator it_cluster_index    = cluster_index.begin();
  std::vector<int>::iterator it_convert_index   = convert_index.begin();
  std::vector<float>::iterator it_cluster_size  = cluster_size.begin();
  std::vector<Node>::iterator it_nodes          = tree.nodes.begin();
  int count = 0;
  for(; it_cluster_index != cluster_index.end();){
    *it_cluster_index = count;  //every leaf is in its own cluster
    *it_convert_index = count;
    *it_cluster_size  = 1.0;
    (*it_nodes).label = count;

    it_cluster_index++;
    it_convert_index++;
    it_cluster_size++;
    it_nodes++;
    count++;
  }

  std::fill(min_values.begin(), min_values.end(), std::numeric_limits<float>::infinity());  //Stores the min values of each row of the distance matrix
  std::fill(min_values_sym.begin(), min_values_sym.end(), std::numeric_limits<float>::infinity());  //Stores the min values of each row of the distance matrix

  best_candidate.dist = std::numeric_limits<float>::infinity();
  best_sym_candidate.dist = std::numeric_limits<float>::infinity();

  Initialize(d);

  //////////////////////////

  //Now I need to coalesce clusters until there is only one lance cluster left
  int no_candidate = 0;
  int updated_cluster_size = 0;
  bool use_sym = false;
  for(int num_nodes = N; num_nodes < N_total; num_nodes++){

    //assert(best_sym_candidate.dist != std::numeric_limits<float>::infinity());
    if(best_candidate.dist == std::numeric_limits<float>::infinity()){ //This is a backup-solution for when MinMatch fails, usually rare

      no_candidate++;
      if(!use_sym){
        sym_d.resize(N,N); //inefficient 
        InitializeSym(sym_d, d);
        use_sym = true;
      }
      assert(best_sym_candidate.dist != std::numeric_limits<float>::infinity());
      i = best_sym_candidate.lin1;
      j = best_sym_candidate.lin2;

    }else{
      i = best_candidate.lin1;
      j = best_candidate.lin2;
    }
    conv_i = convert_index[i];
    conv_j = convert_index[j];
    assert(i != j);
    //i and j are now the cluster indices that are closest to each other.

    //merge these two into one new cluster
    tree.nodes[conv_i].parent            = &tree.nodes[num_nodes];
    tree.nodes[conv_j].parent            = &tree.nodes[num_nodes];  
    tree.nodes[conv_j].num_events        = 0.0;
    tree.nodes[conv_i].num_events        = 0.0;
    tree.nodes[num_nodes].child_left     = &tree.nodes[conv_i];
    tree.nodes[num_nodes].child_right    = &tree.nodes[conv_j];
    tree.nodes[num_nodes].label          = num_nodes; 

    Coalesce(i, j, d);
    if(use_sym) CoalesceSym(i, j, sym_d);

    //if(best_candidate.dist == std::numeric_limits<float>::infinity() && no_candidate == 0) no_candidate = num_nodes;
    //if(best_candidate.dist == std::numeric_limits<float>::infinity() && cluster_index.size() > 2) no_candidate++;

    cluster_size[j]  = cluster_size[i] + cluster_size[j]; //update size of new cluster
    convert_index[j] = num_nodes; //update index of this cluster to new merged one
    //delete cluster i    
    for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){ //using deques instead of lists, which makes this loop slower but iteration faster
      if(*it == i){
        cluster_index.erase(it); //invalidates iterators
        break;
      }
    }

  }

  //std::cerr << "Warning: no candidates " << no_candidate << std::endl;
  //if(no_candidate > N/5.0) std::cerr << "Warning: no candidates " << no_candidate << std::endl;
  //if(no_candidate > 0) std::cerr << "Warning: no candidates " << no_candidate << std::endl;
}

void
MinMatch::UPGMA(CollapsedMatrix<float>& d, Tree& tree){

  int root = N_total-1;
  tree.nodes.resize(N_total);
  tree.nodes[root].label  = root;

  assert(d.size() > 0);
  assert(d.subVectorSize(0) == d.size());

  sym_d.resize(N,N); 

  //Initialize tree builder

  cluster_index.resize(N); //the cluster index will always stay between 0 - N-1 so that I can access the distance matrix.
  std::deque<int>::iterator it_cluster_index    = cluster_index.begin();
  std::vector<int>::iterator it_convert_index   = convert_index.begin();
  std::vector<float>::iterator it_cluster_size  = cluster_size.begin();
  std::vector<Node>::iterator it_nodes          = tree.nodes.begin();
  int count = 0;
  for(; it_cluster_index != cluster_index.end();){
    *it_cluster_index = count;  //every leaf is in its own cluster
    *it_convert_index = count;
    *it_cluster_size  = 1.0;
    (*it_nodes).label = count;

    it_cluster_index++;
    it_convert_index++;
    it_cluster_size++;
    it_nodes++;
    count++;
  }

  std::fill(min_values.begin(), min_values.end(), std::numeric_limits<float>::infinity());  //Stores the min values of each row of the distance matrix
  std::fill(min_values_sym.begin(), min_values_sym.end(), std::numeric_limits<float>::infinity());  //Stores the min values of each row of the distance matrix

  best_candidate.dist = std::numeric_limits<float>::infinity();
  best_sym_candidate.dist = std::numeric_limits<float>::infinity();

  //Initialize(d);
  InitializeSym(sym_d, d);

  //////////////////////////

  //Now I need to coalesce clusters until there is only one lance cluster left
  int no_candidate = 0;
  int updated_cluster_size = 0;
  bool use_sym = true;
  for(int num_nodes = N; num_nodes < N_total; num_nodes++){

    assert(best_sym_candidate.dist != std::numeric_limits<float>::infinity());
    i = best_sym_candidate.lin1;
    j = best_sym_candidate.lin2;

    conv_i = convert_index[i];
    conv_j = convert_index[j];
    assert(i != j);
    //i and j are now the cluster indices that are closest to each other.

    //merge these two into one new cluster
    tree.nodes[conv_i].parent            = &tree.nodes[num_nodes];
    tree.nodes[conv_j].parent            = &tree.nodes[num_nodes];  
    tree.nodes[conv_j].num_events        = 0.0;
    tree.nodes[conv_i].num_events        = 0.0;
    tree.nodes[num_nodes].child_left     = &tree.nodes[conv_i];
    tree.nodes[num_nodes].child_right    = &tree.nodes[conv_j];
    tree.nodes[num_nodes].label          = num_nodes; 

    CoalesceSym(i, j, sym_d);

    //if(best_candidate.dist == std::numeric_limits<float>::infinity() && no_candidate == 0) no_candidate = num_nodes;
    //if(best_candidate.dist == std::numeric_limits<float>::infinity() && cluster_index.size() > 2) no_candidate++;

    cluster_size[j]  = cluster_size[i] + cluster_size[j]; //update size of new cluster
    convert_index[j] = num_nodes; //update index of this cluster to new merged one
    //delete cluster i    
    for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){ //using deques instead of lists, which makes this loop slower but iteration faster
      if(*it == i){
        cluster_index.erase(it); //invalidates iterators
        break;
      }
    }

  }

  //std::cerr << "Warning: no candidates " << no_candidate << std::endl;
  if(no_candidate > N/5.0) std::cerr << "Warning: no candidates " << no_candidate << std::endl;

}



//////////////////////////////////////////


InferBranchLengths::InferBranchLengths(const Data& data){
  N       = data.N;
  N_total = 2*N - 1;
  L       = data.L;
  Ne      = data.Ne;

  old_branch_length.resize(N_total);
  coordinates.resize(N_total);
  sorted_indices.resize(N_total); //node indices in order of coalescent events
  order.resize(N_total); //order of coalescent events
};


//MCMC
void
InferBranchLengths::InitializeBranchLengths(Tree& tree){

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
InferBranchLengths::InitializeMCMC(const Data& data, Tree& tree){

  mut_rate.resize(N_total);
  for(int i = 0; i < N_total; i++){
    int snp_begin = tree.nodes[i].SNP_begin;
    int snp_end   = tree.nodes[i].SNP_end;

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

float
InferBranchLengths::LikelihoodGivenTimes(Tree& tree){

  float log_likelihood = 0.0;

  std::vector<Node>::iterator it_nodes = tree.nodes.begin();
  for(; it_nodes != std::prev(tree.nodes.end(),1); it_nodes++){

    if((*it_nodes).num_events == 0){
      log_likelihood -= mut_rate[(*it_nodes).label] * (*it_nodes).branch_length;
    }else{

      if((*it_nodes).branch_length == 0.0){
        log_likelihood = -std::numeric_limits<float>::infinity();
        break; 
      }else{

        log_likelihood -= mut_rate[(*it_nodes).label] * (*it_nodes).branch_length;

        if((*it_nodes).num_events == 1){
          log_likelihood += fast_log(mut_rate[(*it_nodes).label] * (*it_nodes).branch_length);
        }else if((*it_nodes).num_events < logF.size()-1){
          log_likelihood += (*it_nodes).num_events * fast_log(mut_rate[(*it_nodes).label] * (*it_nodes).branch_length) - logF[(*it_nodes).num_events];
        }else{
          log_likelihood += (*it_nodes).num_events * fast_log(mut_rate[(*it_nodes).label] * (*it_nodes).branch_length);
          log_likelihood -= (*it_nodes).num_events * fast_log((*it_nodes).num_events) - (*it_nodes).num_events; //Stirling
        }

      }

    }

  }

  //assert(log_likelihood <= 0.0);
  if(log_likelihood > 0.0) log_likelihood = 0.0;

  return log_likelihood;

}

void
InferBranchLengths::UpdateAvg(Tree& tree){


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

    }else{

      for(std::vector<int>::iterator it_node = std::next(sorted_indices.begin(), update_node1); it_node != sorted_indices.end(); it_node++){
        avg[*it_node]             += ((count - last_update[*it_node]) * (last_coordinates[*it_node] - avg[*it_node]) + coordinates[*it_node] - last_coordinates[*it_node])/count;
        last_update[*it_node]      = count;
        last_coordinates[*it_node] = coordinates[*it_node];
      }
      update_node1 = -1;

    }

  }


  /*  
      it_avg = std::next(avg.begin(), N);
      it_coords = std::next(coordinates.begin(), N);
  //update avg
  for(int ell = N; ell < N_total; ell++){
   *it_avg  += (*it_coords - *it_avg)/count;
   last_update[ell] = count;
   it_avg++;
   it_coords++;
   }
   */


}

void 
InferBranchLengths::logFactorial(int max){

  logF.resize(max+1); 
  std::vector<float>::iterator it_logF = logF.begin();
  *it_logF = 0;
  it_logF++;
  std::vector<float>::iterator it_logF_prev = logF.begin();
  for(int k = 1; k <= max; k++){
    *it_logF = *it_logF_prev + std::log(k);
    it_logF++;
    it_logF_prev++;
  }

}


//////////////////////////////////////////

void
InferBranchLengths::RandomSwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif){

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
InferBranchLengths::SwitchOrder(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif){

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
  //TODO: With tips that have nonzero coalescent time, we need to reject any transition proposed with these tips.


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

    //log_likelihood_ratio = parent_order - child_order - 2;
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
      tb_new             = tb - delta_tau;
      tb_child_left      = tree.nodes[child_left_label].branch_length;
      tb_child_left_new  = tb_child_left + delta_tau;
      tb_child_right     = tree.nodes[child_right_label].branch_length;
      tb_child_right_new = tb_child_right + delta_tau;

      //mutation and recombination part
      if(tb == 0.0){
        log_likelihood_ratio  = std::numeric_limits<float>::infinity();
      }else if(tb_new <= 0.0){
        log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
      }else{

        if(tb_child_left == 0.0){
          log_likelihood_ratio  = std::numeric_limits<float>::infinity();
        }else if(tb_child_left_new <= 0.0){
          log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
        }else{

          if(tb_child_right == 0.0){
            log_likelihood_ratio  = std::numeric_limits<float>::infinity();
          }else if(tb_child_right_new <= 0.0){
            log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
          }else{

            log_likelihood_ratio += (mut_rate[node_k] - mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
            log_likelihood_ratio += n_num_events * log(tb_new/tb);
            log_likelihood_ratio += child_right_num_events * log(tb_child_right_new/tb_child_right); 
            log_likelihood_ratio += child_left_num_events * log(tb_child_left_new/tb_child_left); 

            delta_tau  *= -1.0;
            child_left_label  = (*tree.nodes[node_swap_k].child_left).label;
            child_right_label = (*tree.nodes[node_swap_k].child_right).label;

            n_num_events           = tree.nodes[node_swap_k].num_events;
            child_left_num_events  = tree.nodes[child_left_label].num_events;
            child_right_num_events = tree.nodes[child_right_label].num_events;

            tb                 = tree.nodes[node_swap_k].branch_length;
            tb_new             = tb - delta_tau;
            tb_child_left      = tree.nodes[child_left_label].branch_length;
            tb_child_left_new  = tb_child_left + delta_tau;
            tb_child_right     = tree.nodes[child_right_label].branch_length;
            tb_child_right_new = tb_child_right + delta_tau;

            if(tb == 0.0){
              log_likelihood_ratio  = std::numeric_limits<float>::infinity();
            }else if(tb_new <= 0.0){
              log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
            }else{

              if(tb_child_left == 0.0){
                log_likelihood_ratio  = std::numeric_limits<float>::infinity();
              }else if(tb_child_left_new <= 0.0){
                log_likelihood_ratio  = -std::numeric_limits<float>::infinity(); 
              }else{

                if(tb_child_right == 0.0){
                  log_likelihood_ratio  = std::numeric_limits<float>::infinity();
                }else if(tb_child_right_new <= 0.0){
                  log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
                }else{ 
                  log_likelihood_ratio += (mut_rate[node_swap_k] - mut_rate[child_left_label] - mut_rate[child_right_label]) * delta_tau;
                  log_likelihood_ratio += n_num_events * log(tb_new/tb);
                  log_likelihood_ratio += child_right_num_events * log(tb_child_right_new/tb_child_right); 
                  log_likelihood_ratio += child_left_num_events * log(tb_child_left_new/tb_child_left);

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

      //std::cerr << k << " " << new_order << " " << log_likelihood_ratio << " " << accept << std::endl;
      //update coordinates and sorted_indices
      if(accept && new_order != k){
        //count_accept++;
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

  //count_proposal++;
}

void
InferBranchLengths::ChangeTimeWhilekAncestors(Tree& tree, int k, std::uniform_real_distribution<double>& dist_unif){

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
    //log_likelihood_ratio = (2.0 * alpha - 1.0) * log(tau_old/tau_new) + alpha * (tau_new/tau_old - tau_old/tau_new);
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

      if(tb == 0.0){
        log_likelihood_ratio  = std::numeric_limits<float>::infinity();
        break;
      }else if(tb_new <= 0.0){
        log_likelihood_ratio  = -std::numeric_limits<float>::infinity();
        break;
      }else{
        log_likelihood_ratio -= mut_rate[n.label] * delta_tau;
        log_likelihood_ratio += n.num_events * fast_log(tb_new/tb);
      }
    }
    if(order[(*tree.nodes[*it_sorted_indices].child_right).label] < k){
      count_number_of_spanning_branches++;

      n = *tree.nodes[*it_sorted_indices].child_right;
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
        log_likelihood_ratio += n.num_events * fast_log(tb_new/tb);
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
    //count_accept++;
    //calculate new branch lengths
    it_sorted_indices = std::next(sorted_indices.begin(), k);
    update_node1 = k;
    for(; it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
      coordinates[*it_sorted_indices]                 += delta_tau;
      child_left_label                                 = (*tree.nodes[*it_sorted_indices].child_left).label;
      tree.nodes[child_left_label].branch_length       = coordinates[*it_sorted_indices] - coordinates[child_left_label];
      child_right_label                                = (*tree.nodes[*it_sorted_indices].child_right).label;
      tree.nodes[child_right_label].branch_length      = coordinates[*it_sorted_indices] - coordinates[child_right_label];
    }
  }
  //count_proposal++;

}

void
InferBranchLengths::ChangeTimeWhilekAncestorsVP(Tree& tree, int k, const std::vector<double>& epoch, const std::vector<double>& coal_rate, std::uniform_real_distribution<double>& dist_unif){

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

  //-k_choose_2 * tau_old + k_choose_2 + tau_new
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
        //assert(coordinates[sorted_indices[k_tmp-1]] + delta_tau >= epoch[ep]);
        //assert(coordinates[sorted_indices[k_tmp-1]] + delta_tau < epoch[ep + 1]);
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

      //assert(tmp_tau >= 0.0);
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
        //k_tmp++;
        //break;
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
      //k_tmp++;
      //If I am already at the last epoch, other intervals cannot be pushed to older epochs, so I can break
      //break;
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
          //assert(coordinates[sorted_indices[k_tmp-1]] >= epoch[ep]);
          //assert(coordinates[sorted_indices[k_tmp-1]] < epoch[ep + 1]);
          tmp_tau = coordinates[sorted_indices[k_tmp]] - coordinates[sorted_indices[k_tmp-1]];
          delta_tmp_tau = epoch[ep+1] - coordinates[sorted_indices[k_tmp-1]];
          //update k_choose_2
          k_choose_2_tmp *= (num_lineages_tmp - 2.0)/num_lineages_tmp;
          num_lineages_tmp--;
        }else{
          //assert(coordinates[sorted_indices[k_tmp-1]] >= epoch[ep]);
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

    //if(k == 101) std::cerr << k << " " << ep_begin << " " << coordinates[sorted_indices[k - 1]] * 30000 * 28 << " " << coal_rate[ep_begin] << " " << delta_tau << " " << log_likelihood_ratio << std::endl;

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

  //if(k == 101) std::cerr << log_likelihood_ratio << std::endl;

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
    //count_accept++;
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
  //count_proposal++;

}


///////////////////////////////////////////


void
InferBranchLengths::GetCoordinates(Node& n, std::vector<double>& coords){

  if(n.child_left != NULL){
    GetCoordinates(*n.child_left, coords);
    GetCoordinates(*n.child_right, coords);

    coords[n.label] = coords[(*n.child_left).label] + (*n.child_left).branch_length;

  }else{
    coords[n.label] = 0.0;
  }

}

void
InferBranchLengths::MCMC(const Data& data, Tree& tree, const int seed){

  //count_accept = 0;
  //count_proposal = 0;
  int delta = std::max(data.N/10.0, 10.0);
  convergence_threshold = 10.0/Ne; //approximately x generations

  //Random number generators
  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_k(N,N_total-1);
  std::uniform_int_distribution<int> dist_switch(N,N_total-2);

  root = N_total - 1;
  logFactorial(N); //precalculates log(k!) for k = 0,...,N

  ////////// Initialize MCMC ///////////

  //Initialize MCMC using coalescent prior
  InitializeMCMC(data, tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) data.N * data.N; j++){
    RandomSwitchOrder(tree, dist_switch(rng), dist_unif);
  }   

  //Apply EM algorithm to calculate MLE of branch lengths given order of coalescences
  InitializeBranchLengths(tree);
  EM(data, tree);

  //EM may set some branch lengths to 0. To get a better starting point, we set the minimum branch length to min_tau
  double min_tau = 1.0/Ne, tau_new, tau;
  double push = 0.0, k_choose_2;
  int node_i, num_lineages;
  for(int i = N; i < N_total; i++){
    num_lineages = 2*N-i;
    k_choose_2   = num_lineages * (num_lineages-1.0)/2.0;

    node_i = sorted_indices[i];
    tau    = push + coordinates[node_i] - coordinates[sorted_indices[i-1]];
    assert(tau >= 0.0);

    if(tau < min_tau){

      do{
        tau_new   = -fast_log(dist_unif(rng))/k_choose_2;
        assert(tau_new >= 0.0);
      }while( coordinates[node_i] + push + tau_new - tau < coordinates[sorted_indices[i-1]] );
      push     += tau_new - tau;

    }
    coordinates[node_i] += push;

    if(coordinates[node_i] < coordinates[sorted_indices[i-1]]) std::cerr << coordinates[node_i] << " " << coordinates[sorted_indices[i-1]] << std::endl;
    assert(coordinates[node_i] >= coordinates[sorted_indices[i-1]]);
    (*tree.nodes[node_i].child_left).branch_length  = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_left).label];
    (*tree.nodes[node_i].child_right).branch_length = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_right).label];
  }

  /////////////// Start MCMC ///////////////

  //transient
  count = 0;
  for(; count < 100*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng < 0.5){
      SwitchOrder(tree, dist_switch(rng), dist_unif);
    }else{ 
      ChangeTimeWhilekAncestors(tree, dist_k(rng), dist_unif);
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

  int num_iterations = 0;
  bool is_count_threshold = false;
  std::vector<int> count_proposals(N_total-N, 0);
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < 0.8){
        SwitchOrder(tree, dist_switch(rng), dist_unif);
        UpdateAvg(tree);
      }else{ 
        int k_candidate = dist_k(rng);
        count_proposals[k_candidate-N]++;
        ChangeTimeWhilekAncestors(tree, k_candidate, dist_unif);
        count++;
        UpdateAvg(tree);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    if(!is_count_threshold){
      //assert(is_avg_increasing);
      is_avg_increasing = true;
      for(std::vector<int>::iterator it_count = count_proposals.begin(); it_count != count_proposals.end(); it_count++){
        if(*it_count < 20){
          is_avg_increasing = false;
          break;
        }
      }
      if(is_avg_increasing) is_count_threshold = true;
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
      is_avg_increasing = true;
      for(it_avg = std::next(avg.begin(),N); it_avg != avg.end(); it_avg++){
        if(ell < root){
          if(*it_avg > avg[(*tree.nodes[ell].parent).label]) is_avg_increasing = false;
        }
        ell++;
      }
    }

  }

  //////////// Caluclate branch lengths from avg ////////////

  //don't need to update avg again because I am guaranteed to finish with having updated all nodes
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
InferBranchLengths::MCMCVariablePopulationSize(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  //count_accept = 0;
  //count_proposal = 0;

  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_k(N,N_total-1);
  std::uniform_int_distribution<int> dist_switch(N,N_total-2);

  convergence_threshold = 100.0/Ne; //approximately x generations

  int delta = std::max(data.N/10.0, 10.0);
  //int delta = 10;
  root = N_total - 1;
  logFactorial(N); //precalculates log(k!) for k = 0,...,N

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

  //std::cerr << "Number of branches with no mutations: branches_with_no_mutations " << branches_with_no_mutations.size() << std::endl;

  ////////////////// Transient /////////////////

  count = 0;
  for(; count < 100*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng < 0.5){
      SwitchOrder(tree, dist_switch(rng), dist_unif);
    }else{ 
      ChangeTimeWhilekAncestorsVP(tree, dist_k(rng), epoch, coal_rate, dist_unif);
      //ChangeTimeWhilekAncestorsLocal(tree, dist_k(rng), dist_unif);
    }

  }

  avg              = coordinates;
  last_coordinates = coordinates;
  last_update.resize(N_total);
  std::fill(last_update.begin(), last_update.end(), 1);
  count = 1;

  int num_iterations = 0;
  bool is_count_threshold = false;
  std::vector<int> count_proposals(N_total-N, 0);
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < 0.8){
        SwitchOrder(tree, dist_switch(rng), dist_unif);
        UpdateAvg(tree);
      }else{ 
        int k_candidate = dist_k(rng);
        count_proposals[k_candidate-N]++;
        ChangeTimeWhilekAncestorsVP(tree, k_candidate, epoch, coal_rate, dist_unif);
        count++;
        UpdateAvg(tree);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //assert(is_avg_increasing);
    if(!is_count_threshold){
      is_avg_increasing = true;
      for(std::vector<int>::iterator it_count = count_proposals.begin(); it_count != count_proposals.end(); it_count++){
        if(*it_count < 20){
          is_avg_increasing = false;
          break;
        }
      }
      if(is_avg_increasing) is_count_threshold = true;
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

      //calculate max diff to previous avg coalescent times
      int ell = N;
      is_avg_increasing = true;
      for(it_avg = std::next(avg.begin(),N); it_avg != avg.end(); it_avg++){
        if(ell < root){
          if(*it_avg > avg[(*tree.nodes[ell].parent).label]) is_avg_increasing = false;
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
InferBranchLengths::MCMCVariablePopulationSizeForRelate(const Data& data, Tree& tree, const std::vector<double>& epoch, std::vector<double>& coal_rate, const int seed){

  //count_accept = 0;
  //count_proposal = 0;
  int delta = std::max(data.N/10.0, 10.0);
  convergence_threshold = 10.0/Ne; //approximately x generations

  //Random number generators
  float uniform_rng;
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0,1);
  std::uniform_int_distribution<int> dist_k(N,N_total-1);
  std::uniform_int_distribution<int> dist_switch(N,N_total-2);

  root = N_total - 1;
  logFactorial(N); //precalculates log(k!) for k = 0,...,N

  ////////// Initialize MCMC ///////////

  //Initialize MCMC using coalescent prior
  InitializeMCMC(data, tree); 

  //Randomly switch around order of coalescences
  for(int j = 0; j < (int) data.N * data.N; j++){
    RandomSwitchOrder(tree, dist_switch(rng), dist_unif);
  }   

  //Apply EM algorithm to calculate MLE of branch lengths given order of coalescences
  InitializeBranchLengths(tree);
  EM(data, tree);

  //EM may set some branch lengths to 0. To get a better starting point, we set the minimum branch length to min_tau
  double min_tau = 1.0/Ne, tau_new, tau;
  double push = 0.0, k_choose_2;
  int node_i, num_lineages;
  for(int i = N; i < N_total; i++){
    num_lineages = 2*N-i;
    k_choose_2   = num_lineages * (num_lineages-1.0)/2.0;

    node_i = sorted_indices[i];
    tau    = push + coordinates[node_i] - coordinates[sorted_indices[i-1]];
    assert(tau >= 0.0);

    if(tau < min_tau){

      do{
        tau_new   = -fast_log(dist_unif(rng))/k_choose_2;
        assert(tau_new >= 0.0);
      }while( coordinates[node_i] + push + tau_new - tau < coordinates[sorted_indices[i-1]] );
      push     += tau_new - tau;

    }
    coordinates[node_i] += push;

    if(coordinates[node_i] < coordinates[sorted_indices[i-1]]) std::cerr << coordinates[node_i] << " " << coordinates[sorted_indices[i-1]] << std::endl;
    assert(coordinates[node_i] >= coordinates[sorted_indices[i-1]]);
    (*tree.nodes[node_i].child_left).branch_length  = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_left).label];
    (*tree.nodes[node_i].child_right).branch_length = coordinates[node_i] - coordinates[(*tree.nodes[node_i].child_right).label];
  }

  /////////////// Start MCMC ///////////////

  //transient
  count = 0;
  for(; count < 200*delta; count++){

    //Either switch order of coalescent event or extent time while k ancestors 
    uniform_rng = dist_unif(rng);
    if(uniform_rng < 0.5){
      SwitchOrder(tree, dist_switch(rng), dist_unif);
    }else{
      //ChangeTimeWhilekAncestors(tree, dist_k(rng), dist_unif);
      ChangeTimeWhilekAncestorsVP(tree, dist_k(rng), epoch, coal_rate, dist_unif); 
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

  int num_iterations = 0;
  bool is_count_threshold = false;
  std::vector<int> count_proposals(N_total-N, 0);
  is_avg_increasing = false;
  while(!is_avg_increasing){

    do{

      count++;

      //Either switch order of coalescent event or extent time while k ancestors 
      uniform_rng = dist_unif(rng);
      if(uniform_rng < 0.5){
        SwitchOrder(tree, dist_switch(rng), dist_unif);
        UpdateAvg(tree);
      }else{ 
        int k_candidate = dist_k(rng);
        count_proposals[k_candidate-N]++;
        //ChangeTimeWhilekAncestors(tree, dist_k(rng), dist_unif);
        ChangeTimeWhilekAncestorsVP(tree, dist_k(rng), epoch, coal_rate, dist_unif);
        count++;
        UpdateAvg(tree);
      }

    }while(count % delta != 0 );

    num_iterations++;

    //MCMC is allowed to terminate if is_count_threshold == true and is_avg_increasing == true.
    //At first, both are set to false, and is_count_threshold will be set to true first (once conditions are met)
    //then is_avg_increasing will be set to true eventually.

    if(!is_count_threshold){
      //assert(is_avg_increasing);
      is_avg_increasing = true;
      for(std::vector<int>::iterator it_count = count_proposals.begin(); it_count != count_proposals.end(); it_count++){
        if(*it_count < 20){
          is_avg_increasing = false;
          break;
        }
      }
      if(is_avg_increasing) is_count_threshold = true;
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
      is_avg_increasing = true;
      for(it_avg = std::next(avg.begin(),N); it_avg != avg.end(); it_avg++){
        if(ell < root){
          if(*it_avg > avg[(*tree.nodes[ell].parent).label]) is_avg_increasing = false;
        }
        ell++;
      }
    }

  }

  //////////// Caluclate branch lengths from avg ////////////

  //don't need to update avg again because I am guaranteed to finish with having updated all nodes
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
InferBranchLengths::EM(const Data& data, Tree& tree, bool called_as_main){

  //MLEReorderNodes(tree);
  if(called_as_main){
    convergence_threshold = 10.0/Ne; //approximately 10 generations
    InitializeMCMC(data, tree);
    InitializeBranchLengths(tree);
  }

  //////////////
  //initializing

  std::vector<int>::iterator it_n        = sorted_indices.begin();
  std::vector<int>::iterator it_n_plus_N = sorted_indices.begin() + N;

  ///////////////////////////////////////
  //iterate until branch lengths converge

  double total_branch_length_diff, total_branch_length = std::numeric_limits<float>::infinity(), prev_total_branch_length;
  double prev_old_coordinate;
  double deltat, num_events_on_subbranch, event_prob;
  std::vector<int>::iterator n;
  std::vector<double>::iterator it_coords;
  double prev_coordinate;
  Node *tree_it, *tree_n;
  std::vector<double>::iterator it_old_branch_length;

  //initialize old_branch_length
  old_branch_length.resize(tree.nodes.size());
  std::vector<Node>::iterator it_node = tree.nodes.begin();
  for(it_old_branch_length = old_branch_length.begin(); it_old_branch_length != old_branch_length.end(); it_old_branch_length++){
    *it_old_branch_length = (*it_node).branch_length;
    it_node++;
  }

  do{

    spanning_branches.resize(N); //consider different data structure
    for(int n = 0; n < N; n++){
      spanning_branches[n] = n;
    }

    prev_total_branch_length = total_branch_length;
    total_branch_length = 0.0;
    prev_old_coordinate = 0.0; 

    n = sorted_indices.begin() + N;
    prev_coordinate = 0.0;
    for(; n != sorted_indices.end(); n++){

      tree_n = &tree.nodes[*n];

      it_coords = std::next(coordinates.begin(),*n);
      deltat = *it_coords - prev_old_coordinate; 
      assert(deltat >= 0.0); 

      //for all branches that span the current region, calculate the number of events in that region
      num_events_on_subbranch = 0.0;
      event_prob              = 0.0;
      for(std::deque<int>::iterator it = spanning_branches.begin(); it != spanning_branches.end();){
        tree_it = &tree.nodes[*it];
        it_old_branch_length = std::next(old_branch_length.begin(), *it);
        if(order[(*(*tree_it).parent).label] >= order[*n]){   
          if(*it_old_branch_length < deltat) std::cerr << *it_old_branch_length << " " << deltat << std::endl;
          assert(*it_old_branch_length >= deltat); 
          if(*it_old_branch_length == 0.0){
            num_events_on_subbranch += (*tree_it).num_events;
          }else{
            num_events_on_subbranch += deltat/(*it_old_branch_length) * (*tree_it).num_events;
          }
          event_prob += mut_rate[(*tree_it).label]; //precalculate this line
          it++;
        }else{ 
          *it_old_branch_length = (*tree_it).branch_length;
          it                    = spanning_branches.erase(it); //two branches will be erased per iteration 
        }
      }
      assert(num_events_on_subbranch >= 0.0);
      assert(event_prob >= 0.0);

      //infer new coordinate of n 
      prev_old_coordinate = *it_coords;
      *it_coords          = prev_coordinate + num_events_on_subbranch/(event_prob + ((double) spanning_branches.size() * (spanning_branches.size() - 1.0)/2.0)); 
      prev_coordinate     = *it_coords;
      //update branch lengths
      (*(*tree_n).child_left).branch_length  = *it_coords - coordinates[(*(*tree_n).child_left).label];
      (*(*tree_n).child_right).branch_length = *it_coords - coordinates[(*(*tree_n).child_right).label];
      //if((*(*tree_n).child_left).branch_length < 1e-10)  (*(*tree_n).child_left).branch_length  = 0.0;
      //if((*(*tree_n).child_right).branch_length < 1e-10) (*(*tree_n).child_right).branch_length = 0.0;
      total_branch_length                   += (*(*tree_n).child_left).branch_length + (*(*tree_n).child_right).branch_length;
      spanning_branches.push_back(*n);

    }

    total_branch_length_diff = std::fabs(total_branch_length - prev_total_branch_length)/((double) N_total); //difference per branch   

    //at this point, spanning_branches should only have 3 entries left. (the two branches to the root + the root)
    for(std::deque<int>::iterator it = spanning_branches.begin(); it != spanning_branches.end(); it++){
      old_branch_length[*it] = tree.nodes[*it].branch_length;
    }    

  }while(total_branch_length_diff > convergence_threshold); 


  if(called_as_main){
    for(auto it = tree.nodes.begin(); it != tree.nodes.end(); it++){
      (*it).branch_length *= ((double) Ne);
    }
  }


}

