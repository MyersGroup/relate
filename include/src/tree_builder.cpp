#include "tree_builder.hpp"


//dist contains symmetric distance
//dist2 contains random number
//dist3 contains sample age
bool operator >(const Candidate& a, const Candidate& b){

  if(a.replace == true && a.dist3 >= b.dist3){
    if(a.dist3 > b.dist3) return 1;
    if(a.dist > b.dist || (a.dist == b.dist && a.dist2 > b.dist2)){
      return 1;
    }
  }

  if(a.dist > b.dist || (a.dist == b.dist && a.dist2 > b.dist2)){
    return 1;
  }

  return 0;

}

Candidate& Candidate::operator =(const Candidate& a){
  lin1    = a.lin1;
  lin2    = a.lin2;
  dist    = a.dist;
  dist2   = a.dist2;
  dist3   = a.dist3;
  replace = a.replace;
  is_from_tmpl = a.is_from_tmpl;
  return *this;
}


//////////////////////////////////////////

MinMatch::MinMatch(Data& data){
  N       = data.N;
  N_total = 2*N-1;
  L       = data.L;
  Ne      = data.Ne;
  threshold = -0.2 * log(data.theta/(1.0 - data.theta)); //this is 0.1 of a mutation
  threshold_CF = -0.001 * log(data.theta/(1.0 - data.theta));

  convert_index.resize(N); //convert_index converts the indices in cluster_index to their actual indices between 0 and N_total-1.
  cluster_size.resize(N);  //stores size of clusters. Accessed using cluster_index.
  min_values.resize(N);
  min_values_sym.resize(N);
  min_values_CF.resize(N);

  mcandidates.resize(N);
  mcandidates_sym.resize(N);
  updated_cluster.resize(N);

};

void
MinMatch::Initialize(CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, Tree *tmpl_tree){

  //I have to go thourgh every pair of clusters and check if they are matching mins.
  //I am also filling in the vector min_values.
  it_min_values_it = min_values.begin(); //iterator for min_values[i]

  std::deque<int>::iterator jt;  //iterator for j
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    mcandidates[*it].dist = std::numeric_limits<float>::infinity();
    mcandidates[*it].dist2 = std::numeric_limits<float>::infinity();
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

          if(tmpl_tree != NULL){
            //check if I coalesced these two
            if((*tmpl_tree).nodes[*it].parent == (*tmpl_tree).nodes[*jt].parent){
              mcandidates[*it].lin1 = *it;
              mcandidates[*it].lin2 = *jt;
              mcandidates[*it].dist = *it_d_it_jt + d[*jt][*it];
              mcandidates[*it].dist2 = dist_unif(rng);
              mcandidates[*it].replace = false;
              mcandidates[*it].is_from_tmpl = true;
              mcandidates[*jt] = mcandidates[*it];
              mcandidates[*jt].lin1 = *jt;
              mcandidates[*jt].lin2 = *it;
            }
          }

          //sym_dist = *it_d_it_jt + d[*jt][*it] + std::max(sample_ages[*it], sample_ages[*jt]);
          //sym_dist = *it_d_it_jt + d[*jt][*it] + (sample_ages[*it] * cluster_size[*it] + sample_ages[*jt] * cluster_size[*jt])/(cluster_size[*it] + cluster_size[*jt]);
          sym_dist = *it_d_it_jt + d[*jt][*it];
          dist_random = dist_unif(rng);
          if(!mcandidates[*it].is_from_tmpl && (mcandidates[*it].dist > sym_dist || (mcandidates[*it].dist == sym_dist && mcandidates[*it].dist2 > dist_random))){
            mcandidates[*it].lin1 = *it;
            mcandidates[*it].lin2 = *jt;
            mcandidates[*it].dist = sym_dist;
            mcandidates[*it].dist2 = dist_random;
          }
          if(!mcandidates[*jt].is_from_tmpl && (mcandidates[*jt].dist > sym_dist || (mcandidates[*jt].dist == sym_dist && mcandidates[*jt].dist2 > dist_random))){
            mcandidates[*jt].lin1 = *it;
            mcandidates[*jt].lin2 = *jt;
            mcandidates[*jt].dist = sym_dist;
            mcandidates[*jt].dist2 = dist_random;
          }

          if(best_candidate.dist > mcandidates[*jt].dist || (best_candidate.dist == mcandidates[*jt].dist && best_candidate.dist2 > mcandidates[*jt].dist2)){ //only need to check for *jt because if it is the absolute minimum, it would also be *it's candidate
            best_candidate.lin1 = *it;
            best_candidate.lin2 = *jt;
            best_candidate.dist = sym_dist;
            best_candidate.dist2 = mcandidates[*jt].dist2;
            best_candidate.is_from_tmpl = mcandidates[*jt].is_from_tmpl;
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
MinMatch::Initialize(CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages, Tree *tmpl_tree){

  //
  //I have to go thourgh every pair of clusters and check if they are matching mins.
  //I am also filling in the vector min_values.
  it_min_values_it = min_values.begin(); //iterator for min_values[i]

  std::deque<int>::iterator jt;  //iterator for j
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    mcandidates[*it].dist = std::numeric_limits<float>::infinity();
    mcandidates[*it].dist2 = std::numeric_limits<float>::infinity();
    mcandidates[*it].dist3 = std::numeric_limits<float>::infinity();
    mcandidates[*it].replace = false;
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

          if(tmpl_tree != NULL){
            //check if I coalesced these two
            if((*tmpl_tree).nodes[*it].parent == (*tmpl_tree).nodes[*jt].parent){
              mcandidates[*it].lin1    = *it;
              mcandidates[*it].lin2    = *jt;
              mcandidates[*it].dist    = *it_d_it_jt + d[*jt][*it];
              mcandidates[*it].dist2   = dist_unif(rng);
              mcandidates[*it].dist3   = std::max(sample_ages[*it], sample_ages[*jt]);
              mcandidates[*it].replace = false;
              mcandidates[*it].is_from_tmpl = true;
              mcandidates[*jt]         = mcandidates[*it];
              mcandidates[*jt].lin1    = *jt;
              mcandidates[*jt].lin2    = *it;
            }
          }

          cand.dist = *it_d_it_jt + d[*jt][*it];
          //cand.dist3 = (sample_ages[*it] * cluster_size[*it] + sample_ages[*jt] * cluster_size[*jt])/(cluster_size[*it] + cluster_size[*jt]);
          cand.dist3 = std::max(sample_ages[*it], sample_ages[*jt]);
          cand.dist2 = dist_unif(rng);
          if( !mcandidates[*it].is_from_tmpl && ((mcandidates[*it].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*it] > cand) ){
            if(cand.dist3 > age){
              cand.replace = true;
            }else{
              cand.replace = false;
            }
            mcandidates[*it] = cand;
            mcandidates[*it].lin1 = *it;
            mcandidates[*it].lin2 = *jt;
          }
          if( !mcandidates[*jt].is_from_tmpl && ((mcandidates[*jt].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*jt] > cand) ){
            if(cand.dist3 > age){  
              cand.replace = true;
            }else{
              cand.replace = false;
            }
            mcandidates[*jt] = cand;
            mcandidates[*jt].lin1 = *it;
            mcandidates[*jt].lin2 = *jt;
          }

          if( (best_candidate.dist == std::numeric_limits<float>::infinity() || mcandidates[*jt].dist3 <= age) && best_candidate > mcandidates[*jt]){ //only need to check for *jt because if it is the absolute minimum, it would also be *it's candidate
            best_candidate = mcandidates[*jt];
            if(best_candidate.dist3 > age){
              best_candidate.replace = true;
            }else{
              best_candidate.replace = false;
            }
            best_candidate.is_from_tmpl = mcandidates[*jt].is_from_tmpl;
          }
          //std::cerr << best_candidate.dist3 << " " << best_candidate.dist << " " << best_candidate.replace << " " << age << " " << mcandidates[*jt].dist3 << " " << best_candidate.lin1 << " " << best_candidate.lin2 << std::endl;

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
      sym_d[*it][*jt] = d[*it][*jt] + d[*jt][*it];// + (sample_ages[*it] + sample_ages[*jt])/2.0; 
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
MinMatch::Coalesce(const int i, const int j, CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, Tree *tmpl_tree){

  /////////////
  //now I have to update the distance matrix. I am simultaneously updating min_values[j], where j is the new cluster     
  float added_cluster_size = cluster_size[i] + cluster_size[j];
  float min_value_k, min_value_j = std::numeric_limits<float>::infinity();
  int updated_cluster_size = 0;

  dj_it = d.rowbegin(j);
  di_it = d.rowbegin(i);
  best_candidate.dist = std::numeric_limits<float>::infinity();
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
  best_candidate.is_from_tmpl = false;
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

      if( dkj != dki || djk != dik || mcandidates[*k].lin1 == j || mcandidates[*k].lin2 == j || mcandidates[*k].lin1 == i || mcandidates[*k].lin2 == i ){ //If *std::next(dj_it,*k) and *std::next(dk_it,j) have not changed, the next candidate for *k will be unchanged

        if(min_value_changed || mcandidates[*k].lin1 == j || mcandidates[*k].lin2 == j || mcandidates[*k].lin1 == i || mcandidates[*k].lin2 == i){

          updated_cluster[updated_cluster_size] = *k; //This is to keep track of which candidates have been changed in current iteration. Needed when I decide that candidates for *k don't change (but candidates of *l might change and could be *k) 
          updated_cluster_size++;

          //count2++;
          mcandidates[*k].dist = std::numeric_limits<float>::infinity();
          mcandidates[*k].dist2 = std::numeric_limits<float>::infinity();
          mcandidates[*k].is_from_tmpl = false;

          for(std::deque<int>::iterator l = cluster_index.begin(); l != k; l++){ //only need *l < *k to avoid going through the same pair twice
            if(*std::next(dk_it,*l) <= min_value_k){//this might be a new candidate
              const float min_value_l = min_values[*l];
              if(*l != j && *l != i){

                //add potential candidates
                if(d[*l][*k] <= min_value_l){ 

                  //feasible candidate
                  if(tmpl_tree != NULL){
                    int n1 = conv_to_tmpl[convert_index[*l]];
                    int n2 = conv_to_tmpl[convert_index[*k]];
                    if(n1 >= 0 && n2 >= 0){
                      if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                        mcandidates[*k].lin1 = *k;
                        mcandidates[*k].lin2 = *l;
                        mcandidates[*k].dist = d[*l][*k] + d[*k][*l];
                        mcandidates[*k].dist2 = dist_unif(rng);
                        mcandidates[*k].is_from_tmpl = true;
                        mcandidates[*l] = mcandidates[*k];
                        mcandidates[*l].lin1 = *l;
                        mcandidates[*l].lin2 = *k;
                      } 
                    }
                  }

                  //sym_dist = d[*l][*k] + d[*k][*l] + std::max(sample_ages[*k], sample_ages[*l]);
                  //sym_dist = d[*l][*k] + d[*k][*l] + (sample_ages[*k] * cluster_size[*k] + sample_ages[*l] * cluster_size[*l])/(cluster_size[*k] + cluster_size[*l]);
                  sym_dist = d[*l][*k] + d[*k][*l];
                  //sym_dist = sample_ages[*k] + sample_ages[*l];
                  dist_random = dist_unif(rng);
                  if( !mcandidates[*k].is_from_tmpl && (mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random))){
                    mcandidates[*k].lin1 = *k;
                    mcandidates[*k].lin2 = *l;
                    mcandidates[*k].dist = sym_dist;
                    mcandidates[*k].dist2 = dist_random;
                    mcandidates[*k].is_from_tmpl = false;
                  }
                  if( !mcandidates[*l].is_from_tmpl && (mcandidates[*l].dist > sym_dist || (mcandidates[*l].dist == sym_dist && mcandidates[*l].dist2 > dist_random))){
                    mcandidates[*l].lin1 = *k;
                    mcandidates[*l].lin2 = *l;
                    mcandidates[*l].dist = sym_dist;
                    mcandidates[*l].dist2 = dist_random;
                    mcandidates[*l].is_from_tmpl = false;
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

                //feasible candidate
                if(tmpl_tree != NULL){
                  int n1 = conv_to_tmpl[convert_index[*l]];
                  int n2 = conv_to_tmpl[convert_index[*k]];
                  if(n1 >= 0 && n2 >= 0){
                    if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                      mcandidates[*k].lin1 = *k;
                      mcandidates[*k].lin2 = *l;
                      mcandidates[*k].dist = d[*l][*k] + d[*k][*l];
                      mcandidates[*k].dist2 = dist_unif(rng);
                      mcandidates[*k].is_from_tmpl = true;
                      mcandidates[*l] = mcandidates[*k];
                      mcandidates[*l].lin1 = *l;
                      mcandidates[*l].lin2 = *k;
                    } 
                  }
                }

                //sym_dist = d[*l][*k] + d[*k][*l] + std::max(sample_ages[*l], sample_ages[*k]);
                //sym_dist = d[*l][*k] + d[*k][*l] + (sample_ages[*l] * cluster_size[*l] + sample_ages[*k] * cluster_size[*k])/(cluster_size[*l] + cluster_size[*k]);
                sym_dist = d[*l][*k] + d[*k][*l]; 
                //sym_dist = sample_ages[*l] + sample_ages[*k]; 
                dist_random = dist_unif(rng);
                if( !mcandidates[*l].is_from_tmpl && (mcandidates[*l].dist > sym_dist || (mcandidates[*l].dist == sym_dist && mcandidates[*l].dist2 > dist_random))){
                  mcandidates[*l].lin1 = *k;
                  mcandidates[*l].lin2 = *l;
                  mcandidates[*l].dist = sym_dist;
                  mcandidates[*l].dist2 = dist_random;
                  mcandidates[*l].is_from_tmpl = false;
                }
                if( !mcandidates[*k].is_from_tmpl && (mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random))){
                  mcandidates[*k].lin1 = *k;
                  mcandidates[*k].lin2 = *l;
                  mcandidates[*k].dist = sym_dist;
                  mcandidates[*k].dist2 = dist_random;
                  mcandidates[*k].is_from_tmpl = false;
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

              //feasible candidate
              if(tmpl_tree != NULL){
                int n1 = conv_to_tmpl[convert_index[*l]];
                int n2 = conv_to_tmpl[convert_index[*k]];
                if(n1 >= 0 && n2 >= 0){
                  if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                    mcandidates[*k].lin1 = *k;
                    mcandidates[*k].lin2 = *l;
                    mcandidates[*k].dist = d[*l][*k] + d[*k][*l];
                    mcandidates[*k].dist2 = dist_unif(rng);
                    mcandidates[*k].is_from_tmpl = true;
                    mcandidates[*l] = mcandidates[*k];
                    mcandidates[*l].lin1 = *l;
                    mcandidates[*l].lin2 = *k;
                  }
                }
              }

              //sym_dist = d[*l][*k] + d[*k][*l] + std::max(sample_ages[*l], sample_ages[*k]);
              //sym_dist = d[*l][*k] + d[*k][*l] + (sample_ages[*l] * cluster_size[*l] + sample_ages[*k] * cluster_size[*k])/(cluster_size[*k] + cluster_size[*l]);
              sym_dist = d[*l][*k] + d[*k][*l];
              //sym_dist = sample_ages[*l] + sample_ages[*k]; 
              dist_random = dist_unif(rng);
              if( !mcandidates[*l].is_from_tmpl && (mcandidates[*l].dist > sym_dist || (mcandidates[*l].dist == sym_dist && mcandidates[*l].dist2 > dist_random))){
                mcandidates[*l].lin1 = *k;
                mcandidates[*l].lin2 = *l;
                mcandidates[*l].dist = sym_dist;
                mcandidates[*l].dist2 = dist_random;
                mcandidates[*l].is_from_tmpl = false;
              }
              if( !mcandidates[*k].is_from_tmpl && (mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random))){
                mcandidates[*k].lin1 = *k;
                mcandidates[*k].lin2 = *l;
                mcandidates[*k].dist = sym_dist;
                mcandidates[*k].dist2 = dist_random;
                mcandidates[*k].is_from_tmpl = false;
              }
            }

          }
        }

      }

      //I might have updated mcandidates[*l] for *l < *k, but if it was the absolute minimum, it has also updated mcandidates[*k]
      if(best_candidate.dist > mcandidates[*k].dist || (best_candidate.dist == mcandidates[*k].dist && best_candidate.dist2 > mcandidates[*k].dist2)){ 
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
  mcandidates[j].dist2 = std::numeric_limits<float>::infinity();
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){ 
    if(*std::next(dj_it,*k) <= min_value_j){ 
      if(d[*k][j] <= min_values[*k]){
        if(*k != i && *k != j){

          //feasible candidate
          if(tmpl_tree != NULL){
            int n1 = conv_to_tmpl[convert_index[j]];
            int n2 = conv_to_tmpl[convert_index[*k]];
            if(n1 >= 0 && n2 >= 0){
              if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                mcandidates[*k].lin1 = *k;
                mcandidates[*k].lin2 = j;
                mcandidates[*k].dist = d[j][*k] + d[*k][j];
                mcandidates[*k].dist2 = dist_unif(rng);
                mcandidates[*k].is_from_tmpl = true;
                mcandidates[j] = mcandidates[*k];
                mcandidates[j].lin1 = j;
                mcandidates[j].lin2 = *k;
              }
            }
          }

          //sym_dist = d[j][*k] + d[*k][j] + std::max(sample_ages[j], sample_ages[*k]);
          //sym_dist = d[j][*k] + d[*k][j] + (sample_ages[j] + cluster_size[j] + sample_ages[*k] * cluster_size[*k])/(cluster_size[j] + cluster_size[*k]);
          sym_dist = d[j][*k] + d[*k][j];
          //sym_dist = sample_ages[j] + sample_ages[*k];
          dist_random = dist_unif(rng);
          if( !mcandidates[*k].is_from_tmpl && (mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random))){
            mcandidates[*k].lin1 = *k;
            mcandidates[*k].lin2 = j;
            mcandidates[*k].dist = sym_dist;
            mcandidates[*k].dist2 = dist_random;
            mcandidates[*k].is_from_tmpl = false;
          }
          if( !mcandidates[j].is_from_tmpl && (mcandidates[j].dist > sym_dist || (mcandidates[j].dist == sym_dist && mcandidates[j].dist2 > dist_random))){
            mcandidates[j].lin1 = *k;
            mcandidates[j].lin2 = j;
            mcandidates[j].dist = sym_dist;
            mcandidates[j].dist2 = dist_random;
            mcandidates[j].is_from_tmpl = false;
          }

        }
      }
    }
  }

  if(best_candidate.dist > mcandidates[j].dist || (best_candidate.dist == mcandidates[j].dist && best_candidate.dist2 > mcandidates[j].dist2)){
    best_candidate   = mcandidates[j];
  }

}

void
MinMatch::Coalesce(const int i, const int j, CollapsedMatrix<float>& d, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages, Tree *tmpl_tree){

  /////////////
  //now I have to update the distance matrix. I am simultaneously updating min_values[j], where j is the new cluster     
  float added_cluster_size = cluster_size[i] + cluster_size[j];
  float min_value_k, min_value_j = std::numeric_limits<float>::infinity();
  int updated_cluster_size = 0;

  dj_it = d.rowbegin(j);
  di_it = d.rowbegin(i);
  best_candidate.dist = std::numeric_limits<float>::infinity();
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
  best_candidate.dist3 = std::numeric_limits<float>::infinity();
  best_candidate.replace = false;
  best_candidate.is_from_tmpl = false;
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
    if(j != *k && i != *k){
      dk_it = d.rowbegin(*k);
      dkj   = *std::next(dk_it,j);
      dki   = *std::next(dk_it,i);
      dik   = *std::next(di_it,*k);
      djk   = *std::next(dj_it,*k);
      min_value_k = min_values[*k];
      if(mcandidates[*k].dist3 <= age) mcandidates[*k].replace = false;

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

      if( dkj != dki || djk != dik || mcandidates[*k].lin1 == j || mcandidates[*k].lin2 == j || mcandidates[*k].lin1 == i || mcandidates[*k].lin2 == i ){ //If *std::next(dj_it,*k) and *std::next(dk_it,j) have not changed, the next candidate for *k will be unchanged

        if(min_value_changed || mcandidates[*k].lin1 == j || mcandidates[*k].lin2 == j || mcandidates[*k].lin1 == i || mcandidates[*k].lin2 == i){

          updated_cluster[updated_cluster_size] = *k; //This is to keep track of which candidates have been changed in current iteration. Needed when I decide that candidates for *k don't change (but candidates of *l might change and could be *k) 
          updated_cluster_size++;

          //count2++;
          mcandidates[*k].dist = std::numeric_limits<float>::infinity();
          mcandidates[*k].dist2 = std::numeric_limits<float>::infinity();
          mcandidates[*k].dist3 = std::numeric_limits<float>::infinity();
          mcandidates[*k].replace = false;
          for(std::deque<int>::iterator l = cluster_index.begin(); l != k; l++){ //only need *l < *k to avoid going through the same pair twice
            if(*std::next(dk_it,*l) <= min_value_k){//this might be a new candidate
              const float min_value_l = min_values[*l];
              if(*l != j && *l != i){

                //add potential candidates
                if(d[*l][*k] <= min_value_l){

                  //feasible candidate
                  if(tmpl_tree != NULL){
                    int n1 = conv_to_tmpl[convert_index[*l]];
                    int n2 = conv_to_tmpl[convert_index[*k]];
                    if(n1 >= 0 && n2 >= 0){
                      if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                        mcandidates[*k].lin1 = *k;
                        mcandidates[*k].lin2 = *l;
                        mcandidates[*k].dist = d[*l][*k] + d[*k][*l];
                        mcandidates[*k].dist2 = dist_unif(rng);
                        mcandidates[*k].dist3 = std::max(sample_ages[*k], sample_ages[*l]);
                        mcandidates[*k].is_from_tmpl = true;
                        mcandidates[*l] = mcandidates[*k];
                        mcandidates[*l].lin1 = *l;
                        mcandidates[*l].lin2 = *k;
                      } 
                    }
                  }

                  //sym_dist = d[*l][*k] + d[*k][*l] + std::max(sample_ages[*k], sample_ages[*l]);
                  cand.dist = d[*l][*k] + d[*k][*l];
                  //cand.dist3 = (sample_ages[*k] * cluster_size[*k] + sample_ages[*l] * cluster_size[*l])/(cluster_size[*k] + cluster_size[*l]);
                  cand.dist3 = std::max(sample_ages[*k], sample_ages[*l]);
                  //sym_dist = d[*l][*k] + d[*k][*l];
                  //sym_dist = sample_ages[*k] + sample_ages[*l];
                  cand.dist2 = dist_unif(rng);
                  if( !mcandidates[*k].is_from_tmpl && ((mcandidates[*k].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*k] > cand)){
                    if(cand.dist3 > age){
                      cand.replace = true;
                    }else{
                      cand.replace = false;
                    }
                    mcandidates[*k] = cand;
                    mcandidates[*k].lin1 = *k;
                    mcandidates[*k].lin2 = *l;
                    mcandidates[*k].is_from_tmpl = false;
                  }
                  if( !mcandidates[*l].is_from_tmpl && ((mcandidates[*l].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*l] > cand)){
                    if(cand.dist3 > age){
                      cand.replace = true;  
                    }else{
                      cand.replace = false;
                    }
                    mcandidates[*l] = cand;
                    mcandidates[*l].lin1 = *k;
                    mcandidates[*l].lin2 = *l;
                    mcandidates[*l].is_from_tmpl = false;
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

                //feasible candidate
                if(tmpl_tree != NULL){
                  int n1 = conv_to_tmpl[convert_index[*l]];
                  int n2 = conv_to_tmpl[convert_index[*k]];
                  if(n1 >= 0 && n2 >= 0){
                    if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                      mcandidates[*k].lin1 = *k;
                      mcandidates[*k].lin2 = *l;
                      mcandidates[*k].dist = d[*l][*k] + d[*k][*l];
                      mcandidates[*k].dist2 = dist_unif(rng);
                      mcandidates[*k].dist3 = std::max(sample_ages[*k], sample_ages[*l]);
                      mcandidates[*k].is_from_tmpl = true;
                      mcandidates[*l] = mcandidates[*k];
                      mcandidates[*l].lin1 = *l;
                      mcandidates[*l].lin2 = *k;
                    } 
                  }
                }


                //sym_dist = d[*l][*k] + d[*k][*l] + std::max(sample_ages[*l], sample_ages[*k]);
                cand.dist = d[*l][*k] + d[*k][*l];
                //cand.dist3 = (sample_ages[*l] * cluster_size[*l] + sample_ages[*k] * cluster_size[*k])/(cluster_size[*l] + cluster_size[*k]);
                cand.dist3 = std::max(sample_ages[*k], sample_ages[*l]);
                //sym_dist = d[*l][*k] + d[*k][*l]; 
                //sym_dist = sample_ages[*l] + sample_ages[*k]; 
                cand.dist2 = dist_unif(rng);
                if( !mcandidates[*l].is_from_tmpl && ( (mcandidates[*l].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*l] > cand)){
                  if(cand.dist3 > age){
                    cand.replace = true;     
                  }else{
                    cand.replace = false;
                  }
                  mcandidates[*l] = cand;
                  mcandidates[*l].lin1 = *k;
                  mcandidates[*l].lin2 = *l;
                  mcandidates[*l].is_from_tmpl = false;
                }
                if( !mcandidates[*k].is_from_tmpl && ((mcandidates[*k].dist == std::numeric_limits<float>::infinity() ||cand.dist3 <= age) && mcandidates[*k] > cand)){
                  if(cand.dist3 > age){
                    cand.replace = true; 
                  }else{
                    cand.replace = false;
                  }
                  mcandidates[*k] = cand;
                  mcandidates[*k].lin1 = *k;
                  mcandidates[*k].lin2 = *l;
                  mcandidates[*k].is_from_tmpl = false;
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

              //feasible candidate
              if(tmpl_tree != NULL){
                int n1 = conv_to_tmpl[convert_index[*l]];
                int n2 = conv_to_tmpl[convert_index[*k]];
                if(n1 >= 0 && n2 >= 0){
                  if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                    mcandidates[*k].lin1 = *k;
                    mcandidates[*k].lin2 = *l;
                    mcandidates[*k].dist = d[*l][*k] + d[*k][*l];
                    mcandidates[*k].dist2 = dist_unif(rng);
                    mcandidates[*k].dist3 = std::max(sample_ages[*k], sample_ages[*l]);
                    mcandidates[*k].is_from_tmpl = true;
                    mcandidates[*l] = mcandidates[*k];
                    mcandidates[*l].lin1 = *l;
                    mcandidates[*l].lin2 = *k;
                  } 
                }
              }


              //sym_dist = d[*l][*k] + d[*k][*l] + std::max(sample_ages[*l], sample_ages[*k]);
              cand.dist = d[*l][*k] + d[*k][*l];
              //cand.dist3 = (sample_ages[*l] * cluster_size[*l] + sample_ages[*k] * cluster_size[*k])/(cluster_size[*k] + cluster_size[*l]);
              cand.dist3 = std::max(sample_ages[*l], sample_ages[*k]);
              //sym_dist = d[*l][*k] + d[*k][*l];
              //sym_dist = sample_ages[*l] + sample_ages[*k]; 
              cand.dist2 = dist_unif(rng);
              if( !mcandidates[*l].is_from_tmpl && ((mcandidates[*l].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*l] > cand)){
                if(cand.dist3 > age){
                  cand.replace = true; 
                }else{
                  cand.replace = false;
                }
                mcandidates[*l] = cand;
                mcandidates[*l].lin1 = *k;
                mcandidates[*l].lin2 = *l;
                mcandidates[*l].is_from_tmpl = false;
              }
              if( !mcandidates[*k].is_from_tmpl && ((mcandidates[*k].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*k] > cand)){
                if(cand.dist3 > age){
                  cand.replace = true; 
                }else{
                  cand.replace = false;
                }
                mcandidates[*k] = cand;
                mcandidates[*k].lin1 = *k;
                mcandidates[*k].lin2 = *l;
                mcandidates[*k].is_from_tmpl = false;
              }
            }

          }
        }

      }

      //I might have updated mcandidates[*l] for *l < *k, but if it was the absolute minimum, it has also updated mcandidates[*k]
      if((best_candidate.dist == std::numeric_limits<float>::infinity() || mcandidates[*k].dist3 <= age) && best_candidate > mcandidates[*k]){ 
        best_candidate   = mcandidates[*k];
        if(best_candidate.dist3 > age){
          best_candidate.replace = true;
        }else{
          best_candidate.replace = false;
        }
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
  mcandidates[j].dist2 = std::numeric_limits<float>::infinity();
  mcandidates[j].dist3 = std::numeric_limits<float>::infinity();
  mcandidates[j].replace = false;
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){ 
    if(*std::next(dj_it,*k) <= min_value_j){ 
      if(d[*k][j] <= min_values[*k]){
        if(*k != i && *k != j){

          //feasible candidate
          if(tmpl_tree != NULL){
            int n1 = conv_to_tmpl[convert_index[j]];
            int n2 = conv_to_tmpl[convert_index[*k]];
            if(n1 >= 0 && n2 >= 0){
              if( (*tmpl_tree).nodes[n1].parent == (*tmpl_tree).nodes[n2].parent ){
                mcandidates[*k].lin1 = *k;
                mcandidates[*k].lin2 = j;
                mcandidates[*k].dist = d[j][*k] + d[*k][j];
                mcandidates[*k].dist2 = dist_unif(rng);
                mcandidates[*k].dist3 = std::max(sample_ages[*k], sample_ages[j]);
                mcandidates[*k].is_from_tmpl = true;
                mcandidates[j] = mcandidates[*k];
                mcandidates[j].lin1 = j;
                mcandidates[j].lin2 = *k;
              } 
            }
          }

          //sym_dist = d[j][*k] + d[*k][j] + std::max(sample_ages[j], sample_ages[*k]);
          cand.dist = d[j][*k] + d[*k][j];
          //cand.dist3 = (sample_ages[j] + cluster_size[j] + sample_ages[*k] * cluster_size[*k])/(cluster_size[j] + cluster_size[*k]);
          cand.dist3 = std::max(sample_ages[j], sample_ages[*k]);

          //sym_dist = d[j][*k] + d[*k][j];
          //sym_dist = sample_ages[j] + sample_ages[*k];
          cand.dist2 = dist_unif(rng);
          if( !mcandidates[*k].is_from_tmpl && ((mcandidates[*k].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*k] > cand)){
            if(cand.dist3 > age){
              cand.replace = true;
            }else{
              cand.replace = false;
            }
            mcandidates[*k] = cand;
            mcandidates[*k].lin1 = *k;
            mcandidates[*k].lin2 = j;
            mcandidates[*k].is_from_tmpl = false;
          }
          if( !mcandidates[j].is_from_tmpl && ((mcandidates[j].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[j] > cand)){
            if(cand.dist3 > age){
              cand.replace = true; 
            }else{
              cand.replace = false;
            }
            mcandidates[j] = cand;
            mcandidates[j].lin1 = *k;
            mcandidates[j].lin2 = j;
            mcandidates[j].is_from_tmpl = false;
          }

        }
      }
    }
  }

  if( (best_candidate.dist == std::numeric_limits<float>::infinity() || mcandidates[j].dist3 <= age) && best_candidate > mcandidates[j]){
    best_candidate   = mcandidates[j];
    if(best_candidate.dist3 > age){
      best_candidate.replace = true;
    }else{
      best_candidate.replace = false;
    }
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
MinMatch::QuickBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& i_sample_ages, Tree *tmpl_tree){

  //store pairs that I coalesced. Try these pairs. 
  //for each tip i, store order of coalescences with others
  //mcandidates is a vector of lengths num_remaining lineages. 
  //I store another vector of size n that converts current lineages to tmpl lineages
  //While building tree, if it's feasible I check if I chose the event in the prev tree, and if so assign it to mcandidates

  if(tmpl_tree != NULL){
    conv_to_tmpl.resize(N_total);
    for(int i = 0; i < N; i++){
      conv_to_tmpl[i] = i;
    }
    for(int i = N; i < N_total; i++){
      conv_to_tmpl[i] = -1;
    }
  }

  rng.seed(1);
  std::uniform_real_distribution<double> dist_unif(0,1);

  std::vector<double> sample_ages = i_sample_ages;

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
  //std::random_shuffle( cluster_index.begin(), cluster_index.end() );

  std::fill(min_values.begin(), min_values.end(), std::numeric_limits<float>::infinity());  //Stores the min values of each row of the distance matrix
  std::fill(min_values_sym.begin(), min_values_sym.end(), std::numeric_limits<float>::infinity());  //Stores the min values of each row of the distance matrix

  best_candidate.dist = std::numeric_limits<float>::infinity();
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
  best_candidate.dist3 = std::numeric_limits<float>::infinity(); 
  best_candidate.replace = false;
  best_candidate.is_from_tmpl = false;
  best_sym_candidate.dist = std::numeric_limits<float>::infinity();

  if(sample_ages.size() == N){

    if(unique_sample_ages.size() == 0){
      std::vector<double> foo = sample_ages;
      std::sort(foo.begin(), foo.end());
      //get vector with unique sample ages and counts
      age = foo[0];
      int i = 0;

      unique_sample_ages.resize(foo.size());
      sample_ages_count.resize(foo.size());
      unique_sample_ages[0] = age;
      sample_ages_count[0]  = 0;
      for(std::vector<double>::iterator it_sam = foo.begin(); it_sam != foo.end(); it_sam++){
        if(*it_sam == age){
          sample_ages_count[i]++;
        }else{
          age = *it_sam;
          i++;
          unique_sample_ages[i] = age;
          sample_ages_count[i]++;
        }
      }
      i++;
      unique_sample_ages.resize(i);
      sample_ages_count.resize(i);
      //for(int k = 0; k < unique_sample_ages.size(); k++){
      //  std::cerr << unique_sample_ages[k] << " " << sample_ages_count[k] << std::endl;
      //}
    }
    int level    = 0;
    int num_lins = sample_ages_count[level];
    age      = unique_sample_ages[level] + 2.0/((double) num_lins * (num_lins - 1.0)) * Ne;

    Initialize(d, dist_unif, sample_ages, tmpl_tree);
    //////////////////////////

    //Now I need to coalesce clusters until there is only one lance cluster left
    int no_candidate = 0;
    int updated_cluster_size = 0;
    bool use_sym = false;
    for(num_nodes = N; num_nodes < N_total; num_nodes++){

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

      if(best_candidate.is_from_tmpl){
        assert(conv_to_tmpl[conv_i] >= 0);
        conv_to_tmpl[num_nodes] = (*(*tmpl_tree).nodes[conv_to_tmpl[conv_i]].parent).label;
      }
      convert_index[j] = num_nodes; //update index of this cluster to new merged one

      Coalesce(i, j, d, dist_unif, sample_ages, tmpl_tree);
      if(use_sym) CoalesceSym(i, j, sym_d);

      //sample_ages[j]   = (sample_ages[i] * cluster_size[i] + sample_ages[j] * cluster_size[j])/(cluster_size[i] + cluster_size[j]);
      sample_ages[j] = std::max(sample_ages[i], sample_ages[j]);
      //std::cerr << num_nodes << " " << age << " " << Ne << " " << sample_ages[j] << " " << num_lins << " " << level << std::endl;
      num_lins--;
      if(unique_sample_ages[level] < sample_ages[j]){
        while(unique_sample_ages[level] < sample_ages[j]){
          level++;
          num_lins += sample_ages_count[level];
        }
      }

      age           += 2.0/((double) num_lins * (num_lins - 1.0)) * Ne;
      //if(age < sample_ages[j]){
      //  age = sample_ages[j];
      //}
      assert(num_lins >= 1); 

      cluster_size[j]  = cluster_size[i] + cluster_size[j]; //update size of new cluster 
      //delete cluster i    
      for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){ //using deques instead of lists, which makes this loop slower but iteration faster
        if(*it == i){
          cluster_index.erase(it); //invalidates iterators
          break;
        }
      }

    }

  }else{

    Initialize(d, dist_unif, tmpl_tree);

    //////////////////////////

    //Now I need to coalesce clusters until there is only one lance cluster left
    int no_candidate = 0;
    int updated_cluster_size = 0;
    bool use_sym = false;
    for(num_nodes = N; num_nodes < N_total; num_nodes++){

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

      if(best_candidate.is_from_tmpl){
        assert(conv_to_tmpl[conv_i] >= 0);
        conv_to_tmpl[num_nodes] = (*(*tmpl_tree).nodes[conv_to_tmpl[conv_i]].parent).label;
      }

      convert_index[j] = num_nodes; //update index of this cluster to new merged one
      Coalesce(i, j, d, dist_unif, tmpl_tree);
      if(use_sym) CoalesceSym(i, j, sym_d);

      cluster_size[j]  = cluster_size[i] + cluster_size[j]; //update size of new cluster
      //delete cluster i    
      for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){ //using deques instead of lists, which makes this loop slower but iteration faster
        if(*it == i){
          cluster_index.erase(it); //invalidates iterators
          break;
        }
      }

    } 

  }

  //std::cerr << "Warning: no candidates " << no_candidate << std::endl;
  //if(no_candidate > N/5.0) std::cerr << "Warning: no candidates " << no_candidate << std::endl;
  //if(no_candidate > 0) std::cerr << "Warning: no candidates " << no_candidate << std::endl;

}

void
MinMatch::TestBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& i_sample_ages){

  int root = N_total-1;
  tree.nodes.resize(N_total);
  tree.nodes[root].label  = root;

  std::vector<double> sample_ages = i_sample_ages;
  if(sample_ages.size() != N){
    sample_ages.resize(N);
    std::fill(sample_ages.begin(), sample_ages.end(), 0.0);
  }

  assert(d.size() > 0);
  assert(d.subVectorSize(0) == d.size());

  rng.seed(1);
  std::uniform_real_distribution<double> dist_unif(0,1);


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
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
  best_sym_candidate.dist = std::numeric_limits<float>::infinity();

  Initialize(d, dist_unif, sample_ages);




}

void
MinMatch::SlowBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& i_sample_ages){

  rng.seed(1);
  std::uniform_real_distribution<double> dist_unif(0,1);

  std::vector<double> sample_ages = i_sample_ages;
  if(sample_ages.size() != N){
    sample_ages.resize(N);
    std::fill(sample_ages.begin(), sample_ages.end(), 0.0);
  }

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
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
  best_sym_candidate.dist = std::numeric_limits<float>::infinity();

  Initialize(d, dist_unif, sample_ages);

  //////////////////////////

  //Now I need to coalesce clusters until there is only one lance cluster left
  int no_candidate = 0;
  int updated_cluster_size = 0;
  bool use_sym = false;
  for(num_nodes = N; num_nodes < N_total; num_nodes++){

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


    //coalesce in CF distance matrix
    float added_cluster_size = cluster_size[i] + cluster_size[j];
    dj_it = d.rowbegin(j);
    di_it = d.rowbegin(i);

    if(0){
    float max = 0;
    for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
      if(i != *k){
        d[i][*k] -= d[i][j];
        if(d[i][*k] < 0) d[i][*k] = 0;
        if(d[i][*k] > max) max = d[i][*k];
      }
    }
    float max2 = 0;
    for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
      if(j != *k){
        d[j][*k] -= d[j][i];
        if(d[j][*k] < 0) d[j][*k] = 0;
        if(d[j][*k] > max) max = d[j][*k];
      }
    }
    if(max2 > 0){
      for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
        d[j][*k] *= max/max2;
      }
    }
    }

    for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
      if(j != *k && i != *k){
        dk_it = d.rowbegin(*k);
        dkj   = *std::next(dk_it,j);
        dki   = *std::next(dk_it,i);
        dik   = *std::next(di_it,*k);
        djk   = *std::next(dj_it,*k);

        if(dik != djk){ //if dik == djk, the distance of j to k does not change
          *std::next(dj_it,*k) = (cluster_size[i] * dik + cluster_size[j] * djk)/added_cluster_size;
        }
        if(dki != dkj){ //if dki == dkj, the distance of k to j does not change
          *std::next(dk_it,j)  = (cluster_size[i] * dki + cluster_size[j] * dkj)/added_cluster_size;
        }

      }
    }

    std::fill(min_values.begin(), min_values.end(), std::numeric_limits<float>::infinity());
    for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){
      for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
        if(min_values[*it] > d[*it][*l] && *l != *it && *l != i){
          min_values[*it] = d[*it][*l];
        }
      }
      min_values[*it] += threshold;
    }

    best_candidate.dist = std::numeric_limits<float>::infinity();
    best_candidate.dist2 = std::numeric_limits<float>::infinity();
    for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){
      if(*it != i){
        for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
          if(min_values[*it] >= d[*it][*l] && *l != *it && *l != i){
            if(min_values[*l] >= d[*l][*it]){
              //double sym_dist = d[*l][*it] + d[*it][*l];
              sym_dist = std::max(sample_ages[*it], sample_ages[*l]) + d[*l][*it] + d[*it][*l];
              //sym_dist = d[*l][*it] + d[*it][*l];
              dist_random = dist_unif(rng);
              if( best_candidate.dist > sym_dist || (best_candidate.dist == sym_dist && dist_random < best_candidate.dist2) ){
                best_candidate.lin1 = *it;
                best_candidate.lin2 = *l; 
                best_candidate.dist = sym_dist;
                best_candidate.dist2 = dist_random;
              }
            }	
          }
        }
      }
    }

    //Coalesce(i, j, d, dist_unif, sample_ages);
    if(use_sym) CoalesceSym(i, j, sym_d);

    //if(best_candidate.dist == std::numeric_limits<float>::infinity() && no_candidate == 0) no_candidate = num_nodes;
    //if(best_candidate.dist == std::numeric_limits<float>::infinity() && cluster_index.size() > 2) no_candidate++;

    cluster_size[j]  = cluster_size[i] + cluster_size[j]; //update size of new cluster
    convert_index[j] = num_nodes; //update index of this cluster to new merged one
    sample_ages[j]   = (sample_ages[i] + sample_ages[j])/2.0;
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

  InitializeSym(sym_d, d);

  //////////////////////////

  //Now I need to coalesce clusters until there is only one lance cluster left
  int no_candidate = 0;
  int updated_cluster_size = 0;
  bool use_sym = true;
  for(num_nodes = N; num_nodes < N_total; num_nodes++){

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

//bacRelate functions

void
MinMatch::Initialize(CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif){

  //I have to go thourgh every pair of clusters and check if they are matching mins.
  //I am also filling in the vector min_values.
  it_min_values_it = min_values.begin(); //iterator for min_values[i]

  std::deque<int>::iterator jt;  //iterator for j
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    mcandidates[*it].dist = std::numeric_limits<float>::infinity();
    mcandidates[*it].dist2 = std::numeric_limits<float>::infinity();
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

  it_min_values_it = min_values_CF.begin();
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    d_it = d_CF.rowbegin(*it);
    for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
      if(*it_min_values_it > *d_it && *l != *it){
        *it_min_values_it = *d_it;
      }
      d_it++;
    }
    *it_min_values_it += threshold_CF;
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

          //sym_dist = *it_d_it_jt + d[*jt][*it];
          sym_dist = 1 - (d_CF[*it][*jt] <= min_values_CF[*it]) * (d_CF[*jt][*it] <= min_values_CF[*jt]);
          if(sym_dist > 0){
            sym_dist = *it_d_it_jt + d[*jt][*it];
          }

          dist_random = dist_unif(rng);
          if(mcandidates[*it].dist > sym_dist || (mcandidates[*it].dist == sym_dist && mcandidates[*it].dist2 > dist_random)){
            mcandidates[*it].lin1 = *it;
            mcandidates[*it].lin2 = *jt;
            mcandidates[*it].dist = sym_dist;
            mcandidates[*it].dist2 = dist_random;
          }
          if(mcandidates[*jt].dist > sym_dist || (mcandidates[*jt].dist == sym_dist && mcandidates[*jt].dist2 > dist_random)){
            mcandidates[*jt].lin1 = *it;
            mcandidates[*jt].lin2 = *jt;
            mcandidates[*jt].dist = sym_dist;
            mcandidates[*jt].dist2 = dist_random;
          }
          if(best_candidate.dist > mcandidates[*jt].dist || (best_candidate.dist == mcandidates[*jt].dist && best_candidate.dist2 > mcandidates[*jt].dist2)){ //only need to check for *jt because if it is the absolute minimum, it would also be *it's candidate
            best_candidate.lin1 = *it;
            best_candidate.lin2 = *jt;
            best_candidate.dist = sym_dist;
            best_candidate.dist2 = mcandidates[*jt].dist2;
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
MinMatch::Initialize(CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages){

  //I have to go thourgh every pair of clusters and check if they are matching mins.
  //I am also filling in the vector min_values.
  it_min_values_it = min_values.begin(); //iterator for min_values[i]

  std::deque<int>::iterator jt;  //iterator for j
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    mcandidates[*it].dist = std::numeric_limits<float>::infinity();
    mcandidates[*it].dist2 = std::numeric_limits<float>::infinity();
    mcandidates[*it].dist3 = std::numeric_limits<float>::infinity();
    mcandidates[*it].replace = false;
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

  it_min_values_it = min_values_CF.begin();
  for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){

    d_it = d_CF.rowbegin(*it);
    for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
      if(*it_min_values_it > *d_it && *l != *it){
        *it_min_values_it = *d_it;
      }
      d_it++;
    }
    *it_min_values_it += threshold_CF;
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

          //sym_dist = *it_d_it_jt + d[*jt][*it];
          sym_dist = 1 - (d_CF[*it][*jt] <= min_values_CF[*it]) * (d_CF[*jt][*it] <= min_values_CF[*jt]);
          if(sym_dist == 0){
            sym_dist = *it_d_it_jt + d[*jt][*it];
          }else{
            sym_dist = std::numeric_limits<float>::infinity();
          }
          cand.dist = sym_dist;
          cand.dist3 = std::max(sample_ages[*it], sample_ages[*jt]);
          cand.dist2 = dist_unif(rng);
          if( (mcandidates[*it].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*it] > cand){
            if(cand.dist3 > age){
              cand.replace = true;
            }else{
              cand.replace = false;
            }
            mcandidates[*it] = cand;
            mcandidates[*it].lin1 = *it;
            mcandidates[*it].lin2 = *jt;
          }
          if( (mcandidates[*jt].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*jt] > cand){
            if(cand.dist3 > age){  
              cand.replace = true;
            }else{
              cand.replace = false;
            }
            mcandidates[*jt] = cand;
            mcandidates[*jt].lin1 = *it;
            mcandidates[*jt].lin2 = *jt;
          }
          if( (best_candidate.dist == std::numeric_limits<float>::infinity() || mcandidates[*jt].dist3 <= age) && best_candidate > mcandidates[*jt]){ //only need to check for *jt because if it is the absolute minimum, it would also be *it's candidate
            best_candidate = mcandidates[*jt];
            if(best_candidate.dist3 > age){
              best_candidate.replace = true;
            }else{
              best_candidate.replace = false;
            }
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
MinMatch::Coalesce(const int i, const int j, CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif){

  /////////////
  //now I have to update the distance matrix. I am simultaneously updating min_values[j], where j is the new cluster     
  float added_cluster_size = cluster_size[i] + cluster_size[j];
  float min_value_k, min_value_j = std::numeric_limits<float>::infinity();
  int updated_cluster_size = 0;

  dj_it = d.rowbegin(j);
  di_it = d.rowbegin(i);
  best_candidate.dist = std::numeric_limits<float>::infinity();
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
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

      if( dkj != dki || djk != dik || mcandidates[*k].lin1 == j || mcandidates[*k].lin2 == j || mcandidates[*k].lin1 == i || mcandidates[*k].lin2 == i){ //If *std::next(dj_it,*k) and *std::next(dk_it,j) have not changed, the next candidate for *k will be unchanged

        if(min_value_changed || mcandidates[*k].lin1 == j || mcandidates[*k].lin2 == j || mcandidates[*k].lin1 == i || mcandidates[*k].lin2 == i){

          updated_cluster[updated_cluster_size] = *k; //This is to keep track of which candidates have been changed in current iteration. Needed when I decide that candidates for *k don't change (but candidates of *l might change and could be *k) 
          updated_cluster_size++;

          //count2++;
          mcandidates[*k].dist = std::numeric_limits<float>::infinity();
          mcandidates[*k].dist2 = std::numeric_limits<float>::infinity();

          for(std::deque<int>::iterator l = cluster_index.begin(); l != k; l++){ //only need *l < *k to avoid going through the same pair twice
            if(*std::next(dk_it,*l) <= min_value_k){//this might be a new candidate
              const float min_value_l = min_values[*l];
              if(*l != j && *l != i){

                //add potential candidates
                if(d[*l][*k] <= min_value_l){ 
                  //sym_dist = d[*l][*k] + d[*k][*l];
                  //sym_dist = d_CF[*l][*k] + d_CF[*k][*l];
                  sym_dist = 1 - (d_CF[*l][*k] <= min_values_CF[*l]) * (d_CF[*k][*l] <= min_values_CF[*k]);
                  if(sym_dist > 0){
                    sym_dist = d[*l][*k] + d[*k][*l];
                  }

                  dist_random = dist_unif(rng);
                  if(mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random)){
                    mcandidates[*k].lin1 = *k;
                    mcandidates[*k].lin2 = *l;
                    mcandidates[*k].dist = sym_dist;
                    mcandidates[*k].dist2 = dist_random;
                  }
                  if(mcandidates[*l].dist > sym_dist || (mcandidates[*l].dist == sym_dist && mcandidates[*l].dist2 > dist_random)){
                    mcandidates[*l].lin1 = *k;
                    mcandidates[*l].lin2 = *l;
                    mcandidates[*l].dist = sym_dist;
                    mcandidates[*l].dist2 = dist_random;
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
                //sym_dist = d[*l][*k] + d[*k][*l];
                //sym_dist = d_CF[*l][*k] + d_CF[*k][*l];
                sym_dist = 1 - (d_CF[*l][*k] <= min_values_CF[*l]) * (d_CF[*k][*l] <= min_values_CF[*k]);
                if(sym_dist > 0){
                  sym_dist = d[*l][*k] + d[*k][*l];
                }
                dist_random = dist_unif(rng);
                if(mcandidates[*l].dist > sym_dist || (mcandidates[*l].dist == sym_dist && mcandidates[*l].dist2 > dist_random)){
                  mcandidates[*l].lin1 = *k;
                  mcandidates[*l].lin2 = *l;
                  mcandidates[*l].dist = sym_dist;
                  mcandidates[*l].dist2 = dist_random;
                }
                if(mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random)){
                  mcandidates[*k].lin1 = *k;
                  mcandidates[*k].lin2 = *l;
                  mcandidates[*k].dist = sym_dist;
                  mcandidates[*k].dist2 = dist_random;
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
              //sym_dist = d[*l][*k] + d[*k][*l];
              //sym_dist = d_CF[*l][*k] + d_CF[*k][*l];
              sym_dist = 1 - (d_CF[*l][*k] <= min_values_CF[*l]) * (d_CF[*k][*l] <= min_values_CF[*k]);
              if(sym_dist > 0){
                sym_dist = d[*l][*k] + d[*k][*l];
              }
              dist_random = dist_unif(rng);
              if(mcandidates[*l].dist > sym_dist || (mcandidates[*l].dist == sym_dist && mcandidates[*l].dist2 > dist_random)){
                mcandidates[*l].lin1 = *k;
                mcandidates[*l].lin2 = *l;
                mcandidates[*l].dist = sym_dist;
                mcandidates[*l].dist2 = dist_random;
              }
              if(mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random)){
                mcandidates[*k].lin1 = *k;
                mcandidates[*k].lin2 = *l;
                mcandidates[*k].dist = sym_dist;
                mcandidates[*k].dist2 = dist_random;
              }
            }

          }
        }

      }

      //I might have updated mcandidates[*l] for *l < *k, but if it was the absolute minimum, it has also updated mcandidates[*k]
      if(best_candidate.dist > mcandidates[*k].dist || (best_candidate.dist == mcandidates[*k].dist && best_candidate.dist2 > mcandidates[*k].dist2)){ 
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
  mcandidates[j].dist2 = std::numeric_limits<float>::infinity();
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){ 
    if(*std::next(dj_it,*k) <= min_value_j){ 
      if(d[*k][j] <= min_values[*k]){
        if(*k != i && *k != j){

          //sym_dist = d[j][*k] + d[*k][j];
          //sym_dist = d_CF[j][*k] + d_CF[*k][j];
          sym_dist = 1 - (d_CF[j][*k] <= min_values_CF[j]) * (d_CF[*k][j] <= min_values_CF[*k]);
          if(sym_dist > 0){
            sym_dist = d[j][*k] + d[*k][j];
          }
          dist_random = dist_unif(rng);
          if(mcandidates[*k].dist > sym_dist || (mcandidates[*k].dist == sym_dist && mcandidates[*k].dist2 > dist_random)){
            mcandidates[*k].lin1 = *k;
            mcandidates[*k].lin2 = j;
            mcandidates[*k].dist = sym_dist;
            mcandidates[*k].dist2 = dist_random;
          }
          if(mcandidates[j].dist > sym_dist || (mcandidates[j].dist == sym_dist && mcandidates[j].dist2 > dist_random)){
            mcandidates[j].lin1 = *k;
            mcandidates[j].lin2 = j;
            mcandidates[j].dist = sym_dist;
            mcandidates[j].dist2 = dist_random;
          }

        }
      }
    }
  }

  if(best_candidate.dist > mcandidates[j].dist || (best_candidate.dist == mcandidates[j].dist && best_candidate.dist2 > mcandidates[j].dist2)){
    best_candidate   = mcandidates[j];
  }

}

void
MinMatch::Coalesce(const int i, const int j, CollapsedMatrix<float>& d, CollapsedMatrix<float>& d_CF, std::uniform_real_distribution<double>& dist_unif, std::vector<double>& sample_ages){

  /////////////
  //now I have to update the distance matrix. I am simultaneously updating min_values[j], where j is the new cluster     
  float added_cluster_size = cluster_size[i] + cluster_size[j];
  float min_value_k, min_value_j = std::numeric_limits<float>::infinity();
  int updated_cluster_size = 0;

  dj_it = d.rowbegin(j);
  di_it = d.rowbegin(i);
  best_candidate.dist = std::numeric_limits<float>::infinity();
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
  best_candidate.dist3 = std::numeric_limits<float>::infinity();
  best_candidate.replace = false;
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
    if(j != *k && i != *k){
      dk_it = d.rowbegin(*k);
      dkj   = *std::next(dk_it,j);
      dki   = *std::next(dk_it,i);
      dik   = *std::next(di_it,*k);
      djk   = *std::next(dj_it,*k);
      min_value_k = min_values[*k];
      if(mcandidates[*k].dist3 <= age) mcandidates[*k].replace = false;

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
          mcandidates[*k].dist2 = std::numeric_limits<float>::infinity();
          mcandidates[*k].dist3 = std::numeric_limits<float>::infinity();
          mcandidates[*k].replace = false;

          for(std::deque<int>::iterator l = cluster_index.begin(); l != k; l++){ //only need *l < *k to avoid going through the same pair twice
            if(*std::next(dk_it,*l) <= min_value_k){//this might be a new candidate
              const float min_value_l = min_values[*l];
              if(*l != j && *l != i){

                //add potential candidates
                if(d[*l][*k] <= min_value_l){ 
                  sym_dist = 1 - (d_CF[*l][*k] <= min_values_CF[*l]) * (d_CF[*k][*l] <= min_values_CF[*k]);
                  if(sym_dist > 0){
                    sym_dist = d[*l][*k] + d[*k][*l];
                  }

                  cand.dist = sym_dist;
                  cand.dist3 = std::max(sample_ages[*k], sample_ages[*l]);
                  cand.dist2 = dist_unif(rng);
                  if( (mcandidates[*k].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*k] > cand){
                    if(cand.dist3 > age){
                      cand.replace = true;
                    }else{
                      cand.replace = false;
                    }
                    mcandidates[*k] = cand;
                    mcandidates[*k].lin1 = *k;
                    mcandidates[*k].lin2 = *l;
                  }
                  if( (mcandidates[*l].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*l] > cand){
                    if(cand.dist3 > age){
                      cand.replace = true;  
                    }else{
                      cand.replace = false;
                    }
                    mcandidates[*l] = cand;
                    mcandidates[*l].lin1 = *k;
                    mcandidates[*l].lin2 = *l;
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

                sym_dist = 1 - (d_CF[*l][*k] <= min_values_CF[*l]) * (d_CF[*k][*l] <= min_values_CF[*k]);
                if(sym_dist > 0){
                  sym_dist = d[*l][*k] + d[*k][*l];
                }
                cand.dist = sym_dist;
                cand.dist3 = std::max(sample_ages[*k], sample_ages[*l]);
                cand.dist2 = dist_unif(rng);
                if( (mcandidates[*l].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*l] > cand){
                  if(cand.dist3 > age){
                    cand.replace = true;     
                  }else{
                    cand.replace = false;
                  }
                  mcandidates[*l] = cand;
                  mcandidates[*l].lin1 = *k;
                  mcandidates[*l].lin2 = *l;
                }
                if( (mcandidates[*k].dist == std::numeric_limits<float>::infinity() ||cand.dist3 <= age) && mcandidates[*k] > cand){
                  if(cand.dist3 > age){
                    cand.replace = true; 
                  }else{
                    cand.replace = false;
                  }
                  mcandidates[*k] = cand;
                  mcandidates[*k].lin1 = *k;
                  mcandidates[*k].lin2 = *l;
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
              
              sym_dist = 1 - (d_CF[*l][*k] <= min_values_CF[*l]) * (d_CF[*k][*l] <= min_values_CF[*k]);
              if(sym_dist > 0){
                sym_dist = d[*l][*k] + d[*k][*l];
              }
              cand.dist = sym_dist;
              cand.dist3 = std::max(sample_ages[*l], sample_ages[*k]);
              cand.dist2 = dist_unif(rng);
              if( (mcandidates[*l].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*l] > cand){
                if(cand.dist3 > age){
                  cand.replace = true; 
                }else{
                  cand.replace = false;
                }
                mcandidates[*l] = cand;
                mcandidates[*l].lin1 = *k;
                mcandidates[*l].lin2 = *l;
              }
              if( (mcandidates[*k].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*k] > cand){
                if(cand.dist3 > age){
                  cand.replace = true; 
                }else{
                  cand.replace = false;
                }
                mcandidates[*k] = cand;
                mcandidates[*k].lin1 = *k;
                mcandidates[*k].lin2 = *l;
              }
            }

          }
        }

      }

      //I might have updated mcandidates[*l] for *l < *k, but if it was the absolute minimum, it has also updated mcandidates[*k]
      if( (best_candidate.dist == std::numeric_limits<float>::infinity() || mcandidates[*k].dist3 <= age) && best_candidate > mcandidates[*k]){ 
        best_candidate   = mcandidates[*k];
        if(best_candidate.dist3 > age){
          best_candidate.replace = true;
        }else{
          best_candidate.replace = false;
        }
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
  mcandidates[j].dist2 = std::numeric_limits<float>::infinity();
  mcandidates[j].dist3 = std::numeric_limits<float>::infinity();
  mcandidates[j].replace = false;
  for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){ 
    if(*std::next(dj_it,*k) <= min_value_j){ 
      if(d[*k][j] <= min_values[*k]){
        if(*k != i && *k != j){

          sym_dist = 1 - (d_CF[j][*k] <= min_values_CF[j]) * (d_CF[*k][j] <= min_values_CF[*k]);
          if(sym_dist > 0){
            sym_dist = d[j][*k] + d[*k][j];
          }
          cand.dist = sym_dist;
          //cand.dist3 = (sample_ages[j] + cluster_size[j] + sample_ages[*k] * cluster_size[*k])/(cluster_size[j] + cluster_size[*k]);
          cand.dist3 = std::max(sample_ages[j], sample_ages[*k]);

          //sym_dist = d[j][*k] + d[*k][j];
          //sym_dist = sample_ages[j] + sample_ages[*k];
          cand.dist2 = dist_unif(rng);
          if( (mcandidates[*k].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[*k] > cand){
            if(cand.dist3 > age){
              cand.replace = true;
            }else{
              cand.replace = false;
            }
            mcandidates[*k] = cand;
            mcandidates[*k].lin1 = *k;
            mcandidates[*k].lin2 = j;
          }
          if( (mcandidates[j].dist == std::numeric_limits<float>::infinity() || cand.dist3 <= age) && mcandidates[j] > cand){
            if(cand.dist3 > age){
              cand.replace = true; 
            }else{
              cand.replace = false;
            }
            mcandidates[j] = cand;
            mcandidates[j].lin1 = *k;
            mcandidates[j].lin2 = j;
          }

        }
      }
    }
  }

  if( (best_candidate.dist == std::numeric_limits<float>::infinity() || mcandidates[j].dist3 <= age) && best_candidate > mcandidates[j]){
    best_candidate   = mcandidates[j];
    if(best_candidate.dist3 > age){
      best_candidate.replace = true;
    }else{
      best_candidate.replace = false;
    }
  }

}

void
MinMatch::QuickBuild(CollapsedMatrix<float>& d, Tree& tree, std::vector<double>& i_sample_ages, const CollapsedMatrix<float>& d_prior){

  //threshold = 0.0;

  rng.seed(1);
  std::uniform_real_distribution<double> dist_unif(0,1);

  std::vector<double> sample_ages = i_sample_ages;

  CollapsedMatrix<float> d_CF = d_prior;
  if(d_CF.size() != d.size()){
    d_CF = d;
  }

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
  best_candidate.dist2 = std::numeric_limits<float>::infinity();
  best_candidate.dist3 = std::numeric_limits<float>::infinity(); 
  best_sym_candidate.dist = std::numeric_limits<float>::infinity();

  if(sample_ages.size() == N){

    if(unique_sample_ages.size() == 0){
      std::vector<double> foo = sample_ages;
      std::sort(foo.begin(), foo.end());
      //get vector with unique sample ages and counts
      age = foo[0];
      int i = 0;

      unique_sample_ages.resize(foo.size());
      sample_ages_count.resize(foo.size());
      unique_sample_ages[0] = age;
      sample_ages_count[0]  = 0;
      for(std::vector<double>::iterator it_sam = foo.begin(); it_sam != foo.end(); it_sam++){
        if(*it_sam == age){
          sample_ages_count[i]++;
        }else{
          age = *it_sam;
          i++;
          unique_sample_ages[i] = age;
          sample_ages_count[i]++;
        }
      }
      i++;
      unique_sample_ages.resize(i);
      sample_ages_count.resize(i);
      //for(int k = 0; k < unique_sample_ages.size(); k++){
      //  std::cerr << unique_sample_ages[k] << " " << sample_ages_count[k] << std::endl;
      //}
    }
    int level    = 0;
    int num_lins = sample_ages_count[level];
    age      = unique_sample_ages[level];

    Initialize(d, d_CF, dist_unif, sample_ages);

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

      //coalesce in CF distance matrix
      //min_values changes for cluster j and also for any cluster that had chosen cluster i or j
      //however, if I assume that d_CF is consistent with a tree, then cluster i or j is only chosen if j is still the best lineage and min_values is unchanged
      min_values_CF[j] = std::numeric_limits<float>::infinity();

      float added_cluster_size = cluster_size[i] + cluster_size[j];
      dj_it = d_CF.rowbegin(j);
      di_it = d_CF.rowbegin(i);
      for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
        if(j != *k && i != *k){
          dk_it = d_CF.rowbegin(*k);
          dkj   = *std::next(dk_it,j);
          dki   = *std::next(dk_it,i);
          dik   = *std::next(di_it,*k);
          djk   = *std::next(dj_it,*k);

          if(dik != djk){ //if dik == djk, the distance of j to k does not change
            *std::next(dj_it,*k) = (cluster_size[i] * dik + cluster_size[j] * djk)/added_cluster_size;
          }
          if(dki != dkj){ //if dki == dkj, the distance of k to j does not change
            *std::next(dk_it,j)  = (cluster_size[i] * dki + cluster_size[j] * dkj)/added_cluster_size;
          }

          if(min_values_CF[j] > d_CF[j][*k]){
            min_values_CF[j] = d_CF[j][*k];
          }
        }
      }
      min_values_CF[j] += threshold_CF;
   
      Coalesce(i, j, d, d_CF, dist_unif, sample_ages);
      if(use_sym) CoalesceSym(i, j, sym_d);

      sample_ages[j] = std::max(sample_ages[i], sample_ages[j]);
      age           += 2.0/((double) num_lins * (num_lins - 1.0)) * Ne;
      //std::cerr << num_nodes << " " << age << " " << Ne << " " << sample_ages[j] << " " << num_lins << " " << level << std::endl;
      num_lins--;
      if(unique_sample_ages[level] < sample_ages[j]){
        while(unique_sample_ages[level] < sample_ages[j]){
          level++;
          num_lins += sample_ages_count[level];
        }
      }
      assert(num_lins >= 1); 

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

  }else{

    Initialize(d, d_CF, dist_unif);
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

      //coalesce in CF distance matrix
      //min_values changes for cluster j and also for any cluster that had chosen cluster i or j
      //however, if I assume that d_CF is consistent with a tree, then cluster i or j is only chosen if j is still the best lineage and min_values is unchanged
      min_values_CF[j] = std::numeric_limits<float>::infinity();
      float added_cluster_size = cluster_size[i] + cluster_size[j];
      dj_it = d_CF.rowbegin(j);
      di_it = d_CF.rowbegin(i);
      for(std::deque<int>::iterator k = cluster_index.begin(); k != cluster_index.end(); k++){
        if(j != *k && i != *k){
          dk_it = d_CF.rowbegin(*k);
          dkj   = *std::next(dk_it,j);
          dki   = *std::next(dk_it,i);
          dik   = *std::next(di_it,*k);
          djk   = *std::next(dj_it,*k);

          if(dik != djk){ //if dik == djk, the distance of j to k does not change
            *std::next(dj_it,*k) = (cluster_size[i] * dik + cluster_size[j] * djk)/added_cluster_size;
          }
          if(dki != dkj){ //if dki == dkj, the distance of k to j does not change
            *std::next(dk_it,j)  = (cluster_size[i] * dki + cluster_size[j] * dkj)/added_cluster_size;
          }

          if(min_values_CF[j] > d_CF[j][*k]){
            min_values_CF[j] = d_CF[j][*k];
          }
        }
      }
      min_values_CF[j] += threshold_CF;
        
      if(0){  
      //inefficient exhaustive min_values update
      std::fill(min_values_CF.begin(), min_values_CF.end(), std::numeric_limits<float>::infinity());
      for(std::deque<int>::iterator it = cluster_index.begin(); it != cluster_index.end(); it++){
        for(std::deque<int>::iterator l = cluster_index.begin(); l != cluster_index.end(); l++){
          if(min_values_CF[*it] > d_CF[*it][*l] && *l != *it && *l != i){
            min_values_CF[*it] = d_CF[*it][*l];
          }
          d_it++;
        }
      }
      }

      Coalesce(i, j, d, d_CF, dist_unif);
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

  }

  //std::cerr << "Warning: no candidates " << no_candidate << std::endl;
  //if(no_candidate > N/5.0) std::cerr << "Warning: no candidates " << no_candidate << std::endl;
  //if(no_candidate > 0) std::cerr << "Warning: no candidates " << no_candidate << std::endl;
}


