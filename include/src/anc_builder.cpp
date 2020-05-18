#include "anc_builder.hpp"

//////////////////////////////////////////

float 
Correlation::Pearson(const Leaves& set1, const Leaves& set2){

  if(set1.num_leaves == N || set2.num_leaves == N){
    if(set1.num_leaves == set2.num_leaves) return 1;
    return 0;
  }

  float prod = 0.0; 
  std::vector<int>::const_iterator it_set1_member = set1.member.begin();
  std::vector<int>::const_iterator it_set2_member = set2.member.begin();

  const std::vector<int>::const_iterator it_set1_member_end = set1.member.end();
  const std::vector<int>::const_iterator it_set2_member_end = set2.member.end();

  while(it_set1_member != it_set1_member_end && it_set2_member != it_set2_member_end){
    if(*it_set1_member == *it_set2_member){
      prod += 1.0;
      it_set1_member++;
      it_set2_member++;
    }else if(*it_set1_member < *it_set2_member){
      it_set1_member++;
    }else{
      it_set2_member++;
    }
  }

  if(prod == set1.num_leaves && prod == set2.num_leaves) return 1.0;
  float r = prod - set1.num_leaves * (((float) set2.num_leaves)/N_float);
  if(r <= 0.0) return 0.0;

  r /= sqrt( ((((float)set1.num_leaves)/N_float) * (N_float - set1.num_leaves)) * ((((float)set2.num_leaves)/N_float) * (N_float - set2.num_leaves)) );
  assert(!std::isnan(r));

  return r;

}

//////////////////////// DistanceMeasure //////////////////

void 
DistanceMeasure::GetTopologyWithRepaint(const int snp){


  CollapsedMatrix<float> alpha_begin, beta_end;
  float logscale_alpha, logscale_beta;
  int boundarySNP_begin, boundarySNP_end;

  FastPainting painter(*data);

  char filename[1024];
  snprintf(filename, sizeof(char) * 1024, "%s_%i.bin", (*data).name.c_str(), section);
  FILE* pFile = fopen(filename, "rb");
  assert(pFile != NULL);
  for(int n = 0; n < N; n++){
    fread(&section_startpos, sizeof(int), 1, pFile);
    fread(&section_endpos, sizeof(int), 1, pFile);
    alpha_begin.ReadFromFile(pFile, boundarySNP_begin, logscale_alpha);
    beta_end.ReadFromFile(pFile, boundarySNP_end, logscale_beta);

    assert(boundarySNP_begin <= section_startpos);
    //std::cerr << boundarySNP_begin << " " << section_startpos << " " << snp << std::endl;
    assert(boundarySNP_end >= section_endpos);
    //Repaint into top and log
    painter.RePaintSection(*data, top[n], log[n], alpha_begin, beta_end, boundarySNP_begin, boundarySNP_end, logscale_alpha, logscale_beta, n);
  }

  topology  = &top; 
  logscales = &log;

  //Reset v_snp_prev. I need to know which derivedSNP 'snp' is.
  //Iterate backwards until I get to section_startpos. 
  //I am always including one SNP prior to section_startpos, so the number of derived SNPs to section startpos gives me the previous derived SNP.
  //This is because when accessing matrices, it is 0 based.
  std::fill(v_snp_prev.begin(), v_snp_prev.end(), 0);
  if(snp > 0){
    int tsnp = snp;
    while(tsnp >= section_startpos){
      for(int n = 0; n < N; n++){
        if((*data).sequence[tsnp][n] == '1'){
          v_snp_prev[n]++;
          assert(v_snp_prev[n] < (int)top[n].size());
        }
      }
      tsnp--;
    } 
  }

  //get v_rpos_prev[n] to the last snp with derived mutation
  for(int n = 0; n < N; n++){
    int tsnp = snp;
    while((*data).sequence[tsnp][n] != '1' && tsnp > 0) tsnp--;  //tsnp > 0 is correct, because I want the previous SNP before section_startpos if there is no derived mutation after that
    v_rpos_prev[n] = (*data).rpos[tsnp];
    v_rpos_next[n] = v_rpos_prev[n];
  }
  assert(section_startpos <= snp);
  assert(section_endpos >= snp);
  section++;

}

void
DistanceMeasure::GetMatrix(const int snp){

  if(snp > section_endpos){
    GetTopologyWithRepaint(snp); //load topology matrix 
    //GetTopology(snp);
  }
  //create distance matrix
  for(int n = 0; n < N; n++){
    float min = std::numeric_limits<float>::infinity();
    if((*data).sequence[snp][n] == '1' || snp == 0 || snp == L-1){

      //no need to average
      std::vector<float>::iterator it_matrix = matrix.rowbegin(n);
      std::vector<float>::iterator it_top    = (*topology)[n].rowbegin(v_snp_prev[n]);
      float logscale_prev                    = (*logscales)[n][v_snp_prev[n]];

      //loop vectorizes
      for(; it_matrix != matrix.rowend(n);){

        //*it_matrix = (std::log(*it_top) + logscale_prev) * scale;
        *it_matrix = (fast_log(*it_top) + logscale_prev) * scale;
        if(*it_matrix < min) min = *it_matrix;

        assert(!std::isnan(*it_matrix));
        //assert(*it_matrix < std::numeric_limits<float>::infinity());
        it_matrix++;
        it_top++;
      }
      matrix[n][n] = 0.0;

    }else{
      if(v_rpos_next[n] <= v_rpos_prev[n]){
        for(int l = snp; l < L; l++){
          if((*data).sequence[l][n] == '1' || l == L-1){
            v_rpos_next[n] = (*data).rpos[l];
            break;
          }
        }
      }
      double rpos_prev = v_rpos_prev[n], rpos_next = v_rpos_next[n];
      double weight_denom, weight_left, weight_right;

      if(rpos_prev == rpos_next){
        weight_left  = 0.5;
        weight_right = 0.5;
      }else{
        weight_denom = (rpos_next - rpos_prev);
        weight_left  = (rpos_next - (*data).rpos[snp])/weight_denom;
        weight_right = ((*data).rpos[snp] - rpos_prev)/weight_denom;
      }
      assert(weight_left  > -std::numeric_limits<double>::infinity());
      assert(weight_right > -std::numeric_limits<double>::infinity());

      std::vector<float>::iterator it_matrix   = matrix.rowbegin(n);
      std::vector<float>::iterator it_top_prev = (*topology)[n].rowbegin(v_snp_prev[n]);
      std::vector<float>::iterator it_top_next = (*topology)[n].rowbegin(v_snp_prev[n] + 1);
      float logscale_prev = (*logscales)[n][v_snp_prev[n]];
      float logscale_next = (*logscales)[n][v_snp_prev[n] + 1];
      float exp_logscale_prev_next = exp(logscale_prev - logscale_next);
      float exp_logscale_next_prev = exp(logscale_next - logscale_prev);

      for(; it_matrix != matrix.rowend(n);){

        if(logscale_prev <= logscale_next){
          //*it_matrix = ( std::log( weight_left * (*it_top_prev) * exp_logscale_prev_next + weight_right * (*it_top_next) ) + logscale_next ) * scale; 
          *it_matrix = ( fast_log( weight_left * (*it_top_prev) * exp_logscale_prev_next + weight_right * (*it_top_next) ) + logscale_next ) * scale; 
        }else{
          //*it_matrix = ( std::log( weight_left * (*it_top_prev) + weight_right * (*it_top_next) * exp_logscale_next_prev ) + logscale_prev ) * scale;
          *it_matrix = ( fast_log( weight_left * (*it_top_prev) + weight_right * (*it_top_next) * exp_logscale_next_prev ) + logscale_prev ) * scale;
        }

        assert(!std::isnan(*it_matrix));
        if(*it_matrix < min) min = *it_matrix;
        it_matrix++;
        it_top_next++;
        it_top_prev++; 

      }
      matrix[n][n] = 0.0;

    }
    for(int j = 0; j < (*data).N; j++){ 
      if(j != n) matrix[n][j] -= min;
    } 

  }


  /* 
  std::cout << snp << std::endl; 
  for(int i = 0; i < (*data).N; i++){
    for(int j = 0; j < (*data).N; j++){
      std::cout << matrix[i][j] << " " ;
    }
    std::cout << std::endl;
  }
  */
  

}


//////////////////////// AncesTreeBuilder //////////////////////

//////////// Constructor

AncesTreeBuilder::AncesTreeBuilder(Data& data){
  N       = data.N;
  N_total = 2*N-1;
  root    = N_total - 1;
  L       = data.L;

  mutations.Init(data); 

  threshold_brancheq = 0.95;
  thr = (int) (0.03 * N) + 1;
 
}


//////////// Members

void 
AncesTreeBuilder::BuildTopology(const int section, const int section_startpos, const int section_endpos, Data& data, AncesTree& anc, const int seed, const bool ancestral_state){

  /////////////////////////////////////////////
  //Tree Building
  //input:  Data and distance matrix
  //output: AncesTree (tree sequence) 

  rng.seed(seed);  
  std::uniform_real_distribution<double> dist_unif(0,1);

  Leaves sequences_carrying_mutation;
  sequences_carrying_mutation.member.resize(N);

  MinMatch tb(data);
  DistanceMeasure d(data, section); //this will calculate the distance measure. Needed because we only painted derived sites, so need to recover d by averaging entries of topology

  float min_value, min_value_alt;
  int is_mapping, is_mapping_alt;

  //build tree topology for snp = section_startpos
  anc.seq.emplace_back();
  CorrTrees::iterator it_seq = anc.seq.begin();
  d.GetMatrix(section_startpos); //calculate d

  if(!ancestral_state){ 
    //naive implementation for symmetrising matrix
    for(int i_row = 0; i_row < data.N; i_row++){
      for(int i_col = i_row+1; i_col < data.N; i_col++){
        d.matrix[i_row][i_col] = (d.matrix[i_row][i_col] + d.matrix[i_col][i_row])/2.0;
        d.matrix[i_col][i_row] = d.matrix[i_row][i_col];  
      }
    }  
  }

  tb.QuickBuild(d.matrix, (*it_seq).tree); //build tree topology and store in (*it_seq).tree
  (*it_seq).pos = section_startpos; //record position for this tree along the genome
  UpdateBranchSNPbegin((*it_seq).tree, section_startpos);

  sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
  for(int i = 0; i < N; i++){
    if(data.sequence[section_startpos][i] == '1'){
      sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
      sequences_carrying_mutation.num_leaves++;
    }else{
      sequences_carrying_mutation.member[i] = 0;
    }
  }

  mutations.info[section_startpos].tree = 0;
  
  if(ancestral_state){
    is_mapping = MapMutation((*it_seq).tree, sequences_carrying_mutation, section_startpos, min_value); 
  }else{
    is_mapping = MapMutation((*it_seq).tree, sequences_carrying_mutation, dist_unif, section_startpos, min_value);
  } 
  if(is_mapping > 2){
    ForceMapMutation((*it_seq).tree, sequences_carrying_mutation, section_startpos, true);
  } 

  int num_tree = 1;
  //build tree topology for snp > 0
  for(int snp = section_startpos+1; snp <= section_endpos; snp++){

    //Check if mutations on snp falls on current tree
    sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
    for(int i = 0; i < N; i++){
      if(data.sequence[snp][i] == '1'){
        d.v_snp_prev[i]++; //This is a help vector to keep track of the previous site with a mutation for individual i. (Needed for calculating d, has nothing to do with checking if mutation falls on tree)
        d.v_rpos_prev[i] = data.rpos[snp]; //similar to v_snp_prev
        sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
        sequences_carrying_mutation.num_leaves++; //not needed anymore
      }else{
        sequences_carrying_mutation.member[i] = 0;
      }
    }

    mutations.info[snp].tree = num_tree-1;
    //if mutation does not fall on current tree, I need to build new tree
   
    if(ancestral_state){
      is_mapping = MapMutation((*it_seq).tree, sequences_carrying_mutation, snp, min_value);
    }else{ 
      is_mapping = MapMutation((*it_seq).tree, sequences_carrying_mutation, dist_unif, snp, min_value);
    }
    
    if(is_mapping > 1){

      //if SNP was flipped, it is mapped to current tree (branch prev_branch)
      int prev_branch;
      if(is_mapping == 2){
        assert(mutations.info[snp].branch.size() == 1.0);
        assert((*it_seq).tree.nodes[*mutations.info[snp].branch.begin()].num_events >= 1.0);
        prev_branch = *mutations.info[snp].branch.begin();
      }

      anc.seq.emplace_back();
      it_seq++;
      d.GetMatrix(snp); //calculates distance matrix d at snp 

      if(!ancestral_state){ 
        //naive implementation for symmetrising matrix
        for(int i_row = 0; i_row < data.N; i_row++){
          for(int i_col = i_row+1; i_col < data.N; i_col++){
            d.matrix[i_row][i_col] = (d.matrix[i_row][i_col] + d.matrix[i_col][i_row])/2.0;
            d.matrix[i_col][i_row] = d.matrix[i_row][i_col];  
          }
        }  
      }

      tb.QuickBuild(d.matrix, (*it_seq).tree); //uses distance matrix d to build tree
      (*it_seq).pos = snp; //store position

      if(ancestral_state){
        is_mapping_alt = MapMutation((*it_seq).tree, sequences_carrying_mutation, snp, min_value_alt);
      }else{
        is_mapping_alt = MapMutation((*it_seq).tree, sequences_carrying_mutation, dist_unif, snp, min_value_alt);
      } 
      if(is_mapping_alt > 1 && min_value_alt >= min_value){ 
        //mutation not mapping to a unique branch and new tree is worse than old tree
        it_seq--;
        anc.seq.pop_back();
        if(is_mapping > 2) ForceMapMutation((*it_seq).tree, sequences_carrying_mutation, snp, true); //otherwise, it was flipped and is already mapped to branch

      }else{

        //is_mapping == 1 or (is_mapping > 1 and min_value_alt < min_value)
        //in other words: mutation is mapping, or new tree is better than old tree
        if(is_mapping == 2){
          //need to subtract from prev tree
          (*std::prev(it_seq)).tree.nodes[prev_branch].num_events -= 1.0;
        }
        if(is_mapping_alt > 2){
          ForceMapMutation((*it_seq).tree, sequences_carrying_mutation, snp, true);
        }

        mutations.info[snp].tree = num_tree; 
        UpdateBranchSNPend((*std::prev(it_seq,1)).tree, snp);
        UpdateBranchSNPbegin((*it_seq).tree, snp);
        num_tree++;
      }

    }
  }
  UpdateBranchSNPend((*it_seq).tree, section_endpos);
  //UpdateBranchLifespan((*it_seq).tree, DeltaLifespan);

  anc.N = N;
  anc.L = num_tree;

}

void
AncesTreeBuilder::AssociateTrees(std::vector<AncesTree>& v_anc, const std::string& dirname){

  ///////////////////////////////////////// TMRCA Inference /////////////////////////
  //At this point we inferred tree topology.
  //Now we need to infer times to the MRCA (i.e., branch lengths) 

  int v_anc_size = v_anc.size() - 1;

  CorrTrees::iterator it_seq_prev;
  CorrTrees::iterator it_seq; 

  std::vector<std::vector<int>> equivalent_branches;
  std::vector<std::vector<int>>::iterator it_equivalent_branches;
  std::vector<std::vector<int>>::reverse_iterator rit_equivalent_branches;
  
  int size, total_size = 0;
  for(int k = 0; k <= v_anc_size; k++){
    std::string output_filename = dirname + "equivalent_branches_" + std::to_string(k) + ".bin";
    FILE* pf = fopen(output_filename.c_str(), "rb");
    assert(pf != NULL);

    fread(&size, sizeof(int), 1, pf); 
    total_size += size;

    equivalent_branches.resize(total_size);
    for(int i = total_size - size; i < total_size; i++){
      equivalent_branches[i].resize(N_total);
      fread(&equivalent_branches[i][0], sizeof(int), N_total, pf);
    }
    fclose(pf); 
  }


  ///////////////////////////////////////////
  //Now carry over information on branches, starting from the first tree.

  it_equivalent_branches = equivalent_branches.begin();
  std::vector<Node>::iterator it_nodes;

  for(int i = 0; i < v_anc_size; i++){
    it_seq_prev = v_anc[i].seq.begin();
    it_seq      = std::next(it_seq_prev,1); 

    for(; it_seq != v_anc[i].seq.end();){
      it_nodes = (*it_seq).tree.nodes.begin();
      for(std::vector<int>::iterator it = (*it_equivalent_branches).begin(); it != (*it_equivalent_branches).end(); it++){
        if(*it != -1){
          (*it_nodes).num_events += (*it_seq_prev).tree.nodes[*it].num_events;
          (*it_nodes).SNP_begin   = (*it_seq_prev).tree.nodes[*it].SNP_begin;
        }
        it_nodes++;
      }
      it_equivalent_branches++;

      it_seq++;
      it_seq_prev++;
    }
    it_seq = v_anc[i+1].seq.begin();

    it_nodes = (*it_seq).tree.nodes.begin();
    for(std::vector<int>::iterator it = (*it_equivalent_branches).begin(); it != (*it_equivalent_branches).end(); it++){
      if(*it != -1){
        (*it_nodes).num_events += (*it_seq_prev).tree.nodes[*it].num_events;
        (*it_nodes).SNP_begin   = (*it_seq_prev).tree.nodes[*it].SNP_begin;
      }
      it_nodes++;
    }
    it_equivalent_branches++;

  }
  it_seq_prev = v_anc[v_anc_size].seq.begin();
  it_seq      = std::next(it_seq_prev, 1); 
  for(; it_seq != v_anc[v_anc_size].seq.end();){

    it_nodes = (*it_seq).tree.nodes.begin();
    for(std::vector<int>::iterator it = (*it_equivalent_branches).begin(); it != (*it_equivalent_branches).end(); it++){
      if(*it != -1){
        (*it_nodes).num_events += (*it_seq_prev).tree.nodes[*it].num_events;
        (*it_nodes).SNP_begin   = (*it_seq_prev).tree.nodes[*it].SNP_begin;
      }
      it_nodes++;
    }
    it_equivalent_branches++;

    it_seq++;
    it_seq_prev++;
  }

  assert(it_equivalent_branches == equivalent_branches.end());

  ///////////////////////////////////////////
  //Now go from the last tree to the first

  rit_equivalent_branches = equivalent_branches.rbegin();

  CorrTrees::reverse_iterator rit_seq_next;
  CorrTrees::reverse_iterator rit_seq; 
  for(int i = v_anc_size; i > 0; i--){
    rit_seq_next = v_anc[i].seq.rbegin();
    rit_seq      = std::next(rit_seq_next,1); 
    for(; rit_seq != v_anc[i].seq.rend();){
      it_nodes = (*rit_seq_next).tree.nodes.begin();
      for(std::vector<int>::iterator it = (*rit_equivalent_branches).begin(); it != (*rit_equivalent_branches).end(); it++){
        if(*it != -1){
          (*rit_seq).tree.nodes[*it].num_events = (*it_nodes).num_events;
          (*rit_seq).tree.nodes[*it].SNP_end    = (*it_nodes).SNP_end;
        }
        it_nodes++;
      }
      rit_equivalent_branches++;

      rit_seq++;
      rit_seq_next++;
    }
    rit_seq = v_anc[i-1].seq.rbegin();

    it_nodes = (*rit_seq_next).tree.nodes.begin();
    for(std::vector<int>::iterator it = (*rit_equivalent_branches).begin(); it != (*rit_equivalent_branches).end(); it++){
      if(*it != -1){
        (*rit_seq).tree.nodes[*it].num_events = (*it_nodes).num_events;
        (*rit_seq).tree.nodes[*it].SNP_end    = (*it_nodes).SNP_end;
      }
      it_nodes++;
    }
    rit_equivalent_branches++;

  } 
  rit_seq_next = v_anc[0].seq.rbegin();
  rit_seq      = std::next(rit_seq_next, 1); 
  for(; rit_seq != v_anc[0].seq.rend();){

    it_nodes = (*rit_seq_next).tree.nodes.begin();
    for(std::vector<int>::iterator it = (*rit_equivalent_branches).begin(); it != (*rit_equivalent_branches).end(); it++){
      if(*it != -1){
        (*rit_seq).tree.nodes[*it].num_events = (*it_nodes).num_events;
        (*rit_seq).tree.nodes[*it].SNP_end    = (*it_nodes).SNP_end;
      }
      it_nodes++;
    }
    rit_equivalent_branches++;

    rit_seq++;
    rit_seq_next++;
  }

  /* 
  //add two recombination events to every branch
  for(int i = v_anc_size; i >= 0; i--){
    rit_seq = v_anc[i].seq.rbegin();
    for(; rit_seq != v_anc[i].seq.rend(); rit_seq++){
      for(std::vector<Node>::iterator it_nodes = (*rit_seq).tree.nodes.begin(); it_nodes != (*rit_seq).tree.nodes.end(); it_nodes++){
        (*it_nodes).num_events += 2;
      }
    }
  }
  */
  
  assert(rit_equivalent_branches == equivalent_branches.rend());

}

/*
std::pair<int, int>
AncesTreeBuilder::OptimizeParameters(const int section_startpos, const int section_endpos, Data& data){

  int window_length = 400;

  int num_nonmapping_snps = 0;
  int num_nonmapping_rare_snps = 0;

  //Paint
  DistanceMeasure d(data); //this will calculate the distance measure.

  std::vector<CollapsedMatrix<float>> topology;
  std::vector<std::vector<float>> logscales;
  CollapsedMatrix<float> alpha_begin, beta_end;

  topology.resize(data.N);
  logscales.resize(data.N);
  alpha_begin.resize(1,data.N);
  beta_end.resize(1,data.N);
  float logscale_alpha = 0.0, logscale_beta = 0.0;
  int boundarySNP_begin = section_startpos, boundarySNP_end = section_endpos;

  // populate variables that are independent of k 
  float prior_theta  = data.theta/(data.N-1.0) - data.ntheta/(data.N-1.0);
  float prior_ntheta = data.ntheta/(data.N-1.0);
  int derived;
  std::fill(beta_end.vbegin(), beta_end.vend(), 1.0);

  // paint
  FastPainting painter(data);
  float rescale = fast_log(data.theta/(1.0-data.theta));

  for(int k = 0; k < data.N; k++){
    for(int n = 0; n < data.N; n++){
      derived                   = (double) (data.sequence[0][k] > data.sequence[0][n]);
      alpha_begin[0][n]         = derived * prior_theta + prior_ntheta;
    }
    painter.RePaintSection(data, topology[k], logscales[k], alpha_begin, beta_end, boundarySNP_begin, boundarySNP_end, logscale_alpha, logscale_beta, k);

    float normalizing_constant = fast_log(N-1.0) - topology.size() * fast_log(data.ntheta); //topology.size() needs to be derived mutations of section
    for(int l = 0; l < (int) logscales[k].size(); l++){
      logscales[k][l] += normalizing_constant;
    }

  }
  int snp = section_startpos + window_length;
  d.AssignTopology(topology, logscales, section_startpos, section_endpos, snp);

  /////////////////////////////////////////////
  //Tree Building
  //input:  Data and distance matrix
  //output: AncesTree (tree sequence)

  Leaves sequences_carrying_mutation;
  sequences_carrying_mutation.member.resize(data.N);

  mutations.Init(data);
  MinMatch tb(data);

  Tree tree;
  float min_value;
  const float log_theta = log(data.theta), log_ntheta = log(data.ntheta);
  float min;

  sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
  for(int i = 0; i < data.N; i++){
    if(data.sequence[snp][i] == '1'){
      sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
      sequences_carrying_mutation.num_leaves++;
    }else{
      sequences_carrying_mutation.member[i]     = 0;
    }
  }

  if(sequences_carrying_mutation.num_leaves > 1){
    d.GetMatrix(snp); //calculate d
    //modify distance matrix so that current SNP is cancelled
    for(int i = 0; i < data.N; i++){
      if(data.sequence[snp][i] == '1'){
        min = std::numeric_limits<float>::infinity();
        for(int j = 0; j < data.N; j++){
          if(data.sequence[snp][j] == '0') d.matrix[i][j] += log_theta; //adding because d.matrix is multiplied by -1
          if(data.sequence[snp][j] == '1') d.matrix[i][j] += log_ntheta;
          if(min > d.matrix[i][j]) min = d.matrix[i][j];
          assert(d.matrix[i][j] < std::numeric_limits<float>::infinity());
        }
        for(int j = 0; j < data.N; j++){
          d.matrix[i][j] -= min;
        }
      }
    }

    tb.QuickBuild(d.matrix, tree); //build tree topology and store in (*it_seq).tree

    if(MapMutation(tree, sequences_carrying_mutation, snp, min_value) > 1){
      if(sequences_carrying_mutation.num_leaves == 2 || sequences_carrying_mutation.num_leaves == 3) num_nonmapping_rare_snps++;
      num_nonmapping_snps++;
    }
  }

  //build tree topology for snp > start_section
  snp++;
  for(; snp <= section_startpos + 2*window_length; snp++){

    //Check if mutations on snp falls on current tree
    sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
    for(int i = 0; i < data.N; i++){
      if(data.sequence[snp][i] == '1'){
        d.v_snp_prev[i]++; //This is a help vector to keep track of the previous site with a mutation for individual i.
        d.v_rpos_prev[i] = data.rpos[snp]; //similar to v_snp_prev
        sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
        sequences_carrying_mutation.num_leaves++; //not needed anymore
      }else{
        sequences_carrying_mutation.member[i] = 0;
      }
    }

    if(sequences_carrying_mutation.num_leaves > 1){
      d.GetMatrix(snp); //calculate d
      //modify distance matrix so that current SNP is cancelled
      for(int i = 0; i < data.N; i++){
        if(data.sequence[snp][i] == '1'){
          min = std::numeric_limits<float>::infinity();
          for(int j = 0; j < data.N; j++){
            if(data.sequence[snp][j] == '0') d.matrix[i][j] += log_theta;
            if(data.sequence[snp][j] == '1') d.matrix[i][j] += log_ntheta;
            if(min > d.matrix[i][j]) min = d.matrix[i][j];
            assert(d.matrix[i][j] < std::numeric_limits<float>::infinity());
          }
          for(int j = 0; j < data.N; j++){
            d.matrix[i][j] -= min;
          }
        }
      }
      tb.QuickBuild(d.matrix, tree); //build tree topology and store in (*it_seq).tree

      if(MapMutation(tree, sequences_carrying_mutation, snp, min_value) > 1){
        if(sequences_carrying_mutation.num_leaves == 2 || sequences_carrying_mutation.num_leaves == 3) num_nonmapping_rare_snps++;
        num_nonmapping_snps++;
      } 
    }

  }

  return std::make_pair(num_nonmapping_snps, num_nonmapping_rare_snps);

}
*/

////////////////////////////////
//PRIVATE HELPER

//check if a mutations lies on the tree. Requires leaves calculated using FindAllLeaves as input

//////////
int 
AncesTreeBuilder::MapMutation(Tree& tree, Leaves& sequences_carrying_mutations, std::uniform_real_distribution<double>& dist_unif, const int snp, float& min_value){

  if(sequences_carrying_mutations.num_leaves == 0 || sequences_carrying_mutations.num_leaves == N) return 1;

  //I want to place the mutation on all branches necessary for no loss of information
  //start with all leaves
  //propagate up and count number of nodes needed.
  //choose flipped or non-flipped depending on which is less.

  PropagateStructGlobal report;
  PropagateMutationGlobal(tree.nodes[root], sequences_carrying_mutations, report);
 
  if(report.min == report.flipped_min && report.min <= thr){
  
    //bool flag = true; //default flag
    bool flag = (dist_unif(rng) < 0.5); 

    if(flag){// not flipped
      min_value = report.min;
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_branch; 
      mutations.info[snp].flipped = false;
      tree.nodes[report.best_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 1;
    }else{
      min_value = report.flipped_min;
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_flipped_branch;
      mutations.info[snp].flipped = true;
      tree.nodes[report.best_flipped_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 2;
    }

  }else if( report.min <= report.flipped_min ){

    min_value = report.min;
    if( report.min <= thr ){
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_branch; 
      mutations.info[snp].flipped = false;
      tree.nodes[report.best_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 1;
    }
    return 3;

  }else{

    min_value = report.flipped_min;
    if( report.flipped_min <= thr ){
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_flipped_branch;
      mutations.info[snp].flipped = true;
      tree.nodes[report.best_flipped_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 2;
    }
    return 3;

  }

}

//////////
//version without random flipping
int 
AncesTreeBuilder::MapMutation(Tree& tree, Leaves& sequences_carrying_mutations, const int snp, float& min_value){

  if(sequences_carrying_mutations.num_leaves == 0 || sequences_carrying_mutations.num_leaves == N) return 1;

  //I want to place the mutation on all branches necessary for no loss of information
  //start with all leaves
  //propagate up and count number of nodes needed.
  //choose flipped or non-flipped depending on which is less.

  PropagateStructGlobal report;
  PropagateMutationGlobal(tree.nodes[root], sequences_carrying_mutations, report);
 
  if(report.min == report.flipped_min && report.min <= thr){
  
    bool flag = true; //default flag
    if(flag){// not flipped
      min_value = report.min;
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_branch; 
      mutations.info[snp].flipped = false;
      tree.nodes[report.best_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 1;
    }else{
      min_value = report.flipped_min;
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_flipped_branch;
      mutations.info[snp].flipped = true;
      tree.nodes[report.best_flipped_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 2;
    }

  }else if( report.min <= report.flipped_min ){

    min_value = report.min;
    if( report.min <= thr ){
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_branch; 
      mutations.info[snp].flipped = false;
      tree.nodes[report.best_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 1;
    }
    return 3;

  }else{

    min_value = report.flipped_min;
    if( report.flipped_min <= thr ){
      mutations.info[snp].branch.resize(1);
      mutations.info[snp].branch[0] = report.best_flipped_branch;
      mutations.info[snp].flipped = true;
      tree.nodes[report.best_flipped_branch].num_events += 1.0;
      assert(mutations.info[snp].branch.size() == 1);
      return 2;
    }
    return 3;

  }

}


int 
AncesTreeBuilder::ForceMapMutation(Tree& tree, Leaves& sequences_carrying_mutations, const int snp, const bool force){

  if(sequences_carrying_mutations.num_leaves == 0 || sequences_carrying_mutations.num_leaves == N) return 1;

  //I want to place the mutation on all branches necessary for no loss of information
  //start with all leaves
  //propagate up and count number of nodes needed.
  //choose flipped or non-flipped depending on which is less.

  std::deque<int> branches;
  std::deque<int> branches_flipped; 

  PropagateStructLocal report;
  PropagateMutationLocal(tree.nodes[root], branches, branches_flipped, sequences_carrying_mutations, report);

  if(branches_flipped.size() == 0){

    assert(branches.size() > 0);

    if(branches.size() == 1 || force){
      mutations.info[snp].branch  = branches;
      float weight = 1.0/((float) branches.size());
      for(std::deque<int>::iterator it = branches.begin(); it != branches.end(); it++){
        tree.nodes[*it].num_events += weight;
      }
    }

    return branches.size();

  }else{

    if( branches.size() <= branches_flipped.size() && branches.size() > 0 ){ 
      if(branches.size() == 1 || force){
        mutations.info[snp].branch  = branches;
        float weight = 1.0/((float) branches.size());
        for(std::deque<int>::iterator it = branches.begin(); it != branches.end(); it++){
          tree.nodes[*it].num_events += weight;
        }
      }

      return branches.size();
    }else{
      if(branches_flipped.size() == 1 || force){
        mutations.info[snp].flipped = true;
        mutations.info[snp].branch  = branches_flipped;
        float weight = 1.0/((float) branches_flipped.size());
        for(std::deque<int>::iterator it = branches_flipped.begin(); it != branches_flipped.end(); it++){
          tree.nodes[*it].num_events += weight;
        }
      }

      return branches_flipped.size();
    }

  }
}

int
AncesTreeBuilder::PropagateMutationExact(Node& node, std::deque<int>& branches, std::deque<int>& branches_flipped, Leaves& sequences_carrying_mutations){

  if(node.child_left != NULL){

    int p1 = PropagateMutationExact(*node.child_left, branches, branches_flipped, sequences_carrying_mutations);
    int p2 = PropagateMutationExact(*node.child_right, branches, branches_flipped, sequences_carrying_mutations);

    if(p1 == p2) return p1;

    //otherwise one of the child_branches need to be included
    if(p1 == 1) branches.push_back((*node.child_left).label);
    if(p2 == 1) branches.push_back((*node.child_right).label);
    if(p1 == -1) branches_flipped.push_back((*node.child_left).label);
    if(p2 == -1) branches_flipped.push_back((*node.child_right).label);

    return 0;

  }else{

    if(sequences_carrying_mutations.member[node.label] == 1){
      return 1;
    }else{
      return -1;
    }

  }

}

void
AncesTreeBuilder::PropagateMutationGlobal(Node& node, Leaves& sequences_carrying_mutations, PropagateStructGlobal& report){

  float total_carriers    = sequences_carrying_mutations.num_leaves; 
  float total_noncarriers = N - total_carriers;  //slightly inefficient 

  if(node.child_left != NULL){

    PropagateStructGlobal report2;

    PropagateMutationGlobal(*node.child_left, sequences_carrying_mutations, report);
    PropagateMutationGlobal(*node.child_right, sequences_carrying_mutations, report2);

    report.num_correct_carriers      += report2.num_correct_carriers;
    report.num_incorrect_noncarriers += report2.num_incorrect_noncarriers;
    report.num_incorrect_carriers     = total_carriers - report.num_correct_carriers;
    report.num_correct_noncarriers    = total_noncarriers - report.num_incorrect_noncarriers;

    int sum = report.num_incorrect_carriers + report.num_incorrect_noncarriers;

    //std::cerr << report.num_incorrect_carriers/total_carriers << std::endl; 
    bool necessary_condition = (((float) report.num_incorrect_carriers)/total_carriers < 0.3);
    necessary_condition *= (((float) report.num_incorrect_noncarriers)/total_noncarriers < 0.3);
    if(report.num_correct_carriers + report.num_incorrect_noncarriers > 0.0){
      necessary_condition *= (((float) report.num_correct_carriers)/(report.num_correct_carriers + report.num_incorrect_noncarriers) > 0.7);
    }
    if(report.num_incorrect_carriers + report.num_correct_noncarriers > 0.0){
      necessary_condition *= (((float) report.num_correct_noncarriers)/(report.num_incorrect_carriers + report.num_correct_noncarriers) > 0.7);
    }
    if( necessary_condition && report.min > sum && report2.min > sum ){
      report.min         = sum;
      report.best_branch = node.label;
    }else{
      if( report.min > report2.min ){
        report.min         = report2.min;
        report.best_branch = report2.best_branch;
      }//else report is correct
    } 

    sum = report.num_correct_carriers + report.num_correct_noncarriers;
    
    necessary_condition = (((float) report.num_correct_carriers)/total_carriers < 0.3);
    necessary_condition *= (((float) report.num_correct_noncarriers)/total_noncarriers < 0.3);
    if(report.num_incorrect_carriers + report.num_correct_noncarriers > 0.0){
      necessary_condition *= (((float) report.num_incorrect_carriers)/(report.num_incorrect_carriers + report.num_correct_noncarriers) > 0.7);
    }
    if(report.num_correct_carriers + report.num_incorrect_noncarriers > 0.0){
      necessary_condition *= (((float) report.num_incorrect_noncarriers)/(report.num_correct_carriers + report.num_incorrect_noncarriers) > 0.7);
    }
    if( necessary_condition && report.flipped_min > sum && report2.flipped_min > sum ){
      report.flipped_min         = sum;
      report.best_flipped_branch = node.label;
    }else{
      if( report.flipped_min > report2.flipped_min ){
        report.flipped_min         = report2.flipped_min;
        report.best_flipped_branch = report2.best_flipped_branch;
      }//else report is correct
    }

  }else{

    if(sequences_carrying_mutations.member[node.label] == 1){
      report.num_correct_carriers      = 1;
      report.num_incorrect_carriers    = total_carriers - 1;
      report.num_correct_noncarriers   = total_noncarriers;
      report.num_incorrect_noncarriers = 0;

      if( report.num_incorrect_carriers/total_carriers < 0.3 ){
        report.min                     = report.num_incorrect_carriers;
        report.best_branch             = node.label;
      }else{
        report.min                     = std::numeric_limits<int>::max();
        report.best_branch             = -1;
      }
      if( report.num_correct_carriers/total_carriers < 0.3 && report.num_correct_noncarriers/total_noncarriers < 0.3 ){
        report.flipped_min             = report.num_correct_noncarriers + report.num_correct_carriers; 
        report.best_flipped_branch     = node.label;
      }else{
        report.flipped_min             = std::numeric_limits<int>::max();
        report.best_flipped_branch     = -1;
      } 
    }else{
      report.num_correct_carriers      = 0;
      report.num_incorrect_carriers    = total_carriers;
      report.num_correct_noncarriers   = total_noncarriers - 1;
      report.num_incorrect_noncarriers = 1;

      if( report.num_incorrect_carriers/total_carriers < 0.3 && report.num_incorrect_noncarriers/total_noncarriers < 0.3 ){
        report.min                     = report.num_incorrect_carriers + report.num_incorrect_noncarriers;
        report.best_branch             = node.label;
      }else{
        report.min                     = std::numeric_limits<int>::max();
        report.best_branch             = -1;
      }
      if( report.num_correct_noncarriers/total_noncarriers < 0.3 ){
        report.flipped_min             = report.num_correct_noncarriers;
        report.best_flipped_branch     = node.label;
      }else{
        report.flipped_min             = std::numeric_limits<int>::max();
        report.best_flipped_branch     = -1;
      } 
    }

  }

}

void
AncesTreeBuilder::PropagateMutationLocal(Node& node, std::deque<int>& branches, std::deque<int>& branches_flipped, Leaves& sequences_carrying_mutations, PropagateStructLocal& report){

  if(node.child_left != NULL){

    PropagateStructLocal report_c1, report_c2;

    PropagateMutationLocal(*node.child_left, branches, branches_flipped, sequences_carrying_mutations, report_c1);
    PropagateMutationLocal(*node.child_right, branches, branches_flipped, sequences_carrying_mutations, report_c2);

    report.num_carriers = report_c1.num_carriers + report_c2.num_carriers;
    report.num_flipped_carriers = report_c1.num_flipped_carriers + report_c2.num_flipped_carriers; 
    float num_leaves = report.num_carriers + report.num_flipped_carriers;

    if(report.num_flipped_carriers/num_leaves < 0.03 && report_c1.best_branch != -1 && report_c2.best_branch != -1){
      if(report_c1.num_carriers > 0 && report_c2.num_carriers > 0){
        report.best_branch             = node.label;
      }else if(report_c1.num_carriers > 0){
        report.best_branch             = report_c1.best_branch;
      }else if(report_c2.num_carriers > 0){
        report.best_branch             = report_c2.best_branch;
      }else{
        assert(false);
      }
    }else{
      if(report_c1.best_branch != -1){
        branches.push_back(report_c1.best_branch);
      }
      if(report_c2.best_branch != -1){
        branches.push_back(report_c2.best_branch);
      }
      report.best_branch = -1;
    }

    if(report.num_carriers/num_leaves < 0.03 && report_c1.best_flipped_branch != -1 && report_c2.best_flipped_branch != -1){
      if(report_c1.num_flipped_carriers > 0 && report_c2.num_flipped_carriers > 0){
        report.best_flipped_branch             = node.label;
      }else if(report_c1.num_flipped_carriers > 0){
        report.best_flipped_branch             = report_c1.best_flipped_branch;
      }else if(report_c2.num_flipped_carriers > 0){
        report.best_flipped_branch             = report_c2.best_flipped_branch;
      }else{
        assert(false);
      }
    }else{
      if(report_c1.best_flipped_branch != -1){
        branches_flipped.push_back(report_c1.best_flipped_branch);
      }
      if(report_c2.best_flipped_branch != -1){
        branches_flipped.push_back(report_c2.best_flipped_branch);
      }
      report.best_flipped_branch = -1;
    }

  }else{

    if(sequences_carrying_mutations.member[node.label] == 1){
      report.num_carriers         = 1;
      report.num_flipped_carriers = 0;
      report.best_branch          = node.label;
      report.best_flipped_branch  = -1;
    }else{
      report.num_carriers         = 0;
      report.num_flipped_carriers = 1;
      report.best_flipped_branch  = node.label;
      report.best_branch          = -1;
    }

  }

}

void  
AncesTreeBuilder::UpdateBranchSNPbegin(Tree& tree, int snp){
  std::vector<Node>::iterator it_node = tree.nodes.begin();
  for(; it_node != tree.nodes.end(); it_node++){
    (*it_node).SNP_begin = snp;
  }
}

void  
AncesTreeBuilder::UpdateBranchSNPend(Tree& tree, int snp){
  std::vector<Node>::iterator it_node = tree.nodes.begin();
  for(; it_node != tree.nodes.end(); it_node++){
    (*it_node).SNP_end = snp;
  }
}


//Find equivalent branches and carry over branch properties (such as number of events, SNP_begin, SNP_end)
void
AncesTreeBuilder::PreCalcPotentialBranches(){

  //the number of leaves a branch needs to be equivalent
  potential_branches.resize(N);
  float threshold_inv = 1/(threshold_brancheq * threshold_brancheq);
  float N_float = N;
  for(int i = 1; i <= N; i++){

    potential_branches[i-1].push_back(i);
    //for branches with i+1 leaves, list the number of leaves a potential equivalent branch needs
    for(int j = i+1; j <= N; j++){
      if(threshold_inv >= j/(N_float-j) * ((N_float-i)/i) ){
        potential_branches[i-1].push_back(j);
        potential_branches[j-1].push_back(i);
      }
    }
  }

}

void
AncesTreeBuilder::BranchAssociation(Tree& ref_tree, Tree& tree, std::vector<int>& equivalent_branches){

  ////
  equivalent_branches.resize(N_total);
  std::fill(equivalent_branches.begin(), equivalent_branches.end(), -1);
  std::vector<int> equivalent_branches_ref(N_total, -1);

  Correlation cor(N);

  //1. calculate leaves
  std::vector<Leaves> tr_leaves;
  std::vector<Leaves> rtr_leaves;

  tree.FindAllLeaves(tr_leaves);
  ref_tree.FindAllLeaves(rtr_leaves);

  //2. Sort nodes according to number of leaves and store in std::vector of size(N_total)

  std::vector<int> sorted_branches(N_total);
  std::size_t n(0);
  std::generate(std::begin(sorted_branches), std::end(sorted_branches), [&]{ return n++; });
  std::sort(std::begin(sorted_branches), std::end(sorted_branches), [&](int i1, int i2) { return rtr_leaves[i1].num_leaves < rtr_leaves[i2].num_leaves; } );

  std::vector<int> index_sorted_branches(N,0); //this vector is for accessing sorted_branches
  for(std::vector<Leaves>::iterator it = rtr_leaves.begin(); it != std::prev(rtr_leaves.end(),1); it++){
    index_sorted_branches[(*it).num_leaves]++;
  }
  int cum = 0;
  for(std::vector<int>::iterator it = index_sorted_branches.begin(); it != index_sorted_branches.end(); it++){
    *it += cum;
    cum = *it;
  }

  //3. For each branch, search for exact equivalent branch in ref_tree
  //   record which branches are unpaired

  std::vector<int> unpaired_branches;

  Node parent, ref_parent;

  //treat leaves separately (this is faster)
  for(int i = 0; i < N; i++){

    if(equivalent_branches[i] == -1){

      parent = *tree.nodes[i].parent;
      ref_parent = *ref_tree.nodes[i].parent;

      if((*parent.child_left).label == i){

        int child_right = (*parent.child_right).label;
        if(child_right < N){
          if( child_right == (*ref_parent.child_right).label ){ 
            equivalent_branches[i]     = i;
            equivalent_branches_ref[i] = i;
            equivalent_branches[child_right]     = child_right;
            equivalent_branches_ref[child_right] = child_right;
          }else if( child_right == (*ref_parent.child_left).label ){ 
            equivalent_branches[i]     = i;
            equivalent_branches_ref[i] = i;
            equivalent_branches[child_right]     = child_right;
            equivalent_branches_ref[child_right] = child_right;
          } 
        }else{
          if(cor.Pearson(tr_leaves[parent.label], rtr_leaves[ref_parent.label]) >= threshold_brancheq){
            equivalent_branches[i]     = i;
            equivalent_branches_ref[i] = i;
          }
        }

      }else{

        int child_left = (*parent.child_left).label;
        if(child_left < N){
          if( child_left == (*ref_parent.child_left).label ){ 
            equivalent_branches[i]     = i;
            equivalent_branches_ref[i] = i;
            equivalent_branches[child_left]     = child_left;
            equivalent_branches_ref[child_left] = child_left;
          }else if( child_left == (*ref_parent.child_right).label ){ 
            equivalent_branches[i]     = i;
            equivalent_branches_ref[i] = i;
            equivalent_branches[child_left]     = child_left;
            equivalent_branches_ref[child_left] = child_left;
          } 
        }else{
          if(cor.Pearson(tr_leaves[parent.label], rtr_leaves[ref_parent.label]) >= threshold_brancheq){
            equivalent_branches[i]     = i;
            equivalent_branches_ref[i] = i;
          }
        }

      }

    }

  }

  for(int i = N; i < N_total-1; i++){ 

    if(cor.Pearson(tr_leaves[i], rtr_leaves[i]) >= 0.9999 && cor.Pearson(tr_leaves[(*tree.nodes[i].parent).label], rtr_leaves[(*ref_tree.nodes[i].parent).label]) >= 0.9999){     
      //branches i and i are equivalent
      equivalent_branches[i] = i;
      equivalent_branches_ref[i] = i;
    }

    if(equivalent_branches[i] == -1){
      int num_leaves = tr_leaves[i].num_leaves;
      for(std::vector<int>::iterator it = std::next(sorted_branches.begin(),index_sorted_branches[num_leaves-1]); it != std::next(sorted_branches.begin(), index_sorted_branches[num_leaves]); it++ ){
        if(cor.Pearson(tr_leaves[i], rtr_leaves[*it]) >= 0.9999 && cor.Pearson(tr_leaves[(*tree.nodes[i].parent).label], rtr_leaves[(*ref_tree.nodes[*it].parent).label]) >= 0.9999){     
          //branches i and *it are equivalent
          equivalent_branches[i]       = *it;
          equivalent_branches_ref[*it] = i;

          break;  
        }
      }
    }

    if(equivalent_branches[i] == -1){
      unpaired_branches.push_back(i);
    }

  }

  //4. For the nodes that are not paired up, find all possible pairs above threshold 
  std::vector<EquivalentNode> possible_pairs;

  for(std::vector<int>::iterator it_unpaired = unpaired_branches.begin(); it_unpaired != unpaired_branches.end(); it_unpaired++){

    int num_leaves = tr_leaves[*it_unpaired].num_leaves - 1;
    for(std::vector<int>::iterator it_k = potential_branches[num_leaves].begin(); it_k != potential_branches[num_leaves].end(); it_k++){

      for(std::vector<int>::iterator it = std::next(sorted_branches.begin(),index_sorted_branches[*it_k-1]); it != std::next(sorted_branches.begin(), index_sorted_branches[*it_k]); it++ ){
        if(equivalent_branches_ref[*it] == -1){
          float score = cor.Pearson(tr_leaves[*it_unpaired], rtr_leaves[*it]);
          if(score >= threshold_brancheq && cor.Pearson(tr_leaves[(*tree.nodes[*it_unpaired].parent).label], rtr_leaves[(*ref_tree.nodes[*it].parent).label]) >= threshold_brancheq){     
            //branches i and *it are equivalent
            possible_pairs.push_back(EquivalentNode(*it_unpaired, *it, score));
          }
        }
      }

    }

  }


  //5. Sort possible_pairs by score
  std::sort(std::begin(possible_pairs), std::end(possible_pairs), std::greater<EquivalentNode>());

  //6. Pair up approximate matches
  for(std::vector<EquivalentNode>::iterator it = possible_pairs.begin(); it != possible_pairs.end(); it++){
    if(equivalent_branches[(*it).node1] == -1 && equivalent_branches_ref[(*it).node2] == -1){
      equivalent_branches[(*it).node1]     = (*it).node2;
      equivalent_branches_ref[(*it).node2] = (*it).node1;
    }
  }

}


