#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.hpp"
#include "collapsed_matrix.hpp"
#include "plot.hpp"
#include "sample.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

void
GetCoordsAndLineages(MarginalTree& mtr, std::vector<float>& coordinates_tree, std::vector<int>& num_lineages){

  /*
  ///
  std::vector<float> vec1(4, 0.0);
  vec1[0] = 10;
  vec1[1] = 5;
  vec1[2] = 6;
  vec1[3] = 7;
  std::vector<int> vec2(4,0); 

  std::vector<int> ind(4);
  std::size_t m2(0);
  std::generate(std::begin(ind), std::end(ind), [&]{ return m2++; });
  for(int i = 0; i < 4; i++){
  std::cerr << ind[i] << std::endl;
  }
  std::sort(std::begin(ind), std::end(ind), [&](int i1, int i2) {
  return std::tie(vec1[i1],i1) < std::tie(vec1[i2],i2); } ); 

  std::size_t m3(0); 
  std::generate(std::begin(vec2), std::end(vec2), [&]{ return m3++; });
  vec2[0] = 0;
  vec2[1] = 0;
  vec2[2] = 0;
  vec2[3] = 0;

  std::sort(std::begin(vec2), std::end(vec2), [&](int i1, int i2) {
  return std::tie(vec1[i1],i1) < std::tie(vec1[i2],i2); } ); 

  for(int i = 0; i < 4; i++){
  std::cerr << vec1[i] << " " << ind[i] << " " << vec2[i] << std::endl;
  }

  exit(1);
  ///
  */

  mtr.tree.GetCoordinates(coordinates_tree); 
  std::vector<int> sorted_indices(coordinates_tree.size());
  int N = (sorted_indices.size() + 1.0)/2.0; 

  std::size_t m1(0);
  std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
  std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
      return std::tie(coordinates_tree[i1],i1) < std::tie(coordinates_tree[i2],i2); } ); 

  int num_lins = 0;
  double age = coordinates_tree[*sorted_indices.begin()];
  std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
  for(std::vector<int>::iterator it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    if(coordinates_tree[*it_sorted_indices] > age){
      for(; it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        num_lineages[*it_sorted_indices_start] = num_lins;          
      }
      age = coordinates_tree[*it_sorted_indices_start];
    }
    if(*it_sorted_indices < N){
      num_lins++;
    }else{
      num_lins--;
    }
    assert(num_lins >= 1);
  }
  assert(num_lins == 1);

  //jointly sort coordinates and num_lineages
  std::vector<int> num_lineages_tmp = num_lineages;
  std::vector<int>::iterator it_num_lin = num_lineages.begin();
  for(std::vector<int>::iterator it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    *it_num_lin = num_lineages_tmp[*it_sorted_indices];
    it_num_lin++;
  }
  std::sort(coordinates_tree.begin(), coordinates_tree.end());

}

void
GetCoordsAndLineagesForPop(MarginalTree& mtr, Sample& samples, std::vector<int>& exclude_groups, std::vector<Leaves>& descendants, std::vector<float>& coordinates_tree, std::vector<int>& num_lineages){

  mtr.tree.GetCoordinates(coordinates_tree); 
  std::vector<int> sorted_indices(coordinates_tree.size());
  int N = (sorted_indices.size() + 1.0)/2.0; 

  std::size_t m1(0);
  std::generate(std::begin(sorted_indices), std::end(sorted_indices), [&]{ return m1++; });
  std::sort(std::begin(sorted_indices), std::end(sorted_indices), [&](int i1, int i2) {
      return std::tie(coordinates_tree[i1],i1) < std::tie(coordinates_tree[i2],i2); } ); 

  int num_lins = 0;
  int num_terminal = 0;
  int num_exclude = 0;
  std::vector<int> exclude_lineages(mtr.tree.nodes.size(), 0);
  double age = coordinates_tree[*sorted_indices.begin()];
  std::vector<int>::iterator it_sorted_indices_start = sorted_indices.begin();
  for(std::vector<int>::iterator it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){

    if(coordinates_tree[*it_sorted_indices] > age){
      for(; it_sorted_indices_start != it_sorted_indices; it_sorted_indices_start++){
        assert(num_lins >= num_terminal + num_exclude);
        num_lineages[*it_sorted_indices_start] = num_lins - num_terminal - num_exclude;          
      }
      age = coordinates_tree[*it_sorted_indices_start];
    }

    bool ignore = true, ignore2 = true, exclude = false;
    if(*it_sorted_indices < N){

      for(int i = 0; i < samples.group_of_interest.size(); i++){
        std::vector<int>::iterator it_desc = descendants[*it_sorted_indices].member.begin();
        for(; it_desc != descendants[*it_sorted_indices].member.end(); it_desc++ ){
          if(samples.group_of_haplotype[*it_desc] == samples.group_of_interest[i]){
            ignore = false;
            break;
          }
        }
      }
      if(!ignore){
        num_lins++;
        num_terminal++;
      }

      for(int i = 0; i < exclude_groups.size(); i++){
        std::vector<int>::iterator it_desc = descendants[*it_sorted_indices].member.begin();
        for(; it_desc != descendants[*it_sorted_indices].member.end(); it_desc++ ){
          if(samples.group_of_haplotype[*it_desc] == exclude_groups[i]){
            exclude_lineages[*it_desc] = 1;
            break;
          }
        }
      }

    }else{

      int child1 = (*mtr.tree.nodes[*it_sorted_indices].child_left).label;
      for(int i = 0; i < samples.group_of_interest.size(); i++){
        std::vector<int>::iterator it_desc = descendants[child1].member.begin();
        for(; it_desc != descendants[child1].member.end(); it_desc++ ){
          if(samples.group_of_haplotype[*it_desc] == samples.group_of_interest[i]){
            ignore = false;
            break;
          }
        }
      }

      if(!ignore && child1 < N){
        num_terminal--;
      }

      ignore2 = true;
      int child2 = (*mtr.tree.nodes[*it_sorted_indices].child_right).label;
      for(int i = 0; i < samples.group_of_interest.size(); i++){
        std::vector<int>::iterator it_desc = descendants[child2].member.begin();
        for(; it_desc != descendants[child2].member.end(); it_desc++){
          if(samples.group_of_haplotype[*it_desc] == samples.group_of_interest[i]){
            ignore2 = false;
            break;
          }
        }
      }
      if(!ignore2 && child2 < N){
        num_terminal--;
      }

      if(!ignore && !ignore2) num_lins--;

      if(exclude_lineages[child1] == 1 || exclude_lineages[child2] == 1){
        exclude_lineages[*it_sorted_indices] = 1;
      }
      if(ignore && !ignore2){
        if(exclude_lineages[child1] == 1 && exclude_lineages[child2] == 0){
          num_exclude++;
        }
      }
      if(ignore2 && !ignore){
        if(exclude_lineages[child2] == 1 && exclude_lineages[child1] == 0){
          num_exclude++;
        }
      }
      if(!ignore2 && !ignore){
        if(exclude_lineages[child1] == 1 && exclude_lineages[child2] == 1){
          num_exclude--;
        }
      }
 

    }
    assert(num_lins >= 0);
    assert(num_terminal >= 0);

  }
  assert(num_lins >= 0);

  //jointly sort coordinates and num_lineages
  std::vector<int> num_lineages_tmp = num_lineages;
  std::vector<int>::iterator it_num_lin = num_lineages.begin();
  for(std::vector<int>::iterator it_sorted_indices = sorted_indices.begin(); it_sorted_indices != sorted_indices.end(); it_sorted_indices++){
    *it_num_lin = num_lineages_tmp[*it_sorted_indices];
    it_num_lin++;
  }
  std::sort(coordinates_tree.begin(), coordinates_tree.end());

}


void
GetBranchLengthsInEpoch(Data& data, std::vector<double>& epoch, std::vector<float>& coordinates, std::vector<int>& num_lineages, std::vector<double>& branch_lengths_in_epoch){

  branch_lengths_in_epoch.resize(epoch.size()-1);
  std::fill(branch_lengths_in_epoch.begin(), branch_lengths_in_epoch.end(), 0.0);

  int ep = 0;

  for(ep = 0; ep < epoch.size(); ep++){
    if(coordinates[0] < epoch[ep]) break;
  }
  ep--;

  branch_lengths_in_epoch[ep] = 0;
  for(int i = 1; i < 2*data.N-1; i++){

    if(coordinates[i] > coordinates[i-1]){

      if(coordinates[i] < epoch[ep+1]){

        if(coordinates[i-1] >= epoch[ep]){
          // epoch[ep], coordinates[i-1], coordinates[i], epoch[ep+1]
          branch_lengths_in_epoch[ep] += num_lineages[i-1] * (coordinates[i] - coordinates[i-1]); 
        }else{
          // coordinates[i-1], epoch[ep], coordinates[i], epoch[ep+1]
          branch_lengths_in_epoch[ep]  = num_lineages[i-1] * (coordinates[i] - epoch[ep]);  
        }

      }else{

        if(coordinates[i-1] >= epoch[ep]){
          // epoch[ep], coordinates[i-1], epoch[ep+1], coordinates[i]
          branch_lengths_in_epoch[ep] += num_lineages[i-1] * (epoch[ep+1] - coordinates[i-1]);
          ep++;
        }else{
          assert(0);
          // coordinates[i-1], epoch[ep], epoch[ep+1], coordinates[i]
          branch_lengths_in_epoch[ep]  = num_lineages[i-1] * (epoch[ep+1] - epoch[ep]);
          ep++;
        }

        if(ep == epoch.size()-1) break;

        while(epoch[ep+1] < coordinates[i] && ep < epoch.size()-1){
          // coordinates[i-1], epoch[ep], epoch[ep+1], coordinates[i]
          branch_lengths_in_epoch[ep]  = num_lineages[i-1] * (epoch[ep+1] - epoch[ep]);
          ep++;
        }

        if(ep < epoch.size()-1){
          assert(coordinates[i] >= epoch[ep]);
          assert(coordinates[i] <= epoch[ep+1]);
          branch_lengths_in_epoch[ep]    = num_lineages[i-1] * (coordinates[i] - epoch[ep]); 
        }else{
          break;
        }

      }

    }

  }

}

////////////// Avg mutation rate //////////////

void 
CalculateAvgMutationRateForChromosome(cxxopts::Options& options, std::vector<double>& mutation_by_epoch, std::vector<double>& opportunity_by_epoch, std::string chr = "NA"){

  std::string line, read;

  AncMutIterators ancmut;
  if(chr == "NA"){
    ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  }else{
    ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + chr + ".anc", options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }   

  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0; 

  Data data(N,L);
  int N_total = 2*data.N-1;

  ///////// read mutations file ///////////

  Mutations mutations(data);
  if(chr == "NA"){
    mutations.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mutations.Read(options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }

  std::vector<int> pos, dist;
  if(options.count("dist")){

    int L_allsnps = 0;

    std::string filename_dist = options["dist"].as<std::string>();
    if(chr != "NA"){
      filename_dist += "_chr" + chr + ".dist";
    }

    igzstream is_L(filename_dist);
    std::string unused;
    std::getline(is_L, unused); 
    while ( std::getline(is_L, unused) ){
      ++L_allsnps;
    }
    is_L.close();

    pos.resize(L_allsnps);
    dist.resize(L_allsnps);
    igzstream is_dist(filename_dist);
    if(is_dist.fail()) is_dist.open(filename_dist + ".gz");
    if(is_dist.fail()){
      std::cerr << "Error while opening " << filename_dist << std::endl;
      exit(1);
    }
    getline(is_dist, line); 
    int snp = 0;
    while(std::getline(is_dist, line)){
      sscanf(line.c_str(), "%d %d", &pos[snp], &dist[snp]);
      snp++;
    }
    is_dist.close();

  }else{

    pos.resize(data.L);
    dist.resize(data.L);
    int snp = 0;
    for(std::vector<SNPInfo>::iterator it_mut = mutations.info.begin(); it_mut != mutations.info.end(); it_mut++){
      pos[snp]   = (*it_mut).pos;
      dist[snp]  = (*it_mut).dist;
      snp++;
    }

  }


  ///////// EPOCHES /////////
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  } 

  int num_epochs;
  std::vector<double> epochs;
  float log_10 = std::log(10);
  if(options.count("bins")){

    double log_age = std::log(0);
    double age = 0;

    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      ep++;
    }
    epoch_boundary = epoch_lower;
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0) epochs.push_back(age);
        epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
    num_epochs = epochs.size();

  }else{

    num_epochs = 31;
    epochs.resize(num_epochs);
    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }

  ////////// Count number of bases by type ///////////
  //double total_num_bases = (*std::prev(mutations.info.end(),1)).pos - (*mutations.info.begin()).pos;
  double total_num_bases = 1e9;
  std::vector<double> count_bases(mutations.info.size(), 0.0);
  std::vector<double>::iterator it_count_bases = count_bases.begin();

  it_mut = mutations.info.begin();
  std::vector<int>::iterator it_dist = dist.begin(), it_pos = pos.begin();

  //if first snp is included
  if((*it_mut).pos == (*it_pos)){
    *it_count_bases = 0.5 * (*it_dist)/total_num_bases;
    //std::cout << *it_count_bases * total_num_bases << std::endl;
    it_mut++;
    it_count_bases++;
  }
  it_dist++;
  it_pos++;
  while(it_count_bases != count_bases.end()){

    if((*it_mut).pos == (*it_pos)){
      it_dist--;
      *it_count_bases  = 0.5 * (*it_dist)/total_num_bases;
      if(it_dist != dist.end()){
        it_dist++;
        *it_count_bases += 0.5 * (*it_dist)/total_num_bases;
      }
      //std::cout << *it_count_bases * total_num_bases << std::endl;
      it_mut++;
      it_count_bases++;
    }

    it_dist++;
    it_pos++;

  }
  assert(it_mut == mutations.info.end());



  ////////////////////
  // Estimate mutation rate through time

  std::vector<double> branch_lengths_in_epoch(num_epochs);
  std::vector<float> coordinates_tree(N_total);
  std::vector<int> num_lineages(N_total);
  int root = N_total-1;
  int i = 0;

  double total_num_mutations = 0;
  double total_branch_length = 0.0;

  //read tree
  ancmut.FirstSNP(mtr, it_mut);

  GetCoordsAndLineages(mtr, coordinates_tree, num_lineages);
  GetBranchLengthsInEpoch(data, epochs, coordinates_tree, num_lineages, branch_lengths_in_epoch);

  SNPInfo snp_info;
  int current_tree = (*it_mut).tree;
  int count_snps = 0;
  for(int snp = 0; snp < data.L; snp++){

    snp_info = (*it_mut);
    if(snp_info.branch.size() == 1){

      //std::cerr << num_tree << " " << snp_info.tree << std::endl;
      if((*it_mut).tree != current_tree){
        current_tree = (*it_mut).tree;
        mtr.tree.GetCoordinates(coordinates_tree);
        GetCoordsAndLineages(mtr, coordinates_tree, num_lineages);
        GetBranchLengthsInEpoch(data, epochs, coordinates_tree, num_lineages, branch_lengths_in_epoch);
      }
      assert(current_tree == snp_info.tree);

      // identify epoch and add to number of mutations per lineage in epoch
      int ep = 0;
      while(epochs[ep] <= snp_info.age_begin){
        ep++;
        if(ep == epochs.size()) break;
      }
      ep--;

      assert(ep >= 0);

      count_snps++; 

      int ep_begin = ep;
      float age_end = snp_info.age_end;
      double branch_length = age_end - snp_info.age_begin;
      if(ep < num_epochs-1){
        if(age_end <= epochs[ep+1]){

          mutation_by_epoch[ep] += 1.0;

        }else{

          mutation_by_epoch[ep] += (epochs[ep+1] - snp_info.age_begin)/branch_length;
          ep++;
          while(epochs[ep+1] <= age_end && ep < num_epochs-1){
            mutation_by_epoch[ep] += (epochs[ep+1]-epochs[ep])/branch_length;
            ep++;
          }
          if(ep + 1 == num_epochs){
            assert(epochs[ep] <= age_end);
          }else{
            mutation_by_epoch[ep] += (age_end-epochs[ep])/branch_length;
            assert(epochs[ep+1] > age_end);
          }

        }
      }

      double foo1 = 0.0, foo2 = 0.0;
      for(int ep_tmp = 0; ep_tmp < num_epochs; ep_tmp++){
        double bl = branch_lengths_in_epoch[ep_tmp];
        opportunity_by_epoch[ep_tmp] += (bl * count_bases[snp]);
        foo1 += bl;
      }
      for(int i = data.N; i < 2*data.N-1; i++){
        total_branch_length += (2*data.N-i) * (coordinates_tree[i] - coordinates_tree[i-1]) * count_bases[snp];
        foo2 += (2*data.N-i) * (coordinates_tree[i] - coordinates_tree[i-1]);
      }

    }

    ancmut.NextSNP(mtr, it_mut);

  }

  ancmut.CloseFiles(); 

}

void AvgMutationRate(cxxopts::Options& options){

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output.  Optional: chr, first_chr, last_chr, years_per_gen, bins, dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate avg mutation rate through time." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  if(options.count("chr")){
    std::cerr << "Calculating average mutation rate for " << options["input"].as<std::string>() << "* ..." << std::endl;
  }else if(options.count("first_chr") && options.count("last_chr")){
    std::cerr << "Calculating average mutation rate for " << options["input"].as<std::string>() << "_chr" << options["first_chr"].as<int>();
    std::cerr << " - " << options["input"].as<std::string>() << "_chr" << options["last_chr"].as<int>() << " ..." << std::endl;
  }else{
    std::cerr << "Calculating average mutation rate for " << options["input"].as<std::string>() << " ..." << std::endl;
  }

  ///////// EPOCHES /////////
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  int num_epochs;
  std::vector<float> epochs;
  float log_10 = std::log(10);
  if(options.count("bins")){

    double log_age = std::log(0);
    double age = 0;

    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      ep++;
    }
    epoch_boundary = epoch_lower;
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0) epochs.push_back(age);
        epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
    num_epochs = epochs.size();

  }else{

    num_epochs = 31;
    epochs.resize(num_epochs);
    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }

  ///////////////////////////////

  std::vector<double> mutation_by_epoch;
  std::vector<double> opportunity_by_epoch;
  mutation_by_epoch.resize(num_epochs);
  opportunity_by_epoch.resize(num_epochs);
  std::fill(mutation_by_epoch.begin(), mutation_by_epoch.end(), 0.0);
  std::fill(opportunity_by_epoch.begin(), opportunity_by_epoch.end(), 0.0);

  if(options.count("chr")){
    igzstream is_chr(options["chr"].as<std::string>());
    if(is_chr.fail()){
      std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
    }
    std::string line;
    while(getline(is_chr, line)){
      CalculateAvgMutationRateForChromosome(options, mutation_by_epoch, opportunity_by_epoch, line);
    }
    is_chr.close();
  }else if(options.count("first_chr") && options.count("last_chr")){
    for(int chr = options["first_chr"].as<int>(); chr <= options["last_chr"].as<int>(); chr++){
      CalculateAvgMutationRateForChromosome(options, mutation_by_epoch, opportunity_by_epoch, std::to_string(chr));
    } 
  }else{
    CalculateAvgMutationRateForChromosome(options, mutation_by_epoch, opportunity_by_epoch);
  }


  // Summarize

  std::ofstream os(options["output"].as<std::string>() + "_avg.rate");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  std::vector<double> rate(num_epochs);

  //divide every entry by epoch delta time and length of genome
  std::vector<double>::iterator it_mutation = mutation_by_epoch.begin();
  std::vector<double>::iterator it_opportunity = opportunity_by_epoch.begin();
  std::vector<float>::iterator it_epoch  = epochs.begin();
  double total_num_bases = 1e9;
  int e = 0;
  while(it_epoch != epochs.end()){
    rate[e] = (*it_mutation/(*it_opportunity))/total_num_bases;
    os << *it_epoch << " " << rate[e] << "\n";
    it_epoch++;
    it_mutation++;
    it_opportunity++;
    e++;
  }
  os.close();

  plot p(60,10);
  p.draw(epochs, rate);

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6);
#ifdef __APPLE__
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000000.0 << "Mb." << std::endl;
#else
  std::cerr << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
#endif
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

}

