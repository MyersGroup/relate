#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.h"
#include "collapsed_matrix.hpp"
#include "plot.hpp"
#include "sample.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"


void
GetBranchLengthsInEpoche(Data& data, std::vector<float>& epoch, std::vector<float>& coordinates, std::vector<double>& branch_lengths_in_epoch){

  branch_lengths_in_epoch.resize(epoch.size()-1);
  std::fill(branch_lengths_in_epoch.begin(), branch_lengths_in_epoch.end(), 0.0);

  int ep = 0;
  double num_lineages;
  branch_lengths_in_epoch[ep] = 0;
  for(int i = data.N; i < 2*data.N-1; i++){

    num_lineages = 2*data.N - i;
    if(coordinates[i] < epoch[ep+1]){

      if(coordinates[i-1] > epoch[ep]){
        // epoch[ep], coordinates[i-1], coordinates[i], epoch[ep+1]
        branch_lengths_in_epoch[ep] += num_lineages * (coordinates[i] - coordinates[i-1]); 
      }else{
        // coordinates[i-1], epoch[ep], coordinates[i], epoch[ep+1]
        branch_lengths_in_epoch[ep]  = num_lineages * (coordinates[i] - epoch[ep]);  
      }

    }else{

      if(coordinates[i-1] >= epoch[ep]){
        // epoch[ep], coordinates[i-1], epoch[ep+1], coordinates[i]
        branch_lengths_in_epoch[ep] += num_lineages * (epoch[ep+1] - coordinates[i-1]);
        ep++;
      }else{
        assert(0);
        // coordinates[i-1], epoch[ep], epoch[ep+1], coordinates[i]
        branch_lengths_in_epoch[ep]  = num_lineages * (epoch[ep+1] - epoch[ep]);
        ep++;
      }

      if(ep == epoch.size()-1) break;

      while(epoch[ep+1] < coordinates[i] && ep < epoch.size()-1){
        // coordinates[i-1], epoch[ep], epoch[ep+1], coordinates[i]
        branch_lengths_in_epoch[ep]  = num_lineages * (epoch[ep+1] - epoch[ep]);
        ep++;
      }

      if(ep < epoch.size()-1){
        assert(coordinates[i] >= epoch[ep]);
        assert(coordinates[i] <= epoch[ep+1]);
        branch_lengths_in_epoch[ep]    = num_lineages * (coordinates[i] - epoch[ep]); 
      }else{
        break;
      }

    }

  }
  //assert(num_lineages == 2.0);

}

////////////// Avg mutation rate //////////////

void 
CalculateAvgMutationRateForChromosome(cxxopts::Options& options, std::vector<double>& mutation_by_epoch, std::vector<double>& opportunity_by_epoch, int chr = -1){

  std::string line, read;

  //parse data
  int N;
  igzstream is_N;
  if(chr == -1){
    is_N.open(options["input"].as<std::string>() + ".anc");
    if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
    if(is_N.fail()) std::cerr << "Error while opening " << options["input"].as<std::string>() << ".anc(.gz)" << std::endl; 
  }else{
    is_N.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc");
    if(is_N.fail()) is_N.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc.gz");
    if(is_N.fail()) std::cerr << "Error while opening " << options["input"].as<std::string>() + "_chr" + std::to_string(chr) << ".anc(.gz)" << std::endl;
  }
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  //make this more efficient
  int L = 0;
  igzstream is_L;
  if(chr == -1){
    is_L.open(options["input"].as<std::string>() + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz");
    if(is_L.fail()) std::cerr << "Error while opening " << options["input"].as<std::string>() << ".mut(.gz)" << std::endl;
  }else{
    is_L.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut.gz");
    if(is_L.fail()) std::cerr << "Error while opening " << options["input"].as<std::string>() + "_chr" + std::to_string(chr) << ".mut(.gz)" << std::endl;
  } 
  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }
  is_L.close();
  L--;

  Data data(N,L);
  int N_total = 2*data.N-1;

  ///////// read mutations file ///////////

  Mutations mutations(data);
  if(chr == -1){
    mutations.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mutations.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
  }

  std::vector<int> pos, dist;
  if(options.count("dist")){

    int L_allsnps = 0;

    std::string filename_dist = options["dist"].as<std::string>();
    if(chr != -1){
      filename_dist += "_chr" + std::to_string(chr) + ".dist";
    }

    is_L.open(filename_dist);
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

  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epoch(num_epochs);
  epoch[0] = 0.0;
  epoch[1] = 1e3/28.0;
  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){
    epoch[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/28.0; 
  }
  epoch[num_epochs-1] = 5e7;


  ////////// Count number of bases by type ///////////
  //double total_num_bases = (*std::prev(mutations.info.end(),1)).pos - (*mutations.info.begin()).pos;
  double total_num_bases = 1e9;
  std::vector<double> count_bases(mutations.info.size(), 0.0);
  std::vector<double>::iterator it_count_bases = count_bases.begin();

  std::vector<SNPInfo>::iterator it_mut = mutations.info.begin();
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

  MarginalTree mtr;
  std::vector<float> coordinates_tree(N_total);
  int root = N_total-1;
  int i = 0;

  double total_num_mutations = 0;
  double total_branch_length = 0.0;

  igzstream is_anc;
  if(chr == -1){
    is_anc.open(options["input"].as<std::string>() + ".anc");
    if(is_anc.fail()) is_anc.open(options["input"].as<std::string>() + ".anc.gz");
  }else{
    is_anc.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc");
    if(is_anc.fail()) is_anc.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc.gz");
  }
  if(is_anc.fail()){
    std::cerr << "Error while opening .anc file." << std::endl;
    exit(1);
  }


  getline(is_anc,line);
  getline(is_anc,line);
  getline(is_anc,line);

  //read tree
  mtr.Read(line, N);
  mtr.tree.GetCoordinates(coordinates_tree);
  std::sort(coordinates_tree.begin(), coordinates_tree.end());
  GetBranchLengthsInEpoche(data, epoch, coordinates_tree, branch_lengths_in_epoch);

  SNPInfo snp_info;
  int num_tree = 0;
  int count_snps = 0;
  for(int snp = 0; snp < data.L; snp++){

    //std::cerr << snp << std::endl;
    snp_info = mutations.info[snp];
    if(snp_info.branch.size() == 1){

      //std::cerr << num_tree << " " << snp_info.tree << std::endl;
      if(num_tree < snp_info.tree){
        while(num_tree < snp_info.tree){
          if(!getline(is_anc,line)){
            break; 
          };
          num_tree++;
        }

        //read tree
        mtr.Read(line, N);
        mtr.tree.GetCoordinates(coordinates_tree);
        std::sort(coordinates_tree.begin(), coordinates_tree.end());
        GetBranchLengthsInEpoche(data, epoch, coordinates_tree, branch_lengths_in_epoch);

      }
      assert(num_tree == snp_info.tree);

      // identify epoch and add to number of mutations per lineage in epoch
      int ep = 0;
      while(epoch[ep] <= snp_info.age_begin){
        ep++;
        if(ep == epoch.size()) break;
      }
      ep--;

      assert(ep >= 0);

      count_snps++; 

      int ep_begin = ep;
      float age_end = snp_info.age_end;
      double branch_length = age_end - snp_info.age_begin;
      if(ep < num_epochs-1){
        if(age_end <= epoch[ep+1]){

          mutation_by_epoch[ep] += 1.0;

        }else{

          mutation_by_epoch[ep] += (epoch[ep+1] - snp_info.age_begin)/branch_length;
          ep++;
          while(epoch[ep+1] <= age_end && ep < num_epochs-1){
            mutation_by_epoch[ep] += (epoch[ep+1]-epoch[ep])/branch_length;
            ep++;
          }
          if(ep + 1 == num_epochs){
            assert(epoch[ep] <= age_end);
          }else{
            mutation_by_epoch[ep] += (age_end-epoch[ep])/branch_length;
            assert(epoch[ep+1] > age_end);
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

  }

  //double foo;
  //for(int ep = 0; ep < epoch.size(); ep++){
  //  total_num_mutations += mutation_by_epoch[ep];
  //  foo += opportunity_by_epoch[ep];
  //}
  //std::cerr << total_num_mutations << " " << data.L << " " << count_snps << std::endl;
  //std::cerr << total_branch_length << " " << foo << std::endl;

  is_anc.close();

}

void AvgMutationRate(cxxopts::Options& options){

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output.  Optional: first_chr, last_chr, num_bins, dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate avg mutation rate through time." << std::endl;
    exit(0);
  }  

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating average mutation rate..." << std::endl;

  ///////// EPOCHES /////////

  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epoch(num_epochs);
  epoch[0] = 0.0;
  epoch[1] = 1e3/28.0;
  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){
    epoch[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/28.0; 
  }
  epoch[num_epochs-1] = 5e7;

  ///////////////////////////////

  std::vector<double> mutation_by_epoch;
  std::vector<double> opportunity_by_epoch;
  mutation_by_epoch.resize(num_epochs);
  opportunity_by_epoch.resize(num_epochs);
  std::fill(mutation_by_epoch.begin(), mutation_by_epoch.end(), 0.0);
  std::fill(opportunity_by_epoch.begin(), opportunity_by_epoch.end(), 0.0);

  if(options.count("first_chr") && options.count("last_chr")){
    for(int chr = options["first_chr"].as<int>(); chr <= options["last_chr"].as<int>(); chr++){
      CalculateAvgMutationRateForChromosome(options, mutation_by_epoch, opportunity_by_epoch, chr);
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
  std::vector<float>::iterator it_epoch  = epoch.begin();
  double total_num_bases = 1e9;
  int e = 0;
  while(it_epoch != epoch.end()){
    rate[e] = (*it_mutation/(*it_opportunity))/total_num_bases;
    os << *it_epoch << " " << rate[e] << "\n";
    it_epoch++;
    it_mutation++;
    it_opportunity++;
    e++;
  }
  os.close();

  plot p(60,10);
  p.draw(epoch, rate);



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

