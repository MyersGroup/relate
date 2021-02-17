#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#include "gzstream.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

#include "mutations.hpp"
#include "coal_tree.hpp"

//////// functions for estimating pairwise coalescence rate /////////

float
GetCoalescentRate(Node n, float factor, std::vector<float>& epoch, std::vector<CollapsedMatrix<float>>& coalescence_rate_data, std::vector<int>& leaves){

  //go through tree and fill in the matrix
  if(n.child_left != NULL){

    std::vector<int> leaves_child_left, leaves_child_right;

    Node child_left  = *n.child_left;
    Node child_right = *n.child_right;

    int num_children_left, num_children_right;

    float coalescent_time = GetCoalescentRate(child_left, factor, epoch, coalescence_rate_data, leaves_child_left) + child_left.branch_length;
    GetCoalescentRate(child_right, factor, epoch, coalescence_rate_data, leaves_child_right);

    leaves.resize(leaves_child_left.size() + leaves_child_right.size());
    std::vector<int>::iterator it_leaves = leaves.begin();

    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      *it_leaves = *it_leaves_child_left;
      it_leaves++;
    }
    for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){
      *it_leaves = *it_leaves_child_right;
      it_leaves++;
    }

    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){

        assert(*it_leaves_child_left != *it_leaves_child_right);

        if(*it_leaves_child_left < *it_leaves_child_right){
          for(int e = 0; e < (int) epoch.size() - 1; e++){
            if(coalescent_time < epoch[e+1]){
              assert(coalescent_time >= epoch[e]);
              coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
              coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - epoch[e]);
              break;
            }else{
              assert(coalescent_time >= epoch[e+1]);
              coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - epoch[e]);
            }
          }
        }else{
          for(int e = 0; e < (int) epoch.size() - 1; e++){
            if(coalescent_time < epoch[e+1]){
              assert(coalescent_time >= epoch[e]);
              coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
              coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - epoch[e]);
              break;
            }else{
              assert(coalescent_time >= epoch[e+1]);
              coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - epoch[e]);
            }
          } 
        }

      }
    }

    return coalescent_time;

  }else{

    leaves.push_back(n.label);
    return 0.0;

  }


}

float
GetCoalescentRate(Node n, float factor, std::vector<float>& epoch, std::vector<double>& sample_ages, std::vector<CollapsedMatrix<float>>& coalescence_rate_data, std::vector<int>& leaves){

  //go through tree and fill in the matrix
  if(n.child_left != NULL){

    std::vector<int> leaves_child_left, leaves_child_right;

    Node child_left  = *n.child_left;
    Node child_right = *n.child_right;

    int num_children_left, num_children_right;

    float coalescent_time = GetCoalescentRate(child_left, factor, epoch, sample_ages, coalescence_rate_data, leaves_child_left) + child_left.branch_length;
    GetCoalescentRate(child_right, factor, epoch, sample_ages, coalescence_rate_data, leaves_child_right);

    leaves.resize(leaves_child_left.size() + leaves_child_right.size());
    std::vector<int>::iterator it_leaves = leaves.begin();

    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      *it_leaves = *it_leaves_child_left;
      it_leaves++;
    }
    for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){
      *it_leaves = *it_leaves_child_right;
      it_leaves++;
    }

    double max_sample_age = 0.0;
    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){

        assert(*it_leaves_child_left != *it_leaves_child_right);

        max_sample_age = sample_ages[*it_leaves_child_left];
        if(max_sample_age < sample_ages[*it_leaves_child_right]){
          max_sample_age = sample_ages[*it_leaves_child_right];
        }

        if(max_sample_age == 0.0){
          if(*it_leaves_child_left < *it_leaves_child_right){
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(coalescent_time < epoch[e+1]){
                assert(coalescent_time >= epoch[e]);
                coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
                coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - epoch[e]);
                break;
              }else{
                assert(coalescent_time >= epoch[e+1]);
                coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - epoch[e]);
              }
            }
          }else{
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(coalescent_time < epoch[e+1]){
                assert(coalescent_time >= epoch[e]);
                coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
                coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - epoch[e]);
                break;
              }else{
                assert(coalescent_time >= epoch[e+1]);
                coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - epoch[e]);
              }
            } 
          }
        }else{        
          if(*it_leaves_child_left < *it_leaves_child_right){
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(max_sample_age < epoch[e+1]){
                if(max_sample_age >= epoch[e]){                
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - max_sample_age);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - max_sample_age);
                  }
                }else{
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - epoch[e]);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - epoch[e]);
                  }
                }
              }
            }
          }else{
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(max_sample_age < epoch[e+1]){
                if(max_sample_age >= epoch[e]){                
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - max_sample_age);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - max_sample_age);
                  }
                }else{
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - epoch[e]);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - epoch[e]);
                  }
                }
              }
            }
          }
        }

      }
    }

    return coalescent_time;

  }else{

    leaves.push_back(n.label);
    return sample_ages[n.label];

  }


}

int 
CoalescentRateForSection(cxxopts::Options& options, std::string chr = "NA"){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: years_per_gen, dist, bins." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .anc file and calculates pairwise coalescence rate. Output is bin file. Use SummarizeCoalesecntRate to obtain coalescence rates." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  if(chr == "NA"){
    std::cerr << "Calculating coalescence rate for " << options["input"].as<std::string>() << " ..." << std::endl;
  }else{
    std::cerr << "Calculating coalescence rate for " << options["input"].as<std::string>() << "_chr" << chr << " ..." << std::endl;
  }

  ////////////////////////
  //read in anc file

  fasta mask;
  double cutoff = 0.9;
  if(options.count("mask")){
    if(chr == "NA"){
      mask.Read(options["mask"].as<std::string>());
    }else{
      mask.Read(options["mask"].as<std::string>() + "_chr" + chr + ".fa");
    }
  }

  AncMutIterators ancmut;
  if(options.count("dist")){
    if(chr == "NA"){
      ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut", options["dist"].as<std::string>());
    }else{
      ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + chr + ".anc", options["input"].as<std::string>() + "_chr" + chr + ".mut", options["dist"].as<std::string>() + "_chr" + chr + ".dist");
    } 
  }else{  
    if(chr == "NA"){
      ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
    }else{
      ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + chr + ".anc", options["input"].as<std::string>() + "_chr" + chr + ".mut");
    } 
  } 

  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0; 
  Data data(N,L);

  //Read mut
  Mutations mut;

  if(chr == "NA"){
    mut.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mut.Read(options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }

  //calculate epoch times

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

  std::vector<CollapsedMatrix<float>> coalescence_rate_data(num_epochs);
  coalescence_rate_data[0].resize(data.N, data.N);
  std::fill(coalescence_rate_data[0].vbegin(), coalescence_rate_data[0].vend(), 0.0);
  for(int e = 1; e < num_epochs; e++){
    coalescence_rate_data[e].resize(data.N, data.N);
    std::fill(coalescence_rate_data[e].vbegin(), coalescence_rate_data[e].vend(), 0.0);
  }

  ////////////////////////////////
  //Pairwise coalescence rate
  //In each tree, find the coalescent time. Then update count_per_epoch and coalescent_time_in_epoch. 

  if(ancmut.sample_ages.size() > 0){
    float factor = 0.0;
    num_bases_tree_persists = 0.0;
    double num_passing = 1.0;
    while(num_bases_tree_persists >= 0.0){
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

      num_passing = 1.0;
      if(options.count("mask") > 0){

        int tree_index = (*it_mut).tree;
        int pos_start = (*it_mut).pos;
        int pos_end = pos_start;
        if(it_mut != ancmut.mut_end()){
          while((*it_mut).tree == tree_index){
            pos_end = (*it_mut).pos;
            it_mut++;
            if(it_mut == ancmut.mut_end()) break;
          }
        }
        num_passing = 0.0;
        if(pos_start < mask.seq.size() && pos_end < mask.seq.size()){
          for(int bp = pos_start; bp < pos_end; bp++){
            if(mask.seq[bp-1] == 'P') num_passing++; 
          }
        }

        if(pos_end - pos_start + 1 <= 0){ 
          num_passing = 0.0;
        }else{
          num_passing /= (pos_end - pos_start + 1);
        }

      }

      if(num_passing >= cutoff){
        std::vector<int> leaves;
        factor = num_bases_tree_persists;
        GetCoalescentRate(*std::prev(mtr.tree.nodes.end(),1), factor, epochs, ancmut.sample_ages, coalescence_rate_data, leaves);
      }

    }  
  }else{
    float factor = 0.0;
    num_bases_tree_persists = 0.0;
    double num_passing = 1.0;
    while(num_bases_tree_persists >= 0.0){
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

      num_passing = 1.0;
      num_passing = 1.0;
      if(options.count("mask") > 0){

        int tree_index = (*it_mut).tree;
        int pos_start = (*it_mut).pos;
        int pos_end = pos_start;
        if(it_mut != ancmut.mut_end()){
          while((*it_mut).tree == tree_index){
            pos_end = (*it_mut).pos;
            it_mut++;
            if(it_mut == ancmut.mut_end()) break;
          }
        }
        num_passing = 0.0;
        if(pos_start < mask.seq.size() && pos_end < mask.seq.size()){
          for(int bp = pos_start; bp < pos_end; bp++){
            if(mask.seq[bp-1] == 'P') num_passing++; 
          }
        }

        if(pos_end - pos_start + 1 <= 0){ 
          num_passing = 0.0;
        }else{
          num_passing /= (pos_end - pos_start + 1);
        }

      }

      if(num_passing >= cutoff){
        std::vector<int> leaves;
        factor = num_bases_tree_persists;
        GetCoalescentRate(*std::prev(mtr.tree.nodes.end(),1), factor, epochs, coalescence_rate_data, leaves);
      }
    }  
  }

  if(ancmut.sample_ages.size() > 0){
    std::vector<float> epochs_new;
    std::vector<double> all_sample_ages;
    std::vector<int> epoch_old_index;
    all_sample_ages = ancmut.sample_ages;
    std::sort(all_sample_ages.begin(), all_sample_ages.end());


    double ages = *all_sample_ages.begin();
    int ep = 0;
    if(ages == 0.0){
      epochs_new.push_back(ages);
      epoch_old_index.push_back(ep);
      ep++;
    }else{
      while(epochs[ep] < ages){
        epochs_new.push_back(epochs[ep]);
        epoch_old_index.push_back(ep);
        ep++;
        if(ep == epochs.size()) break;
      }
      if(ages != epochs[ep]){
        epochs_new.push_back(ages);
        epoch_old_index.push_back(ep-1);
      }
    }
    for(std::vector<double>::iterator it_sample_ages = all_sample_ages.begin(); it_sample_ages != all_sample_ages.end(); it_sample_ages++){
      if(ages < *it_sample_ages){
        ages = *it_sample_ages;
        while(epochs[ep] < ages){
          epochs_new.push_back(epochs[ep]);
          epoch_old_index.push_back(ep);
          ep++;
          if(ep == epochs.size()) break;
        }
        if(ep == epochs.size()) break;
        if(ages != epochs[ep]){
          epochs_new.push_back(ages);
          epoch_old_index.push_back(ep-1);
        }
      } 
    }
    for(; ep < epochs.size(); ep++){
      epochs_new.push_back(epochs[ep]);
      epoch_old_index.push_back(ep);
    }
    num_epochs = epochs_new.size();

    std::vector<CollapsedMatrix<float>> coalescence_rate_data_new(num_epochs);
    for(ep = 0; ep < num_epochs-1; ep++){
      coalescence_rate_data_new[ep] = coalescence_rate_data[epoch_old_index[ep]];
      //need to know which samples have sample_ages < epochs[ep]
      for(int i = 0; i < data.N; i++){
        if(ancmut.sample_ages[i] >= epochs_new[ep+1]){
          for(int j = 0; j < data.N; j++){
            coalescence_rate_data_new[ep][i][j] = 0.0;
            coalescence_rate_data_new[ep][j][i] = 0.0;
          } 
        }
      }
    }
    ep = num_epochs-1;
    coalescence_rate_data_new[ep] = coalescence_rate_data[epoch_old_index[ep]];

    //output as bin
    FILE* fp;
    if(chr == "NA"){
      fp = fopen((options["output"].as<std::string>() + ".bin" ).c_str(), "wb");  
    }else{
      fp = fopen((options["output"].as<std::string>() + "_chr" + chr + ".bin" ).c_str(), "wb");  
    }

    fwrite(&num_epochs, sizeof(int), 1, fp);
    fwrite(&epochs_new[0], sizeof(float), epochs_new.size(), fp);
    for(std::vector<CollapsedMatrix<float>>::iterator it_coalescence_rate_data = coalescence_rate_data_new.begin(); it_coalescence_rate_data != coalescence_rate_data_new.end();){
      (*it_coalescence_rate_data).DumpToFile(fp);
      it_coalescence_rate_data++;
    }

    fclose(fp);


  }else{

    //output as bin
    FILE* fp;
    if(chr == "NA"){
      fp = fopen((options["output"].as<std::string>() + ".bin" ).c_str(), "wb");  
    }else{
      fp = fopen((options["output"].as<std::string>() + "_chr" + chr + ".bin" ).c_str(), "wb");  
    }

    fwrite(&num_epochs, sizeof(int), 1, fp);
    fwrite(&epochs[0], sizeof(float), epochs.size(), fp);
    for(std::vector<CollapsedMatrix<float>>::iterator it_coalescence_rate_data = coalescence_rate_data.begin(); it_coalescence_rate_data != coalescence_rate_data.end();){
      (*it_coalescence_rate_data).DumpToFile(fp);
      it_coalescence_rate_data++;
    }

    fclose(fp);

  }

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

  return 0;

}

////// function for estimating coalescence rates in tree /////////
void
CoalescenceRateForTree(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: bins, dist, chr, coal." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Calculating coalescence rates for " << options["input"].as<std::string>() << "..." << std::endl;

  /////////////////////////////////
  //get TMRCA at each SNP

  ////////////////////////////////////////

  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  //decide on epochs 
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

  int num_bootstrap = 1;
  int block_size = 1000;

  coal_tree ct(epochs, num_bootstrap, block_size);

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file

  std::string line;
  std::vector<std::string> chromosomes;
  std::vector<std::string> filename_mut, filename_anc, filename_dist;
  if(options.count("chr") > 0){

    igzstream is_chr(options["chr"].as<std::string>());
    if(is_chr.fail()){
      std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
    }
    while(getline(is_chr, line)){
      chromosomes.push_back(line);
      filename_anc.push_back(options["input"].as<std::string>() + "_chr" + line + ".anc");
      filename_mut.push_back(options["input"].as<std::string>() + "_chr" + line + ".mut");
      if(options.count("dist")) filename_dist.push_back(options["dist"].as<std::string>() + "_chr" + line + ".dist");
    }
    is_chr.close();

  }else{
    chromosomes.resize(1);
    chromosomes[0] = options["input"].as<std::string>(); 
    filename_anc.push_back(options["input"].as<std::string>() + ".anc");
    filename_mut.push_back(options["input"].as<std::string>() + ".mut");
    if(options.count("dist")) filename_dist.push_back(options["dist"].as<std::string>());
  }

  if(!options.count("dist")){
    for(int chr = 0; chr < chromosomes.size(); chr++){

      AncMutIterators ancmut(filename_anc[chr], filename_mut[chr]);
      float num_bases_tree_persists = 0.0;

      ct.update_ancmut(ancmut);

      int tree_count = 0, perc = -1;
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      while(num_bases_tree_persists >= 0.0){
        if( (int) (((double)tree_count)/ancmut.NumTrees() * 100.0) > perc ){
          perc = (int) (((double)tree_count)/ancmut.NumTrees() * 100.0);
          std::cerr << "[" << perc << "%]\r";
        }
        tree_count++;
        ct.populate(mtr.tree, num_bases_tree_persists);	
        num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      }
      std::cerr << "[100%]\r";
      std::cerr << std::endl;

    }
  }else{

    for(int chr = 0; chr < chromosomes.size(); chr++){

      AncMutIterators ancmut(filename_anc[chr], filename_mut[chr], filename_dist[chr]);
      float num_bases_tree_persists = 0.0;

      ct.update_ancmut(ancmut);

      int tree_count = 0, perc = -1;
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      while(num_bases_tree_persists >= 0.0){
        if( (int) (((double)tree_count)/ancmut.NumTrees() * 100.0) > perc ){
          perc = (int) (((double)tree_count)/ancmut.NumTrees() * 100.0);
          std::cerr << "[" << perc << "%]\r";
        }
        tree_count++;
        ct.populate(mtr.tree, num_bases_tree_persists);	
        num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      }
      std::cerr << "[100%]\r";
      std::cerr << std::endl;

    }
  }

  igzstream is_coal(options["input"].as<std::string>() + ".coal");
  bool is_coal_fail = is_coal.fail();
  is_coal_fail = true;

  if(!is_coal_fail){

    float tmp;
    std::vector<float> coal_epochs;
    std::vector<double> coal_values;
    getline(is_coal, line); //group assignment
    getline(is_coal, line); //epoch boundaries
    std::istringstream is_coal_epoch(line);
    int e = 0;
    while(is_coal_epoch >> tmp){
      coal_epochs.push_back(tmp);
      e++;
    }
    assert(coal_epochs.size() == epochs.size());
    getline(is_coal, line);
    std::istringstream is_coal_values(line);
    is_coal_values >> tmp >> tmp;
    while(is_coal_values >> tmp){
      if(tmp == 0.0 && coal_values.size() > 0){
        if(*std::prev(coal_values.end(),1) > 0.0){
          coal_values.push_back(*std::prev(coal_values.end(),1));
        }
      }else{
        coal_values.push_back(tmp);
      }
    }

    for(int i = (int)coal_values.size()-1; i > 0; i--){
      if(coal_values[i-1] == 0){
        if(coal_values[i] > 0.0){
          coal_values[i-1] = coal_values[i];
        }else{
          coal_values[i-1] = 1.0;
        }
      } 
    }
    is_coal.close(); 
    ct.Dump(options["output"].as<std::string>() + ".coal", coal_values);

  }else{
    ct.Dump(options["output"].as<std::string>() + ".coal");
  }

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

void
GenerateConstCoal(cxxopts::Options& options){

  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: bins, years_per_gen." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate coalescence rates for sample." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Generating coalescence rates file with const rate " << options["input"].as<std::string>() << "..." << std::endl;

  /////////////////////////////////
  //get TMRCA at each SNP

  ////////////////////////////////////////

  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  //decide on epochs 
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


  std::ofstream os(options["output"].as<std::string>() + ".coal");
  os << "group1\n";

  for(int e = 0; e < num_epochs; e++){
    os << epochs[e]  << " "; 
  }
  os << "\n";

  std::vector<double> coal(epochs.size());

  double Ne = std::stof(options["input"].as<std::string>());
  os << 0 << " " << 0 << " ";
  for(int e = 0; e < num_epochs; e++){
    os << 1.0/Ne << " ";
  }
  os << "\n"; 

  os.close();



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


