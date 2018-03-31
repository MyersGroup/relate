#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#include "gzstream.h"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"




float
GetCoalescentRate(Node n, float factor, std::vector<float>& epoch, std::vector<CollapsedMatrix<float>>& coalescent_rate_data, std::vector<int>& leaves){

  //go through tree and fill in the matrix
  if(n.child_left != NULL){

    std::vector<int> leaves_child_left, leaves_child_right;

    Node child_left  = *n.child_left;
    Node child_right = *n.child_right;

    int num_children_left, num_children_right;

    float coalescent_time = GetCoalescentRate(child_left, factor, epoch, coalescent_rate_data, leaves_child_left) + child_left.branch_length;
    GetCoalescentRate(child_right, factor, epoch, coalescent_rate_data, leaves_child_right);
   
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
              coalescent_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
              coalescent_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - epoch[e]);
              break;
            }else{
              assert(coalescent_time >= epoch[e+1]);
              coalescent_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - epoch[e]);
            }
          }
        }else{
          for(int e = 0; e < (int) epoch.size() - 1; e++){
            if(coalescent_time < epoch[e+1]){
              assert(coalescent_time >= epoch[e]);
              coalescent_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
              coalescent_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - epoch[e]);
              break;
            }else{
              assert(coalescent_time >= epoch[e+1]);
              coalescent_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - epoch[e]);
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

int CoalescentRateForSection(cxxopts::Options& options, int chr = -1){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: num_bins." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .anc file and calculates pairwise coalescent rate. Output is bin file. Use SummarizeCoalesecntRate to obtain coalescent rates." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  if(chr == -1){
    std::cerr << "Calculating coalescent rate for " << options["input"].as<std::string>() << " ..." << std::endl;
  }else{
    std::cerr << "Calculating coalescent rate for " << options["input"].as<std::string>() << "_chr" << chr << " ..." << std::endl;
  }

  ////////////////////////
  //read in anc file

  int N;
  igzstream is_N;
  if(chr == -1){
    is_N.open(options["input"].as<std::string>() + ".anc");
    if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
  }else{
    is_N.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc");
    if(is_N.fail()) is_N.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc.gz");
  }
  if(is_N.fail()){
    std::cerr << "Error while opening .anc file." << std::endl;
    exit(1);
  } 
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  int L = 0;
  igzstream is_L;
  if(chr == -1){
    is_L.open(options["input"].as<std::string>() + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz");
  }else{
    is_L.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut.gz");
  }
  if(is_L.fail()){
    std::cerr << "Error while opening .mut file." << std::endl;
    exit(1);
  } 
  
  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }
  is_L.close();
  Data data(N,L);
  
  //Read anc
  AncesTree anc;
  Mutations mut;

  if(chr == -1){
    anc.Read(options["input"].as<std::string>() + ".anc");
    mut.Read(options["input"].as<std::string>() + ".mut");
  }else{
    anc.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc");
    mut.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
  }

  //calculate epoch times
  
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }
  
  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epoch(num_epochs);
  std::vector<CollapsedMatrix<float>> coalescent_rate_data(num_epochs);
  epoch[0] = 0.0;
  coalescent_rate_data[0].resize(data.N, data.N);
  std::fill(coalescent_rate_data[0].vbegin(), coalescent_rate_data[0].vend(), 0.0);
  epoch[1] = 1e3/years_per_gen;
  coalescent_rate_data[1].resize(data.N, data.N);
  std::fill(coalescent_rate_data[1].vbegin(), coalescent_rate_data[1].vend(), 0.0);
  
  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){

    epoch[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    coalescent_rate_data[e].resize(data.N, data.N);
    std::fill(coalescent_rate_data[e].vbegin(), coalescent_rate_data[e].vend(), 0.0);
  
  }
  epoch[num_epochs-1] = 5e7;
  coalescent_rate_data[num_epochs-1].resize(data.N, data.N);
  std::fill(coalescent_rate_data[num_epochs-1].vbegin(), coalescent_rate_data[num_epochs-1].vend(), 0.0);
  

  ////////////////////////////////
  //Pairwise coalescent rate
  //In each tree, find the coalescent time. Then update count_per_epoch and coalescent_time_in_epoch. 

  float factor = 0.0;
  CorrTrees::iterator it_anc = anc.seq.begin();
  for(; it_anc != std::prev(anc.seq.end(),1); it_anc++){

    std::vector<int> leaves;
    factor = 1.0;
    GetCoalescentRate(*std::prev((*it_anc).tree.nodes.end(),1), factor, epoch, coalescent_rate_data, leaves);
    
  }


  //output as bin
  FILE* fp;
  if(chr == -1){
    fp = fopen((options["output"].as<std::string>() + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + std::to_string(chr) + ".bin" ).c_str(), "wb");  
  }

  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epoch[0], sizeof(float), epoch.size(), fp);
  for(std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_data = coalescent_rate_data.begin(); it_coalescent_rate_data != coalescent_rate_data.end();){
    (*it_coalescent_rate_data).DumpToFile(fp);
    it_coalescent_rate_data++;
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


