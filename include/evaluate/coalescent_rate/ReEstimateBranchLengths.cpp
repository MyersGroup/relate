//Inference of branch lengths using EM and MCMC.

#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include "gzstream.h"
#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "branch_length_estimator.hpp"
#include "anc_builder.hpp"

void ShowProgress(int progress){

  std::cerr << "[" << progress << "%]\r";
  std::cerr.flush();
  /*
  int p = 0;
  std::cerr << "[";
  for(; p < progress; p++){
    std::cerr << "*";
  }
  for(; p < 100; p++){
    std::cerr << " ";
  }
  std::cerr << "]" << progress << "%\r";
  std::cerr.flush();
  */

}

int ReEstimateBranchLengths(cxxopts::Options& options){

  int seed;
  if(!options.count("seed")){
    seed = std::time(0) + getpid();
  }else{
    seed = options["seed"].as<int>();
  }

  int Ne = 3e4;
  double mutation_rate = options["mutation_rate"].as<float>();
  std::string line;
  double tmp;

  //parse data
  int N;
  igzstream is_N(options["input"].as<std::string>() + ".anc");
  if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
  if(is_N.fail()){
    std::cerr << "Error while opening .anc file." << std::endl;
    exit(1);
  } 
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  //make this more efficient
  int L = 0;
  igzstream is_L;
  
  if(options.count("dist")){
    is_L.open(options["dist"].as<std::string>());
    if(is_L.fail()){
      std::cerr << "Error while opening .dist file." << std::endl;
      exit(1);
    } 
  }else{
    is_L.open(options["input"].as<std::string>() + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".gz");
    if(is_L.fail()){
      std::cerr << "Error while opening .mut file." << std::endl;
      exit(1);
    } 
  }
  
  while(std::getline(is_L, line)){
    ++L;
  }
  L--;
  is_L.close();

  Data data(N, L, Ne, mutation_rate);

  Mutations mut(data);
  mut.Read(options["input"].as<std::string>() + ".mut");

  data.pos.resize(L);

  if(options.count("dist")){
    igzstream is_dist(options["dist"].as<std::string>());
    if(is_dist.fail()){
      std::cerr << "Error while opening " << options["dist"].as<std::string>() << std::endl;
      exit(1);
    }
    getline(is_dist, line); 
    int dtmp, snp = 0;
    while(std::getline(is_dist, line)){
      sscanf(line.c_str(), "%d %d", &dtmp, &data.pos[snp]);
      snp++;
    }
    is_dist.close();
  }else{
    std::vector<int>::iterator it_pos = data.pos.begin();
    for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
      *it_pos = (*it_mut).dist;
      it_pos++;
    }
  }

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Reinferring branch lengths for " << options["input"].as<std::string>() << " ..." << std::endl;

  
  // read epochs and population size 
  igzstream is(options["coal"].as<std::string>()); 
  std::vector<double> epoch, coalescent_rate;
  getline(is, line);
  getline(is, line);
  std::istringstream is_epoch(line);
  while(is_epoch){
    is_epoch >> tmp;
    epoch.push_back(tmp/data.Ne);
  }
  getline(is, line);
  is.close();

  std::istringstream is_pop_size(line);
  is_pop_size >> tmp >> tmp;
  while(is_pop_size){
    is_pop_size >> tmp;
    //tmp = 1.0/data.Ne; 
    if(tmp == 0.0 && coalescent_rate.size() > 0){
      if(*std::prev(coalescent_rate.end(),1) > 0.0){
        coalescent_rate.push_back(*std::prev(coalescent_rate.end(),1));
      }
      //coalescent_rate.push_back(1);
    }else{
      coalescent_rate.push_back(tmp * data.Ne);
    }
  }

  for(int i = (int)coalescent_rate.size()-1; i > 0; i--){
    if(coalescent_rate[i-1] == 0){
      if(coalescent_rate[i] > 0.0){
        coalescent_rate[i-1] = coalescent_rate[i];
      }else{
        coalescent_rate[i-1] = 1.0;
      }
    } 
  } 

  //multiply by mutation rate
  is.open(options["mrate"].as<std::string>());
  double mepoch, mrate;
  int e = 0;
  while(getline(is, line)){
    sscanf(line.c_str(), "%lf %lf", &mepoch, &mrate);
    assert(mepoch/data.Ne == epoch[e]);
    if(mrate > 0){
      coalescent_rate[e] *= data.mu/mrate;
    }
    e++;
  }

  /* 
  for(int i = 0; i < (int)coalescent_rate.size(); i++){
    std::cerr << coalescent_rate[i] << " ";
  }
  std::cerr << std::endl;
  */

  ///////////////////////////////////////// TMRCA Inference /////////////////////////
  //Infer Branchlengths

  AncesTree anc;
  anc.Read(options["input"].as<std::string>() + ".anc");

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  //Infer branch lengths
  InferBranchLengths bl(data);
  //EstimateBranchLengths bl2(data);

  int num_trees = anc.seq.size();
  int progress_interval = (int)(num_trees/100.0) + 1;
  int count_trees = 0, progress = 0;
  CorrTrees::iterator it_seq   = anc.seq.begin();
  for(; it_seq != anc.seq.end(); it_seq++){

    if(count_trees % progress_interval == 0){
      progress++;
      ShowProgress(progress); 
    }
    count_trees++; 
    //for(std::vector<Node>::iterator it_node = (*it_seq).tree.nodes.begin(); it_node != (*it_seq).tree.nodes.end(); it_node++){
    //  (*it_node).num_events = 0.0;
    //}
    bl.MCMCVariablePopulationSize(data, (*it_seq).tree, epoch, coalescent_rate, seed); //this is estimating times
    //bl2.MCMCVariablePopulationSize(data, (*it_seq).tree, epoch, coalescent_rate, seed); //this is estimating times
  }
  ShowProgress(100);
  std::cerr << std::endl;
  //Dump to file
  anc.Dump(options["output"].as<std::string>() + ".anc");

  ////////////////////////// Update mutation file
  
  CorrTrees::iterator it_anc = anc.seq.begin();
  std::vector<float> coordinates(2*data.N-1);
  int num_tree = mut.info[0].tree;
  int root = 2*data.N-2;

  (*it_anc).tree.GetCoordinates(coordinates);

  std::vector<SNPInfo>::iterator it_mut = mut.info.begin();
  for(; it_mut != mut.info.end(); it_mut++){
    //need the tree such that snp_of_next_tree > (*it_mut).snp_id
    //and snp_of_current_tree <= (*it_mut).snp_id
    if((*it_mut).tree > num_tree){
      while((*it_mut).tree > num_tree){
        it_anc++;
        if(it_anc == anc.seq.end()){
          it_anc--;
          break;
        }
        num_tree++;
      }
      (*it_anc).tree.GetCoordinates(coordinates);
    }
    if((*it_mut).tree != num_tree) std::cerr << (*it_mut).tree << " " << num_tree << std::endl;
    if((*it_mut).branch.size() == 1){
      int branch = *(*it_mut).branch.begin();
      (*it_mut).age_begin = coordinates[branch];
      (*it_mut).age_end   = coordinates[(*(*it_anc).tree.nodes[branch].parent).label]; 
    }
  }
  mut.Dump(options["output"].as<std::string>() + ".mut"); 

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
