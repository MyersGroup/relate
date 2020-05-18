//Inference of branch lengths using EM and MCMC.

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "branch_length_estimator.hpp"
#include "anc_builder.hpp"

int GetBranchLengths(cxxopts::Options& options, int chunk_index, int first_section, int last_section){

  bool popsize = false;
  if(!options.count("effectiveN") && !options.count("coal")) popsize = true;
  bool help = false;
  if(popsize || !options.count("mutation_rate") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: effectiveN, mutation_rate, first_section, last_section, output. Optional: coal, sample_ages." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use after FindEquivalentBranches to infer branch lengths." << std::endl;
    exit(0);
  }

  int seed;
  if(!options.count("seed")){
    seed = std::time(0) + getpid();
  }else{
    seed = options["seed"].as<int>();
  }

  std::string file_out = options["output"].as<std::string>() + "/";

  int N, L, num_windows;
  FILE* fp = fopen((file_out + "parameters_c" + std::to_string(chunk_index) + ".bin").c_str(), "r");
  assert(fp != NULL);
  fread(&N, sizeof(int), 1, fp);
  fread(&L, sizeof(int), 1, fp);
  fread(&num_windows, sizeof(int), 1, fp);
  fclose(fp);
  num_windows--;

  int Ne = 30000;
  if(options.count("effectiveN")) Ne = (int) options["effectiveN"].as<float>();
  double mutation_rate = options["mutation_rate"].as<float>();
  Data data((file_out + "chunk_" + std::to_string(chunk_index) + ".bp").c_str(), (file_out + "parameters_c" + std::to_string(chunk_index) + ".bin").c_str(), Ne, mutation_rate); //struct data is defined in data.hpp 
  const std::string dirname = file_out + "chunk_" + std::to_string(chunk_index) + "/";

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Inferring branch lengths of AncesTrees in sections " << first_section << "-" << last_section << "..." << std::endl;
 
  bool is_coal = false;
  std::vector<double> epoch, coalescent_rate;
  if(options.count("coal")){

    std::string line;
    double tmp;
    is_coal = true;
    // read epochs and population size 
    std::ifstream is(options["coal"].as<std::string>()); 
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

  }


  //////////////////////////////////
  //Parse Data
  if(first_section >= num_windows) return 1;

  std::vector<double> sample_ages(N);
  if(options.count("sample_ages")){
    igzstream is_ages(options["sample_ages"].as<std::string>());
    int i = 0; 
    while(is_ages >> sample_ages[i]){
      i++;
      //sample_ages[i] = sample_ages[i-1];
      //i++;
      if(i == N) break;
    }
    if(i < N) sample_ages.clear();
  }else{
    sample_ages.clear(); 
  }

  ///////////////////////////////////////// TMRCA Inference /////////////////////////
  //Infer Branchlengths

  if(sample_ages.size() == 0){

    last_section = std::min(num_windows-1, last_section);
    for(int section = first_section; section <= last_section; section++){

      std::cerr << "[" << section << "/" << last_section << "]\r";
      std::cerr.flush(); 

      //Read anc
      AncesTree anc;
      std::string filename; 
      filename = dirname + options["output"].as<std::string>() + "_" + std::to_string(section) + ".anc";
      anc.ReadBin(filename);

      //Infer branch lengths
      InferBranchLengths bl(data);
      //EstimateBranchLengths bl2(data);
      //EstimateBranchLengthsWithSampleAge bl2(data, sample_ages);

      int num_sec = (int) anc.seq.size()/100.0 + 1;

      if(is_coal){
        for(CorrTrees::iterator it_seq = anc.seq.begin(); it_seq != anc.seq.end(); it_seq++){
          //bl.MCMCVariablePopulationSizeForRelate(data, (*it_seq).tree, epoch, coalescent_rate, seed); //this is estimating times
          bl.MCMCVariablePopulationSizeForRelate(data, (*it_seq).tree, epoch, coalescent_rate, seed); //this is estimating times
        }
      }else{
        int count = 0;
        CorrTrees::iterator it_seq = anc.seq.begin();
        for(; it_seq != anc.seq.end(); it_seq++){
          if(count % num_sec == 0){
            std::cerr << "[" << section << "/" << last_section << "] " << "[" << count << "/" << anc.seq.size() << "]\r";
            std::cerr.flush(); 
          }
          count++;
          //bl.MCMC(data, (*it_seq).tree, seed); //this is estimating times
          bl.MCMC(data, (*it_seq).tree, seed); //this is estimating times
        }
        std::cerr << "[" << section << "/" << last_section << "] " << "[" << count << "/" << anc.seq.size() << "]\r";
      }

      //Dump to file
      anc.DumpBin(filename);

    }

  }else{

    last_section = std::min(num_windows-1, last_section);
    for(int section = first_section; section <= last_section; section++){

      std::cerr << "[" << section << "/" << last_section << "]\r";
      std::cerr.flush(); 

      //Read anc
      AncesTree anc;
      std::string filename; 
      filename = dirname + options["output"].as<std::string>() + "_" + std::to_string(section) + ".anc";
      anc.ReadBin(filename);
      anc.sample_ages = sample_ages;

      //Infer branch lengths
      EstimateBranchLengthsWithSampleAge bl(data, anc.sample_ages);

      int num_sec = (int) anc.seq.size()/100.0 + 1;

      if(is_coal){
        for(CorrTrees::iterator it_seq = anc.seq.begin(); it_seq != anc.seq.end(); it_seq++){
          bl.MCMCVariablePopulationSizeForRelate(data, (*it_seq).tree, epoch, coalescent_rate, seed); //this is estimating times
        }
      }else{
        int count = 0;
        CorrTrees::iterator it_seq = anc.seq.begin();
        for(; it_seq != anc.seq.end(); it_seq++){
          if(count % num_sec == 0){
            std::cerr << "[" << section << "/" << last_section << "] " << "[" << count << "/" << anc.seq.size() << "]\r";
            std::cerr.flush(); 
          }
          count++;
          bl.MCMC(data, (*it_seq).tree, seed); //this is estimating times
        }
        std::cerr << "[" << section << "/" << last_section << "] " << "[" << count << "/" << anc.seq.size() << "]\r";
      }

      //Dump to file
      anc.DumpBin(filename);

    }

  
  
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
