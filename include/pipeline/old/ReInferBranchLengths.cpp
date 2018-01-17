//Inference of branch lengths using EM and MCMC.

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include "cxxopts.hpp"
#include "data.hpp"
#include "arg.hpp"
#include "arg_builder.hpp"

int ReInferBranchLengths(cxxopts::Options& options){

  bool help = false;
  if(!options.count("effectiveN") || !options.count("mutation_rate") || !options.count("output") || !options.count("first_section")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: effectiveN, mutation_rate, first_section, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use after AssociateTrees to infer branch lengths." << std::endl;
    exit(0);
  }

  int seed;
  if(!options.count("seed")){
    seed = std::time(0) + getpid();
  }

  int Ne = (int) options["effectiveN"].as<float>();
  double mutation_rate = options["mutation_rate"].as<float>();
  int first_section = options["first_section"].as<int>();
  Data data("sequences.bin", "pos.bin", "recombination_rate.bin", "rpos.bin", Ne, mutation_rate); //struct data is defined in data.hpp 
  int num_windows = (int) ((data.L+data.window_length-1)/data.window_length);

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Reinferring branch lengths in tree " << first_section << "..." << std::endl;

  ///////////////////////////////////////// TMRCA Inference /////////////////////////
  //Infer Branchlengths

  Arg arg;
  arg.seq.emplace_back();

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  std::ifstream is_arg(options["output"].as<std::string>() + ".arg");
  if(is_arg.fail()){
    std::cerr << "Error while reading arg." << std::endl;
    exit(1);
  }
  std::string line, read;

  getline(is_arg,line);
  getline(is_arg,line);

  for(int i = 0; i <= first_section; i++){
    assert(getline(is_arg,line));
  }

  int i = 0;
  read.clear();
  while(line[i] != ':'){
    read += line[i];
    i++;
  }
  i += 2;
  (*arg.seq.begin()).tree.ReadTree(line.substr(i,line.size()-1), data.N);


  //Infer branch lengths
  InferBranchLengths bl(data);
  for(CorrTrees::iterator it_seq = arg.seq.begin(); it_seq != arg.seq.end(); it_seq++){
    //bl.EM(data,(*it_seq).tree, true); //this is estimating times
    double start_time = time(NULL);
    bl.MCMC(data, (*arg.seq.begin()).tree, seed); //this is estimating times
    double end_time = time(NULL);
    std::cerr << end_time - start_time << " secs." << std::endl;
    
    //bl.GradientEM(data, (*it_seq).tree);
  }

  /*

  std::vector<std::string> poplabels = {"ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU","JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI"};
  std::vector<std::string> pop_of_interest = {"CEU"}; //has to be a subset of poplabels
  std::sort(pop_of_interest.begin(), pop_of_interest.end());

  //read sample
  std::ifstream is("../data/1000GP_Phase3_sub.sample");
  getline(is, line);
  int ind = 0;
  std::vector<int> subpop;
  while(getline(is, line)){

    int i = 0;
    while(line[i] != ' ') i++;
    i++;
    read.clear();
    while(line[i] != ' '){
      read += line[i];
      i++;
    }
    for(std::vector<std::string>::iterator it_pop = pop_of_interest.begin(); it_pop != pop_of_interest.end(); it_pop++){
      if(!read.compare(*it_pop)){
        subpop.push_back(2*ind);
        subpop.push_back(2*ind+1);
        break;
      }
    }
    ind++;

  }
  int N_total = 2*subpop.size()-1;

  Tree subtr;
  (*arg.seq.begin()).tree.GetSubTree(subpop, subtr);
  (*arg.seq.begin()).tree = subtr;

  */
  int N_total = 2*data.N-1;
  //Dump to file
  //arg.Dump(filename);
  (*arg.seq.begin()).tree.WriteNewick(options["output"].as<std::string>() + ".newick");

  for(int i = 0; i < N_total; i++){
    (*arg.seq.begin()).tree.nodes[i].branch_length = (*arg.seq.begin()).tree.nodes[i].SNP_end - (*arg.seq.begin()).tree.nodes[i].SNP_begin;
    assert(!std::isnan((*arg.seq.begin()).tree.nodes[i].branch_length));
  }
  (*arg.seq.begin()).tree.WriteNewick(options["output"].as<std::string>() + "_lifespan.newick");

  for(int i = 0; i < N_total; i++){
    (*arg.seq.begin()).tree.nodes[i].branch_length = (*arg.seq.begin()).tree.nodes[i].num_events;
    assert(!std::isnan((*arg.seq.begin()).tree.nodes[i].branch_length));
  }
  (*arg.seq.begin()).tree.WriteNewick(options["output"].as<std::string>() + "_events.newick"); 

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6) << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
  std::cerr << "------------------------------------------------------" << std::endl << std::endl;

  return 0;
}
