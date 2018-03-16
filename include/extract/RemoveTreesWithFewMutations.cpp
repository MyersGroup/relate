#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>

#include "collapsed_matrix.hpp"
#include "anc.hpp"
#include "cxxopts.hpp"

void
GetDistFromMut(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("mut") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mut, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Extract dist file from mut." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Extracting dist file from " << options["mut"].as<std::string>() << " ... " << std::endl;

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());

 
  FILE* fp_dist = fopen((options["output"].as<std::string>() + ".dist").c_str(), "w");
  fprintf(fp_dist, "#pos dist\n");
  for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end();){
    fprintf(fp_dist, "%d %d\n", (*it_mut).pos, (*it_mut).dist);
    it_mut++;
  }
  fclose(fp_dist);


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
RemoveTreesWithFewMutations(cxxopts::Options& options){

  //////////////////////////////////
  //Program options
  
  bool help = false;
  if(!options.count("threshold") || !options.count("anc") || !options.count("mut") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: threshold, anc, mut, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Estimate population size using coalescent rate." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Removing trees with few mutations from " << options["mut"].as<std::string>() << " ... " << std::endl;
    
  int N, num_trees;

  igzstream is_N(options["anc"].as<std::string>());
  if(is_N.fail()) is_N.open(options["anc"].as<std::string>() + ".gz");
  if(is_N.fail()){
    std::cerr << "Error opening .anc file" << std::endl;
    exit(1);
  }
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  int L = 0;
  igzstream is_L(options["mut"].as<std::string>());
  if(is_L.fail()) is_L.open(options["mut"].as<std::string>() + ".gz");
  if(is_L.fail()){
    std::cerr << "Error opening .mut file" << std::endl;
    exit(1);
  }
  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }
  is_L.close();
 
  Data data(N,L);

  AncesTree anc;
  anc.Read(options["anc"].as<std::string>());

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());

  Mutations mut_subset;
  mut_subset.header = mut.header;

  int num_mutations_threshold = options["threshold"].as<int>(), num_mutations;
  //count the number of mutatinos on each tree.
  
  CorrTrees::iterator it_anc = anc.seq.begin();
  std::vector<Node>::iterator it_node;
  int num_tree_before = mut.info[0].tree, num_tree_after = 0, snp = 0;
  for(; it_anc != anc.seq.end();){ 
    num_mutations = 0;
    it_node = (*it_anc).tree.nodes.begin();
    for(; it_node != (*it_anc).tree.nodes.end(); it_node++){
      num_mutations += (*it_node).num_events;
    }

    if(num_mutations >= num_mutations_threshold){ 
      //include SNPs to mut_subset
      while(mut.info[snp].tree < num_tree_before){
        snp++;
        if(snp == (int)mut.info.size()) break;
      }
      if(snp == (int)mut.info.size()) break;
      if(mut.info[snp].tree != num_tree_before) std::cerr << mut.info[snp].tree << " " << num_tree_before << std::endl;
      assert(mut.info[snp].tree == num_tree_before);
      while(mut.info[snp].tree == num_tree_before){
        mut.info[snp].tree = num_tree_after;
        mut_subset.info.push_back(mut.info[snp]);
        snp++;
        if(snp == (int)mut.info.size()) break;
      }
      if(snp == (int)mut.info.size()) break;
      num_tree_after++;
      it_anc++;
    }else{
      it_anc = anc.seq.erase(it_anc);
    }
    
    num_tree_before++;
  }


  // dump anc and mut
  anc.Dump(options["output"].as<std::string>() + ".anc");
  mut_subset.Dump(options["output"].as<std::string>() + ".mut");

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
