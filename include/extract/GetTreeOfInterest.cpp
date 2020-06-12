#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"
#include <ctime>


void
GetTreeOfInterest(cxxopts::Options& options){
  
  //////////////////////////////////
  //Program options
  
  bool help = false;
  if(!options.count("anc") || !options.count("mut") || (!options.count("bp_of_interest") && (!options.count("first_bp") || !options.count("last_bp")) ) || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, output, bp_of_interest or (first_bp and last_bp). Optional: years_per_gen (Default: 28)" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Outputs tree of interest as .newick file." << std::endl;
    exit(0);
  }  
 
	int first_bp, last_bp;
	if(options.count("bp_of_interest")){
    first_bp = options["bp_of_interest"].as<int>();
		last_bp  = first_bp;
	}else if(options.count("first_bp") && options.count("last_bp")){
    first_bp = options["first_bp"].as<int>();
		last_bp  = options["last_bp"].as<int>();
	}else{
    std::cerr << "Error: need either --bp_of_interest or both --first_bp and --last_bp" << std::endl;
		exit(1);
	}

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Get tree from " << options["anc"].as<std::string>() << " in region [" << first_bp << "," << last_bp  << "]..." << std::endl;

	double years_per_gen = 28;
	if(options.count("years_per_gen")) years_per_gen = options["years_per_gen"].as<float>();

  //////////////////////////////////
  //Parse Data
  AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0;
  Data data(N,L);

  int i; 
  std::string line, read;

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());
  int index_of_first_bp = -1;
  for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
		index_of_first_bp++;
		if((*it_mut).pos >= first_bp) break;
  }
  if(index_of_first_bp == -1){
    std::cerr << "BP not covered by anc/mut" << std::endl;	
    exit(1);
  }
  int tree_index_start = mut.info[index_of_first_bp].tree;

	int index_of_last_bp = index_of_first_bp;
	if(last_bp > first_bp && it_mut != mut.info.end()){
		if((*it_mut).pos < last_bp){
			for(; it_mut != mut.info.end(); it_mut++){
				index_of_last_bp++;
				if((*it_mut).pos >= last_bp) break;
			}
			if(it_mut == mut.info.end()) index_of_last_bp = mut.info.size() - 1;
		}
	}
  int tree_index_end = mut.info[index_of_last_bp].tree;

	std::string filename = options["output"].as<std::string>() + ".newick";
	std::ofstream os(filename);
	std::ofstream os_pos(options["output"].as<std::string>() + ".pos");

  int count_trees = 0;
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  while(num_bases_tree_persists >= 0.0){
 
    if(count_trees >= tree_index_start && count_trees <= tree_index_end){
      //tree is in line
      os_pos << mut.info[mtr.pos].pos << "\n";
			mtr.tree.WriteNewick(os, years_per_gen);
    }
		if(count_trees == tree_index_end) break;
  
    count_trees++;
    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  }
  os.close();
	os_pos.close();

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

