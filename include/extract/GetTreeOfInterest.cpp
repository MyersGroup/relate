#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.h"
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
  if(!options.count("anc") || !options.count("mut") || !options.count("snp_of_interest")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, snp_of_interest." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Outputs tree of interest as .newick file." << std::endl;
    exit(0);
  }  
  
  int snp_of_interest = options["snp_of_interest"].as<int>();
  int index_of_snp_of_interest;

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Get tree from " << options["anc"].as<std::string>() << " at BP " << snp_of_interest << "..." << std::endl;


  //////////////////////////////////
  //Parse Data
  int N;
  igzstream is_N(options["anc"].as<std::string>());
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  int i; 
  std::string line, read;

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());
  index_of_snp_of_interest = 0;
  for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
    if((*it_mut).pos == snp_of_interest) break;
    index_of_snp_of_interest++;
  }
  int tree_index_of_interest = mut.info[index_of_snp_of_interest].tree;

  igzstream is_anc(options["anc"].as<std::string>());
  if(is_anc.fail()){
    std::cerr << "Error while reading anc." << std::endl;
    exit(1);
  }

  getline(is_anc,line);
  getline(is_anc,line);

  MarginalTree mtr;

  int count_trees = 0;
  while(getline(is_anc,line)){
 
    if(count_trees == tree_index_of_interest){
      
      //tree is in line
      mtr.Read(line, N);

      for(int k = 0; k < (int) mtr.tree.nodes.size(); k++){
        mtr.tree.nodes[k].branch_length *= 28.0;
      }
      mtr.tree.WriteNewick("tree_at_" + std::to_string(mtr.pos) + ".newick");

      for(int i = 0; i < (int) mtr.tree.nodes.size(); i++){
        mtr.tree.nodes[i].branch_length = mtr.tree.nodes[i].SNP_end - mtr.tree.nodes[i].SNP_begin;
        assert(!std::isnan(mtr.tree.nodes[i].branch_length));
      }
      mtr.tree.WriteNewick("tree_at_" + std::to_string(mtr.pos) + "_lifespan.newick");
      
      for(int i = 0; i < (int) mtr.tree.nodes.size(); i++){
        mtr.tree.nodes[i].branch_length = mtr.tree.nodes[i].num_events;
        assert(!std::isnan(mtr.tree.nodes[i].branch_length));
      }
      mtr.tree.WriteNewick("tree_at_" + std::to_string(mtr.pos) + "_events.newick"); 

      break; 

    }
  
    count_trees++;
  
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

