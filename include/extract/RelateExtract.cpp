#include "GetTreeOfInterest.cpp"
#include "CreateAncesTreeFileForSubpopulation.cpp"
#include "RemoveTreesWithFewMutations.cpp"

#include "cxxopts.hpp"
#include <string>

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("Relate");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("poplabels", "Filename of file containing population labels.", cxxopts::value<std::string>()) 
    ("anc", "Filename of file containing trees.", cxxopts::value<std::string>())
    ("mut", "Filename of file containing mut.", cxxopts::value<std::string>())
    ("pop_of_interest", "Population label. If not specified, use all haplotypes.", cxxopts::value<std::string>())
    ("snp_of_interest", "BP of SNP of interest.", cxxopts::value<int>())
    ("threshold", "Threshold used in RemoveTreesWithFewMutations.", cxxopts::value<int>())
    ("o,output", "Filename of output (excl file extension).", cxxopts::value<std::string>());

  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("TreeAtSNPAsNewick")){
  
    GetTreeOfInterest(options);

  }else if(!mode.compare("SubTreesForSubpopulation")){
  
    CreateAncesTreeFileForSubpopulation(options);

  }else if(!mode.compare("RemoveTreesWithFewMutations")){
  
    RemoveTreesWithFewMutations(options);
    GetDistFromMut(options);

  }else if(!mode.compare("ExtractDistFromMut")){
  
    GetDistFromMut(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "TreeAtSNPAsNewick, SubTreesForSubpopulation, RemoveTreesWithFewMutations, ExtractDistFromMut." << std::endl;
  
  }

  bool help = false;
  if(!options.count("mode")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

  return 0;

}

