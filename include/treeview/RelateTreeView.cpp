#include "Treeview.cpp"

#include "cxxopts.hpp"
#include <string>

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("RelateTreeview");
  options.add_options()
    ("help", "Print help.")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("anc", "Filename of file containing trees.", cxxopts::value<std::string>())
    ("mut", "Filename of file containing mut.", cxxopts::value<std::string>())
    ("haps", "Filename of haps file (Output file format of Shapeit).", cxxopts::value<std::string>())
    ("sample", "Filename of sample file (Output file format of Shapeit).", cxxopts::value<std::string>())
    ("dist", "Filename of dist file", cxxopts::value<std::string>())
    ("mask", "Filename of mask file", cxxopts::value<std::string>())
    ("snp_of_interest", "BP of SNP of interest.", cxxopts::value<int>())
    ("i,input", "Filename of anc and mut files without file extension.", cxxopts::value<std::string>())
    ("o,output", "Filename for updated anc and mut files without file extension.", cxxopts::value<std::string>())
    ("seed", "Seed for MCMC in branch lengths estimation.", cxxopts::value<int>());
  
  options.parse(argc, argv);

  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("TreeView")){

    TreeView(options);

  }else if(!mode.compare("MutationsOnBranches")){

    MutationsOnBranches(options);

  }else if(!mode.compare("BranchesBelowMutation")){

    BranchesBelowMutation(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "TreeView, MutationsOnBranches, BranchesBelowMutation." << std::endl;
  
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

