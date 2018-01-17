#include <iostream>
#include "collapsed_matrix.hpp"
#include "cxxopts.hpp"

#include <ctime>
#include <tgmath.h>

int SummarizeCoalescentRateForChr(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("index_first_tree") || !options.count("index_last_tree") || !options.count("trees_per_section") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: index_first_tree, index_last_tree, trees_per_section, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .anc file and calculates pairwise coalescent rate for trees with index [index_first_tree, index_last_tree]. Output is bin file. Use SummarizeCoalesecntRate to obtain coalescent rates." << std::endl;
    exit(0);
  }  

  //calculate epoche times
  int num_epoches = 40;

  //get filenames
  int index_first_tree     = options["index_first_tree"].as<int>(); 
  int trees_per_section     = options["trees_per_section"].as<int>();
  int index_last_tree       = options["index_last_tree"].as<int>();

  std::vector<std::string> filenames;
  std::string filename_base = options["output"].as<std::string>();

  int index_index_first_tree = index_first_tree, index_index_last_tree = index_first_tree;
  for(; index_index_first_tree < index_last_tree;){
    index_index_last_tree    = index_index_first_tree + trees_per_section - 1;
    if(index_last_tree - index_index_last_tree < trees_per_section) index_index_last_tree = index_last_tree;
    filenames.push_back(filename_base + "_" + std::to_string(index_index_first_tree) + "_" + std::to_string(index_index_last_tree) + ".bin");
    index_index_first_tree  = index_index_last_tree + 1;
  }

  //populate coalescent_rate_data
  std::vector<CollapsedMatrix<float>> coalescent_rate_data(num_epoches);
  std::vector<float>::iterator it_coalescent_rate_data;

  FILE* fp = fopen(filenames[0].c_str(),"rb");
  assert(fp != NULL);
  for(int e = 0; e < num_epoches; e++){
    coalescent_rate_data[e].ReadFromFile(fp); 
  }
  fclose(fp);

  CollapsedMatrix<float> coalescent_rate_data_section;
  for(int i = 1; i < (int) filenames.size(); i++){
    fp = fopen(filenames[i].c_str(),"rb");
    assert(fp != NULL);
    for(int e = 0; e < num_epoches; e++){    
      coalescent_rate_data_section.ReadFromFile(fp);
      it_coalescent_rate_data     = coalescent_rate_data[e].vbegin();
      for(std::vector<float>::iterator it_coalescent_rate_data_section = coalescent_rate_data_section.vbegin(); it_coalescent_rate_data_section != coalescent_rate_data_section.vend();){
        *it_coalescent_rate_data += *it_coalescent_rate_data_section;
        it_coalescent_rate_data_section++;
        it_coalescent_rate_data++;
      }
    }
    fclose(fp);
  }

  //output as bin
  fp = fopen((options["output"].as<std::string>() + ".bin" ).c_str(), "wb");  
  for(std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_data = coalescent_rate_data.begin(); it_coalescent_rate_data != coalescent_rate_data.end();){
    (*it_coalescent_rate_data).DumpToFile(fp);
    it_coalescent_rate_data++;
  }
  fclose(fp);

  return 0;

}
