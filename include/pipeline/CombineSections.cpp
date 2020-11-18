//Combining all .anc and .mut files into one.

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "mutations.hpp"

#include <sys/types.h> // required for stat.h
#include <sys/stat.h> // no clue why required -- man pages say so

int CombineSections(cxxopts::Options& options, int chunk_index = 0){

  bool help = false;
  if(!options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: output. Optional: effectiveN (should be consistent across Relate run)" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use after InferBranchLengths to combine files containing trees in small chunks to one file for a section." << std::endl;
    exit(0);
  }

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Combining AncesTrees in Sections..." << std::endl;

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

  //////////////////////////////////
  //Parse Data

  Data data(N, L, Ne); //struct data is defined in data.hpp
  const std::string dirname = file_out + "chunk_" + std::to_string(chunk_index) + "/";
  const std::string output_file = dirname + options["output"].as<std::string>();

  ///////////////////////////////////////// Combine AncesTrees /////////////////////////

  AncesTree anc, anc_tmp;

  int num_tree = 0;
  std::string filename;
 
  //read ancs and splice together
  CorrTrees::iterator it_anc = anc.seq.end();
  for(int i = 0; i < num_windows; i++){
    it_anc = anc.seq.end();
    filename = output_file + "_" + std::to_string(i) + ".anc";
    anc_tmp.ReadBin(filename);
    anc.seq.splice(it_anc, anc_tmp.seq);
    num_tree += anc_tmp.seq.size();
  }
  anc.N = data.N;
  anc.L = num_tree; 

  ///////////////////////////////////////// Create Mutations File /////////////////////////

  std::vector<std::string> mut_filenames(num_windows);
  for(int i = 0; i < num_windows; i++){
    mut_filenames[i] = output_file + "_" + std::to_string(i) + ".mut";
  }
  Mutations mutations(data);

  mutations.ReadShortFormat(mut_filenames);
  mutations.GetAge(anc);
 
  //////////////////////////////////////// Output //////////////////////////////////
  //Dump AncesTree to File

  anc.DumpBin(output_file + "_c" + std::to_string(chunk_index) + ".anc");
  mutations.DumpShortFormat(output_file + "_c" + std::to_string(chunk_index) + ".mut");

  //////////////////////////////////////// Delete tmp files //////////////////////////////////

  for(int i = 0; i < num_windows; i++){
    std::remove((output_file + "_" + std::to_string(i) + ".anc").c_str());
    std::remove((output_file + "_" + std::to_string(i) + ".mut").c_str());
  }
  std::remove((file_out + "chunk_" + std::to_string(chunk_index) + ".bp").c_str());
  std::remove((file_out + "parameters_c" + std::to_string(chunk_index) + ".bin").c_str());

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
