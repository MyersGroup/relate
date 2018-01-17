//Build topology for trees in regions of 500 SNPs.
// Input: Forward and Backward probabilities calculating in Painting
// Output: .anc and .mut file 

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"

int BuildTopology(cxxopts::Options& options,int chunk_index, int first_section, int last_section){

  //////////////////////////////////
  //Parse Data

  int N, L, num_windows;
  std::vector<int> window_boundaries;
  FILE* fp = fopen(("parameters_c" + std::to_string(chunk_index) + ".bin").c_str(), "r");
  assert(fp != NULL);
  fread(&N, sizeof(int), 1, fp);
  fread(&L, sizeof(int), 1, fp);
  fread(&num_windows, sizeof(int), 1, fp);
  window_boundaries.resize(num_windows);
  fread(&window_boundaries[0], sizeof(int), num_windows, fp);
  fclose(fp);
  num_windows--;

  Data data(("chunk_" + std::to_string(chunk_index) + ".hap").c_str(), ("chunk_" + std::to_string(chunk_index) + ".bp").c_str(), ("chunk_" + std::to_string(chunk_index) + ".r").c_str(), ("chunk_" + std::to_string(chunk_index) + ".rpos").c_str()); //struct data is defined in data.hpp
  data.name = ("chunk_" + std::to_string(chunk_index) + "/paint/relate");
  const std::string dirname = "chunk_" + std::to_string(chunk_index) + "/";
  if(first_section >= num_windows) return 1;

  ///////////////////////////////////////////// Build AncesTree //////////////////////////
  //input:  Data and distance matrix
  //output: AncesTree (tree sequence)

  last_section = std::min(num_windows-1, last_section);

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Estimating topologies of AncesTrees in sections " << first_section << "-" << last_section << "..." << std::endl;

  for(int section = first_section; section <= last_section; section++){

    std::cerr << "[" << section << "/" << last_section << "]\r";
    std::cerr.flush(); 

    AncesTree anc;
    AncesTreeBuilder ancbuilder(data);

    int section_startpos = window_boundaries[section];
    int section_endpos   = window_boundaries[section+1]-1;
    if(section_endpos >= data.L) section_endpos = data.L-1;

    ancbuilder.BuildTopology(section ,section_startpos, section_endpos, data, anc);

    /////////////////////////////////////////// Dump AncesTree to File //////////////////////

    anc.DumpBin(dirname + options["output"].as<std::string>() + "_" + std::to_string(section) + ".anc");
    ancbuilder.mutations.DumpShortFormat(dirname + options["output"].as<std::string>() + "_" + std::to_string(section) + ".mut", section_startpos, section_endpos);

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
