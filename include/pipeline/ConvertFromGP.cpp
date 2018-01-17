/*! \file ConvertFromGP.cpp 
 *  \brief Convert Data from GP format to custom format
 */

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "cxxopts.hpp"
#include "data.hpp"

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options  
  cxxopts::Options options("ConvertFromGP");
  options.add_options()
    ("help", "Print help")
    ("h,haps", "Filename of file containing haplotype data", cxxopts::value<std::string>()) //!option
    ("l,legend", "Filename of file containing legend", cxxopts::value<std::string>())
    ("m,map", "Filename of file containing recombination map", cxxopts::value<std::string>())
    ("f,fasta", "Filename of file containing ancestral genome as fasta", cxxopts::value<std::string>())
    ("s,samples", "Filename of file containing sample labels", cxxopts::value<std::string>())
    ("x,excluded_samples", "Filename of file containing excluded samples", cxxopts::value<std::string>())
    ("a,ancestral_state", "Filename of file containing ancestral states as fasta", cxxopts::value<std::string>())
    ("c,mask", "Filename of file containing mask as fasta", cxxopts::value<std::string>());
  
  options.parse(argc, argv);
  bool help = false;
  if(!options.count("haps") || !options.count("legend")  || !options.count("map") || !options.count("fasta") || !options.count("samples") || !options.count("excluded_samples") || !options.count("ancestral_state") || !options.count("mask")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, legend, map, fasta, samples, excluded_samples, ancestral_state, mask." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use to convert 1000 GP data to file format needed by Relate." << std::endl;
    exit(0);
  }

  std::cerr << "############" << std::endl;
  std::cerr << "Preparing Data..." << std::endl;

  //////////////////////////////////
  //Parse Data

  GPData data; //struct data is defined in data.hpp
  data.ReadGP(options["haps"].as<std::string>(),options["legend"].as<std::string>(),options["map"].as<std::string>(),options["ancestral_state"].as<std::string>(), options["mask"].as<std::string>(), options["excluded_samples"].as<std::string>());
  data.PrepareMutationsFile(options["fasta"].as<std::string>(), options["samples"].as<std::string>(), options["excluded_samples"].as<std::string>());

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
