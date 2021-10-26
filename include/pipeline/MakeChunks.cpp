//Make overlapping chunks with fixed number of SNPs from data set
#ifndef MAKE_CHUNKS
#define MAKE_CHUNKS

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "filesystem.hpp"
#include "cxxopts.hpp"
#include "data.hpp"

int MakeChunks(cxxopts::Options& options, int chunk_size = 0){

  bool help = false;
  if(!options.count("haps") || !options.count("sample") || !options.count("map") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, sample, map, output. Optional: memory, dist, transversion." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use to make smaller chunks from the data." << std::endl;
    exit(0);
  }

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Parsing data.." << std::endl;

  bool use_transitions = true;
  if(options.count("transversion")) use_transitions = false;
 
  //////////////////////////////////
  //Parse Data

  struct stat info;
  //check if directory exists
  if( stat( (options["output"].as<std::string>() + "/").c_str(), &info ) == 0 ){
    std::cerr << "Error: Directory " << options["output"].as<std::string>() << " already exists. Relate will use this directory to store temporary files." << std::endl;
    exit(1);
  }
  filesys f;
  f.MakeDir((options["output"].as<std::string>() + "/").c_str());

  Data data;
  /*
  if(options.count("chunk_size")){

    if(options.count("dist")){
      data.MakeChunks(options["haps"].as<std::string>(), options["sample"].as<std::string>(), options["map"].as<std::string>(), options["dist"].as<std::string>(), options["chunk_size"].as<int>());
    }else{
      data.MakeChunks(options["haps"].as<std::string>(), options["sample"].as<std::string>(), options["map"].as<std::string>(), "unspecified", options["chunk_size"].as<int>());
    }

  }else{

  }
  */

  if(options.count("memory")){
    if(options.count("dist")){
      data.MakeChunks(options["haps"].as<std::string>(), options["sample"].as<std::string>(), options["map"].as<std::string>(), options["dist"].as<std::string>(), options["output"].as<std::string>(), use_transitions, options["memory"].as<float>());
    }else{
      data.MakeChunks(options["haps"].as<std::string>(), options["sample"].as<std::string>(), options["map"].as<std::string>(), "unspecified", options["output"].as<std::string>(), use_transitions, options["memory"].as<float>());
    }
  }else{
    if(options.count("dist")){
      data.MakeChunks(options["haps"].as<std::string>(), options["sample"].as<std::string>(), options["map"].as<std::string>(), options["dist"].as<std::string>(), options["output"].as<std::string>(), use_transitions);
    }else{
      data.MakeChunks(options["haps"].as<std::string>(), options["sample"].as<std::string>(), options["map"].as<std::string>(), "unspecified", options["output"].as<std::string>(), use_transitions);
    }
  }
  
  std::vector<double> sample_ages(data.N);
  if(options.count("sample_ages")){
    igzstream is_ages(options["sample_ages"].as<std::string>());
		if(is_ages.fail()){
			std::cerr << "Warning: unable to open sample ages file" << std::endl;
		}else{
			int i = 0; 
			while(is_ages >> sample_ages[i]){
				i++;
				//sample_ages[i] = sample_ages[i-1];
				//i++;
				if(i == data.N) break;
			}
			if(i != data.N){
				std::cerr << "Warning: sample ages specified but number of samples does not agree with other input files. Ignoring sample ages." << std::endl;
				sample_ages.clear();
			}
		}
  }else{

    sample_ages.clear(); 
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

#endif //MAKE_CHUNKS
