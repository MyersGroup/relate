#include <iostream>
#include "collapsed_matrix.hpp"
#include "cxxopts.hpp"

#include <ctime>
#include <tgmath.h>

int SummarizeCoalescentRateForGenome(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( ((!options.count("first_chr") || !options.count("last_chr")) && !options.count("chr")) || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: chr or (first_chr, last_chr), output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Summarizes .bin files for chromosomes." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Summarizing over chromosomes..." << std::endl;  

  //calculate epoch times
	std::vector<std::string> filenames;
	std::string filename_base = options["output"].as<std::string>();
	if(options.count("chr")){

		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		std::string line;
		while(getline(is_chr, line)){
			filenames.push_back(filename_base + "_chr" + line + ".bin");
		}
		is_chr.close();

	}else{

    int start     = options["first_chr"].as<int>(); 
    int end       = options["last_chr"].as<int>();
		for(int chr = start; chr <= end; chr++){
			filenames.push_back(filename_base + "_chr" + std::to_string(chr) + ".bin");
		}

	}

  //populate coalescent_rate_data
  FILE* fp = fopen(filenames[0].c_str(),"rb");
  assert(fp != NULL);

  int num_epochs;
  std::vector<float> epochs;
  fread(&num_epochs, sizeof(int), 1, fp);
  epochs.resize(num_epochs);
  fread(&epochs[0], sizeof(float), num_epochs, fp);

  std::vector<CollapsedMatrix<float>> coalescent_rate_data(num_epochs);
  for(int e = 0; e < num_epochs; e++){
    coalescent_rate_data[e].ReadFromFile(fp); 
  }
  fclose(fp);
  std::remove(filenames[0].c_str());

  std::vector<float>::iterator it_coalescent_rate_data;
  CollapsedMatrix<float> coalescent_rate_data_section;
  for(int i = 1; i < (int) filenames.size(); i++){
    fp = fopen(filenames[i].c_str(),"rb");
    assert(fp != NULL);
    fread(&num_epochs, sizeof(int), 1, fp); 
    fread(&epochs[0], sizeof(float), num_epochs, fp);

    for(int e = 0; e < num_epochs; e++){    
      coalescent_rate_data_section.ReadFromFile(fp);
      it_coalescent_rate_data     = coalescent_rate_data[e].vbegin();
      for(std::vector<float>::iterator it_coalescent_rate_data_section = coalescent_rate_data_section.vbegin(); it_coalescent_rate_data_section != coalescent_rate_data_section.vend();){
        *it_coalescent_rate_data += *it_coalescent_rate_data_section;
        it_coalescent_rate_data_section++;
        it_coalescent_rate_data++;
      }
    }
    fclose(fp);
    std::remove(filenames[i].c_str());
  }

  //output as bin
  fp = fopen((options["output"].as<std::string>() + ".bin").c_str(), "wb");  

  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(float), epochs.size(), fp);
  for(std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_data = coalescent_rate_data.begin(); it_coalescent_rate_data != coalescent_rate_data.end();){
    (*it_coalescent_rate_data).DumpToFile(fp);
    it_coalescent_rate_data++;
  }
  fclose(fp);


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
