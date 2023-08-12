//All against all painting.
#ifndef PAINTING
#define PAINTING

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <limits.h>

#include <string> // required for std::string

#include "cxxopts.hpp"
#include "filesystem.hpp"
#include "data.hpp"
#include "fast_painting.hpp"

int Paint(cxxopts::Options& options, int chunk_index){

  std::string file_out = options["output"].as<std::string>() + "/"; 

  int N, L, num_windows;
  std::vector<int> window_boundaries;
  FILE* fp = fopen((file_out + "parameters_c" + std::to_string(chunk_index) + ".bin").c_str(), "r");
  assert(fp != NULL);
  fread(&N, sizeof(int), 1, fp);
  fread(&L, sizeof(int), 1, fp);
  fread(&num_windows, sizeof(int), 1, fp);
  window_boundaries.resize(num_windows);
  fread(&window_boundaries[0], sizeof(int), num_windows, fp);
  fclose(fp);
  num_windows--;

  bool use_transitions = true;

	Data data((file_out + "chunk_" + std::to_string(chunk_index) + ".hap").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".bp").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".dist").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".r").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".rpos").c_str(),  (file_out + "chunk_" + std::to_string(chunk_index) + ".state").c_str()); //struct data is defined in data.hpp 
  data.name = (file_out + "chunk_" + std::to_string(chunk_index) + "/paint/relate");

  if(options.count("painting")){

	  std::string val;
		std::string painting = options["painting"].as<std::string>();
		int i = 0;
		for(;i < painting.size(); i++){
      if(painting[i] == ',') break;
			val += painting[i];
		}
		data.theta = std::stof(val);
		data.ntheta = 1.0 - data.theta;

		i++;
		val.clear();
		for(;i < painting.size(); i++){
			if(painting[i] == ',') break;
			val += painting[i];
		}
    double rho = std::stof(val);
		for(int l = 0; l < (int)data.r.size(); l++){
			data.r[l] *= rho;
		}

	}

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Painting sequences..." << std::endl;

  //create directory called paint/ if not existent
  filesys f;
  f.MakeDir((file_out + "chunk_" + std::to_string(chunk_index) + "/").c_str());
  f.MakeDir((file_out + "chunk_" + std::to_string(chunk_index) + "/paint/").c_str());

  //////////////////////////////////////////// Paint sequence ////////////////////////////

  char filename[1024];
  std::vector<FILE*> pfiles(num_windows);
  for(int w = 0; w < num_windows; w++){
    snprintf(filename, sizeof(char) * 1024, "%s_%i.bin", data.name.c_str(), w);
    pfiles[w] = fopen(filename, "wb");
    assert(pfiles[w] != NULL);
  }

  for(int hap = 0; hap < data.N; hap++){
    //std::cerr << hap << std::endl;
    FastPainting painter(data);
    painter.PaintSteppingStones(data, window_boundaries, pfiles, hap);
    //painter.PaintSteppingStonesNormal(data, window_boundaries, pfiles, hap);

  }

  for(int w = 0; w < num_windows; w++){  
    fclose(pfiles[w]);
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

#endif //PAINTING
