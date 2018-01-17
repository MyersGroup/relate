//All against all painting.

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

  /*
  int L = 10;
  int N = 10;  
  for(int i = 0; i < L; i++){
    for(int j = 0; j < N; j++){
      std::cerr << data.sequence[i][j] << " ";
    }
    std::cerr << std::endl;
  }
  for(int i = 0; i < L; i++){
    std::cerr << data.pos[i] << " " << data.rpos[i] << " " << data.r[i] << std::endl;
  }

  std::cerr << data.sequence.size() << " " << data.sequence.subVectorSize(0) << std::endl;
  std::cerr << data.pos.size() << std::endl;
  std::cerr << data.rpos.size() << std::endl;
  std::cerr << data.r.size() << std::endl;
  */

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Painting sequences..." << std::endl;

  //create directory called paint/ if not existent
  filesys f;
  f.MakeDir(("chunk_" + std::to_string(chunk_index) + "/").c_str());
  f.MakeDir(("chunk_" + std::to_string(chunk_index) + "/paint/").c_str());

  //////////////////////////////////////////// Paint sequence ////////////////////////////

  char filename[32];
  std::vector<FILE*> pfiles(num_windows);
  for(int w = 0; w < num_windows; w++){
    snprintf(filename, sizeof(char) * 32, "%s_%i.bin", data.name.c_str(), w);
    pfiles[w] = fopen(filename, "wb");
    assert(pfiles[w] != NULL);
  }

  for(int hap = 0; hap < data.N; hap++){
    //std::cerr << hap << std::endl;
    FastPainting painter(data);
    painter.PaintSteppingStones(data, window_boundaries, pfiles, hap);
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
