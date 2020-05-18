//Finalizing

#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>
#include <string>

#include "cxxopts.hpp"
#include "filesystem.hpp"
#include "collapsed_matrix.hpp"
#include "data.hpp"

int Clean(cxxopts::Options& options){

  bool help = false;
  if(!options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "This function will attempt to delete all temporary files created by Relate. Use when Relate crashes." << std::endl;
    exit(0);
  }


  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Cleaning directory..." << std::endl;

  std::string file_out = options["output"].as<std::string>() + "/";

  int N, L, num_chunks;
  FILE* fp = fopen((file_out + "parameters.bin").c_str(), "r");
  if(fp == NULL){
    std::cerr << "Cannot delete files. Please delete temporary files manually." << std::endl;
    exit(1);
  }
  fread(&N, sizeof(int), 1, fp);
  fread(&L, sizeof(int), 1, fp);
  fread(&num_chunks, sizeof(int), 1, fp);
  fclose(fp); 
 
  /////////////////////////////
  //delete painting and data binaries
  
  std::string filename, output_filename; 
  for(int chunk_index = 0; chunk_index < num_chunks; chunk_index++){
  
    struct stat info;

    std::string dirname = file_out + "chunk_" + std::to_string(chunk_index) + "/";
    //check if directory exists
    if( stat( dirname.c_str(), &info ) == 0 ){

      int N, L, num_windows;
      fp = fopen((file_out + "parameters_c" + std::to_string(chunk_index) + ".bin").c_str(), "r");
      if(fp != NULL){

        fread(&N, sizeof(int), 1, fp);
        fread(&L, sizeof(int), 1, fp);
        fread(&num_windows, sizeof(int), 1, fp);
        fclose(fp);
        num_windows--;

        for(int i = 0; i < num_windows; i++){
          filename        = dirname + "equivalent_branches_" + std::to_string(i) + ".bin";
          output_filename = dirname + options["output"].as<std::string>();
          std::remove(filename.c_str());
          std::remove((output_filename + "_" + std::to_string(i) + ".anc").c_str());
          std::remove((output_filename + "_" + std::to_string(i) + ".mut").c_str());
        }


        //check if directory exists
        if( stat( (file_out + "chunk_" + std::to_string(chunk_index) + "/paint/").c_str(), &info ) == 0 ){
          //paint/ exists so delete it.  
          char painting_filename[32];
          for(int w = 0; w < num_windows; w++){
            snprintf(painting_filename, sizeof(char) * 32, "%s_%i.bin", (file_out + "chunk_" + std::to_string(chunk_index) + "/paint/relate").c_str(), w);
            std::remove(painting_filename);
          }
        }

      }

      std::string file_prefix = dirname + options["output"].as<std::string>();
      std::remove((file_prefix + "_c" + std::to_string(chunk_index) + ".mut").c_str());
      std::remove((file_prefix + "_c" + std::to_string(chunk_index) + ".anc").c_str());

    }
 
    std::remove((file_out + "chunk_" + std::to_string(chunk_index) + ".hap").c_str());
    std::remove((file_out + "chunk_" + std::to_string(chunk_index) + ".r").c_str());
    std::remove((file_out + "chunk_" + std::to_string(chunk_index) + ".rpos").c_str());
    std::remove((file_out + "chunk_" + std::to_string(chunk_index) + ".bp").c_str());
    std::remove((file_out + "parameters_c" + std::to_string(chunk_index) + ".bin").c_str());

  }

  filesys f;
  for(int c = 0; c < num_chunks; c++){

    struct stat info;

    //now delete directories
    if( stat( (file_out + "chunk_" + std::to_string(c) + "/paint/").c_str() , &info ) == 0 ) f.RmDir( (file_out + "chunk_" + std::to_string(c) + "/paint/").c_str() );
    if( stat( (file_out + "chunk_" + std::to_string(c) + "/").c_str() , &info ) == 0 ) f.RmDir( (file_out + "chunk_" + std::to_string(c) + "/").c_str() );
  }

  std::remove((file_out + "parameters.bin").c_str());
  std::remove((file_out + "props.bin").c_str());

  struct stat info;
  if( stat( (file_out).c_str() , &info ) == 0 ) f.RmDir( (file_out).c_str() );

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
