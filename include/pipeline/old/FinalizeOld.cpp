//Finalizing

#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>
#include <string>

#include "cxxopts.hpp"
#include "collapsed_matrix.hpp"
#include "data.hpp"
#include "arg.hpp"
#include "arg_builder.hpp"

int Finalize(cxxopts::Options& options){

  bool help = false;
  if(!options.count("seq") || !options.count("prep") || !options.count("chunk_size") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use at the end to finalize results. This will summarize all sections into one file." << std::endl;
    exit(0);
  }

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing..." << std::endl;

  int N, L, num_chunks;
  std::ifstream is_param(options["seq"].as<std::string>());
  if(is_param.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  } 
  is_param >> N >> L >> num_chunks;
  is_param.close();
  
  int chunk_size = options["chunk_size"].as<int>();

  int i, j; 
  std::string line, line2, read;

  ///////////////////////////////////////// Combine Args /////////////////////////

  //Merge files to one. 

  std::ifstream is;
  std::ofstream os;

  int num_flips = 0, num_non_mapping = 0;
  int delta_chunk = chunk_size-10000;
  int num_trees_cum = 0, first_tree = 0;

  ///////////////////////////////////////// Combine Mutation Files /////////////////////////

  std::ifstream is_prep(options["prep"].as<std::string>());
  if(is_prep.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  } 
  getline(is_prep,line2);
  
  os.open("result/" + options["output"].as<std::string>() + ".mut");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  } 

  //Header-Legend
  i = 0;
  while(line2[i] != ';') i++;
  i++;
  while(line2[i] != ';') i++;
  i++;
  while(line2[i] != ';') i++;
  os << line2.substr(0,i) << ";";
  os << "tree_index;branch_indices;is_mapping;is_flipped;(age_begin,age_end);";
  os << line2.substr(i+1,line2.size()-1) << "\n";

  //Rest
  for(int c = 0; c < num_chunks; c++){

    is.open("results_c" + std::to_string(c) + "/" + options["output"].as<std::string>() + ".mut");
    if(is.fail()){
      std::cerr << "Error while opening file." << std::endl;
      exit(1);
    }  
    getline(is,line);

    if(c > 0){
      for(int snp = 0; snp < 10000; snp++){
        getline(is,line);
      }
    }

    int num_trees_chunk = 0;

    //read in end_chunk lines
    //need to record num_flips and num_non_mapping (3rd,4th entry)
    int snp = 0;
    while(snp < delta_chunk && getline(is, line)){
      
      read.clear();
      i = 0;
      while(line[i] != ';'){
        read += line[i];
        i++;
      }
      j = i;
      i++;
      if(snp == 0){
        num_trees_chunk = std::stoi(read);
        first_tree = num_trees_chunk;
      }else if(std::stoi(read) > num_trees_chunk){
        num_trees_chunk++;
      }
      read.clear();

      while(line[i] != ';') i++;
      i++;
      if(line[i] == '1') num_non_mapping++;
      i+=2;
      
      
      if(line[i] == '1') num_flips++;
      line.erase(0,j);

      getline(is_prep, line2);

      i = 0;
      while(line2[i] != ';') i++;
      i++;
      while(line2[i] != ';') i++;
      i++;
      while(line2[i] != ';') i++;

      os << line2.substr(0,i) << ";";
      os << num_trees_chunk + num_trees_cum - first_tree << line;
      os << line2.substr(i+1,line2.size()-1) << "\n";

      snp++;
    }

    num_trees_cum += num_trees_chunk - first_tree + 1;
    delta_chunk = chunk_size;
    is.close();

  }

  os.close();

  std::cerr << "num_non_mapping " << num_non_mapping << "; num_flips " << num_flips << std::endl;

  //////////////////////////////////////////////////

  os.open("result/" + options["output"].as<std::string>() + ".arg");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  } 

  os << "NUM_HAPLOTYPES " << N << "\n";
  os << "NUM_TREES " << num_trees_cum << "\n";

  int num_trees = 0;
  int start_chunk = 0;
  int end_chunk = chunk_size-10000;
  int overlap_chunk_size = 10000;
  for(int c = 0; c < num_chunks; c++){
 
    is.open("results_c" + std::to_string(c) + "/" + options["output"].as<std::string>() + ".arg");
    if(is.fail()){
      std::cerr << "Error while opening file." << std::endl;
      exit(1);
    } 
    
    getline(is,line);
    getline(is,line);

    int position, prev_position;
    if(c == 0){
      getline(is,line);
      //first tree
      i = 0;
      read.clear();
      while(line[i] != ':'){
        read += line[i];
        i++;
      }
      line.replace(0,i,std::to_string(start_chunk));
      os << line << "\n";
      num_trees++;
    }else{
      while(getline(is,line)){
        //check startpos of tree
        //if it is less than end_chunk, accept
        i = 0;
        read.clear();
        while(line[i] != ':'){
          read += line[i];
          i++;
        }
        position = std::stoi(line);
        if(position > overlap_chunk_size) break;
        line2 = line;
      }
      //first tree
      i = 0;
      read.clear();
      while(line2[i] != ':'){
        read += line2[i];
        i++;
      }
      line2.replace(0,i,std::to_string(overlap_chunk_size + start_chunk));
      os << line2 << "\n";
      num_trees++;
      position = position + start_chunk;
      //second tree
      if(position < end_chunk){
        line.replace(0,i,std::to_string(position));
        os << line << "\n";
        num_trees++;
      }
    }

    //rest
    while(getline(is,line)){
      //check startpos of tree
      //if it is less than end_chunk, accept
      i = 0;
      read.clear();
      while(line[i] != ':'){
        read += line[i];
        i++;
      }
      position = std::stoi(read) + start_chunk;
      
      if(position < end_chunk){
        line.replace(0,i,std::to_string(position));
        os << line << "\n";
        num_trees++;
      }else{ 
        break;
      }
    
    }

    start_chunk = end_chunk - overlap_chunk_size;
    end_chunk  += chunk_size;
    //if(end_chunk > L) end_chunk = L;
    is.close();

  }


  os.close();

  assert(num_trees == num_trees_cum);

  /////////////////////////////////////////////
  //Resource Usage

  rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6) << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
  std::cerr << "Output written to result/" << options["output"].as<std::string>() + ".arg and result/" << options["output"].as<std::string>() + ".mut." << std::endl;
  std::cerr << "------------------------------------------------------" << std::endl << std::endl;

  return 0;
}
