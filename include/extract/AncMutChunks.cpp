#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include "gzstream.hpp"
#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"

void
DivideAncMut(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("threads") || !options.count("anc") || !options.count("mut") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: threads, anc, mut, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Dividing .anc/.mut files into smaller files for parallelization." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Dividing .anc/.mut files into smaller files for parallelization..." << std::endl;

  int N, num_trees;

  igzstream is_N(options["anc"].as<std::string>());
  if(is_N.fail()) is_N.open(options["anc"].as<std::string>() + ".gz");
  if(is_N.fail()){
    std::cerr << "Error opening .anc file" << std::endl;
    exit(1);
  }
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.ignore(256, ' ');
  is_N >> num_trees;
  is_N.close();
  const int num_trees_check = num_trees;
    
  int L = 0;
  igzstream is_L(options["mut"].as<std::string>());
  if(is_L.fail()) is_L.open(options["mut"].as<std::string>() + ".gz");
  if(is_L.fail()){
    std::cerr << "Error opening .mut file" << std::endl;
    exit(1);
  }
  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }
  is_L.close();

  Data data(N,L);

  int num_threads = options["threads"].as<int>();
  //divide anc/mut into chunks with ((int) num_trees/(100.0 + num_threads)) + 1 trees each
  int num_trees_per_chunk = ((int) num_trees/(5.0 * num_threads)) + 1;
  if(num_trees_per_chunk < 10) num_trees_per_chunk = 10;

  std::string line;
  igzstream is(options["anc"].as<std::string>());
  if(is.fail()) is.open(options["anc"].as<std::string>() + ".gz");
  if(is.fail()){
    std::cerr << "Error opening .anc file" << std::endl;
    exit(1);
  }
  assert(getline(is,line));
  assert(getline(is,line));
  igzstream is_mut(options["mut"].as<std::string>());
  if(is_mut.fail()) is_mut.open(options["mut"].as<std::string>() + ".gz");
  if(is_mut.fail()){
    std::cerr << "Error opening .mut file" << std::endl;
    exit(1);
  }
  std::string header;
  assert(getline(is_mut, header));

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());

  int i = 0, snp = 0, tree_index = mut.info[snp].tree;

  while(num_trees > num_trees_per_chunk + 10){
    ogzstream os(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".anc.gz");
    ogzstream os_mut(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".mut.gz");

    os << "NUM_HAPLOTYPES " << data.N << "\n";
    os << "NUM_TREES " << num_trees_per_chunk << "\n";
    os_mut << header << "\n";

    for(int k = 0; k < num_trees_per_chunk; k++){
      assert(getline(is, line));
      os << line << "\n";

      if(snp < data.L){
        while(mut.info[snp].tree == tree_index){
          assert(getline(is_mut,line));
          os_mut << line << "\n";
          snp++;
          if(snp >= data.L) break;
        }
      }else{
        std::cerr << "Mutation file does not seem to contain all SNPs.\n";
        exit(1);
      }
      tree_index++;
    }

    os.close();
    os_mut.close();
    num_trees -= num_trees_per_chunk;
    i++;
  }

  //last chunk
  ogzstream os(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".anc.gz");
  ogzstream os_mut(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".mut.gz");

  os << "NUM_HAPLOTYPES " << data.N << "\n";
  os << "NUM_TREES " << num_trees << "\n";
  os_mut << header << "\n";

  while(getline(is, line)){
    os << line << "\n";
    if(snp < data.L){
      while(mut.info[snp].tree == tree_index){
        assert(getline(is_mut,line));
        os_mut << line << "\n";
        snp++;
        if(snp >= data.L) break;
      }
    }else{
      std::cerr << "Mutation file does not seem to contain all SNPs.\n";
      exit(1);
    }
    tree_index++;
  }

  assert(tree_index - mut.info[0].tree == num_trees_check);
  os.close();
  os_mut.close();

  std::ofstream os_param(options["output"].as<std::string>() + ".param");
  if(os_param.fail()){
    std::cerr << "Error while opening " << options["output"].as<std::string>() + ".param" << std::endl;
    exit(1);
  }
  os_param << "NUM_HAPLOTYPES NUM_SNPS NUM_TREES NUM_CHUNKS\n";
  os_param << data.N << " " << data.L << " " << num_trees_check << " " << i + 1 << std::endl;
  os_param.close();

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

}

void
CombineAncMut(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Combining .anc/.mut files into one file." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Combining .anc/.mut files into one file..." << std::endl;

  std::string line;
  Data data;
  int num_trees, num_chunks;

  std::ifstream is_param(options["output"].as<std::string>() + ".param");
  if(is_param.fail()){
    std::cerr << "Unable to combine anc/mut files because file " << options["output"].as<std::string>() + ".param" << " has been lost." << std::endl;
    std::cerr << "You can manually create this file which has the following header followed by the corresponding entries on the next line" << std::endl;
    std::cerr << "NUM_HAPLOTYPES NUM_SNPS NUM_TREES NUM_CHUNKS" << std::endl;
    std::cerr << std::endl;
  }
  assert(getline(is_param, line));
  assert(getline(is_param, line));
  sscanf(line.c_str(), "%d %d %d %d", &data.N, &data.L, &num_trees, &num_chunks);
  is_param.close();

  ogzstream os(options["output"].as<std::string>() + ".anc.gz");
  ogzstream os_mut(options["output"].as<std::string>() + ".mut.gz");

  os << "NUM_HAPLOTYPES " << data.N << "\n";
  os << "NUM_TREES " << num_trees << "\n";

  for(int i = 0; i < num_chunks; i++){
    bool is_gzipped = false;
    igzstream is(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".anc");
    if(is.fail()){
      is.open(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".anc.gz");
      is_gzipped = true;
    }
    if(is.fail()){
      std::cerr << "Error opening .anc file" << std::endl;
      exit(1);
    }
    assert(getline(is,line));
    assert(getline(is,line));
    while(getline(is, line)){
      os << line << "\n";
    }
    is.close();
    if(is_gzipped){
      std::remove((options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".anc.gz").c_str());
    }else{
      std::remove((options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".anc").c_str());
    } 

    is_gzipped = false;
    igzstream is_mut(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".mut");
    if(is_mut.fail()){
      is_mut.open(options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".mut.gz");
      is_gzipped = true;
    }
    if(is_mut.fail()){
      std::cerr << "Error opening .mut file" << std::endl;
      exit(1);
    }
    assert(getline(is_mut, line));
    if(i == 0) os_mut << line << "\n";
    while(getline(is_mut, line)){
      os_mut << line << "\n";
    }
    is_mut.close();

    if(is_gzipped){
      std::remove((options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".mut.gz").c_str());
    }else{
      std::remove((options["output"].as<std::string>() + "_chr" + std::to_string(i) + ".mut").c_str());
    }
  }
  os.close();
  os_mut.close();

  std::remove((options["output"].as<std::string>() + ".param").c_str());

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

}

void
AncMutForSubregion(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("anc") || !options.count("mut") || !options.count("first_bp") || !options.count("last_bp") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: anc, mut, first_bp, last_bp, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Extracting .anc/.mut files for subregion." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Extracting .anc/.mut files for subregion " << options["first_bp"].as<int>() << " - " << options["last_bp"].as<int>() << "..." << std::endl;

  int N, num_trees;

  igzstream is_N(options["anc"].as<std::string>());
  if(is_N.fail()) is_N.open(options["anc"].as<std::string>() + ".gz");
  if(is_N.fail()){
    std::cerr << "Error opening .anc file" << std::endl;
    exit(1);
  }
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.ignore(256, ' ');
  is_N >> num_trees;
  is_N.close();
  const int num_trees_check = num_trees;
    
  int L = 0;
  igzstream is_L(options["mut"].as<std::string>());
  if(is_L.fail()) is_L.open(options["mut"].as<std::string>() + ".gz");
  if(is_L.fail()){
    std::cerr << "Error opening .mut file" << std::endl;
    exit(1);
  }
  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }
  is_L.close();

  Data data(N,L);

  int first_bp = options["first_bp"].as<int>();
  int last_bp  = options["last_bp"].as<int>();

  std::string line;
  igzstream is(options["anc"].as<std::string>());
  if(is.fail()) is.open(options["anc"].as<std::string>() + ".gz");
  if(is.fail()){
    std::cerr << "Error opening .anc file" << std::endl;
    exit(1);
  }
  assert(getline(is,line));
  assert(getline(is,line));
  igzstream is_mut(options["mut"].as<std::string>());
  if(is_mut.fail()) is_mut.open(options["mut"].as<std::string>() + ".gz");
  if(is_mut.fail()){
    std::cerr << "Error opening .mut file" << std::endl;
    exit(1);
  }
  std::string header;
  assert(getline(is_mut, header));
  is_mut.close();

  //output
  std::ofstream os(options["output"].as<std::string>() + ".anc");
  if(os.fail()){
    std::cerr << "Error opening " << options["output"].as<std::string>() + ".anc" << std::endl;
    exit(1);
  }

  Mutations mut;
  mut.Read(options["mut"].as<std::string>());
  Mutations mut_subregion;
  mut_subregion.header = header;

  if(last_bp < mut.info[0].pos || first_bp > mut.info[mut.info.size()-1].pos){
    std::cerr << "Error: Region is outside of anc/mut files." << std::endl;
    exit(1);
  }

  int tree_index_begin = -1;
  int tree_index_end   = -1;
  for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
    if( (*it_mut).pos >= first_bp && tree_index_begin == -1 ){
      tree_index_begin = (*it_mut).tree;
    }
    
    if((*it_mut).pos >= first_bp && (*it_mut).pos <= last_bp){
      (*it_mut).tree -= tree_index_begin;
      mut_subregion.info.push_back(*it_mut);
      (*it_mut).tree += tree_index_begin;
    }

    if(std::next(it_mut,1) == mut.info.end()){
      tree_index_end = (*it_mut).tree;
      break;
    }else if( (*it_mut).pos <= last_bp && (*std::next(it_mut,1)).pos > last_bp ){
      tree_index_end = (*it_mut).tree;
      break;
    }
  }

  mut_subregion.Dump(options["output"].as<std::string>() + ".mut");

  int tree_index = mut.info[0].tree;
  assert(tree_index == 0);

  os << "NUM_HAPLOTYPES " << data.N << "\n";
  os << "NUM_TREES " << tree_index_end - tree_index_begin + 1 << "\n";
  while(getline(is, line)){
    if(tree_index >= tree_index_begin && tree_index <= tree_index_end){
      os << line << "\n";
    } 
    tree_index++;
  }
  os.close();

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

}

