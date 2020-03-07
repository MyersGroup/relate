#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>
#include <string>

#include "cxxopts.hpp"
#include "filesystem.hpp"
#include "collapsed_matrix.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"

int Finalize(cxxopts::Options& options){

  bool help = false;
  if(!options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: output. Optional: annot, sample_ages." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use at the end to finalize results. This will summarize all sections into one file." << std::endl;
    exit(0);
  }


  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing..." << std::endl;

  int N, L, num_chunks;
  double tmp;
  int overlap_chunk_size = 10000;
  int delta_chunk;
 
  int i, j; 
  std::string line, line2, read;

  std::vector<int> section_boundary_start, section_boundary_end;
  FILE* fp = fopen("parameters.bin", "r");
  assert(fp != NULL);
  fread(&N, sizeof(int), 1, fp);
  fread(&L, sizeof(int), 1, fp);
  fread(&num_chunks, sizeof(int), 1, fp);
  section_boundary_start.resize(num_chunks);
  section_boundary_end.resize(num_chunks);
  fread(&tmp, sizeof(double), 1, fp);
  fread(&section_boundary_start[0], sizeof(int), num_chunks, fp);
  fread(&section_boundary_end[0], sizeof(int), num_chunks, fp);
  fclose(fp); 

  std::vector<double> sample_ages(N);
  if(options.count("sample_ages")){
    igzstream is_ages(options["sample_ages"].as<std::string>());
    int i = 0; 
    while(is_ages >> sample_ages[i]){
      i++;
      sample_ages[i] = sample_ages[i-1];
      i++;
      if(i == N) break;
    }
    if(i < N) sample_ages.clear();
  }else{
    sample_ages.clear(); 
  }

  ///////////////////////////////////////// Combine AncesTrees /////////////////////////

  //Merge files to one. 

  std::ifstream is;
  std::ofstream os;

  int num_flips = 0, num_non_mapping = 0;
  int num_trees_cum = 0, first_tree = 0;

  ///////////////////////////////////////// Combine Mutation Files /////////////////////////

  FILE* fp_props = fopen("props.bin", "rb");
  std::ifstream is_annot;
  bool exists_annot = false;
  if(options.count("annot")){
    is_annot.open(options["annot"].as<std::string>());
    getline(is_annot, line2);
    exists_annot = true;
  }

  os.open(options["output"].as<std::string>() + ".mut");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  } 

  //Header-Legend
  os << "snp;pos_of_snp;dist;rs-id;tree_index;branch_indices;is_not_mapping;is_flipped;age_begin;age_end;ancestral_allele/alternative_allele;";
  if(exists_annot) os << line2;
  os << "\n";

  int pos, dist, snp_ind;
  char rsid[1024];
  char ancestral[1024], alternative[1024];

  //Rest
  for(int c = 0; c < num_chunks; c++){

    std::string file_prefix = "chunk_" + std::to_string(c) + "/" + options["output"].as<std::string>();

    is.open(file_prefix + "_c" + std::to_string(c) + ".mut");
    std::remove((file_prefix + "_c" + std::to_string(c) + ".mut").c_str());
    if(is.fail()){
      std::cerr << "Error while opening file." << std::endl;
      exit(1);
    }  
    getline(is,line);

    delta_chunk = section_boundary_end[c] - section_boundary_start[c];
    if(c > 0){
      for(int snp = 0; snp < overlap_chunk_size; snp++){
        getline(is,line);
      }
      if(c + 1 != num_chunks) delta_chunk -= overlap_chunk_size;
    }
    if(num_chunks > 1) delta_chunk -= overlap_chunk_size;
    
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
      }else if(std::stoi(read) > num_trees_chunk){ //new tree
        num_trees_chunk++;
      }
      read.clear();

      while(line[i] != ';') i++;
      i++;
      if(line[i] == '1') num_non_mapping++;
      i+=2;
      
      if(line[i] == '1') num_flips++;
      line.erase(0,j);

      fread(&snp_ind, sizeof(int), 1, fp_props);
      fread(&pos, sizeof(int), 1, fp_props);
      fread(&dist, sizeof(int), 1, fp_props);
      fread(&rsid[0], sizeof(char), 1024, fp_props);
      fread(&ancestral[0], sizeof(char), 1024, fp_props);
      fread(&alternative[0], sizeof(char), 1024, fp_props);
      
      os << snp_ind << ";" << pos << ";" << dist << ";" << rsid << ";";
      os << num_trees_chunk + num_trees_cum - first_tree << line;
      os << ancestral << "/" << alternative << ";";
      if(exists_annot){ 
        getline(is_annot, line2);
        os << line2;
      }
      os << "\n";

      snp++;
    }

    num_trees_cum += num_trees_chunk - first_tree + 1;
    is.close();

  }
  if(exists_annot) is_annot.close();
  fclose(fp_props);
  os.close();

  std::cerr << "Number of not mapping SNPs: " << num_non_mapping << "\n";
  std::cerr << "Number of flipped SNPs    : " << num_flips << std::endl;

  //////////////////////////////////////////////////

  std::string filename_os = options["output"].as<std::string>() + ".anc";
  FILE *pfile = std::fopen(filename_os.c_str(), "w");

  if(pfile == NULL){

    std::cerr << "Error while writing to " << filename_os << "." << std::endl;

  }else{

    if(sample_ages.size() == 0){
      fprintf(pfile, "NUM_HAPLOTYPES %d\n", N);
    }else{
      fprintf(pfile, "NUM_HAPLOTYPES %d ", N);
      for(std::vector<double>::iterator it_ages = sample_ages.begin(); it_ages != sample_ages.end(); it_ages++){
        fprintf(pfile, "%f ", *it_ages);
      }     
      fprintf(pfile, "\n");
    }
    fprintf(pfile, "NUM_TREES %d\n", num_trees_cum);

  }


  int num_trees = 0; 
  for(int c = 0; c < num_chunks; c++){

    int start_chunk = section_boundary_start[c];
    int end_chunk   = section_boundary_end[c];
    if(num_chunks > 1 && c+1 != num_chunks){
      end_chunk -= overlap_chunk_size;
    }

    std::string file_prefix = "chunk_" + std::to_string(c) + "/" + options["output"].as<std::string>();

    int position;
    AncesTree anc;
    anc.ReadBin(file_prefix + "_c" + std::to_string(c) + ".anc");  
    std::remove((file_prefix + "_c" + std::to_string(c) + ".anc").c_str());

    //first tree
    if(c == 0){      
      (*anc.seq.begin()).pos = start_chunk;
      for(std::vector<Node>::iterator it_node = (*anc.seq.begin()).tree.nodes.begin(); it_node != (*anc.seq.begin()).tree.nodes.end(); it_node++){
        (*it_node).SNP_begin += start_chunk;
        (*it_node).SNP_end   += start_chunk;
      }
      num_trees++;
    }else{

      CorrTrees::iterator it_anc = anc.seq.begin();
      while(std::next(it_anc,1) != anc.seq.end()){
        if((*std::next(it_anc,1)).pos <= overlap_chunk_size){
          it_anc = anc.seq.erase(it_anc);
        }else{
          break;
        }        
      }

      (*anc.seq.begin()).pos = overlap_chunk_size + start_chunk;
      for(std::vector<Node>::iterator it_node = (*anc.seq.begin()).tree.nodes.begin(); it_node != (*anc.seq.begin()).tree.nodes.end(); it_node++){
        (*it_node).SNP_begin += start_chunk;
        (*it_node).SNP_end   += start_chunk;
      }
      num_trees++;
    
    }
    for(CorrTrees::iterator it_anc = std::next(anc.seq.begin(),1); it_anc != anc.seq.end(); ){
      position = (*it_anc).pos + start_chunk;

      if(position < end_chunk){
        (*it_anc).pos = position;
        for(std::vector<Node>::iterator it_node = (*it_anc).tree.nodes.begin(); it_node != (*it_anc).tree.nodes.end(); it_node++){
          (*it_node).SNP_begin += start_chunk;
          (*it_node).SNP_end   += start_chunk;
        }
        num_trees++;
        it_anc++;
      }else{ 
        it_anc = anc.seq.erase(it_anc);
      }
    }

    //Dump to file.
    anc.Dump(pfile);

  }


  fclose(pfile);

  assert(num_trees == num_trees_cum);
  std::remove("parameters.bin");
  std::remove("props.bin");

  filesys f;
  for(int c = 0; c < num_chunks; c++){
    //now delete directories
    f.RmDir( ("chunk_" + std::to_string(c) + "/paint/").c_str() );
    f.RmDir( ("chunk_" + std::to_string(c) + "/").c_str() );
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
