#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>

#include "collapsed_matrix.hpp"
#include "sample.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

#include "AvgMutationRate.cpp"

char
IsCharNucl(char c){

  std::string nucl = "ACGT";
  std::string nucl_small = "acgt";

  for(int i = 0; i < 4; i++){
    if(nucl[i] == c){
      return c;
    }else if(nucl_small[i] == c){
      return nucl[i];
    }
  }

  return 'N';

}

void
CountBasesByType(Data& data, const std::string& filename_mask, const std::string& filename_ancestor, CollapsedMatrix<double>& count_bases_by_type, std::map<std::string, int>& dict_mutation_pattern, Mutations& mutations, std::vector<int>& pos){

  count_bases_by_type.resize(mutations.info.size(),dict_mutation_pattern.size());

  std::ifstream is;
  std::string line;

  int i = 0;

  //read ancestral type
  std::string ancestor;
  is.open(filename_ancestor);
  if(is.fail()){
    std::cerr << "Error while opening ancestor." << std::endl;
    exit(1);
  }
  getline(is,line);
  while(getline(is,line)){
    ancestor += line;
  }
  is.close();
  std::cerr << "Done with reading ancestor." << std::endl;

  //read mask
  std::string mask;
  is.open(filename_mask);
  if(is.fail()){
    std::cerr << "Error while opening mask." << std::endl;
    exit(1);
  }
  getline(is,line);
  while(getline(is,line)){
    mask += line;
  }
  is.close();
  std::cerr << "Done with reading mask." << std::endl;

  //count bases by type

  std::string pattern;
  std::string nucl = "ACGT";
  int mask_threshold = 1801;
  std::string::iterator it_mask, it_start, it_end, it_ancestor; 

  std::cerr << mask.size() << " " << ancestor.size() << std::endl;
  assert(mask.size() == ancestor.size());

  it_start = mask.begin();
  it_end   = std::next(mask.begin(), std::min((int)mask.size(),1001));
  int d_num_nonpass_vincity = 0;
  for(it_mask = it_start; it_mask != it_end; it_mask++){
    if(*it_mask != 'P'){
      d_num_nonpass_vincity++;
    } 
  }
  it_end--;

  it_mask = mask.begin();
  it_ancestor = ancestor.begin();
  std::vector<SNPInfo>::iterator it_info = mutations.info.begin();
  std::vector<int>::iterator it_pos = pos.begin();
  int p = 0;
  int snp = 0;
  while(it_end != mask.end() && it_mask != std::next(mask.begin(), 1001) && p < (*it_info).pos){
    it_end++;
    if(*it_end != 'P') d_num_nonpass_vincity++;

    it_mask++;
    it_ancestor++;
    p++;
  }

  if(it_mask == std::next(mask.begin(), 1001)){

    while(it_end != mask.end() && it_mask != std::next(mask.begin(), 1001)){
      it_end++;
      if(*it_end != 'P') d_num_nonpass_vincity++;

      //only add if its between the previous and the next snp (in the mutations file with all SNPs)
      if( p >= 0.5*((*it_pos) + (*std::prev(it_pos,1))) && p < 0.5 * ((*it_pos) + (*std::next(it_pos,1))) ){

        if(*it_mask == 'P' && d_num_nonpass_vincity < 0.5 * mask_threshold && (*it_info).branch.size() == 1){
          //add to count 
          if(IsCharNucl(*std::prev(it_ancestor,1)) != 'N' && IsCharNucl(*std::next(it_ancestor,1)) != 'N' && IsCharNucl(*it_ancestor) != 'N'){
            pattern  = toupper(*std::prev(it_ancestor,1));
            pattern += toupper(*std::next(it_ancestor,1));
            pattern += toupper(*it_ancestor);  

            for(std::string::iterator it_nucl = nucl.begin(); it_nucl != nucl.end(); it_nucl++){
              if(*it_nucl != IsCharNucl(*it_ancestor)){
                count_bases_by_type[snp][dict_mutation_pattern[pattern + "/" + *it_nucl]] += 1.0;
              }
            }
          }
        }

      }

      if(p >= 0.5 * ((*std::next(it_pos,1)) + (*it_pos))){
        it_info++;
        snp++;
        if(it_info == mutations.info.end()) break;
      }
      while((*it_pos) < (*it_info).pos){
        it_pos++;
      }

      it_mask++;
      it_ancestor++;
      p++;
    }

  }else{
    while(it_end != mask.end() && p < (*it_info).pos){
      if(*it_start != 'P') d_num_nonpass_vincity--;
      it_start++;
      it_end++;
      if(*it_end != 'P') d_num_nonpass_vincity++;
      p++;
      it_mask++;
      it_ancestor++;
    } 
  }

  assert(p <= (*it_info).pos);

  while(it_end != std::prev(mask.end(),1) && it_info != std::prev(mutations.info.end(),1)){
    if(*it_start != 'P') d_num_nonpass_vincity--;
    it_start++;
    it_end++;
    if(*it_end != 'P') d_num_nonpass_vincity++;

    assert(d_num_nonpass_vincity >= 0);

    //only add if its between the previous and the next snp (in the mutations file with all SNPs)
    if( p >= 0.5*((*it_pos) + (*std::prev(it_pos,1))) && p < 0.5 * ((*it_pos) + (*std::next(it_pos,1))) ){

      if(*it_mask == 'P' && d_num_nonpass_vincity < mask_threshold && (*it_info).branch.size() == 1){
        //add to count
        if(IsCharNucl(*std::prev(it_ancestor,1)) != 'N' && IsCharNucl(*std::next(it_ancestor,1)) != 'N' && IsCharNucl(*it_ancestor) != 'N'){
          pattern  = toupper(*std::prev(it_ancestor,1));
          pattern += toupper(*std::next(it_ancestor,1));
          pattern += toupper(*it_ancestor);  

          int count = 0;
          for(std::string::iterator it_nucl = nucl.begin(); it_nucl != nucl.end(); it_nucl++){
            if(*it_nucl != IsCharNucl(*it_ancestor)){
              count++;
              count_bases_by_type[snp][dict_mutation_pattern[pattern + "/" + *it_nucl]] += 1.0;
            }
          }
          assert(count == 3);
        }
      }  

    }

    if(p >= 0.5 * ((*std::next(it_pos,1)) + (*it_pos))){
      it_info++;
      snp++;
      if(it_info == mutations.info.end()) break;
    }
    while((*it_pos) < (*it_info).pos){
      it_pos++;
    }
    if(it_info == std::prev(mutations.info.end(),1)) break;

    it_mask++;
    it_ancestor++;
    p++;

  }

  while(it_mask != std::prev(mask.end(),1) && it_info != std::prev(mutations.info.end(),1)){
    if(*it_start != 'P') d_num_nonpass_vincity--;
    it_start++;

    //if(snp != data.L-1){
    //  int begin = std::max(0,snp-10);
    //  int end   = std::min(data.L-1, snp+10);
    //  rec = (data.rpos[end] - data.rpos[begin])/(mutations.info[end].pos - mutations.info[begin].pos);
    //}

    assert(d_num_nonpass_vincity >= 0);
    //only add if its between the previous and the next snp (in the mutations file with all SNPs)
    if( p >= 0.5*((*it_pos) + (*std::prev(it_pos,1))) && p < 0.5 * ((*it_pos) + (*std::next(it_pos,1))) ){

      if(*it_mask == 'P' && d_num_nonpass_vincity < 0.5 * mask_threshold && (*it_info).branch.size() == 1){
        //add to count 
        if(IsCharNucl(*std::prev(it_ancestor,1)) != 'N' && IsCharNucl(*std::next(it_ancestor,1)) != 'N' && IsCharNucl(*it_ancestor) != 'N'){
          pattern  = toupper(*std::prev(it_ancestor,1));
          pattern += toupper(*std::next(it_ancestor,1));
          pattern += toupper(*it_ancestor);  

          for(std::string::iterator it_nucl = nucl.begin(); it_nucl != nucl.end(); it_nucl++){
            if(*it_nucl != IsCharNucl(*it_ancestor)){
              count_bases_by_type[snp][dict_mutation_pattern[pattern + "/" + *it_nucl]] += 1.0;
            }
          }
        }
      }  

    }

    if(p >= 0.5 * ((*std::next(it_pos,1)) + (*it_pos))){
      it_info++;
      snp++;
      if(it_info == mutations.info.end()) break;
    }
    while((*it_pos) < (*it_info).pos){
      it_pos++;
    }

    it_mask++;
    it_ancestor++;
    p++;
  }

}

/////////////// Estimate mutation rate ////////////

void FinalizeAvg(cxxopts::Options& options){


  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Estimate population size using coalescent rate." << std::endl;
    exit(0);
  } 

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing mutation rate..." << std::endl;


  //epoch
  /*
  int num_epochs = 40;
  std::vector<float> epoch(num_epochs);
  epoch[0] = 0.0;
  for(int e = 1; e < num_epochs; e++){
    epoch[e] = std::exp(5.0*(e+9)/15.0);
  }
  */

  int num_epochs;
  std::vector<float> epoch;

  CollapsedMatrix<double> mutation_by_type_and_epoch, opportunity_by_type_and_epoch;
  FILE* fp; 
  fp = fopen((options["input"].as<std::string>() + "_mut" + ".bin" ).c_str(), "r");  
  assert(fp != NULL);

  fread(&num_epochs, sizeof(int), 1, fp);
  epoch.resize(num_epochs);
  fread(&epoch[0], sizeof(float), num_epochs, fp);

  mutation_by_type_and_epoch.ReadFromFile(fp);
  fclose(fp);
  fp = fopen((options["input"].as<std::string>() + "_opp" + ".bin" ).c_str(), "r");  
  opportunity_by_type_and_epoch.ReadFromFile(fp); 


  std::ofstream os(options["output"].as<std::string>() + ".rate");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  //output
  for(int ep = 0; ep < num_epochs-1; ep++){
    std::vector<double>::iterator it_row_mut = mutation_by_type_and_epoch.rowbegin(ep);
    std::vector<double>::iterator it_row_opp = opportunity_by_type_and_epoch.rowbegin(ep);
    os << epoch[ep] << " ";
    float mut = 0.0, opp = 0.0;
    for(; it_row_mut != mutation_by_type_and_epoch.rowend(ep);){
      mut += *it_row_mut;
      opp += *it_row_opp;
      //std::cerr << epoch[ep] << " " << *it_row_opp << std::endl;
      it_row_mut++;
      it_row_opp++;
    }
    os << mut/opp*3.0 << "\n";
    //std::cerr << mut << " " << opp << " " << mut/opp*3.0 << std::endl; 
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

void FinalizeMutationRate(cxxopts::Options& options){


  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Estimate population size using coalescent rate." << std::endl;
    exit(0);
  } 

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing mutation rate..." << std::endl;


  //epoch
  /*
  int num_epochs = 40;
  std::vector<float> epoch(num_epochs);
  epoch[0] = 0.0;
  for(int e = 1; e < num_epochs; e++){
    epoch[e] = std::exp(5.0*(e+9)/15.0);
  }
  */

  int num_epochs;
  std::vector<float> epoch;

  CollapsedMatrix<double> mutation_by_type_and_epoch, opportunity_by_type_and_epoch;
  FILE* fp; 
  fp = fopen((options["input"].as<std::string>() + "_mut" + ".bin" ).c_str(), "r");  
  assert(fp != NULL);

  fread(&num_epochs, sizeof(int), 1, fp);
  epoch.resize(num_epochs);
  fread(&epoch[0], sizeof(float), num_epochs, fp);

  mutation_by_type_and_epoch.ReadFromFile(fp);
  fclose(fp);
  fp = fopen((options["input"].as<std::string>() + "_opp" + ".bin" ).c_str(), "r");  
  opportunity_by_type_and_epoch.ReadFromFile(fp); 


  std::ofstream os(options["output"].as<std::string>() + ".rate");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  std::string alphabet = "ACGT", pattern;
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      pattern  = *it_str1;
      pattern += *it_str2;
      os << pattern + "C/A " << pattern + "C/G " << pattern + "C/T " << pattern + "T/A " << pattern + "T/C " << pattern + "T/G "; 
    }   
  }
  os << "\n";

  //output
  for(int ep = 0; ep < num_epochs-1; ep++){
    std::vector<double>::iterator it_row_mut = mutation_by_type_and_epoch.rowbegin(ep);
    std::vector<double>::iterator it_row_opp = opportunity_by_type_and_epoch.rowbegin(ep);
    os << epoch[ep] << " ";
    for(; it_row_mut != mutation_by_type_and_epoch.rowend(ep);){
      os << *it_row_mut/(*it_row_opp) << " ";
      it_row_mut++;
      it_row_opp++;
    }
    os << "\n"; 
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

void SummarizeWholeGenome(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("first_chr") || !options.count("last_chr") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: first_chr, last_chr, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .bin files and calculates summarizes them." << std::endl;
    exit(0);
  }  

  int start     = options["first_chr"].as<int>(); 
  int end       = options["last_chr"].as<int>();
  std::vector<std::string> filenames;
  std::string filename_base = options["output"].as<std::string>();

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Summarizing mutation rates for chr " << start << " " << end << "..." << std::endl;

  for(int chr = start; chr <= end; chr++){
    filenames.push_back(filename_base + "_chr" + std::to_string(chr) + "_mut" + ".bin");
  }

  //open files and add together

  FILE* fp = fopen(filenames[0].c_str(),"rb");
  assert(fp != NULL);  

  int num_epochs;
  std::vector<float> epochs;
  fread(&num_epochs, sizeof(int), 1, fp);
  epochs.resize(num_epochs);
  fread(&epochs[0], sizeof(float), num_epochs, fp);

  CollapsedMatrix<double> mut_by_type_and_epoch, mut_by_type_and_epoch_tmp;
  mut_by_type_and_epoch.ReadFromFile(fp);
  fclose(fp);

  for(int i = 1; i < (int) filenames.size(); i++){
    fp = fopen(filenames[i].c_str(),"rb");
    
    fread(&num_epochs, sizeof(int), 1, fp);
    fread(&epochs[0], sizeof(float), num_epochs, fp);
    
    mut_by_type_and_epoch_tmp.ReadFromFile(fp);
    fclose(fp);
    std::vector<double>::iterator it_mut = mut_by_type_and_epoch.vbegin();
    for(std::vector<double>::iterator it_mut_tmp = mut_by_type_and_epoch_tmp.vbegin(); it_mut_tmp != mut_by_type_and_epoch_tmp.vend();){
      *it_mut += *it_mut_tmp;
      it_mut++;
      it_mut_tmp++;
    }

  }

  filenames.clear();
  for(int chr = start; chr <= end; chr++){
    filenames.push_back(filename_base + "_chr" + std::to_string(chr) + "_opp" + ".bin");
  }

  CollapsedMatrix<double> opp_by_type_and_epoch, opp_by_type_and_epoch_tmp;
  fp = fopen(filenames[0].c_str(),"rb");
  opp_by_type_and_epoch.ReadFromFile(fp);
  fclose(fp);

  for(int i = 1; i < (int) filenames.size(); i++){

    fp = fopen(filenames[i].c_str(),"rb");
    opp_by_type_and_epoch_tmp.ReadFromFile(fp);
    fclose(fp);
    std::vector<double>::iterator it_opp = opp_by_type_and_epoch.vbegin();
    for(std::vector<double>::iterator it_opp_tmp = opp_by_type_and_epoch_tmp.vbegin(); it_opp_tmp != opp_by_type_and_epoch_tmp.vend();){
      *it_opp += *it_opp_tmp;
      it_opp++;
      it_opp_tmp++;
    }

  }

  for(int chr = start; chr <= end; chr++){ 
    std::remove((options["input"].as<std::string>() + "_chr" + std::to_string(chr) + "_mut" + ".bin").c_str());
    std::remove((options["input"].as<std::string>() + "_chr" + std::to_string(chr) + "_opp" + ".bin").c_str());
  }  

  fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb"); 
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(float), epochs.size(), fp);
  mut_by_type_and_epoch.DumpToFile(fp);
  fclose(fp);
  fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  opp_by_type_and_epoch.DumpToFile(fp);
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

}

void MutationRateWithContext(cxxopts::Options& options, int chr = -1){

  bool help = false;
  if(!options.count("mask") || !options.count("ancestor") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mask, ancestor, input, output. Optional: num_bins, dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate mutation rate for 96 categories (Do not apply after using RemoveTreesWithFewMutations).." << std::endl;
    exit(0);
  }  

  std::string line, line2 , read;

  //////////// PARSE DATA ///////////

  int N;
  std::ifstream is_N(options["input"].as<std::string>() + ".anc");
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  int L = 0;
  std::ifstream is_L(options["input"].as<std::string>() + ".mut");
  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }
  is_L.close();

  Data data(N,L);
  int N_total = 2*data.N-1;

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating mutation rate for 96 categories.." << std::endl;

  ////////// read mutations file ///////////

  Mutations mutations;
  if(chr == -1){
    mutations.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mutations.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
  }

  std::vector<int> pos;
  if(options.count("dist")){

    int L_allsnps = 0;
    is_L.open(options["dist"].as<std::string>());
    std::string unused;
    std::getline(is_L, unused); 
    while ( std::getline(is_L, unused) ){
      ++L_allsnps;
    }
    is_L.close();

    pos.resize(L_allsnps);
    std::ifstream is_dist(options["dist"].as<std::string>());
    if(is_dist.fail()){
      std::cerr << "Error while opening file." << std::endl;
      exit(1);
    }
    getline(is_dist, line); 
    int snp = 0, dist;
    while(std::getline(is_dist, line)){
      sscanf(line.c_str(), "%d %d", &pos[snp], &dist);
      snp++;
    }
    is_dist.close();

  }else{
  
    pos.resize(data.L);
    int snp = 0;
    for(std::vector<SNPInfo>::iterator it_mut = mutations.info.begin(); it_mut != mutations.info.end(); it_mut++){
      pos[snp]  = (*it_mut).pos;
      snp++;
    }
  
  }

  ////////// EPOCHES ///////////

  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epoch(num_epochs);
  epoch[0] = 0.0;
  epoch[1] = 1e3/28.0;
  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){
    epoch[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/28.0; 
  }
  epoch[num_epochs-1] = 5e7; 

  ////////// define mutation types /////////
  std::string alphabet = "ACGT";
  std::string reverse_alphabet = "TGCA";
  std::map<std::string, int> dict_mutation_pattern;

  //A-C, A-G, A-T, C-A, C-G, C-T
  //T-G, T-C, T-A, G-T, G-C, G-A
  //plus 16 flanks per mutation type, i.e. 6*16 = 96

  std::string pattern, reverse_pattern;
  int index = 0;
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      pattern = *it_str1;
      pattern += *it_str2;
      dict_mutation_pattern[pattern + "C/A"] = index;
      index++;
      dict_mutation_pattern[pattern + "C/G"] = index;
      index++;
      dict_mutation_pattern[pattern + "C/T"] = index;
      index++;
      dict_mutation_pattern[pattern + "A/T"] = index;
      index++;
      dict_mutation_pattern[pattern + "A/G"] = index;
      index++;
      dict_mutation_pattern[pattern + "A/C"] = index;
      index++;
    }   
  }
  index = 0;
  for(std::string::iterator it_str1 = reverse_alphabet.begin(); it_str1 != reverse_alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = reverse_alphabet.begin(); it_str2 != reverse_alphabet.end(); it_str2++){
      reverse_pattern = *it_str2;
      reverse_pattern += *it_str1;
      dict_mutation_pattern[reverse_pattern + "G/T"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "G/C"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "G/A"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "T/A"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "T/C"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "T/G"] = index;
      index++;
    }   
  }
  int num_mutation_cathegories = index;

  ///////// Count number of bases by type //////////

  CollapsedMatrix<double> count_bases_by_type;
  //count_bases_by_type.resize(mutations.info.size(),dict_mutation_pattern.size());
  CountBasesByType(data, options["mask"].as<std::string>(), options["ancestor"].as<std::string>(), count_bases_by_type, dict_mutation_pattern, mutations, pos);

  ////////////////////
  // Estimate mutation rate through time

  CollapsedMatrix<double> mutation_by_type_and_epoch;
  CollapsedMatrix<double> opportunity_by_type_and_epoch;
  mutation_by_type_and_epoch.resize(num_epochs, num_mutation_cathegories);
  opportunity_by_type_and_epoch.resize(num_epochs, num_mutation_cathegories);
  std::vector<double> branch_lengths_in_epoch(num_epochs);

  MarginalTree mtr;
  //Tree subtr;
  std::vector<float> coordinates_tree(N_total);
  int root = N_total-1;
  int i = 0;

  std::ifstream is_anc;
  if(chr == -1){
    is_anc.open(options["input"].as<std::string>() + ".anc");
  }else{
    is_anc.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc");
  }
  if(is_anc.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  getline(is_anc,line);
  getline(is_anc,line);
  getline(is_anc,line);

  //read tree
  mtr.Read(line, N);
  mtr.tree.GetCoordinates(coordinates_tree);
  std::sort(coordinates_tree.begin(), coordinates_tree.end());
  GetBranchLengthsInEpoche(data, epoch, coordinates_tree, branch_lengths_in_epoch);

  SNPInfo snp_info;
  float rec;
  int num_tree = 0;
  double total_branch_length = 0.0;
  int count_snps = 0;
  for(int snp = 0; snp < data.L; snp++){

    snp_info = mutations.info[snp];
    if(snp_info.branch.size() == 1){

        if(num_tree < snp_info.tree){
          while(num_tree < snp_info.tree){
            if(!getline(is_anc,line)){
              break; 
            };
            num_tree++;
          }

          //read tree
          mtr.Read(line, N);
          mtr.tree.GetCoordinates(coordinates_tree);
          std::sort(coordinates_tree.begin(), coordinates_tree.end());
          GetBranchLengthsInEpoche(data, epoch, coordinates_tree, branch_lengths_in_epoch);
        }
        assert(num_tree == snp_info.tree);

        //if(snp_info.age_begin <= coordinates_tree[root]){

          // identify cathegory of mutation
          pattern = snp_info.upstream_base + snp_info.downstream_base + snp_info.mutation_type;
          // then identify its age and number of lineages at the time
          int ind = dict_mutation_pattern[pattern];

          // identify epoch and add to number of mutations per lineage in epoch
          int ep = 0;
          while(epoch[ep] <= snp_info.age_begin){
            ep++;
            if(ep == epoch.size()) break;
          }
          ep--;

          //std::cerr << ep << " " << ind << " " << mutation_by_type_and_epoch.size() << " " << mutation_by_type_and_epoch.subVectorSize(ep) << std::endl;
          //std::cerr << ep << " " << ind << " " << age << " " << epoch[0] << " " << epoch[1] << std::endl;
          assert(ep >= 0);
          assert(mutation_by_type_and_epoch.subVectorSize(ep) == num_mutation_cathegories);
          //assert(count_bases_by_type[snp][ind] > 0);

          int ep_begin = ep;
          float age_end = std::min(snp_info.age_end,coordinates_tree[root]);
          assert(age_end < epoch[num_epochs-1]);
          double branch_length = age_end - snp_info.age_begin;

          if(age_end <= epoch[ep+1]){

            mutation_by_type_and_epoch[ep][ind] += 1.0;

          }else{

            mutation_by_type_and_epoch[ep][ind]   += (epoch[ep+1] - snp_info.age_begin)/branch_length;
            ep++;
            while(epoch[ep+1] <= age_end){
              mutation_by_type_and_epoch[ep][ind] += (epoch[ep+1]-epoch[ep])/branch_length;
              ep++;
            }
            mutation_by_type_and_epoch[ep][ind]   += (age_end-epoch[ep])/branch_length;

          }

          for(int ep_tmp = 0; ep_tmp < num_epochs; ep_tmp++){
            double bl = branch_lengths_in_epoch[ep_tmp];
            for(int ind_tmp = 0; ind_tmp < num_mutation_cathegories; ind_tmp++){
              opportunity_by_type_and_epoch[ep_tmp][ind_tmp] += bl * count_bases_by_type[snp][ind_tmp];
            }
          }
 
          //count_snps++;

    }

  }

  is_anc.close();

  //output mutation_by_type_and_epoch
  //       opportunity_by_type_and_epoch 
  FILE* fp;
  if(chr == -1){
    fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + std::to_string(chr) + "_mut" + ".bin" ).c_str(), "wb");  
  }
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epoch[0], sizeof(float), epoch.size(), fp);
  mutation_by_type_and_epoch.DumpToFile(fp);
  fclose(fp);
  if(chr == -1){
    fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + std::to_string(chr) + "_opp" + ".bin" ).c_str(), "wb");  
  }
  opportunity_by_type_and_epoch.DumpToFile(fp);
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


}


//////////////// Other applications ///////////////

void BranchLengthVsMutations(cxxopts::Options& options){

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating number of mutations vs opportunity..." << std::endl;

  bool help = false;
  if(!options.count("pos") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: pos, input, output. Optional: mut." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "..." << std::endl;
    exit(0);
  }  

  std::string line, read;

  //parse data
  int N;
  std::ifstream is_N(options["input"].as<std::string>() + ".anc");
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  //make this more efficient
  int L = 0;
  std::ifstream is_L(options["pos"].as<std::string>());
  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }

  Data data(N,L);
  //data.ReadPosition(options["pos"].as<std::string>());
  int N_total = 2*data.N-1;

  //epoch
  int num_epochs = 40;
  std::vector<float> epoch(num_epochs);
  epoch[0] = 0.0;
  for(int e = 1; e < num_epochs; e++){
    epoch[e] = std::exp(5.0*(e+9)/15.0);
  }

  //read mutations file
  Mutations mutations(data);
  mutations.Read(options["input"].as<std::string>() + ".mut");

  ///////////////////
  // Count number of bases by type
  double total_num_bases = (*std::prev(mutations.info.end(),1)).pos - (*mutations.info.begin()).pos;
  std::vector<double> count_bases(mutations.info.size(), 0.0);
  std::vector<double>::iterator it_count_bases   = count_bases.begin();

  std::vector<SNPInfo>::iterator it_mut = mutations.info.begin();
  *it_count_bases = (0.5*(*it_mut).dist)/total_num_bases;
  it_count_bases++;
  while(it_count_bases != count_bases.end()){
    *it_count_bases = (0.5*(*it_mut).dist)/total_num_bases;
    if(it_mut != mutations.info.end()){
      it_mut++;
      *it_count_bases += (0.5*(*it_mut).dist)/total_num_bases;
    }
    it_count_bases++;
  }
  it_mut++;
  assert(it_mut == mutations.info.end());



  ////////////////////
  // Branch lengths vs opportunity 

  std::vector<double> num_mutations_in_epoch(num_epochs);
  std::vector<double> branch_lengths_in_epoch(num_epochs);

  MarginalTree mtr;
  std::vector<float> coordinates(N_total);
  int root = N_total-1;
  int i = 0;
  int snp_of_next_tree, snp;


  std::ifstream is_anc(options["input"].as<std::string>() + ".anc");
  if(is_anc.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  getline(is_anc,line);
  getline(is_anc,line);
  getline(is_anc,line);

  //read tree
  mtr.Read(line, N);
  mtr.tree.GetCoordinates(coordinates);

  std::ofstream os(options["output"].as<std::string>() + ".xy");

  std::fill(num_mutations_in_epoch.begin(), num_mutations_in_epoch.end(), 0.0);
  std::fill(branch_lengths_in_epoch.begin(), branch_lengths_in_epoch.end(), 0.0);
  for(int i = 0; i < N_total - 1; i++){
    float num_events = mtr.tree.nodes[i].num_events;
    float bl         = mtr.tree.nodes[i].branch_length;
    int parent       = (*mtr.tree.nodes[i].parent).label;
    int snp_begin    = mtr.tree.nodes[i].SNP_begin;
    int snp_end      = mtr.tree.nodes[i].SNP_end;
    int delta_pos    = mutations.info[snp_end].pos - mutations.info[snp_begin].pos;      
    if(delta_pos < 0.0) std::cerr << 0 << " " << snp_begin << " " << snp_end << " " << mutations.info[snp_end].pos << " " << mutations.info[snp_begin].pos << std::endl;
    assert(delta_pos >= 0.0);

    int ep = 0;
    while(epoch[ep] < coordinates[i]) ep++;
    if(epoch[ep] <= coordinates[parent]){
      assert(epoch[ep] >= coordinates[i]);
      num_mutations_in_epoch[ep-1]    += num_events * (epoch[ep] - coordinates[i])/bl;
      branch_lengths_in_epoch[ep-1]   += delta_pos * (epoch[ep] - coordinates[i]);
      ep++;
      while(epoch[ep] < coordinates[parent]){
        num_mutations_in_epoch[ep-1]  += num_events * (epoch[ep] - epoch[ep-1])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (epoch[ep] - epoch[ep-1]);
        ep++;
      }
      assert(coordinates[parent] >= epoch[ep-1]);
      num_mutations_in_epoch[ep-1]  += num_events * (coordinates[parent] - epoch[ep-1])/bl;
      branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - epoch[ep-1]);
    }else{
      num_mutations_in_epoch[ep-1]  += num_events * (coordinates[parent] - coordinates[i])/bl;
      branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - coordinates[i]);
    }
  }

  for(int ep = 0; ep < epoch.size() - 1; ep++){
    os << mtr.pos << " " << (int) 28 * (epoch[ep] + epoch[ep+1])/2.0 << " " << data.mu * branch_lengths_in_epoch[ep] << " " << num_mutations_in_epoch[ep] << "\n"; 
  }

  while(getline(is_anc, line)){

    //read tree
    mtr.Read(line, N);
    mtr.tree.GetCoordinates(coordinates);

    std::fill(num_mutations_in_epoch.begin(), num_mutations_in_epoch.end(), 0.0);
    std::fill(branch_lengths_in_epoch.begin(), branch_lengths_in_epoch.end(), 0.0);
    for(int i = 0; i < N_total - 1; i++){
      float num_events = mtr.tree.nodes[i].num_events;
      float bl         = mtr.tree.nodes[i].branch_length;
      int parent       = (*mtr.tree.nodes[i].parent).label;
      int snp_begin    = mtr.tree.nodes[i].SNP_begin;
      int snp_end      = mtr.tree.nodes[i].SNP_end;
      if(snp_end >= data.L) snp_end = data.L-1;
      int delta_pos    = mutations.info[snp_end].pos - mutations.info[snp_begin].pos;

      if(delta_pos < 0.0) std::cerr << snp << " " << snp_begin << " " << snp_end << " " << mutations.info[snp_end].pos << " " << mutations.info[snp_begin].pos << std::endl;
      assert(delta_pos >= 0.0);
      if(num_events < 0.0)  std::cerr << snp << " " << snp_begin << " " << snp_end << " " << mutations.info[snp_end].pos << " " << mutations.info[snp_begin].pos << std::endl;
      assert(num_events >= 0.0);

      int ep = 0;
      while(epoch[ep] < coordinates[i]) ep++;
      assert(epoch[ep] >= coordinates[i]);

      if(epoch[ep] <= coordinates[parent]){
        num_mutations_in_epoch[ep-1] += num_events * (epoch[ep] - coordinates[i])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (epoch[ep] - coordinates[i]);
        ep++;
        while(epoch[ep] < coordinates[parent]){
          num_mutations_in_epoch[ep-1] += num_events * (epoch[ep] - epoch[ep-1])/bl;
          branch_lengths_in_epoch[ep-1] += delta_pos * (epoch[ep] - epoch[ep-1]);
          ep++;
        }
        assert(coordinates[parent] >= epoch[ep-1]);
        num_mutations_in_epoch[ep-1] += num_events * (coordinates[parent] - epoch[ep-1])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - epoch[ep-1]);
      }else{
        num_mutations_in_epoch[ep-1] += num_events * (coordinates[parent] - coordinates[i])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - coordinates[i]);
      }
    }

    for(int ep = 0; ep < epoch.size()-1; ep++){
      os << mtr.pos << " " << (int) 28 * (epoch[ep] + epoch[ep+1])/2.0 << " " << data.mu * branch_lengths_in_epoch[ep] << " " << num_mutations_in_epoch[ep] << "\n"; 
    }

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

void FinalizeMutationCount(cxxopts::Options& options){


  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Estimate population size using coalescent rate." << std::endl;
    exit(0);
  } 

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing mutation count.." << std::endl;


  int num_epochs;
  std::vector<float> epoch;

  CollapsedMatrix<double> mutation_by_type_and_epoch, opportunity_by_type_and_epoch;
  FILE* fp; 
  fp = fopen((options["input"].as<std::string>() + "_mut" + ".bin" ).c_str(), "r");  
  assert(fp != NULL);

  fread(&num_epochs, sizeof(int), 1, fp);
  epoch.resize(num_epochs);
  fread(&epoch[0], sizeof(float), num_epochs, fp);
  mutation_by_type_and_epoch.ReadFromFile(fp);
  fclose(fp);

  std::ofstream os(options["output"].as<std::string>() + ".mcount");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  std::string alphabet = "ACGT", pattern;
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      pattern  = *it_str1;
      pattern += *it_str2;
      os << pattern + "C/A " << pattern + "C/G " << pattern + "C/T " << pattern + "T/A " << pattern + "T/C " << pattern + "T/G "; 
    }   
  }
  os << "\n";

  //output
  for(int ep = 0; ep < num_epochs-1; ep++){
    std::vector<double>::iterator it_row_mut = mutation_by_type_and_epoch.rowbegin(ep);
    os << epoch[ep] << " ";
    for(; it_row_mut != mutation_by_type_and_epoch.rowend(ep);){
      os << *it_row_mut << " ";
      it_row_mut++;
    }
    os << "\n"; 
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


int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options
  cxxopts::Options options("RelateMutationRate");
  options.add_options()
    ("help", "Print help")
    ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
    ("first_chr", "Index of fist chr", cxxopts::value<int>())
    ("last_chr", "Index of last chr", cxxopts::value<int>())
    ("num_bins", "Number of bins.", cxxopts::value<int>())
     ("dist", "Filename of file containing dist.", cxxopts::value<std::string>())
    ("mask", "Filename of file containing mask", cxxopts::value<std::string>())
    ("ancestor", "Filename of file containing human ancestor genome.", cxxopts::value<std::string>())
    ("i,input", "Filename of .anc and .mut file without file extension", cxxopts::value<std::string>())
    ("o,output", "Output file", cxxopts::value<std::string>());

  options.parse(argc, argv);
  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("WithContext")){

    //////////////////////////////////
    //Program options
    bool help = false;
    if( !options.count("input") || !options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: mask, ancestor, input, output. Optional: first_chr, last_chr, num_bins, dist." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Calculate mutation rate for 96 categories (Do not apply after using RemoveTreesWithFewMutations)." << std::endl;
      exit(0);
    } 
    if(options["input"].as<std::string>() != options["output"].as<std::string>()){
      std::cerr << "Sorry, in this mode input and output need to be the same!." << std::endl;
      exit(1);
    }

    if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      for(int chr = options["first_chr"].as<int>(); chr <= options["last_chr"].as<int>(); chr++){ 
        MutationRateWithContext(options, chr);
      }
      SummarizeWholeGenome(options);      
      FinalizeMutationRate(options);
    }else{
      MutationRateWithContext(options);
      FinalizeMutationRate(options);
    }

  }else if(!mode.compare("WithContextForChromosome")){

    MutationRateWithContext(options);

  }else if(!mode.compare("SummarizeForGenome")){

    SummarizeWholeGenome(options);

  }else if(!mode.compare("Finalize")){

    FinalizeMutationRate(options);

  }else if(!mode.compare("FinalizeMutationCount")){

    FinalizeMutationCount(options);

  }else if(!mode.compare("FinalizeAvg")){

    FinalizeAvg(options);

  }else if(!mode.compare("Avg")){

    AvgMutationRate(options);

  }else if(!mode.compare("XY")){

    BranchLengthVsMutations(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "Avg, WithContext, WithContextForChromosome, SummarizeForGenome, Finalize, FinalizeMutationCount, XY." << std::endl;

  }

  bool help = false;
  if(!options.count("mode")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    exit(0);
  }

}

