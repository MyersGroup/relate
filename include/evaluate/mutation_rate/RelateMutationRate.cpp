#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.hpp"
#include "collapsed_matrix.hpp"
#include "sample.hpp"
#include "anc.hpp"
#include "mutations.hpp"
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

  igzstream is;
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
  fasta mask;
  mask.Read(filename_mask);
  std::cerr << "Done with reading mask." << std::endl;

  //count bases by type

  std::string pattern;
  std::string nucl = "ACGT";
  //int mask_threshold = 1801;
  int mask_threshold = 2000;
  std::string::iterator it_mask, it_start, it_end, it_ancestor; 

  std::cerr << mask.seq.size() << " " << ancestor.size() << std::endl;
  if(mask.seq.size() < ancestor.size()){
    int mask_size = mask.seq.size();
    mask.seq.resize(ancestor.size());
    std::fill(std::next(mask.seq.begin(),mask_size), mask.seq.end(), 'N');
  }else{
    int anc_size = ancestor.size();
    ancestor.resize(mask.seq.size());
    std::fill(std::next(ancestor.begin(),anc_size), ancestor.end(), 'N');
  }
  assert(mask.seq.size() == ancestor.size());


  it_start = mask.seq.begin();
  it_end   = std::next(mask.seq.begin(), std::min((int)mask.seq.size(),1001));
  int d_num_nonpass_vincity = 0;
  for(it_mask = it_start; it_mask != it_end; it_mask++){
    if(*it_mask != 'P'){
      d_num_nonpass_vincity++;
    } 
  }
  it_end--;

  it_mask = mask.seq.begin();
  it_ancestor = ancestor.begin();
  std::vector<SNPInfo>::iterator it_info = mutations.info.begin();
  std::vector<int>::iterator it_pos = pos.begin();
  int p = 0;
  int snp = 0;
  while(it_end != mask.seq.end() && it_mask != std::next(mask.seq.begin(), 1001) && p < (*it_info).pos){
    it_end++;
    if(*it_end != 'P') d_num_nonpass_vincity++;

    it_mask++;
    it_ancestor++;
    p++;
  }

  if(it_mask == std::next(mask.seq.begin(), 1001)){

    while(it_end != mask.seq.end() && it_mask != std::next(mask.seq.begin(), 1001)){
      it_end++;
      if(*it_end != 'P') d_num_nonpass_vincity++;

      //only add if its between the previous and the next snp (in the mutations file with all SNPs)
      if( p >= 0.5*((*it_pos) + (*std::prev(it_pos,1))) && p < 0.5 * ((*it_pos) + (*std::next(it_pos,1))) ){

        if(*it_mask == 'P' && d_num_nonpass_vincity <= 0.5 * mask_threshold && (*it_info).branch.size() == 1){
          //add to count 
          if(IsCharNucl(*std::prev(it_ancestor,1)) != 'N' && IsCharNucl(*std::next(it_ancestor,1)) != 'N' && IsCharNucl(*it_ancestor) != 'N'){
            pattern  = toupper(*std::prev(it_ancestor,1));
            pattern += toupper(*std::next(it_ancestor,1));
            pattern += toupper(*it_ancestor);  

            for(std::string::iterator it_nucl = nucl.begin(); it_nucl != nucl.end(); it_nucl++){
              if(*it_nucl != IsCharNucl(*it_ancestor)){
                count_bases_by_type[snp][dict_mutation_pattern[pattern + *it_nucl]] += 1.0;
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
    while(it_end != mask.seq.end() && p < (*it_info).pos){
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

  while(it_end != std::prev(mask.seq.end(),1) && it_info != std::prev(mutations.info.end(),1)){
    if(*it_start != 'P') d_num_nonpass_vincity--;
    it_start++;
    it_end++;
    if(*it_end != 'P') d_num_nonpass_vincity++;

    assert(d_num_nonpass_vincity >= 0);

    //only add if its between the previous and the next snp (in the mutations file with all SNPs)
    if( p >= 0.5*((*it_pos) + (*std::prev(it_pos,1))) && p < 0.5 * ((*it_pos) + (*std::next(it_pos,1))) ){

      if(*it_mask == 'P' && d_num_nonpass_vincity <= mask_threshold && (*it_info).branch.size() == 1){
        //add to count
        if(IsCharNucl(*std::prev(it_ancestor,1)) != 'N' && IsCharNucl(*std::next(it_ancestor,1)) != 'N' && IsCharNucl(*it_ancestor) != 'N'){
          pattern  = toupper(*std::prev(it_ancestor,1));
          pattern += toupper(*std::next(it_ancestor,1));
          pattern += toupper(*it_ancestor);  

          int count = 0;
          for(std::string::iterator it_nucl = nucl.begin(); it_nucl != nucl.end(); it_nucl++){
            if(*it_nucl != IsCharNucl(*it_ancestor)){
              count++;
              count_bases_by_type[snp][dict_mutation_pattern[pattern + *it_nucl]] += 1.0;
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

  while(it_mask != std::prev(mask.seq.end(),1) && it_info != std::prev(mutations.info.end(),1)){
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

      if(*it_mask == 'P' && d_num_nonpass_vincity <= 0.5 * mask_threshold && (*it_info).branch.size() == 1){
        //add to count 
        if(IsCharNucl(*std::prev(it_ancestor,1)) != 'N' && IsCharNucl(*std::next(it_ancestor,1)) != 'N' && IsCharNucl(*it_ancestor) != 'N'){
          pattern  = toupper(*std::prev(it_ancestor,1));
          pattern += toupper(*std::next(it_ancestor,1));
          pattern += toupper(*it_ancestor);  

          for(std::string::iterator it_nucl = nucl.begin(); it_nucl != nucl.end(); it_nucl++){
            if(*it_nucl != IsCharNucl(*it_ancestor)){
              count_bases_by_type[snp][dict_mutation_pattern[pattern + *it_nucl]] += 1.0;
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
    std::cout << "Needed: input, output. Optional: chr, first_chr, last_chr." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Extract avg mutation rate from .bin file." << std::endl;
    exit(0);
  } 

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing mutation rate..." << std::endl;

  int num_epochs;
  std::vector<double> epoch;

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
    std::cout << "Needed: input, output. Optional: chr, first_chr, last_chr." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Extract mutation rate of 96 categories from .bin file." << std::endl;
    exit(0);
  } 

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing mutation rate..." << std::endl;

  int num_epochs;
  int num_categories = 0;
  std::vector<double> epoch;

  CollapsedMatrix<double> mutation_by_type_and_epoch, opportunity_by_type_and_epoch;
  FILE* fp; 
  fp = fopen((options["input"].as<std::string>() + "_mut" + ".bin" ).c_str(), "r");  
  assert(fp != NULL);

  fread(&num_epochs, sizeof(int), 1, fp);
  epoch.resize(num_epochs);
  fread(&epoch[0], sizeof(double), num_epochs, fp);

  mutation_by_type_and_epoch.ReadFromFile(fp);
  num_categories = mutation_by_type_and_epoch.subVectorSize(0);
  fclose(fp);
  fp = fopen((options["input"].as<std::string>() + "_opp" + ".bin" ).c_str(), "r");  
  opportunity_by_type_and_epoch.ReadFromFile(fp); 


  std::ofstream os(options["output"].as<std::string>() + ".rate");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  os << "epoch.start ";
  //for(int i = 0; i < num_categories; i++){
  //  os << i + 1 << " ";
  //}
  //os << "\n";

  std::string alphabet = "ACGT";
  //A-C, A-G, A-T, C-A, C-G, C-T
  //T-G, T-C, T-A, G-T, G-C, G-A
  //plus 16 flanks per mutation type, i.e. 6*16 = 96
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      os << *it_str1 << "C/A" << *it_str2 << " ";
      os << *it_str1 << "C/G" << *it_str2 << " ";
      os << *it_str1 << "C/T" << *it_str2 << " ";
      os << *it_str1 << "A/T" << *it_str2 << " ";
      os << *it_str1 << "A/G" << *it_str2 << " ";
      os << *it_str1 << "A/C" << *it_str2 << " ";
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
  if( ((!options.count("first_chr") || !options.count("last_chr")) && !options.count("chr")) || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: first_chr, last_chr, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .bin files and summarizes them into one .bin file." << std::endl;
    exit(0);
  }  

	std::vector<std::string> filenames;
	std::string filename_base = options["output"].as<std::string>();
	std::vector<std::string> chromosomes;


	if(options.count("chr")){
		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		std::cerr << "------------------------------------------------------" << std::endl;
		std::cerr << "Summarizing mutation rates for chr in " << options["chr"].as<std::string>() << "..." << std::endl;
		std::string line;
		while(getline(is_chr, line)){
			filenames.push_back(filename_base + "_chr" + line + "_mut" + ".bin");
			chromosomes.push_back(line);
		}
		is_chr.close();
	}else{
		int start     = options["first_chr"].as<int>(); 
		int end       = options["last_chr"].as<int>();
		std::cerr << "------------------------------------------------------" << std::endl;
		std::cerr << "Summarizing mutation rates for chr " << start << " " << end << "..." << std::endl;
		for(int chr = start; chr <= end; chr++){
			filenames.push_back(filename_base + "_chr" + std::to_string(chr) + "_mut" + ".bin");
			chromosomes.push_back(std::to_string(chr));
		}
	}

  //open files and add together

  std::cerr << filenames[0] << std::endl;
  FILE* fp = fopen(filenames[0].c_str(),"rb");
  assert(fp != NULL);  

  int num_epochs;
  std::vector<double> epochs;
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
	for(int chr = 0; chr < chromosomes.size(); chr++){
		filenames.push_back(filename_base + "_chr" + chromosomes[chr] + "_opp" + ".bin");
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

  for(int chr = 0; chr < chromosomes.size(); chr++){
		std::remove((options["input"].as<std::string>() + "_chr" + chromosomes[chr] + "_mut" + ".bin").c_str());
    std::remove((options["input"].as<std::string>() + "_chr" + chromosomes[chr] + "_opp" + ".bin").c_str());
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

void MutationRateWithContext(cxxopts::Options& options, std::string chr = "NA"){

  bool help = false;
  if(!options.count("mask") || !options.count("ancestor") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mask, ancestor, input, output. Optional: years_per_gen, bins, dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate mutation rate for 96 categories (Do not apply after using RemoveTreesWithFewMutations).." << std::endl;
    exit(0);
  }  

  std::string line, line2 , read;

  //////////// PARSE DATA ///////////

  AncMutIterators ancmut;
  if(chr == "NA"){
    ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  }else{
    ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + chr + ".anc", options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }   
  
  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0; 

  Data data(N,L);
  int N_total = 2*data.N-1;

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating mutation rate for 96 categories " << options["input"].as<std::string>() << " ..." << std::endl;

  ////////// read mutations file ///////////

  Mutations mutations;
  if(chr == "NA"){
    mutations.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mutations.Read(options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }

  std::vector<int> pos;
  if(options.count("dist")){

    int L_allsnps = 0;
    igzstream is_L(options["dist"].as<std::string>());
    std::string unused;
    std::getline(is_L, unused); 
    while ( std::getline(is_L, unused) ){
      ++L_allsnps;
    }
    is_L.close();

    pos.resize(L_allsnps);
    igzstream is_dist(options["dist"].as<std::string>());
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
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }
 
  int num_epochs;
  std::vector<double> epochs;
  float log_10 = std::log(10);
  if(options.count("bins")){

    double log_age = std::log(0);
    double age = 0;
 
    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      ep++;
    }
    epoch_boundary = epoch_lower;
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0) epochs.push_back(age);
        epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();

  }else{

    num_epochs = 31;
    epochs.resize(num_epochs);
    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }

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
      dict_mutation_pattern[pattern + "CA"] = index;
      index++;
      dict_mutation_pattern[pattern + "CG"] = index;
      index++;
      dict_mutation_pattern[pattern + "CT"] = index;
      index++;
      dict_mutation_pattern[pattern + "AT"] = index;
      index++;
      dict_mutation_pattern[pattern + "AG"] = index;
      index++;
      dict_mutation_pattern[pattern + "AC"] = index;
      index++;
    }   
  }
  index = 0;
  for(std::string::iterator it_str1 = reverse_alphabet.begin(); it_str1 != reverse_alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = reverse_alphabet.begin(); it_str2 != reverse_alphabet.end(); it_str2++){
      reverse_pattern = *it_str2;
      reverse_pattern += *it_str1;
      dict_mutation_pattern[reverse_pattern + "GT"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "GC"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "GA"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "TA"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "TC"] = index;
      index++;
      dict_mutation_pattern[reverse_pattern + "TG"] = index;
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

  //Tree subtr;
  std::vector<float> coordinates_tree(N_total);
  std::vector<int> num_lineages(N_total);
  int root = N_total-1;
  int i = 0;

  //read tree
  float num_bases_SNP_persists;
  num_bases_SNP_persists = ancmut.FirstSNP(mtr, it_mut);

  GetCoordsAndLineages(mtr, coordinates_tree, num_lineages);
  GetBranchLengthsInEpoch(data, epochs, coordinates_tree, num_lineages, branch_lengths_in_epoch);

  SNPInfo snp_info;
  float rec;
  int current_tree = (*it_mut).tree;
  double total_branch_length = 0.0;
  int count_snps = 0;
  for(int snp = 0; snp < data.L; snp++){

    snp_info = (*it_mut);
    if(snp_info.branch.size() == 1){
     
      if((*it_mut).tree != current_tree){
        current_tree = (*it_mut).tree;
        mtr.tree.GetCoordinates(coordinates_tree);
        GetCoordsAndLineages(mtr, coordinates_tree, num_lineages);
        GetBranchLengthsInEpoch(data, epochs, coordinates_tree, num_lineages, branch_lengths_in_epoch);
      }
      assert(current_tree == snp_info.tree);

      if(snp_info.upstream_base != "NA" && snp_info.downstream_base != "NA" && snp_info.mutation_type[0] != snp_info.mutation_type[2]){

        if(snp_info.mutation_type.size() == 3){

          if(snp_info.mutation_type[0] == 'A' || snp_info.mutation_type[0] == 'C' || snp_info.mutation_type[0] == 'G' || snp_info.mutation_type[0] == 'T'){

            if(snp_info.mutation_type[2] == 'A' || snp_info.mutation_type[2] == 'C' || snp_info.mutation_type[2] == 'G' || snp_info.mutation_type[2] == 'T'){

              // identify cathegory of mutation
              pattern = snp_info.upstream_base + snp_info.downstream_base + snp_info.mutation_type[0] + snp_info.mutation_type[2];
              // then identify its age and number of lineages at the time
              int ind = dict_mutation_pattern[pattern];

              // identify epoch and add to number of mutations per lineage in epoch
              int ep = 0;
              while(epochs[ep] <= snp_info.age_begin){
                ep++;
                if(ep == epochs.size()) break;
              }
              ep--;

              assert(ep >= 0);
              assert(mutation_by_type_and_epoch.subVectorSize(ep) == num_mutation_cathegories);
              //assert(count_bases_by_type[snp][ind] > 0);

              int ep_begin = ep;
              float age_end = std::min(snp_info.age_end,coordinates_tree[root]);
              assert(age_end < epochs[num_epochs-1]);
              double branch_length = age_end - snp_info.age_begin;

              if(age_end <= epochs[ep+1]){

                mutation_by_type_and_epoch[ep][ind] += 1.0;

              }else{

                mutation_by_type_and_epoch[ep][ind]   += (epochs[ep+1] - snp_info.age_begin)/branch_length;
                ep++;
                while(epochs[ep+1] <= age_end){
                  mutation_by_type_and_epoch[ep][ind] += (epochs[ep+1]-epochs[ep])/branch_length;
                  ep++;
                }
                mutation_by_type_and_epoch[ep][ind]   += (age_end-epochs[ep])/branch_length;

              }
 
              for(int ep_tmp = 0; ep_tmp < num_epochs; ep_tmp++){
                double bl = branch_lengths_in_epoch[ep_tmp];
                for(int ind_tmp = 0; ind_tmp < num_mutation_cathegories; ind_tmp++){
                  opportunity_by_type_and_epoch[ep_tmp][ind_tmp] += bl * count_bases_by_type[snp][ind_tmp];
                }
              }

              //count_snps++;
              //
            }

          }

        }

      }

    }
    num_bases_SNP_persists = ancmut.NextSNP(mtr, it_mut);

  }

  for(int i = 0; i < num_epochs; i++){
    std::cerr << mutation_by_type_and_epoch[i][0] << " " << opportunity_by_type_and_epoch[i][0] << std::endl;
  }

  //output mutation_by_type_and_epoch
  //       opportunity_by_type_and_epoch 
  FILE* fp;
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_mut" + ".bin" ).c_str(), "wb");  
  }
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(double), epochs.size(), fp);
  mutation_by_type_and_epoch.DumpToFile(fp);
  fclose(fp);
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_opp" + ".bin" ).c_str(), "wb");  
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


////////////////////
void MutationRateForCategory(cxxopts::Options& options, std::string chr = "NA"){

  bool help = false;
  if(!options.count("mask") || !options.count("ancestor") || !options.count("mutcat") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mask, ancestor, mutcat, input, output. Optional: years_per_gen, bins, dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate mutation rate for categories (Do not apply after using RemoveTreesWithFewMutations).." << std::endl;
    exit(0);
  }  

  std::string line, line2 , read;

  //////////// PARSE DATA ///////////

  MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;

  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.
  AncMutIterators ancmut;

  if(chr == "NA"){
    ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  }else{
    ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + chr + ".anc", options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  int N = (mtr.tree.nodes.size() + 1)/2.0;
  int L = ancmut.NumSnps();
  Data data(N,L);
  int N_total = 2*data.N-1;

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating mutation rate for categories " << options["input"].as<std::string>() << " ..." << std::endl;

  ////////// read mutations file ///////////

  Mutations mutations;
  if(chr == "NA"){
    mutations.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mutations.Read(options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }

  std::vector<int> pos;
  if(options.count("dist")){

    int L_allsnps = 0;
    igzstream is_L;
    is_L.open(options["dist"].as<std::string>());
    if(is_L.fail()){
      std::cerr << "Error opening file " << options["dist"].as<std::string>() << std::endl;
      exit(1);
    }
    std::string unused;
    std::getline(is_L, unused); 
    while ( std::getline(is_L, unused) ){
      ++L_allsnps;
    }
    is_L.close();

    pos.resize(L_allsnps);
    igzstream is_dist(options["dist"].as<std::string>());
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
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }
 
  int num_epochs;
  std::vector<double> epochs;
  float log_10 = std::log(10);
  if(options.count("bins")){

    double log_age = std::log(0);
    double age = 0;
  
    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      ep++;
    }
    epoch_boundary = epoch_lower;
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0) epochs.push_back(age);
        epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();

  }else{

    num_epochs = 31;
    epochs.resize(num_epochs);
    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }


  /////////////////////////////////////////////////////////////////
  //Mutation specific

  ////////// define mutation types /////////
  std::string alphabet = "ACGT";
  std::map<char, std::string> complement;
  complement['A'] = "T";
  complement['C'] = "G";
  complement['G'] = "C";
  complement['T'] = "A";
  std::map<std::string, int> dict_mutation_pattern;

  //A-C, A-G, A-T, C-A, C-G, C-T
  //T-G, T-C, T-A, G-T, G-C, G-A
  //plus 16 flanks per mutation type, i.e. 6*16 = 96

  //parse file storing
  //upstream downstream ancestral derived category
  igzstream is_cat(options["mutcat"].as<std::string>());
  if(is_cat.fail()){
    std::cerr << "Error: unable to open file " << options["mutcat"].as<std::string>() << std::endl;
  }
  getline(is_cat, line);
  //I have to make sure all 96 categories are represented in this file

  char mutation_type[5];
  std::string pattern, reverse_pattern;
  int category, num_categories = 0;
  std::vector<int> check_num_categories;
  while(getline(is_cat, line)){
    sscanf(line.c_str(), "%c %c %c %c %d", &mutation_type[0], &mutation_type[1], &mutation_type[2], &mutation_type[3], &category);
    pattern = mutation_type;
    dict_mutation_pattern[pattern] = category;
    pattern = complement[mutation_type[1]] + complement[mutation_type[0]] + complement[mutation_type[2]] + complement[mutation_type[3]];
    dict_mutation_pattern[pattern] = category;
    if(category >= num_categories){
      check_num_categories.resize(category+1);
      std::fill(std::next(check_num_categories.begin(),num_categories), check_num_categories.end(), 0);
      num_categories = category + 1;      
      check_num_categories[category]++;
    }else{
      check_num_categories[category]++;
    }
  }
  is_cat.close();

  for(std::vector<int>::iterator it_check = check_num_categories.begin(); it_check != check_num_categories.end(); it_check++){
    if(*it_check == 0){
      std::cerr << "Error: category indices not 0-indexed or contiguous." << std::endl;
      exit(1);
    }
  }

  //check that I got all 96 categories and that categories are contiguous integers
  int index = 0;
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      pattern = *it_str1;
      pattern += *it_str2;
      reverse_pattern = complement[*it_str2] + complement[*it_str1];

      if ( dict_mutation_pattern.find(pattern + "CA") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "GT") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "CG") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "GC") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "CT") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "AG") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AT") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TA") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AG") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TC") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AC") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TG") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
    }   
  }

  ///////// Count number of bases by type //////////

  CollapsedMatrix<double> count_bases_by_type;
  //count_bases_by_type.resize(mutations.info.size(),dict_mutation_pattern.size());
  CountBasesByType(data, options["mask"].as<std::string>(), options["ancestor"].as<std::string>(), count_bases_by_type, dict_mutation_pattern, mutations, pos);
  fasta mask;
  mask.Read(options["mask"].as<std::string>());

  ////////////////////////////////////////////////////////////////////////
  // Estimate mutation rate through time

  std::vector<CollapsedMatrix<double>> mutation_by_type_and_epoch(ancmut.NumTrees());
  std::vector<CollapsedMatrix<double>> opportunity_by_type_and_epoch(ancmut.NumTrees());
  for(std::vector<CollapsedMatrix<double>>::iterator it_m = mutation_by_type_and_epoch.begin(); it_m != mutation_by_type_and_epoch.end(); it_m++){
    (*it_m).resize(num_epochs, num_categories);
  }
  for(std::vector<CollapsedMatrix<double>>::iterator it_o = opportunity_by_type_and_epoch.begin(); it_o != opportunity_by_type_and_epoch.end(); it_o++){
    (*it_o).resize(num_epochs, num_categories);
  } 

  std::vector<double> branch_lengths_in_epoch(num_epochs);

  //Tree subtr;
  std::vector<float> coordinates_tree(N_total);
  std::vector<int> num_lineages(N_total);
  int root = N_total-1;
  int i = 0;

  //iterate over trees
  SNPInfo snp_info;
  float rec;
  int num_tree = 0;
  double total_branch_length = 0.0;
  int count_snps = 0;

  int snp = 0;
  while(num_bases_tree_persists >= 0){

    mtr.tree.GetCoordinates(coordinates_tree);
    GetCoordsAndLineages(mtr, coordinates_tree, num_lineages);
    GetBranchLengthsInEpoch(data, epochs, coordinates_tree, num_lineages, branch_lengths_in_epoch);
    num_tree = mutations.info[snp].tree;

    while(num_tree == mutations.info[snp].tree){

      snp_info = mutations.info[snp];

      if(snp_info.branch.size() == 1 && mask.seq[snp_info.pos-1] != 'N'){

        assert(ancmut.get_treecount() == snp_info.tree);

        //if statements to make sure we have a well defined biallelic SNP
        if(snp_info.upstream_base != "NA" && snp_info.downstream_base != "NA" && snp_info.mutation_type[0] != snp_info.mutation_type[2]){

          if(snp_info.mutation_type.size() == 3){

            if(snp_info.mutation_type[0] == 'A' || snp_info.mutation_type[0] == 'C' || snp_info.mutation_type[0] == 'G' || snp_info.mutation_type[0] == 'T'){

              if(snp_info.mutation_type[2] == 'A' || snp_info.mutation_type[2] == 'C' || snp_info.mutation_type[2] == 'G' || snp_info.mutation_type[2] == 'T'){

                // identify category of mutation
                pattern = snp_info.upstream_base + snp_info.downstream_base + snp_info.mutation_type[0] + snp_info.mutation_type[2];
                // then identify its age and number of lineages at the time
                assert(dict_mutation_pattern.find(pattern) != dict_mutation_pattern.end());
                int ind = dict_mutation_pattern[pattern];

                // identify epoch and add to number of mutations per lineage in epoch
                int ep = 0;
                while(epochs[ep] <= snp_info.age_begin){
                  ep++;
                  if(ep == epochs.size()) break;
                }
                ep--;

                assert(ep >= 0);
                assert(mutation_by_type_and_epoch[num_tree].subVectorSize(ep) == num_categories);

                int ep_begin = ep;
                float age_end = std::min(snp_info.age_end,coordinates_tree[root]);
                assert(age_end < epochs[num_epochs-1]);
                double branch_length = age_end - snp_info.age_begin;

                if(age_end <= epochs[ep+1]){

                  mutation_by_type_and_epoch[num_tree][ep][ind] += 1.0;

                }else{

                  mutation_by_type_and_epoch[num_tree][ep][ind]   += (epochs[ep+1] - snp_info.age_begin)/branch_length;
                  ep++;
                  while(epochs[ep+1] <= age_end){
                    mutation_by_type_and_epoch[num_tree][ep][ind] += (epochs[ep+1]-epochs[ep])/branch_length;
                    ep++;
                  }
                  mutation_by_type_and_epoch[num_tree][ep][ind]   += (age_end-epochs[ep])/branch_length;

                }

                for(int ep_tmp = 0; ep_tmp < num_epochs; ep_tmp++){
                  double bl = branch_lengths_in_epoch[ep_tmp];
                  assert(bl >= 0.0);
                  for(int ind_tmp = 0; ind_tmp < num_categories; ind_tmp++){
                    //std::cerr << count_bases_by_type[snp][ind_tmp] << " ";
                    opportunity_by_type_and_epoch[num_tree][ep_tmp][ind_tmp] += bl * count_bases_by_type[snp][ind_tmp];
                  }
                  //std::cerr << std::endl;
                }

              }

            }

          } 

        }

      }

      snp++;

    }

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

  }


  //bootstrap
  //sample indices from 0 too ancmut.NumTrees()-1 at random with replacement
  //dump to files and write summarise and finalise functions for bootstrap

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> sam(0, (ancmut.NumTrees()-1.0)/1000.0);

  int n_boot = 100;
  std::vector<CollapsedMatrix<double>> boot_mutation_by_type_and_epoch(n_boot);
  std::vector<CollapsedMatrix<double>> boot_opportunity_by_type_and_epoch(n_boot);
  for(std::vector<CollapsedMatrix<double>>::iterator it_m = boot_mutation_by_type_and_epoch.begin(); it_m != boot_mutation_by_type_and_epoch.end(); it_m++){
    (*it_m).resize(num_epochs, num_categories);
  }
  for(std::vector<CollapsedMatrix<double>>::iterator it_o = boot_opportunity_by_type_and_epoch.begin(); it_o != boot_opportunity_by_type_and_epoch.end(); it_o++){
    (*it_o).resize(num_epochs, num_categories);
  } 

  for(int n = 0; n < n_boot; n++){
  
    std::vector<int> boot_trees(ancmut.NumTrees());
    if(0){
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end(); it_boot_trees++){
      *it_boot_trees = sam(gen);
    }
    }

    int size = 0;
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end();){
      int start = 1000*sam(gen);
      for(int k = start; k < start + 1000 && size < boot_trees.size() && k < boot_trees.size(); k++){
        *it_boot_trees = k;
        it_boot_trees++;
        size++;
      }
    }
    boot_trees.resize(size);
    //std::sort(boot_trees.begin(), boot_trees.end());

    //fill in boot_mutation_by_type_and_epoch[n] and boot_opportunity_by_type_and_epoch[n] by summing over the trees  
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end(); it_boot_trees++){
      std::vector<double>::iterator it_bmut_by_type_and_epoch = boot_mutation_by_type_and_epoch[n].vbegin();
      std::vector<double>::iterator it_bopp_by_type_and_epoch = boot_opportunity_by_type_and_epoch[n].vbegin();
      std::vector<double>::iterator it_mut_by_type_and_epoch  = mutation_by_type_and_epoch[*it_boot_trees].vbegin();
      std::vector<double>::iterator it_opp_by_type_and_epoch  = opportunity_by_type_and_epoch[*it_boot_trees].vbegin();
      for(; it_mut_by_type_and_epoch != mutation_by_type_and_epoch[*it_boot_trees].vend();){    
        *it_bmut_by_type_and_epoch += *it_mut_by_type_and_epoch;
        *it_bopp_by_type_and_epoch += *it_opp_by_type_and_epoch;
        it_mut_by_type_and_epoch++;
        it_opp_by_type_and_epoch++;      
        it_bmut_by_type_and_epoch++;
        it_bopp_by_type_and_epoch++;
      }
    }
  
  }

  //output mutation_by_type_and_epoch
  //       opportunity_by_type_and_epoch 
  FILE* fp;
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_mut" + ".bin" ).c_str(), "wb");  
  }
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(double), epochs.size(), fp);
  for(int n = 0; n < n_boot; n++){
    boot_mutation_by_type_and_epoch[n].DumpToFile(fp);
  }
  fclose(fp);
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_opp" + ".bin" ).c_str(), "wb");  
  }
  for(int n = 0; n < n_boot; n++){
    boot_opportunity_by_type_and_epoch[n].DumpToFile(fp);
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


}

void MutationRateForCategoryForGroup(cxxopts::Options& options, std::string chr = "NA"){

  bool help = false;
  if(!options.count("mask") || !options.count("ancestor") || !options.count("mutcat") || !options.count("pop_of_interest") || !options.count("poplabels") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mask, ancestor, mutcat, input, output, pop_of_interest, poplabels. Optional: years_per_gen, bins, dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate mutation rate for categories (Do not apply after using RemoveTreesWithFewMutations).." << std::endl;
    exit(0);
  }  

  std::string line, line2 , read;

  //////////// PARSE DATA ///////////

  MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;

  Sample samples;
  samples.Read(options["poplabels"].as<std::string>());

  std::string label;
  if(!options.count("pop_of_interest")){
    label = samples.AssignPopOfInterest("All");
  }else{
    label = samples.AssignPopOfInterest(options["pop_of_interest"].as<std::string>());
  }

  if(label.compare("All")){
    for(std::vector<int>::iterator it_group_of_interest = samples.group_of_interest.begin(); it_group_of_interest != samples.group_of_interest.end(); it_group_of_interest++ ){
      std::cerr << samples.groups[*it_group_of_interest] << " ";
    }
    std::cerr << std::endl;
  }else{
    std::cerr << "All" << std::endl;
  }

  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.
  AncMutIterators ancmut;

  if(chr == "NA"){
    ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  }else{
    ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + chr + ".anc", options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  Data data(N,L);
  int N_total = 2*data.N-1;

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating mutation rate for categories " << options["input"].as<std::string>() << " ..." << std::endl;

  ////////// read mutations file ///////////

  Mutations mutations;
  if(chr == "NA"){
    mutations.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mutations.Read(options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }

  std::vector<int> pos;
  if(options.count("dist")){

    int L_allsnps = 0;
    igzstream is_L;
    is_L.open(options["dist"].as<std::string>());
    if(is_L.fail()){
      std::cerr << "Error opening file " << options["dist"].as<std::string>() << std::endl;
      exit(1);
    }
    std::string unused;
    std::getline(is_L, unused); 
    while ( std::getline(is_L, unused) ){
      ++L_allsnps;
    }
    is_L.close();

    pos.resize(L_allsnps);
    igzstream is_dist(options["dist"].as<std::string>());
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
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }
 
  int num_epochs;
  std::vector<double> epochs;
  float log_10 = std::log(10);
  if(options.count("binsfile")){
  
    std::ifstream is(options["binsfile"].as<std::string>());
    std::string line;
    while(getline(is, line)){
      if(epochs.size() == 0){
        if(std::stof(line) > 0){
          epochs.push_back(0);
        }
      }
      epochs.push_back(std::stof(line));
    }
    if(epochs[epochs.size()-1] < 1e8) epochs.push_back(1e8);
    num_epochs = epochs.size();

    for(int i = 0; i < num_epochs - 1; i++){
      assert(epochs[i] >= 0);
      assert(epochs[i] <= epochs[i+1]);
    }

  }else if(options.count("bins")){

    double age = 0.0;
    if(options.count("sample_age") > 0) age = options["sample_age"].as<float>();
    double count = 0.0;
    double log_age = std::log(age * years_per_gen)/log_10;

    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    epoch_boundary = epoch_lower;
    
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      if(epoch_boundary - log_age < 0.5*epoch_step) epoch_boundary += epoch_step;
      ep++;
    }
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0){
          epochs.push_back(age);
          if(epoch_boundary - log_age < 0.5*epoch_step) epoch_boundary += epoch_step;
        }
        if(std::fabs(log_age - epoch_boundary) > 1e-3){
          epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        }
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10*epochs[epochs.size()-1])/years_per_gen );
    num_epochs = epochs.size();	
    
  }else{

    num_epochs = 31;
    epochs.resize(num_epochs);
    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }

  //for(int i = 0; i < num_epochs; i++){
  //  std::cerr << epochs[i] << " ";
  //}
  //std::cerr << std::endl;

  /////////////////////////////////////////////////////////////////
  //Mutation specific

  ////////// define mutation types /////////
  std::string alphabet = "ACGT";
  std::map<char, std::string> complement;
  complement['A'] = "T";
  complement['C'] = "G";
  complement['G'] = "C";
  complement['T'] = "A";
  std::map<std::string, int> dict_mutation_pattern;

  //A-C, A-G, A-T, C-A, C-G, C-T
  //T-G, T-C, T-A, G-T, G-C, G-A
  //plus 16 flanks per mutation type, i.e. 6*16 = 96

  //parse file storing
  //upstream downstream ancestral derived category
  igzstream is_cat(options["mutcat"].as<std::string>());
  if(is_cat.fail()){
    std::cerr << "Error: unable to open file " << options["mutcat"].as<std::string>() << std::endl;
  }
  getline(is_cat, line);
  //I have to make sure all 96 categories are represented in this file

  char mutation_type[5];
  std::string pattern, reverse_pattern;
  int category, num_categories = 0;
  std::vector<int> check_num_categories;
  while(getline(is_cat, line)){
    sscanf(line.c_str(), "%c %c %c %c %d", &mutation_type[0], &mutation_type[1], &mutation_type[2], &mutation_type[3], &category);
    pattern = mutation_type;
    dict_mutation_pattern[pattern] = category;
    pattern = complement[mutation_type[1]] + complement[mutation_type[0]] + complement[mutation_type[2]] + complement[mutation_type[3]];
    dict_mutation_pattern[pattern] = category;
    if(category >= num_categories){
      check_num_categories.resize(category+1);
      std::fill(std::next(check_num_categories.begin(),num_categories), check_num_categories.end(), 0);
      num_categories = category + 1;      
      check_num_categories[category]++;
    }else{
      check_num_categories[category]++;
    }
  }
  is_cat.close();

  for(std::vector<int>::iterator it_check = check_num_categories.begin(); it_check != check_num_categories.end(); it_check++){
    if(*it_check == 0){
      std::cerr << "Error: category indices not 0-indexed or contiguous." << std::endl;
      exit(1);
    }
  }

  //check that I got all 96 categories and that categories are contiguous integers
  int index = 0;
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      pattern = *it_str1;
      pattern += *it_str2;
      reverse_pattern = complement[*it_str2] + complement[*it_str1];

      if ( dict_mutation_pattern.find(pattern + "CA") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "GT") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "CG") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "GC") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "CT") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "AG") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AT") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TA") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AG") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TC") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AC") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TG") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
    }   
  }

  ///////// Count number of bases by type //////////

  CollapsedMatrix<double> count_bases_by_type;
  //count_bases_by_type.resize(mutations.info.size(),dict_mutation_pattern.size());
  CountBasesByType(data, options["mask"].as<std::string>(), options["ancestor"].as<std::string>(), count_bases_by_type, dict_mutation_pattern, mutations, pos);
  fasta mask;
  mask.Read(options["mask"].as<std::string>());

  ////////////////////////////////////////////////////////////////////////
  // Estimate mutation rate through time

  std::vector<CollapsedMatrix<double>> mutation_by_type_and_epoch(ancmut.NumTrees());
  std::vector<CollapsedMatrix<double>> opportunity_by_type_and_epoch(ancmut.NumTrees());
  for(std::vector<CollapsedMatrix<double>>::iterator it_m = mutation_by_type_and_epoch.begin(); it_m != mutation_by_type_and_epoch.end(); it_m++){
    (*it_m).resize(num_epochs, num_categories);
  }
  for(std::vector<CollapsedMatrix<double>>::iterator it_o = opportunity_by_type_and_epoch.begin(); it_o != opportunity_by_type_and_epoch.end(); it_o++){
    (*it_o).resize(num_epochs, num_categories);
  } 

  std::vector<double> branch_lengths_in_epoch(num_epochs);
  //Tree subtr;
  std::vector<float> coordinates_tree(N_total);
  std::vector<int> num_lineages(N_total);
  std::vector<Leaves> descendants;
  int root = N_total-1;
  int i = 0;

  //iterate over trees
  SNPInfo snp_info;
  float rec;
  int num_tree = 0;
  double total_branch_length = 0.0;
  int count_snps = 0;

  std::vector<int> exclude;
  if(0){
    std::vector<std::string> id = {"LBK", "NE1", "BR2", "JP14", "NG10", "PB675", "atp016"};
    //std::vector<std::string> id = {"KK1", "Loschbour", "sf12", "SRA62"};
    //std::vector<std::string> id = {"WestEurasia"};
    for(int j = 0; j < id.size(); j++){
      i = 0;
      for(; i < samples.groups.size(); i++){
        if(samples.groups[i] == id[j]) break;
      }
      if(i == 0 && samples.groups[0] != id[j]){
        std::cerr << "Group to exclude not found" << std::endl;
        exit(1);
      }
      exclude.push_back(i);
    }
  }

  int snp = 0;
  while(num_bases_tree_persists >= 0){

    //need to get branch lengths in epoch for only the pop of interest
    mtr.tree.GetCoordinates(coordinates_tree);
    mtr.tree.FindAllLeaves(descendants);
    GetCoordsAndLineagesForPop(mtr, samples, exclude, descendants, coordinates_tree, num_lineages);
    GetBranchLengthsInEpoch(data, epochs, coordinates_tree, num_lineages, branch_lengths_in_epoch);
    num_tree = mutations.info[snp].tree;

    //for(int i = 0; i < num_lineages.size(); i++){
    //  std::cerr << num_lineages[i] << " ";
    //}
    //std::cerr << std::endl;

    while(num_tree == mutations.info[snp].tree){

      snp_info = mutations.info[snp];

      if(snp_info.branch.size() == 1 && mask.seq[snp_info.pos-1] != 'N'){

        assert(ancmut.get_treecount() == snp_info.tree);
        //check for mutation whether it is segregating in pop_of_interest
        bool use = false, excl = false;

        if(descendants[snp_info.branch[0]].num_leaves > 1){
          for(int i = 0; i < samples.group_of_interest.size(); i++){
            for(std::vector<int>::iterator it_mem = descendants[snp_info.branch[0]].member.begin(); it_mem != descendants[snp_info.branch[0]].member.end(); it_mem++){
              if(i == 0){
                for(int j = 0; j < exclude.size(); j++){
                  if(samples.group_of_haplotype[*it_mem] == exclude[j]){
                    excl = true;
                    break;
                  }
                }
                if(excl){
                  break;
                }
              }
              if(samples.group_of_haplotype[*it_mem] == samples.group_of_interest[i]){
                use = true;
                if(i != 0) break;
              }
            }
            if(excl){
              use = false;
              break;
            }
          } 
        }

        //if statements to make sure we have a well defined biallelic SNP
        if(use && snp_info.upstream_base != "NA" && snp_info.downstream_base != "NA" && snp_info.mutation_type[0] != snp_info.mutation_type[2]){

          if(snp_info.mutation_type.size() == 3){

            if(snp_info.mutation_type[0] == 'A' || snp_info.mutation_type[0] == 'C' || snp_info.mutation_type[0] == 'G' || snp_info.mutation_type[0] == 'T'){

              if(snp_info.mutation_type[2] == 'A' || snp_info.mutation_type[2] == 'C' || snp_info.mutation_type[2] == 'G' || snp_info.mutation_type[2] == 'T'){

                if(0){
                if(snp_info.upstream_base == "T" && snp_info.downstream_base == "C" && snp_info.mutation_type == "C/T"){
                  std::cerr << snp_info.pos << " " << snp_info.tree << " " << snp_info.branch[0] << std::endl;

                  std::cerr << "Excl: " << exclude[0] << std::endl;
                  for(int i = 0; i < descendants[snp_info.branch[0]].member.size(); i++){
                    std::cerr << samples.group_of_haplotype[descendants[snp_info.branch[0]].member[i]] << " ";
                  }
                  std::cerr << std::endl;


                }
                }

                // identify category of mutation
                pattern = snp_info.upstream_base + snp_info.downstream_base + snp_info.mutation_type[0] + snp_info.mutation_type[2];
                // then identify its age and number of lineages at the time
                assert(dict_mutation_pattern.find(pattern) != dict_mutation_pattern.end());
                int ind = dict_mutation_pattern[pattern];

                // identify epoch and add to number of mutations per lineage in epoch
                int ep = 0;
                while(epochs[ep] <= snp_info.age_begin){
                  ep++;
                  if(ep == epochs.size()) break;
                }
                ep--;

                assert(ep >= 0);
                assert(mutation_by_type_and_epoch[num_tree].subVectorSize(ep) == num_categories);

                int ep_begin = ep;
                double age_end = std::min(snp_info.age_end,coordinates_tree[root]);
                //if(age_end > epochs[num_epochs-1]){
                //  std::cerr << snp << " " << snp_info.age_end << " " << coordinates_tree[root] << std::endl;
                //}
                assert(age_end <= epochs[num_epochs-1]);
                double branch_length = age_end - snp_info.age_begin;

                if(age_end <= epochs[ep+1]){

                  mutation_by_type_and_epoch[num_tree][ep][ind] += 1.0;

                }else{

                  mutation_by_type_and_epoch[num_tree][ep][ind]   += (epochs[ep+1] - snp_info.age_begin)/branch_length;
                  ep++;
                  while(epochs[ep+1] <= age_end){
                    mutation_by_type_and_epoch[num_tree][ep][ind] += (epochs[ep+1]-epochs[ep])/branch_length;
                    ep++;
                  }
                  mutation_by_type_and_epoch[num_tree][ep][ind]   += (age_end-epochs[ep])/branch_length;

                }

                for(int ep_tmp = 0; ep_tmp < num_epochs; ep_tmp++){
                  double bl = branch_lengths_in_epoch[ep_tmp];
                  assert(bl >= 0.0);
                  for(int ind_tmp = 0; ind_tmp < num_categories; ind_tmp++){
                    //std::cerr << opportunity_by_type_and_epoch[num_tree][ep_tmp][ind_tmp] << " ";
                    opportunity_by_type_and_epoch[num_tree][ep_tmp][ind_tmp] += bl * count_bases_by_type[snp][ind_tmp];
                  }
                  //std::cerr << std::endl;
                }

              }

            }

          } 

        }

      }

      snp++;

    }

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

  }


  //bootstrap
  //sample indices from 0 too ancmut.NumTrees()-1 at random with replacement
  //dump to files and write summarise and finalise functions for bootstrap

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  if(options.count("seed") > 0){
    gen.seed(options["seed"].as<int>());
  }
  std::uniform_int_distribution<> sam(0, (ancmut.NumTrees()-1.0)/1000.0);

  int n_boot = 100;
  std::vector<CollapsedMatrix<double>> boot_mutation_by_type_and_epoch(n_boot);
  std::vector<CollapsedMatrix<double>> boot_opportunity_by_type_and_epoch(n_boot);
  for(std::vector<CollapsedMatrix<double>>::iterator it_m = boot_mutation_by_type_and_epoch.begin(); it_m != boot_mutation_by_type_and_epoch.end(); it_m++){
    (*it_m).resize(num_epochs, num_categories);
  }
  for(std::vector<CollapsedMatrix<double>>::iterator it_o = boot_opportunity_by_type_and_epoch.begin(); it_o != boot_opportunity_by_type_and_epoch.end(); it_o++){
    (*it_o).resize(num_epochs, num_categories);
  } 

  for(int n = 0; n < n_boot; n++){
  
    std::vector<int> boot_trees(ancmut.NumTrees());
    if(0){
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end(); it_boot_trees++){
      *it_boot_trees = sam(gen);
    }
    }

    int size = 0;
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end();){
      int start = 1000*sam(gen);
      std::cerr << start << " ";
      for(int k = start; k < start + 1000 && size < boot_trees.size() && k < boot_trees.size(); k++){
        *it_boot_trees = k;
        it_boot_trees++;
        size++;
      }
    }
    std::cerr << std::endl;
    boot_trees.resize(size);
    //std::sort(boot_trees.begin(), boot_trees.end());

    //fill in boot_mutation_by_type_and_epoch[n] and boot_opportunity_by_type_and_epoch[n] by summing over the trees  
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end(); it_boot_trees++){
      std::vector<double>::iterator it_bmut_by_type_and_epoch = boot_mutation_by_type_and_epoch[n].vbegin();
      std::vector<double>::iterator it_bopp_by_type_and_epoch = boot_opportunity_by_type_and_epoch[n].vbegin();
      std::vector<double>::iterator it_mut_by_type_and_epoch  = mutation_by_type_and_epoch[*it_boot_trees].vbegin();
      std::vector<double>::iterator it_opp_by_type_and_epoch  = opportunity_by_type_and_epoch[*it_boot_trees].vbegin();
      for(; it_mut_by_type_and_epoch != mutation_by_type_and_epoch[*it_boot_trees].vend();){    
        *it_bmut_by_type_and_epoch += *it_mut_by_type_and_epoch;
        *it_bopp_by_type_and_epoch += *it_opp_by_type_and_epoch;
        it_mut_by_type_and_epoch++;
        it_opp_by_type_and_epoch++;      
        it_bmut_by_type_and_epoch++;
        it_bopp_by_type_and_epoch++;
      }
    }
  
  }

  //output mutation_by_type_and_epoch
  //       opportunity_by_type_and_epoch 
  FILE* fp;
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_mut" + ".bin" ).c_str(), "wb");  
  }
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(double), epochs.size(), fp);
  for(int n = 0; n < n_boot; n++){
    boot_mutation_by_type_and_epoch[n].DumpToFile(fp);
  }
  fclose(fp);
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_opp" + ".bin" ).c_str(), "wb");  
  }
  for(int n = 0; n < n_boot; n++){
    boot_opportunity_by_type_and_epoch[n].DumpToFile(fp);
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


}

void SummarizeWholeGenomeForCategory(cxxopts::Options& options){

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
    std::cout << "Reads .bin files and summarizes them into one .bin file." << std::endl;
    exit(0);
  }  

	std::vector<std::string> filenames;
	std::vector<std::string> chromosomes;
	std::string filename_base = options["output"].as<std::string>();

	if(options.count("chr")){
		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		std::cerr << "------------------------------------------------------" << std::endl;
		std::cerr << "Summarizing mutation rates for chr in " << options["chr"].as<std::string>() << "..." << std::endl;
		std::string line;
		while(getline(is_chr, line)){
			filenames.push_back(filename_base + "_chr" + line + "_mut" + ".bin");
			chromosomes.push_back(line);
		}
		is_chr.close();
	}else{
		int start     = options["first_chr"].as<int>(); 
		int end       = options["last_chr"].as<int>();
		std::cerr << "------------------------------------------------------" << std::endl;
		std::cerr << "Summarizing mutation rates for chr " << start << " " << end << "..." << std::endl;
		for(int chr = start; chr <= end; chr++){
			filenames.push_back(filename_base + "_chr" + std::to_string(chr) + "_mut" + ".bin");
			chromosomes.push_back(std::to_string(chr));
		}
	}

  //open files and add together

  FILE* fp = fopen(filenames[0].c_str(),"rb");
  assert(fp != NULL);  

  int num_epochs;
  std::vector<double> epochs;
  fread(&num_epochs, sizeof(int), 1, fp);
  epochs.resize(num_epochs);
  fread(&epochs[0], sizeof(double), num_epochs, fp);

  std::cerr << num_epochs << std::endl;
  for(int e = 0; e < num_epochs; e++){
    std::cerr << epochs[e] << " ";
  }
  std::cerr << std::endl;

  int n_boot = 100;
  std::vector<CollapsedMatrix<double>> mut_by_type_and_epoch(n_boot);
  CollapsedMatrix<double> mut_by_type_and_epoch_tmp;
  for(int n = 0; n < n_boot; n++){
    mut_by_type_and_epoch[n].ReadFromFile(fp);
  }
  fclose(fp);
  for(int i = 1; i < (int) filenames.size(); i++){
    fp = fopen(filenames[i].c_str(),"rb");

    fread(&num_epochs, sizeof(int), 1, fp);
    fread(&epochs[0], sizeof(double), num_epochs, fp);

    for(int n = 0; n < n_boot; n++){
      mut_by_type_and_epoch_tmp.ReadFromFile(fp);
      std::vector<double>::iterator it_mut = mut_by_type_and_epoch[n].vbegin();
      for(std::vector<double>::iterator it_mut_tmp = mut_by_type_and_epoch_tmp.vbegin(); it_mut_tmp != mut_by_type_and_epoch_tmp.vend();){
        *it_mut += *it_mut_tmp;
        it_mut++;
        it_mut_tmp++;
      }
    }

    fclose(fp);

  }

	filenames.clear();
	for(int chr = 0; chr < chromosomes.size(); chr++){ 
		filenames.push_back(filename_base + "_chr" + chromosomes[chr] + "_opp" + ".bin");
	}  

  std::vector<CollapsedMatrix<double>> opp_by_type_and_epoch(n_boot);
  CollapsedMatrix<double> opp_by_type_and_epoch_tmp;
  fp = fopen(filenames[0].c_str(),"rb");
  for(int n = 0; n < n_boot; n++){
    opp_by_type_and_epoch[n].ReadFromFile(fp);
  }
  fclose(fp);

  for(int i = 1; i < (int) filenames.size(); i++){

    fp = fopen(filenames[i].c_str(),"rb");
    for(int n = 0; n < n_boot; n++){
      opp_by_type_and_epoch_tmp.ReadFromFile(fp);
      std::vector<double>::iterator it_opp = opp_by_type_and_epoch[n].vbegin();
      for(std::vector<double>::iterator it_opp_tmp = opp_by_type_and_epoch_tmp.vbegin(); it_opp_tmp != opp_by_type_and_epoch_tmp.vend();){
        *it_opp += *it_opp_tmp;
        it_opp++;
        it_opp_tmp++;
      }
    }

    fclose(fp);
  }

	for(int chr = 0; chr < chromosomes.size(); chr++){ 
		std::remove((filename_base + "_chr" + chromosomes[chr] + "_mut" + ".bin").c_str());
		std::remove((filename_base + "_chr" + chromosomes[chr] + "_opp" + ".bin").c_str());
	}  

  fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb"); 
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(double), epochs.size(), fp);
  for(int n = 0; n < n_boot; n++){
    mut_by_type_and_epoch[n].DumpToFile(fp);
  }
  fclose(fp);
  fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  for(int n = 0; n < n_boot; n++){
    opp_by_type_and_epoch[n].DumpToFile(fp);
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

}

void FinalizeMutationRateForCategory(cxxopts::Options& options){


  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: chr, first_chr, last_chr." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Extract mutation rate of 96 categories from .bin file." << std::endl;
    exit(0);
  } 

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing mutation rate..." << std::endl;

  int num_epochs;
  int num_categories = 0;
  std::vector<double> epoch;

  FILE* fp; 
  fp = fopen((options["input"].as<std::string>() + "_mut" + ".bin" ).c_str(), "r");  
  assert(fp != NULL);

  fread(&num_epochs, sizeof(int), 1, fp);
  epoch.resize(num_epochs);
  fread(&epoch[0], sizeof(double), num_epochs, fp);

  int n_boot = 100;
  std::vector<CollapsedMatrix<double>> boot_mutation_by_type_and_epoch(n_boot);
  std::vector<CollapsedMatrix<double>> boot_opportunity_by_type_and_epoch(n_boot);

  for(int n = 0; n < n_boot; n++){
    boot_mutation_by_type_and_epoch[n].ReadFromFile(fp);
    num_categories = boot_mutation_by_type_and_epoch[n].subVectorSize(0);
  }
  fclose(fp);
  fp = fopen((options["input"].as<std::string>() + "_opp" + ".bin" ).c_str(), "r");    
  for(int n = 0; n < n_boot; n++){
    boot_opportunity_by_type_and_epoch[n].ReadFromFile(fp);
  }
  fclose(fp);


  std::ofstream os(options["output"].as<std::string>() + ".rate");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  os << "epoch.start ";
  for(int i = 0; i < num_categories; i++){
    os << i + 1 << " ";
  }
  os << "\n";

  //output
  for(int ep = 0; ep < num_epochs-1; ep++){
    
    for(int n = 0; n < n_boot; n++){
      std::vector<double>::iterator it_row_mut = boot_mutation_by_type_and_epoch[n].rowbegin(ep);
      std::vector<double>::iterator it_row_opp = boot_opportunity_by_type_and_epoch[n].rowbegin(ep);
      os << epoch[ep] << " ";
      for(; it_row_mut != boot_mutation_by_type_and_epoch[n].rowend(ep);){
        os << *it_row_mut/(*it_row_opp) << " ";
        it_row_mut++;
        it_row_opp++;
      }
      os << "\n"; 
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


////////////////////
void MutationRateForPattern(cxxopts::Options& options, std::string chr = "NA"){

  bool help = false;
  if(!options.count("mask") || !options.count("ancestor") || !options.count("mutcat") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: mask, ancestor, mutcat, input, output. Optional: years_per_gen, bins, dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Calculate mutation rate for specific pattern (Do not apply after using RemoveTreesWithFewMutations).." << std::endl;
    exit(0);
  }  

  std::string line, line2 , read;

  //////////// PARSE DATA ///////////

  MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;

  ////////// 1. Read one tree at a time /////////

  //We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
  //The mut file is read once, file is closed after constructor is called.
  AncMutIterators ancmut;

  if(chr == "NA"){
    ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  }else{
    ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + chr + ".anc", options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  int N = (mtr.tree.nodes.size() + 1)/2.0;
  int L = ancmut.NumSnps();
  Data data(N,L);
  int N_total = 2*data.N-1;

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating mutation rate for categories " << options["input"].as<std::string>() << " ..." << std::endl;

  ////////// read mutations file ///////////

  Mutations mutations;
  if(chr == "NA"){
    mutations.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mutations.Read(options["input"].as<std::string>() + "_chr" + chr + ".mut");
  }

  std::vector<int> pos;
  if(options.count("dist")){

    int L_allsnps = 0;
    igzstream is_L;
    is_L.open(options["dist"].as<std::string>());
    if(is_L.fail()){
      std::cerr << "Error opening file " << options["dist"].as<std::string>() << std::endl;
      exit(1);
    }
    std::string unused;
    std::getline(is_L, unused); 
    while ( std::getline(is_L, unused) ){
      ++L_allsnps;
    }
    is_L.close();

    pos.resize(L_allsnps);
    igzstream is_dist(options["dist"].as<std::string>());
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
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

 
  int num_epochs;
  std::vector<double> epochs;
  float log_10 = std::log(10);
  if(options.count("bins")){

    double log_age = std::log(0);
    double age = 0; 
 
    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      ep++;
    }
    epoch_boundary = epoch_lower;
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0) epochs.push_back(age);
        epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();

  }else{

    num_epochs = 31;
    epochs.resize(num_epochs);
    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }


  /////////////////////////////////////////////////////////////////
  //Mutation specific

  ////////// define mutation types /////////
  std::string alphabet = "ACGT";
  std::map<char, std::string> complement;
  complement['A'] = "T";
  complement['C'] = "G";
  complement['G'] = "C";
  complement['T'] = "A";
  std::map<std::string, int> dict_mutation_pattern;

  //A-C, A-G, A-T, C-A, C-G, C-T
  //T-G, T-C, T-A, G-T, G-C, G-A
  //plus 16 flanks per mutation type, i.e. 6*16 = 96

  //parse file storing
  //upstream downstream ancestral derived category
  igzstream is_cat(options["mutcat"].as<std::string>());
  if(is_cat.fail()){
    std::cerr << "Error: unable to open file " << options["mutcat"].as<std::string>() << std::endl;
  }
  getline(is_cat, line);
  //I have to make sure all 96 categories are represented in this file

  char mutation_type[5];
  std::string pattern, reverse_pattern;
  int category, num_categories = 0;
  std::vector<int> check_num_categories;
  while(getline(is_cat, line)){
    sscanf(line.c_str(), "%c %c %c %c %d", &mutation_type[0], &mutation_type[1], &mutation_type[2], &mutation_type[3], &category);
    pattern = mutation_type;
    dict_mutation_pattern[pattern] = category;
    pattern = complement[mutation_type[1]] + complement[mutation_type[0]] + complement[mutation_type[2]] + complement[mutation_type[3]];
    dict_mutation_pattern[pattern] = category;
    if(category >= num_categories){
      check_num_categories.resize(category+1);
      std::fill(std::next(check_num_categories.begin(),num_categories), check_num_categories.end(), 0);
      num_categories = category + 1;      
      check_num_categories[category]++;
    }else{
      check_num_categories[category]++;
    }
  }
  is_cat.close();

  for(std::vector<int>::iterator it_check = check_num_categories.begin(); it_check != check_num_categories.end(); it_check++){
    if(*it_check == 0){
      std::cerr << "Error: category indices not 0-indexed or contiguous." << std::endl;
      exit(1);
    }
  }

  //check that I got all 96 categories and that categories are contiguous integers
  int index = 0;
  for(std::string::iterator it_str1 = alphabet.begin(); it_str1 != alphabet.end(); it_str1++){
    for(std::string::iterator it_str2 = alphabet.begin(); it_str2 != alphabet.end(); it_str2++){
      pattern = *it_str1;
      pattern += *it_str2;
      reverse_pattern = complement[*it_str2] + complement[*it_str1];

      if ( dict_mutation_pattern.find(pattern + "CA") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "GT") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "CG") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "GC") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "CT") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "AG") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AT") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TA") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AG") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TC") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
      if ( dict_mutation_pattern.find(pattern + "AC") == dict_mutation_pattern.end() && dict_mutation_pattern.find(reverse_pattern + "TG") == dict_mutation_pattern.end() ) {
        // not found
        std::cerr << "Error: not all 96 mutation categories provided." << std::endl;
        exit(1);
      }
    }   
  }

  ///////// Count number of bases by type //////////

  CollapsedMatrix<double> count_bases_by_type;
  //count_bases_by_type.resize(mutations.info.size(),dict_mutation_pattern.size());
  CountBasesByType(data, options["mask"].as<std::string>(), options["ancestor"].as<std::string>(), count_bases_by_type, dict_mutation_pattern, mutations, pos);

  ////////////////////////////////////////////////////////////////////////
  // Estimate mutation rate through time

  std::vector<CollapsedMatrix<double>> mutation_by_type_and_epoch(ancmut.NumTrees());
  std::vector<CollapsedMatrix<double>> opportunity_by_type_and_epoch(ancmut.NumTrees());
  for(std::vector<CollapsedMatrix<double>>::iterator it_m = mutation_by_type_and_epoch.begin(); it_m != mutation_by_type_and_epoch.end(); it_m++){
    (*it_m).resize(num_epochs, num_categories);
  }
  for(std::vector<CollapsedMatrix<double>>::iterator it_o = opportunity_by_type_and_epoch.begin(); it_o != opportunity_by_type_and_epoch.end(); it_o++){
    (*it_o).resize(num_epochs, num_categories);
  } 

  std::vector<double> branch_lengths_in_epoch(num_epochs);

  //Tree subtr;
  std::vector<float> coordinates_tree(N_total);
  std::vector<int> num_lineages(N_total);
  int root = N_total-1;
  int i = 0;

  //iterate over trees
  SNPInfo snp_info;
  float rec;
  int num_tree = 0;
  double total_branch_length = 0.0;
  int count_snps = 0;

  int snp = 0;
  while(num_bases_tree_persists >= 0){

    mtr.tree.GetCoordinates(coordinates_tree);
    GetCoordsAndLineages(mtr, coordinates_tree, num_lineages);
    GetBranchLengthsInEpoch(data, epochs, coordinates_tree, num_lineages, branch_lengths_in_epoch);
    num_tree = mutations.info[snp].tree;

    while(num_tree == mutations.info[snp].tree){

      snp_info = mutations.info[snp];
      if(snp_info.branch.size() == 1){

        assert(ancmut.get_treecount() == snp_info.tree);

        //if statements to make sure we have a well defined biallelic SNP
        if(snp_info.upstream_base != "NA" && snp_info.downstream_base != "NA" && snp_info.mutation_type[0] != snp_info.mutation_type[2]){

          if(snp_info.mutation_type.size() == 3){

            if(snp_info.mutation_type[0] == 'A' || snp_info.mutation_type[0] == 'C' || snp_info.mutation_type[0] == 'G' || snp_info.mutation_type[0] == 'T'){

              if(snp_info.mutation_type[2] == 'A' || snp_info.mutation_type[2] == 'C' || snp_info.mutation_type[2] == 'G' || snp_info.mutation_type[2] == 'T'){


                // identify category of mutation
                pattern = snp_info.upstream_base + snp_info.downstream_base + snp_info.mutation_type[0] + snp_info.mutation_type[2];
                // then identify its age and number of lineages at the time
                assert(dict_mutation_pattern.find(pattern) != dict_mutation_pattern.end());
                int ind = dict_mutation_pattern[pattern];

                // identify epoch and add to number of mutations per lineage in epoch
                int ep = 0;
                while(epochs[ep] <= snp_info.age_begin){
                  ep++;
                  if(ep == epochs.size()) break;
                }
                ep--;

                assert(ep >= 0);
                assert(mutation_by_type_and_epoch[num_tree].subVectorSize(ep) == num_categories);

                int ep_begin = ep;
                float age_end = std::min(snp_info.age_end,coordinates_tree[root]);
                assert(age_end < epochs[num_epochs-1]);
                double branch_length = age_end - snp_info.age_begin;

                if(age_end <= epochs[ep+1]){

                  mutation_by_type_and_epoch[num_tree][ep][ind] += 1.0;

                }else{

                  mutation_by_type_and_epoch[num_tree][ep][ind]   += (epochs[ep+1] - snp_info.age_begin)/branch_length;
                  ep++;
                  while(epochs[ep+1] <= age_end){
                    mutation_by_type_and_epoch[num_tree][ep][ind] += (epochs[ep+1]-epochs[ep])/branch_length;
                    ep++;
                  }
                  mutation_by_type_and_epoch[num_tree][ep][ind]   += (age_end-epochs[ep])/branch_length;

                }

                for(int ep_tmp = 0; ep_tmp < num_epochs; ep_tmp++){
                  double bl = branch_lengths_in_epoch[ep_tmp];
                  for(int ind_tmp = 0; ind_tmp < num_categories; ind_tmp++){
                    opportunity_by_type_and_epoch[num_tree][ep_tmp][ind_tmp] += bl * count_bases_by_type[snp][ind_tmp];
                  }
                }

              }

            }

          } 

        }

      }

      snp++;

    }

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

  }


  //bootstrap
  //sample indices from 0 too ancmut.NumTrees()-1 at random with replacement
  //dump to files and write summarise and finalise functions for bootstrap

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  //std::uniform_int_distribution<> sam(0, ancmut.NumTrees()-1);
  std::uniform_int_distribution<> sam(0, (ancmut.NumTrees()-1)/10000.0);

  int n_boot = 1000;
  std::vector<CollapsedMatrix<double>> boot_mutation_by_type_and_epoch(n_boot);
  std::vector<CollapsedMatrix<double>> boot_opportunity_by_type_and_epoch(n_boot);
  for(std::vector<CollapsedMatrix<double>>::iterator it_m = boot_mutation_by_type_and_epoch.begin(); it_m != boot_mutation_by_type_and_epoch.end(); it_m++){
    (*it_m).resize(num_epochs, num_categories);
  }
  for(std::vector<CollapsedMatrix<double>>::iterator it_o = boot_opportunity_by_type_and_epoch.begin(); it_o != boot_opportunity_by_type_and_epoch.end(); it_o++){
    (*it_o).resize(num_epochs, num_categories);
  } 

  for(int n = 0; n < n_boot; n++){
 
    //random trees 
    std::vector<int> boot_trees(ancmut.NumTrees());
    if(0){
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end(); it_boot_trees++){
      *it_boot_trees = sam(gen);
    }
    }

    int size = 0;
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end();){
      int start = 10000*sam(gen);
      for(int k = start; k < start + 10000 || k < ancmut.NumTrees(); k++){
        *it_boot_trees = k;
        it_boot_trees++;
        size++;
      }
    }
    std::cerr << size << " " << ancmut.NumTrees() << std::endl;
    boot_trees.resize(size);
    //std::sort(boot_trees.begin(), boot_trees.end());

    //fill in boot_mutation_by_type_and_epoch[n] and boot_opportunity_by_type_and_epoch[n] by summing over the trees  
    for(std::vector<int>::iterator it_boot_trees = boot_trees.begin(); it_boot_trees != boot_trees.end(); it_boot_trees++){
      std::vector<double>::iterator it_bmut_by_type_and_epoch = boot_mutation_by_type_and_epoch[n].vbegin();
      std::vector<double>::iterator it_bopp_by_type_and_epoch = boot_opportunity_by_type_and_epoch[n].vbegin();
      std::vector<double>::iterator it_mut_by_type_and_epoch  = mutation_by_type_and_epoch[*it_boot_trees].vbegin();
      std::vector<double>::iterator it_opp_by_type_and_epoch  = opportunity_by_type_and_epoch[*it_boot_trees].vbegin();
      for(; it_mut_by_type_and_epoch != mutation_by_type_and_epoch[*it_boot_trees].vend();){    
        *it_bmut_by_type_and_epoch += *it_mut_by_type_and_epoch;
        *it_bopp_by_type_and_epoch += *it_opp_by_type_and_epoch;
        it_mut_by_type_and_epoch++;
        it_opp_by_type_and_epoch++;      
        it_bmut_by_type_and_epoch++;
        it_bopp_by_type_and_epoch++;
      }
    }
  
  }

  //output mutation_by_type_and_epoch
  //       opportunity_by_type_and_epoch 
  FILE* fp;
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_mut" + ".bin" ).c_str(), "wb");  
  }
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(float), epochs.size(), fp);
  for(int n = 0; n < n_boot; n++){
    boot_mutation_by_type_and_epoch[n].DumpToFile(fp);
  }
  fclose(fp);
  if(chr == "NA"){
    fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + chr + "_opp" + ".bin" ).c_str(), "wb");  
  }
  for(int n = 0; n < n_boot; n++){
    boot_opportunity_by_type_and_epoch[n].DumpToFile(fp);
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


}

void SummarizeWholeGenomeForPattern(cxxopts::Options& options){

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
    std::cout << "Reads .bin files and summarizes them into one .bin file." << std::endl;
    exit(0);
  }  

	std::vector<std::string> filenames;
	std::vector<std::string> chromosomes;
	std::string filename_base = options["output"].as<std::string>();

	if(options.count("chr")){
		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		std::cerr << "------------------------------------------------------" << std::endl;
		std::cerr << "Summarizing mutation rates for chr in " << options["chr"].as<std::string>() << "..." << std::endl;
		std::string line;
		while(getline(is_chr, line)){
		  filenames.push_back(filename_base + "_chr" + line + "_mut" + ".bin");
			chromosomes.push_back(line);
		}
		is_chr.close();
	}else{
		int start     = options["first_chr"].as<int>(); 
		int end       = options["last_chr"].as<int>();
		std::cerr << "------------------------------------------------------" << std::endl;
		std::cerr << "Summarizing mutation rates for chr " << start << " " << end << "..." << std::endl;
		for(int chr = start; chr <= end; chr++){
			filenames.push_back(filename_base + "_chr" + std::to_string(chr) + "_mut" + ".bin");
			chromosomes.push_back(std::to_string(chr));
		}
	}

  //open files and add together

  FILE* fp = fopen(filenames[0].c_str(),"rb");
  assert(fp != NULL);  

  int num_epochs;
  std::vector<double> epochs;
  fread(&num_epochs, sizeof(int), 1, fp);
  epochs.resize(num_epochs);
  fread(&epochs[0], sizeof(float), num_epochs, fp);

  int n_boot = 1000;
  std::vector<CollapsedMatrix<double>> mut_by_type_and_epoch(n_boot);
  CollapsedMatrix<double> mut_by_type_and_epoch_tmp;
  for(int n = 0; n < n_boot; n++){
    mut_by_type_and_epoch[n].ReadFromFile(fp);
  }
  fclose(fp);
  for(int i = 1; i < (int) filenames.size(); i++){
    fp = fopen(filenames[i].c_str(),"rb");

    fread(&num_epochs, sizeof(int), 1, fp);
    fread(&epochs[0], sizeof(float), num_epochs, fp);

    for(int n = 0; n < n_boot; n++){
      mut_by_type_and_epoch_tmp.ReadFromFile(fp);
      std::vector<double>::iterator it_mut = mut_by_type_and_epoch[n].vbegin();
      for(std::vector<double>::iterator it_mut_tmp = mut_by_type_and_epoch_tmp.vbegin(); it_mut_tmp != mut_by_type_and_epoch_tmp.vend();){
        *it_mut += *it_mut_tmp;
        it_mut++;
        it_mut_tmp++;
      }
    }

    fclose(fp);

  }

  filenames.clear();
  for(int chr = 0; chr < chromosomes.size(); chr++){
    filenames.push_back(filename_base + "_chr" + chromosomes[chr] + "_opp" + ".bin");
  }

  std::vector<CollapsedMatrix<double>> opp_by_type_and_epoch(n_boot);
  CollapsedMatrix<double> opp_by_type_and_epoch_tmp;
  fp = fopen(filenames[0].c_str(),"rb");
  for(int n = 0; n < n_boot; n++){
    opp_by_type_and_epoch[n].ReadFromFile(fp);
  }
  fclose(fp);

  for(int i = 1; i < (int) filenames.size(); i++){

    fp = fopen(filenames[i].c_str(),"rb");
    for(int n = 0; n < n_boot; n++){
      opp_by_type_and_epoch_tmp.ReadFromFile(fp);
      std::vector<double>::iterator it_opp = opp_by_type_and_epoch[n].vbegin();
      for(std::vector<double>::iterator it_opp_tmp = opp_by_type_and_epoch_tmp.vbegin(); it_opp_tmp != opp_by_type_and_epoch_tmp.vend();){
        *it_opp += *it_opp_tmp;
        it_opp++;
        it_opp_tmp++;
      }
    }

    fclose(fp);
  }

  for(int chr = 0; chr < chromosomes.size(); chr++){ 
    std::remove((options["input"].as<std::string>() + "_chr" + chromosomes[chr] + "_mut" + ".bin").c_str());
    std::remove((options["input"].as<std::string>() + "_chr" + chromosomes[chr] + "_opp" + ".bin").c_str());
  }  

  fp = fopen((options["output"].as<std::string>() + "_mut" + ".bin" ).c_str(), "wb"); 
  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epochs[0], sizeof(float), epochs.size(), fp);
  for(int n = 0; n < n_boot; n++){
    mut_by_type_and_epoch[n].DumpToFile(fp);
  }
  fclose(fp);
  fp = fopen((options["output"].as<std::string>() + "_opp" + ".bin" ).c_str(), "wb");  
  for(int n = 0; n < n_boot; n++){
    opp_by_type_and_epoch[n].DumpToFile(fp);
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

}

void FinalizeMutationRateForPattern(cxxopts::Options& options){


  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: chr, first_chr, last_chr." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Extract mutation rate of 96 categories from .bin file." << std::endl;
    exit(0);
  } 

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Finalizing mutation rate..." << std::endl;

  int num_epochs;
  int num_categories = 0;
  std::vector<double> epoch;


  int n_boot = 1000;
  std::vector<CollapsedMatrix<double>> boot_mutation_by_type_and_epoch(n_boot);
  std::vector<CollapsedMatrix<double>> boot_opportunity_by_type_and_epoch(n_boot);
  for(std::vector<CollapsedMatrix<double>>::iterator it_m = boot_mutation_by_type_and_epoch.begin(); it_m != boot_mutation_by_type_and_epoch.end(); it_m++){
    (*it_m).resize(num_epochs, num_categories);
  }
  for(std::vector<CollapsedMatrix<double>>::iterator it_o = boot_opportunity_by_type_and_epoch.begin(); it_o != boot_opportunity_by_type_and_epoch.end(); it_o++){
    (*it_o).resize(num_epochs, num_categories);
  } 
   
  FILE* fp; 
  fp = fopen((options["input"].as<std::string>() + "_mut" + ".bin" ).c_str(), "r");  
  assert(fp != NULL);

  fread(&num_epochs, sizeof(int), 1, fp);
  epoch.resize(num_epochs);
  fread(&epoch[0], sizeof(float), num_epochs, fp);

  for(int n = 0; n < n_boot; n++){
    boot_mutation_by_type_and_epoch[n].ReadFromFile(fp);
    num_categories = boot_mutation_by_type_and_epoch[n].subVectorSize(0);
  }
  fclose(fp);
  fp = fopen((options["input"].as<std::string>() + "_opp" + ".bin" ).c_str(), "r");    
  for(int n = 0; n < n_boot; n++){
    boot_opportunity_by_type_and_epoch[n].ReadFromFile(fp);
  }
  fclose(fp);


  std::ofstream os(options["output"].as<std::string>() + ".rate");
  if(os.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }

  os << "epoch.start ";
  for(int i = 0; i < num_categories; i++){
    os << i + 1 << " ";
  }
  os << "\n";

  //output
  for(int ep = 0; ep < num_epochs-1; ep++){
    
    for(int n = 0; n < n_boot; n++){
      std::vector<double>::iterator it_row_mut = boot_mutation_by_type_and_epoch[n].rowbegin(ep);
      std::vector<double>::iterator it_row_opp = boot_opportunity_by_type_and_epoch[n].rowbegin(ep);
      os << epoch[ep] << " ";
      for(; it_row_mut != boot_mutation_by_type_and_epoch[n].rowend(ep);){
        os << *it_row_mut/(*it_row_opp) << " ";
        it_row_mut++;
        it_row_opp++;
      }
      os << "\n"; 
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



//////////////// Other applications ///////////////

void BranchLengthVsMutations(cxxopts::Options& options){

  std::cerr << "------------------------------------------------------" << std::endl;
  std::cerr << "Calculating number of mutations vs opportunity..." << std::endl;

  bool help = false;
  if(!options.count("pos") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: pos, input, output. Optional: years_per_gen, mut." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "..." << std::endl;
    exit(0);
  }  

  std::string line, read;

  //parse data
  AncMutIterators ancmut;
  ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  
  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0; 

  Data data(N,L);
  //data.ReadPosition(options["pos"].as<std::string>());
  int N_total = 2*data.N-1;

  //epoch
  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }
 
  int num_epochs;
  std::vector<double> epochs;
  float log_10 = std::log(10);
  if(options.count("bins")){

    double log_age = std::log(0);
    double age = 0; 
 
    double epoch_lower, epoch_upper, epoch_step;
    std::string str_epochs = options["bins"].as<std::string>();
    std::string tmp;
    int i = 0;
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_lower = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_upper = std::stof(tmp);
    i++;
    if(i >= str_epochs.size()){
      std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
      exit(1);
    }
    tmp = "";
    while(str_epochs[i] != ','){
      tmp += str_epochs[i];
      i++;
      if(i == str_epochs.size()) break;
    }
    epoch_step = std::stof(tmp);

    int ep = 0;
    epochs.resize(1);
    epochs[ep] = 0.0;
    ep++; 
    double epoch_boundary = 0.0;
    if(log_age < epoch_lower && age != 0.0){
      epochs.push_back(age);
      ep++;
    }
    epoch_boundary = epoch_lower;
    while(epoch_boundary < epoch_upper){
      if(log_age < epoch_boundary){
        if(ep == 1 && age != 0.0) epochs.push_back(age);
        epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
        ep++;
      }
      epoch_boundary += epoch_step;
    }
    epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
    epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();

  }else{

    num_epochs = 31;
    epochs.resize(num_epochs);
    epochs[0] = 0.0;
    epochs[1] = 1e3/years_per_gen;
    for(int e = 2; e < num_epochs-1; e++){
      epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    }
    epochs[num_epochs-1] = 1e8/years_per_gen;

  }


  //read mutations file
  Mutations mutations(data);
  mutations.Read(options["input"].as<std::string>() + ".mut");

  ///////////////////
  // Count number of bases by type
  double total_num_bases = (*std::prev(mutations.info.end(),1)).pos - (*mutations.info.begin()).pos;
  std::vector<double> count_bases(mutations.info.size(), 0.0);
  std::vector<double>::iterator it_count_bases   = count_bases.begin();

  it_mut = mutations.info.begin();
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

  std::vector<float> coordinates(N_total);
  int root = N_total-1;
  int i = 0;
  int snp_of_next_tree, snp;


  //read tree
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
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
    while(epochs[ep] < coordinates[i]) ep++;
    if(epochs[ep] <= coordinates[parent]){
      assert(epochs[ep] >= coordinates[i]);
      num_mutations_in_epoch[ep-1]    += num_events * (epochs[ep] - coordinates[i])/bl;
      branch_lengths_in_epoch[ep-1]   += delta_pos * (epochs[ep] - coordinates[i]);
      ep++;
      while(epochs[ep] < coordinates[parent]){
        num_mutations_in_epoch[ep-1]  += num_events * (epochs[ep] - epochs[ep-1])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (epochs[ep] - epochs[ep-1]);
        ep++;
      }
      assert(coordinates[parent] >= epochs[ep-1]);
      num_mutations_in_epoch[ep-1]  += num_events * (coordinates[parent] - epochs[ep-1])/bl;
      branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - epochs[ep-1]);
    }else{
      num_mutations_in_epoch[ep-1]  += num_events * (coordinates[parent] - coordinates[i])/bl;
      branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - coordinates[i]);
    }
  }

  for(int ep = 0; ep < epochs.size() - 1; ep++){
    os << mtr.pos << " " << (int) years_per_gen * (epochs[ep] + epochs[ep+1])/2.0 << " " << data.mu * branch_lengths_in_epoch[ep] << " " << num_mutations_in_epoch[ep] << "\n"; 
  }

  while(num_bases_tree_persists >= 0.0){

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
      while(epochs[ep] < coordinates[i]) ep++;
      assert(epochs[ep] >= coordinates[i]);

      if(epochs[ep] <= coordinates[parent]){
        num_mutations_in_epoch[ep-1] += num_events * (epochs[ep] - coordinates[i])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (epochs[ep] - coordinates[i]);
        ep++;
        while(epochs[ep] < coordinates[parent]){
          num_mutations_in_epoch[ep-1] += num_events * (epochs[ep] - epochs[ep-1])/bl;
          branch_lengths_in_epoch[ep-1] += delta_pos * (epochs[ep] - epochs[ep-1]);
          ep++;
        }
        assert(coordinates[parent] >= epochs[ep-1]);
        num_mutations_in_epoch[ep-1] += num_events * (coordinates[parent] - epochs[ep-1])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - epochs[ep-1]);
      }else{
        num_mutations_in_epoch[ep-1] += num_events * (coordinates[parent] - coordinates[i])/bl;
        branch_lengths_in_epoch[ep-1] += delta_pos * (coordinates[parent] - coordinates[i]);
      }
    }

    for(int ep = 0; ep < epochs.size()-1; ep++){
      os << mtr.pos << " " << (int) years_per_gen * (epochs[ep] + epochs[ep+1])/2.0 << " " << data.mu * branch_lengths_in_epoch[ep] << " " << num_mutations_in_epoch[ep] << "\n"; 
    }

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    mtr.tree.GetCoordinates(coordinates);

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
    std::cout << "Needed: input, output. Optional: chr, first_chr, last_chr." << std::endl;
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
  std::vector<double> epoch;

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
		("chr", "Optional: File specifying chromosomes to use. Overrides first_chr, last_chr.", cxxopts::value<std::string>()) 
    ("first_chr", "Index of fist chr", cxxopts::value<int>())
    ("last_chr", "Index of last chr", cxxopts::value<int>())
    ("years_per_gen", "Years per generation (float). Default: 28.", cxxopts::value<float>())
    ("bins", "Specify epoch bins. Format: lower, upper, stepsize for function c(0,10^seq(lower, upper, stepsize)).", cxxopts::value<std::string>())
    ("binsfile", "Filename containing bins. One per line in generations.", cxxopts::value<std::string>())
    ("sample_age", "Sample age in generations.", cxxopts::value<float>())
    ("dist", "Filename of file containing dist.", cxxopts::value<std::string>())
    ("mask", "Filename of file containing mask", cxxopts::value<std::string>())
    ("ancestor", "Filename of file containing human ancestor genome.", cxxopts::value<std::string>())
    ("mutcat", "Filename of file containing mutation categories.", cxxopts::value<std::string>())
    ("poplabels", "Optional: Filename of file containing population labels. If ='hap', each haplotype is in its own group.", cxxopts::value<std::string>()) 
    ("pop_of_interest", "Optional: Name of pop of interest.", cxxopts::value<std::string>()) 
    ("i,input", "Filename of .anc and .mut file without file extension", cxxopts::value<std::string>())
    ("o,output", "Output file", cxxopts::value<std::string>())
    ("seed", "Optional. Random seed.", cxxopts::value<int>());

  options.parse(argc, argv);
  std::string mode = options["mode"].as<std::string>();

  if(!mode.compare("WithContext")){

    //////////////////////////////////
    //Program options
    bool help = false;
    if( !options.count("input") || !options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: mask, ancestor, input, output. Optional: chr, first_chr, last_chr, years_per_gen, bins, dist." << std::endl;
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

		if(options.count("chr")){
			igzstream is_chr(options["chr"].as<std::string>());
			if(is_chr.fail()){
				std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
			}
			std::string line;
			while(getline(is_chr, line)){
				MutationRateWithContext(options, line);
			}
			is_chr.close();
			SummarizeWholeGenome(options);      
			FinalizeMutationRate(options);
		}else if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      for(int chr = options["first_chr"].as<int>(); chr <= options["last_chr"].as<int>(); chr++){ 
        MutationRateWithContext(options, std::to_string(chr));
      }
      SummarizeWholeGenome(options);      
      FinalizeMutationRate(options);
    }else{
      MutationRateWithContext(options);
      FinalizeMutationRate(options);
    }

  }else if(!mode.compare("MutationRateForCategory")){

    //////////////////////////////////
    //Program options
    bool help = false;
    if( !options.count("input") || !options.count("output")){
      std::cout << "Not enough arguments supplied." << std::endl;
      std::cout << "Needed: mask, ancestor, mutcat, input, output. Optional: chr, first_chr, last_chr, years_per_gen, bins, dist." << std::endl;
      help = true;
    }
    if(options.count("help") || help){
      std::cout << options.help({""}) << std::endl;
      std::cout << "Calculate mutation rate for categories specified in input.mutcat (Do not apply after using RemoveTreesWithFewMutations)." << std::endl;
      exit(0);
    } 
    if(options["input"].as<std::string>() != options["output"].as<std::string>()){
      std::cerr << "Sorry, in this mode input and output need to be the same!." << std::endl;
      exit(1);
    }

		if(options.count("chr")){
			igzstream is_chr(options["chr"].as<std::string>());
			if(is_chr.fail()){
				std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
			}
			std::string line;
			while(getline(is_chr, line)){
				MutationRateForCategory(options, line);
			}
			is_chr.close();
			SummarizeWholeGenomeForCategory(options);      
			FinalizeMutationRateForCategory(options);
		}else if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      for(int chr = options["first_chr"].as<int>(); chr <= options["last_chr"].as<int>(); chr++){ 
        MutationRateForCategory(options, std::to_string(chr));
      }
      SummarizeWholeGenomeForCategory(options);      
      FinalizeMutationRateForCategory(options);
    }else{
      MutationRateForCategory(options);
      FinalizeMutationRateForCategory(options);
    }


  }else if(!mode.compare("ForCategoryForChromosome")){

    MutationRateForCategory(options);

  }else if(!mode.compare("ForCategoryForPopForChromosome")){

    MutationRateForCategoryForGroup(options);

  }else if(!mode.compare("WithContextForChromosome")){

    MutationRateWithContext(options);

  }else if(!mode.compare("SummarizeForGenome")){

    SummarizeWholeGenome(options);

  }else if(!mode.compare("SummarizeForGenomeForCategory")){

    SummarizeWholeGenomeForCategory(options);

  }else if(!mode.compare("Finalize")){

		if(options.count("chr")){
			SummarizeWholeGenome(options);
		}else if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      SummarizeWholeGenome(options);
    }
    FinalizeMutationRate(options);

  }else if(!mode.compare("FinalizeForCategory")){

		if(options.count("chr")){
			SummarizeWholeGenomeForCategory(options);
		}else if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      SummarizeWholeGenomeForCategory(options);
    }
    FinalizeMutationRateForCategory(options);

  }else if(!mode.compare("FinalizeMutationCount")){

		if(options.count("chr")){
			SummarizeWholeGenome(options);
		}else if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      SummarizeWholeGenome(options);
    }

    FinalizeMutationCount(options);

  }else if(!mode.compare("FinalizeAvg")){

		if(options.count("chr")){
			SummarizeWholeGenome(options);
		}else if(options.count("first_chr") && options.count("last_chr")){
      if(options["first_chr"].as<int>() < 0 || options["last_chr"].as<int>() < 0){
        std::cerr << "Do not use negative chr indices." << std::endl;
        exit(1);
      }
      SummarizeWholeGenome(options);
    }

    FinalizeAvg(options);

  }else if(!mode.compare("Avg")){

    AvgMutationRate(options);

  }else if(!mode.compare("XY")){

    BranchLengthVsMutations(options);

  }else{

    std::cout << "####### error #######" << std::endl;
    std::cout << "Invalid or missing mode." << std::endl;
    std::cout << "Options for --mode are:" << std::endl;
    std::cout << "Avg, MutationRateForCategory, ForCategoryForChromosome, WithContext, WithContextForChromosome, SummarizeForGenome, SummarizeForGenomeForCategory, Finalize, FinalizeForCategory, FinalizeMutationCount, XY." << std::endl;

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

