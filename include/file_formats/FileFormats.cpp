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
#include "data.hpp"
#include "sample.hpp"
#include "mutations.hpp"
#include "cxxopts.hpp"

void
ConvertFromHapLegendSample(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("haps") || !options.count("sample")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, haps, sample. Optional: chr." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Converts hap/legend/sample file format (Impute2) to haps/sample file format (Shapeit).." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Converting hap/legend/sample file format to haps/sample file format.." << std::endl;

  int chr = 0;
  if(options.count("chr")){
    chr = options["chr"].as<int>();
  }
  std::string line_hap, line_legend;

  //parse hap and legend
  //combine hap and legend to haps

  //legend
  //assume 4 columns with optional 5th column

  igzstream is_hap(options["input"].as<std::string>() + ".hap");
  if(is_hap.fail()) is_hap.open(options["input"].as<std::string>() + ".hap.gz");
  igzstream is_legend(options["input"].as<std::string>() + ".legend");
  if(is_legend.fail()) is_legend.open(options["input"].as<std::string>() + ".legend.gz");
  FILE* fp_haps   = fopen((options["haps"].as<std::string>()).c_str(), "w");

  if(is_hap.fail()){
    std::cerr << "Error while opening file " << options["input"].as<std::string>() + ".hap(.gz)" << "." << std::endl;
    exit(1);
  } 
  if(is_legend.fail()){
    std::cerr << "Error while opening file " << options["input"].as<std::string>() + ".legend(.gz)" << "." << std::endl;
    exit(1);
  } 

  int bp1, bp2, bp3;
  char rsid1[1024], rsid2[1024], rsid3[1024];
  char type1[1024], type2[1024], type3[1024];
  char ancestral1[1024], alternative1[1024], ancestral2[1024], alternative2[1024], ancestral3[1024], alternative3[1024];
  int matches1, matches2, matches3;
  int snp = 1, snp_accepted = 1;

  getline(is_legend, line_legend); //skip header
  assert(getline(is_legend, line_legend));
  matches1 = sscanf(line_legend.c_str(), "%s %d %s %s %s", rsid1, &bp1, ancestral1, alternative1, type1);
  assert(getline(is_legend, line_legend));
  matches2 = sscanf(line_legend.c_str(), "%s %d %s %s %s", rsid2, &bp2, ancestral2, alternative2, type2);

  while(getline(is_legend, line_legend)){

    assert(getline(is_hap, line_hap));
    matches3 = sscanf(line_legend.c_str(), "%s %d %s %s %s", rsid3, &bp3, ancestral3, alternative3, type3);

    if(snp == 1 && bp2 > bp1){ //if first SNP is unique

      if(matches1 == 4){
        fprintf(fp_haps, "%d %s %d %s %s", chr, rsid1, bp1, ancestral1, alternative1);
        fprintf(fp_haps, " %s\n", line_hap.c_str());
        snp_accepted++;
      }else if(matches1 == 5){

        if(strcmp(type1, "Biallelic_SNP") == 0){ 
          fprintf(fp_haps, "%d %s %d %s %s", chr, rsid1, bp1, ancestral1, alternative1);
          fprintf(fp_haps, " %s\n", line_hap.c_str());
          snp_accepted++;
        }

      }else{
        std::cerr << "An error occurred while reading line " << snp << std::endl;
        exit(1);
      }

      assert(getline(is_hap, line_hap));
      snp++;

    }

    if(bp3 > bp2 && bp2 > bp1){ //if SNP is unique

      if(matches2 == 4){
        fprintf(fp_haps, "%d %s %d %s %s", chr, rsid2, bp2, ancestral2, alternative2);
        fprintf(fp_haps, " %s\n", line_hap.c_str());
        snp_accepted++;
      }else if(matches2 == 5){

        if(strcmp(type2, "Biallelic_SNP") == 0){ 
          fprintf(fp_haps, "%d %s %d %s %s", chr, rsid2, bp2, ancestral2, alternative2);
          fprintf(fp_haps, " %s\n", line_hap.c_str());
          snp_accepted++;
        } 

      }else{
        std::cerr << "An error occurred while reading line " << snp << std::endl;
        exit(1);
      }

    }else if(bp2 < bp1){
      std::cerr << "Error: snp are not sorted by bp." << std::endl;
      exit(1);
    }

    bp1          = bp2;
    strcpy(ancestral1, ancestral2);
    strcpy(alternative1, alternative2);
    matches1     = matches2;
    strcpy(type1, type2);
    strcpy(rsid1, rsid2);

    bp2          = bp3;
    strcpy(ancestral2, ancestral3);
    strcpy(alternative2, alternative3);
    matches2     = matches3;
    strcpy(type2, type3);
    strcpy(rsid2, rsid3);

    snp++;

  }

  //last SNP
  if(bp2 > bp1){
    if(matches2 == 4){
      fprintf(fp_haps, "%d %s %d %s %s", chr, rsid2, bp2, ancestral2, alternative2);
      fprintf(fp_haps, " %s\n", line_hap.c_str());
      snp_accepted++;
    }else if(matches2 == 5){

      if(strcmp(type2, "Biallelic_SNP") == 0){ 
        fprintf(fp_haps, "%d %s %d %s %s", chr, rsid2, bp2, ancestral2, alternative2);
        fprintf(fp_haps, " %s\n", line_hap.c_str());
        snp_accepted++;
      } 

    }else{
      std::cerr << "An error occurred while reading line " << snp << std::endl;
      exit(1);
    } 
  }
  snp++;


  fclose(fp_haps);
  is_hap.close();
  is_legend.close();

  //create sample file (not dependent on input)
  igzstream is_sample(options["input"].as<std::string>() + ".sample");
  if(is_sample.fail()) is_sample.open(options["input"].as<std::string>() + ".sample.gz");
  if(is_sample.fail()){
    std::cerr << "Error while opening file " << options["input"].as<std::string>() + ".sample(.gz)" << "." << std::endl;
    exit(1);
  } 
  std::ofstream os_sample(options["sample"].as<std::string>());

  os_sample << "ID_1\tID_2\tmissing\n";
  os_sample << "0\t0\t0\n";

  std::string line;
  char id[1024];
  getline(is_sample,line);
  while(getline(is_sample, line)){
    sscanf(line.c_str(), "%s", id);
    os_sample << id << "\t" << id << "\t0\n";
  }
  is_sample.close();
  os_sample.close();

  std::cerr << "Removed " << snp - snp_accepted << " non-biallelic SNPs." << std::endl;
  std::cerr << "Output written to " << options["haps"].as<std::string>() << " and " << options["sample"].as<std::string>() << "." << std::endl;

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
ConvertFromVcf(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("haps") || !options.count("sample")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, haps, sample. Optional: flag (0: all variants, 1 only SNP)" << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Converts vcf file format to haps/sample file format (Shapeit).." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Converting vcf file format to haps/sample file format.." << std::endl;


  std::string line;

  //parse vcf 

  igzstream is_vcf(options["input"].as<std::string>() + ".vcf");
  if(is_vcf.fail()) is_vcf.open(options["input"].as<std::string>() + ".vcf.gz");
  if(is_vcf.fail()){
    std::cerr << "Error opening file " << options["input"].as<std::string>() + ".vcf(.gz)" << std::endl;
    exit(1);
  }
  FILE* fp_haps = fopen((options["haps"].as<std::string>()).c_str(), "w");
  assert(fp_haps);

  char chr[1024];
  int bp;
  char rsid[1024];
  char type[1024];
  char ancestral[1024], alternative[1024];

  std::string line_id;
  getline(is_vcf, line);
  while(line[0] == '#'){
    line_id = line;
    getline(is_vcf,line);
  }


  //create sample file (not dependent on input)
  std::ofstream os_sample(options["sample"].as<std::string>());

  int c = 0;
  for(int k = 0; k < 9; k++){
    while(line_id[c] != '\t' && line_id[c] != ' '){
      c++;
    }
    c++;
  }

  os_sample << "ID_1\tID_2\tmissing\n";
  os_sample << "0\t0\t0\n";

  char id[1024];
  int N = 0;
  while(c < line_id.size()){
    N++;
    sscanf(&(line_id.c_str())[c], "%s", id);
    os_sample << id << "\t" << id << "\t0\n";
    while(line_id[c] != '\t' && line_id[c] != ' ' && c < line_id.size()){
      c++;
    }
    c++;
  }
  os_sample.close();

  //create haps file
  int N_prev = N;

  bool only_snps = false;
  if(options.count("flag")){
    if(options["flag"].as<int>() == 1) only_snps = true;
  }

  std::vector<char> sequence;
  do{

    sscanf(line.c_str(), "%s %d %s %s %s", chr, &bp, rsid, ancestral, alternative);

    if(only_snps && strlen(ancestral) == 1 && strlen(alternative) == 1){

      int c = 0;
      for(int k = 0; k < 9; k++){
        while(line[c] != '\t' && line[c] != ' '){
          c++;
        }
        c++;
      }

      std::vector<char> seq(2*N_prev);
      N = 0;
      int freq = 0;
      while(c < line.size()-2){
        if(line[c] == '0' && line[c+1] == '|' && line[c+2] == '0'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          c += 2;
        }else if(line[c] == '0' && line[c+1] == '|' && line[c+2] == '1'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          freq++;
          c += 2;
        }else if(line[c] == '1' && line[c+1] == '|' && line[c+2] == '0'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          freq++;
          c += 2;
        }else if(line[c] == '1' && line[c+1] == '|' && line[c+2] == '1'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          freq += 2;
          c += 2;
        }else if(line[c] == '0' && line[c+1] == '/' && line[c+2] == '0'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          c += 2;
        }else if(line[c] == '0' && line[c+1] == '/' && line[c+2] == '1'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          freq++;
          c += 2;
        }else if(line[c] == '1' && line[c+1] == '/' && line[c+2] == '0'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          freq++;
          c += 2;
        }else if(line[c] == '1' && line[c+1] == '/' && line[c+2] == '1'){
          if(N >= N_prev) break;
          seq[2*N] = line[c];
          seq[2*N+1] = line[c+2];
          N++;
          freq += 2;
          c += 2;
        }
        if(c < line.size()){
          while(line[c] != ' ' && line[c] != '\t' && line[c] != '\n' && c < line.size()-1){
            c++;
          }
          c++;
        }
      }

      if(N == N_prev){

        for(int i = 0; i < strlen(rsid); i++){
          if(rsid[i] == ';') rsid[i] = ',';
        }
        fprintf(fp_haps, "%s %s %d %s %s", chr, rsid, bp, ancestral, alternative);
        for(std::vector<char>::iterator it_seq = seq.begin(); it_seq != seq.end(); it_seq++){
          fprintf(fp_haps, " %c", *it_seq);
        }
        fprintf(fp_haps, "\n");
      }

    }

  }while(getline(is_vcf, line));
  fclose(fp_haps);
  is_vcf.close();

  std::cerr << "Output written to " << options["haps"].as<std::string>() << " and " << options["sample"].as<std::string>() << "." << std::endl;

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
RemoveNonBiallelicSNPs(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("haps") | !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Removes SNPs at positions with multiple SNPs." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Removing non-biallelic SNPs.." << std::endl;

  igzstream is_haps(options["haps"].as<std::string>());
  std::ofstream os_haps(options["output"].as<std::string>() + ".haps");

  if(is_haps.fail()){
    std::cerr << "Error while opening file " << options["haps"].as<std::string>() << "." << std::endl;
    exit(1);
  } 

  char chr[1024];
  char rsid[1024];
  int bp1, bp2, bp3;
  std::string line1, line2, line3;

  assert(getline(is_haps,line1));
  sscanf(line1.c_str(), "%s %s %d", chr, rsid, &bp1);
  assert(getline(is_haps,line2));
  sscanf(line2.c_str(), "%s %s %d", chr, rsid, &bp2);

  int snp = 1, snp_accepted = 1;
  while(getline(is_haps, line3)){

    sscanf(line3.c_str(), "%s %s %d", chr, rsid, &bp3);
    if(snp == 1 && bp2 > bp1){
      os_haps << line1 << "\n";
      snp_accepted++;
      snp++;
    }

    if(bp3 > bp2 && bp2 > bp1){
      os_haps << line2 << "\n";
      snp_accepted++;
    }

    if(bp2 < bp1){
      std::cerr << "An error occurred while reading line " << snp << ". Input file might not be sorted by bp." << std::endl;
      exit(1);
    }

    bp1 = bp2;
    bp2 = bp3;
    line1 = line2;
    line2 = line3;

    snp++;
  }
  if(bp2 > bp1){
    os_haps << line2 << "\n";
    snp_accepted++;
  }
  snp++;
  is_haps.close();
  os_haps.close();

  std::cerr << "Removed " << snp - snp_accepted << " non-biallelic SNPs." << std::endl;
  std::cerr << "Output written to " << options["output"].as<std::string>() + ".haps" << "." << std::endl;


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
RemoveSamples(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("haps") || !options.count("sample") || !options.count("input") || !options.count("output") ){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, sample, input, output. Optional: poplabels, flag (0: remove fixed mutations - default, 1: output all mutations)." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Remove samples specified in input.." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Removing samples specified in input.. " << std::endl;

  bool remove_fixed = true;
  if(options.count("flag")){
    if((options["flag"].as<std::string>() != "0") && (options["flag"].as<std::string>() != "1")){
      std::cerr << "Error: flag does not exist." << std::endl;
      std::cerr << "Use --flag 0: remove fixed mutations - default, --flag 1: output all mutations." << std::endl;
      exit(1);
    }
    remove_fixed = (options["flag"].as<std::string>() == "0");
  }

  haps m_hap(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  FILE* fp = fopen((options["output"].as<std::string>() + ".haps").c_str(), "w");
  Data data(m_hap.GetN(), m_hap.GetL());

  std::cerr << m_hap.GetN() << " " << m_hap.GetL() << std::endl;

  // first read input file containing ids of individuals to remove
  char id[1024], id2[1024];
  std::string line, line2;
  std::vector<std::string> id_remove;
  igzstream is_rem(options["input"].as<std::string>());
  if(is_rem.fail()){
    std::cerr << "Error while opening file " << options["input"].as<std::string>() << "." << std::endl;
    exit(1);
  } 
  while(getline(is_rem, line)){
    id_remove.push_back(line);
  }
  is_rem.close();

  igzstream is(options["sample"].as<std::string>());
  std::ofstream os(options["output"].as<std::string>() + ".sample");
  std::vector<int> remaining_haps(data.N - 2*id_remove.size());

  if(is.fail()){
    std::cerr << "Error while opening file " << options["sample"].as<std::string>() << "." << std::endl;
    exit(1);
  }

  igzstream is_pop;
  std::ofstream os_pop;
  bool poplabels = false;
  if(options.count("poplabels")){
    poplabels = true;
    is_pop.open(options["poplabels"].as<std::string>());
    if(is_pop.fail()){
      std::cerr << "Error while opening file " << options["poplabels"].as<std::string>() << "." << std::endl;
    }
    os_pop.open(options["output"].as<std::string>() + ".poplabels");
    getline(is_pop, line2);
    os_pop << line2 << "\n";
  } 

  os << "ID_1\tID_2\tmissing\n";
  os << "0\t0\t0\n";

  getline(is, line);
  getline(is, line); 
  // then update .sample file. store indices of remaining haplotypes in remaining_haps
  int i = 0, j = 0;
  while(getline(is, line)){

    if(poplabels){
      if(!getline(is_pop, line2)){
        std::cerr << "Error: poplabels file has fewer samples than the .sample file.\n";
        exit(1);
      }
    }

    sscanf(line.c_str(), "%s %s", id, id2);

    bool remove = false;
    for(std::vector<std::string>::iterator it_rem = id_remove.begin(); it_rem != id_remove.end(); it_rem++){
      if(strcmp((*it_rem).c_str(), id) == 0){
        remove = true;
        break;
      }
    }

    if(remove == false){
      os << line << "\n";
      if(poplabels) os_pop << line2 << "\n";
      remaining_haps[i] = j;
      i++;
      j++;
      if(strcmp(id, id2) == 0){ //diploid
        remaining_haps[i] = j;
        i++;
        j++;
      }
    }else{
      j++;
      if(strcmp(id, id2) == 0) j++;
    }

  }
  assert(i == remaining_haps.size());
  if(poplabels){
    if(getline(is_pop,line2)){
      std::cerr << "Error: poplabels file has more samples than the .sample file.\n";
      exit(1);
    }
    is_pop.close();
    os_pop.close();
  }
  is.close();
  os.close();

  // update haps file.
  std::vector<char> sequence(data.N), sequence_new(remaining_haps.size());
  int bp, num_carriers;
  int L_new = 0;
  for(int snp = 0; snp < data.L; snp++){

    m_hap.ReadSNP(sequence, bp);  
    num_carriers = 0;
    for(int k = 0; k < remaining_haps.size(); k++){
      sequence_new[k] = sequence[remaining_haps[k]];
      if(sequence_new[k] == '1'){
        num_carriers++;
      }
    }
    if(!remove_fixed || (num_carriers > 0 && num_carriers < remaining_haps.size())){
      m_hap.DumpSNP(sequence_new, bp, fp);
      L_new++;
    }

  }
  fclose(fp);
  m_hap.CloseFile();

  std::cerr << "Removed " << data.L - L_new << " SNPs." << std::endl;
  std::cerr << "Output written to " << options["output"].as<std::string>() + ".haps" << " and " << options["output"].as<std::string>() + ".sample" << std::endl;


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
FilterHapsUsingMask(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("haps") || !options.count("sample") || !options.count("mask") || !options.count("output") ){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, sample, mask, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Filter haps using mask.." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Filtering haps using mask.. " << std::endl;


  //read in mask
  fasta mask;
  mask.Read(options["mask"].as<std::string>());

  //parse haps line by line and delete SNPs appropriately
  haps m_hap(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str()); //use only to get N and L;
  m_hap.CloseFile();
  std::vector<char> sequence(m_hap.GetN());
  std::vector<int> pos(m_hap.GetL(), 0), dist(m_hap.GetL(), 0);
  Data data(m_hap.GetN(), m_hap.GetL());

  igzstream is(options["haps"].as<std::string>());
  if(is.fail()){
    std::cerr << "Error opening file." << std::endl;
    exit(1);
  }
  std::ofstream os(options["output"].as<std::string>() + ".haps");
  if(os.fail()){
    std::cerr << "Error opening file." << std::endl;
    exit(1);
  }

  char chr[1024];
  char rsid[1024];
  int bp;
  int p_next, p_prev = 0;
  std::string line;
  std::string::iterator it_mask, it_start, it_end;

  //int mask_threshold = 1801;
  int mask_threshold = 2000;
  int d_num_nonpass_vincity;
  int passing_snp = 0;
  for(int snp = 0; snp < data.L; snp++){

    assert(getline(is, line));
    sscanf(line.c_str(), "%s %s %d", chr, rsid, &bp);

    d_num_nonpass_vincity = 0; 
    if(mask.seq[bp-1] != 'P'){
      d_num_nonpass_vincity = mask_threshold;
    }else{
      it_start = std::next(mask.seq.begin(), std::max(0, bp - 1000));
      it_end   = std::next(mask.seq.begin(), std::min((int) mask.seq.size(), bp + 1001));
      for(it_mask = it_start; it_mask != it_end; it_mask++){
        if(*it_mask != 'P'){
          d_num_nonpass_vincity++;
        } 
      }
    }

    //only include if we are confident that the snp is correct
    if(d_num_nonpass_vincity < mask_threshold){

      //write SNP
      os << line << "\n";

      pos[passing_snp] = bp;
      //go through mask from p_prev to p and count the number of bases that pass the mask
      int distance = 0; 
      if(passing_snp > 0){
        it_start = std::next(mask.seq.begin(), std::max(0, p_prev - 1000));
        it_end   = std::next(mask.seq.begin(), std::min((int)mask.seq.size(), p_prev + 1001));

        d_num_nonpass_vincity = 0;
        for(it_mask = it_start; it_mask != it_end; it_mask++){
          if(*it_mask != 'P'){
            d_num_nonpass_vincity++;
          } 
        }
        it_end--;

        for(it_mask = std::next(mask.seq.begin(),p_prev); it_mask != std::next(mask.seq.begin(), bp); it_mask++){
          if(*it_start != 'P') d_num_nonpass_vincity--;
          it_start++;
          if(it_end != mask.seq.end()){
            it_end++;
            if(*it_end != 'P') d_num_nonpass_vincity++;
          }
          assert(d_num_nonpass_vincity >= 0);
          if(*it_mask == 'P' && d_num_nonpass_vincity < mask_threshold){
            distance++;
          }  
        }

        if(distance == 0) distance = 1;
        dist[passing_snp - 1] = distance;
      }

      p_prev = bp;
      passing_snp++;

    } 
  }
  is.close();
  os.close();
  dist[passing_snp-1] = 1;

  pos.resize(passing_snp);
  dist.resize(passing_snp);

  FILE* fp_dist = fopen((options["output"].as<std::string>() + ".dist").c_str(), "w");
  std::vector<int>::iterator it_pos = pos.begin();
  fprintf(fp_dist, "#pos dist\n");
  for(std::vector<int>::iterator it_dist = dist.begin(); it_dist != dist.end();){
    fprintf(fp_dist, "%d %d\n", (*it_pos), (*it_dist));
    it_pos++;
    it_dist++;  
  }
  fclose(fp_dist);

  std::cerr << "Removed " << data.L - passing_snp << " SNPs." << std::endl;
  std::cerr << "Output written to " << options["output"].as<std::string>() + ".haps" << " and " << options["output"].as<std::string>() + ".dist" << std::endl;


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
FlipHapsUsingAncestor(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("haps") || !options.count("sample") || !options.count("ancestor") || !options.count("output") ){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, sample, ancestor, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Detemine ancestral allele and flip SNPs." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Detemining ancestral allele and flipping SNPs if necessary... " << std::endl;

  fasta ancestor;
  ancestor.Read(options["ancestor"].as<std::string>());
  haps m_hap(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  m_hap.CloseFile();
  Data data(m_hap.GetN(), m_hap.GetL());

  igzstream is(options["haps"].as<std::string>().c_str());
  if(is.fail()){
    std::cerr << "Error opening file." << std::endl;
    exit(1);
  }
  std::ofstream os(options["output"].as<std::string>() + ".haps");
  if(os.fail()){
    std::cerr << "Error opening file." << std::endl;
    exit(1);
  }
  std::string line;

  int bp;
  char chr[1024];
  char rsid[1024];
  char ancestral[1024], alternative[1024];
  char ancestral_allele;
  std::vector<char> sequence(data.N);
  std::string::iterator it_line;
  int number_flipped = 0;
  int removed_snps = 0;
  for(int snp = 0; snp < data.L; snp++){

    assert(getline(is, line));
    sscanf(line.c_str(), "%s %s %d %s %s", chr, rsid, &bp, ancestral, alternative);
    ancestral_allele = std::toupper(ancestor.seq[bp-1]);

    if(strlen(ancestral) == 1 || strlen(alternative) == 1){
      //if(ancestral_allele == 'A' || ancestral_allele == 'C' || ancestral_allele == 'G' || ancestral_allele == 'T'){

      if(ancestral_allele == ancestral[0] && strlen(ancestral) == 1){

        it_line = line.begin();
        //chr
        while(*it_line != ' ') it_line++;
        it_line++;
        //rsid
        while(*it_line != ' ') it_line++;
        it_line++;
        //bp
        while(*it_line != ' ') it_line++;
        it_line++;
        //ancestral
        while(*it_line != ' ') it_line++;
        it_line++;
        //alternative
        while(*it_line != ' ') it_line++;
        it_line++;

        bool is_snp = false;
        for(; it_line != line.end(); it_line++){
          if(*it_line == '1'){
            is_snp = true;
            break;
          }
        }

        if(is_snp){
          os << line << "\n";
        }else{
          removed_snps++;
        }

      }else if(ancestral_allele == alternative[0] && strlen(alternative) == 1){

        number_flipped++;
        it_line = line.begin();
        //chr
        while(*it_line != ' ') it_line++;
        it_line++;
        //rsid
        while(*it_line != ' ') it_line++;
        it_line++;
        //bp
        while(*it_line != ' ') it_line++;
        it_line++;
        //ancestral
        *it_line = alternative[0];
        it_line++;
        *it_line=' ';
        it_line++;
        //alternative
        for(int c = 0; c < strlen(ancestral); c++){
          *it_line=ancestral[c];
          it_line++;
        }
        assert(*it_line == ' ');
        it_line++;

        bool is_snp = false;
        for(; it_line != line.end(); it_line++){
          if(*it_line == '0'){
            *it_line = '1';
            is_snp = true;
          }else if(*it_line == '1'){
            *it_line = '0';
          }
        }

        if(is_snp){
          os << line << "\n";
        }else{
          removed_snps++;
        }

      }else{
        removed_snps++;
        //std::cerr << ancestral_allele << " " << ancestral << " " << alternative << " " << ancestor.seq[bp-1] << " " << ancestor.seq[bp+1] << std::endl;
      }

    }else{
      removed_snps++;
    }

    //}

  }
  is.close();
  os.close();


  std::cerr << "Had to remove " << removed_snps << " SNPs because of non-matching nucleotides" << std::endl;
  std::cerr << "Number of flipped SNPs is " << number_flipped << "." << std::endl;
  std::cerr << "Output written to " << options["output"].as<std::string>() + ".haps." << std::endl;

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
GenerateSNPAnnotations(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("haps") || !options.count("sample") || !options.count("poplabels") || !options.count("output") ){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, sample, poplabels, output. Optional: ancestor, mut." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Generate additional annotation for SNPs." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Generating additional annotation for SNPs... " << std::endl;

  bool is_ancestor = false;
  if(options.count("ancestor")) is_ancestor = true;

  fasta ancestor;
  if(is_ancestor) ancestor.Read(options["ancestor"].as<std::string>());

  haps m_hap(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  Data data(m_hap.GetN(), m_hap.GetL());

  Mutations mut;
  bool is_mut = false;
  if(options.count("mut")){
    is_mut = true;
    mut.Read(options["mut"].as<std::string>());
  }

  Sample sample;
  sample.Read(options["poplabels"].as<std::string>());

  std::ofstream os(options["output"].as<std::string>() + ".annot");
  if(os.fail()){
    std::cerr << "Error opening file." << std::endl;
    exit(1);
  }

  os << "upstream_allele;downstream_allele;";
  for(int p = 0; p < sample.groups.size(); p++){
    os << sample.groups[p] << ";";
  }
  os << "\n";

  int bp;
  char nucl;
  std::vector<char> sequence(data.N);
  std::vector<int> carriers_by_pop(sample.groups.size());
  int percentage = 0;
  std::cerr << "[" << percentage << "%]\r";
  std::cerr.flush();
  for(int snp = 0; snp < data.L; snp++){

    m_hap.ReadSNP(sequence, bp);

    //std::cerr << (int)(data.L/100) << " " << snp % (int)(data.L/100) << std::endl;
    if((snp % (int)(data.L/100)) == 0){
      std::cerr << "[" << percentage << "%]\r";
      percentage++;
      std::cerr.flush();
    }

    //get nucleotide before and after in sequence
    if(is_ancestor){
      if(bp > 1){
        nucl = std::toupper(ancestor.seq[bp-2]);
        if(nucl == 'A' || nucl == 'C' || nucl == 'G' || nucl == 'T'){
          os << nucl << ";";
          if(is_mut) mut.info[snp].upstream_base = nucl;
        }else{
          os << "NA;";
        }
      }else{
        os << "NA;";
      }
      if(bp < ancestor.seq.size()){
        nucl = std::toupper(ancestor.seq[bp]);
        if(nucl == 'A' || nucl == 'C' || nucl == 'G' || nucl == 'T'){
          os << nucl << ";";
          if(is_mut) mut.info[snp].downstream_base = nucl;
        }else{
          os << "NA;";
        }
      }else{
        os << "NA;";
      }
    }else{
      os << "NA;NA;";
    }

    //calculate number of carriers in each population
    std::fill(carriers_by_pop.begin(), carriers_by_pop.end(), 0);
    for(int i = 0; i < data.N; i++){
      if(sequence[i] == '1'){
        carriers_by_pop[sample.group_of_haplotype[i]]++;
      }
    }
    for(std::vector<int>::iterator it_carriers = carriers_by_pop.begin(); it_carriers != carriers_by_pop.end(); it_carriers++){
      os << *it_carriers << ";";
    }
    os << "\n";

    if(is_mut){
      mut.info[snp].freq.resize(carriers_by_pop.size());
      int pop = 0;
      for(std::vector<int>::iterator it_carriers = carriers_by_pop.begin(); it_carriers != carriers_by_pop.end(); it_carriers++){
        mut.info[snp].freq[pop] = *it_carriers;
        pop++;
      } 
    }

  }

  if(options.count("mut")){

    mut.header  = "snp;pos_of_snp;dist;rs-id;tree_index;branch_indices;is_not_mapping;is_flipped;age_begin;age_end;ancestral_allele/alternative_allele;"; 
    mut.header += "upstream_allele;downstream_allele;";
    for(std::vector<std::string>::iterator it_groups = sample.groups.begin(); it_groups != sample.groups.end(); it_groups++){
      mut.header += *it_groups + ";";
    }
    mut.Dump(options["output"].as<std::string>() + ".mut");
    std::cerr << "Output written to " << options["output"].as<std::string>() + ".annot and " << options["output"].as<std::string>() + ".mut" << std::endl;
  }else{
    std::cerr << "Output written to " << options["output"].as<std::string>() + ".annot" << std::endl;
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

}

