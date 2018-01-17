#include <iostream>
#include <utility>

#include "collapsed_matrix.hpp"
#include "arg.hpp"
#include "arg_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

void
SummarizeAgeLineagesDistribution(const std::vector<std::string>& filenames_bin, const std::vector<std::string>& filenames_mut, CollapsedMatrix<float>& age_lineages){

  std::string line, last_line, read;
  int first_pos, last_pos;

  std::fill(age_lineages.vbegin(), age_lineages.vend(), 0.0);
  CollapsedMatrix<float> age_lineages_for_one_chr;
  for(int i = 0; i < (int) filenames_bin.size(); i++){
    FILE* fp = fopen(filenames_bin[i].c_str(), "rb");
    age_lineages_for_one_chr.ReadFromFile(fp); 
    fclose(fp);

    //read length of chromosome;
    int j = 0;

    std::ifstream is(filenames_mut[i]);
    if(is.fail()){
      std::cerr << "Error while opening file." << std::endl;
      exit(1);
    }   
    getline(is, line);
    getline(is, line);
    getline(is, line);
    read.clear();
    while(line[j] != ';') j++;
    j++; 
    while(line[j] != ';'){
      read += line[j];
      j++;
    }
    first_pos = std::stoi(read);
    while(getline(is, line)) last_line = line;
    j = 0;
    read.clear();
    while(last_line[j] != ';') j++;
    j++; 
    while(last_line[j] != ';'){
      read += last_line[j];
      j++;
    }
    last_pos = std::stoi(read);
    is.close();

    //add to age_lineages
    float weight = (last_pos - first_pos)/1e9;
    std::vector<float>::iterator it_age_lineages_for_one_chr = age_lineages_for_one_chr.vbegin();
    std::vector<float>::iterator it_age_lineages             = age_lineages.vbegin();

    for(; it_age_lineages != age_lineages.vend(); it_age_lineages++){
      *it_age_lineages += weight * (*it_age_lineages_for_one_chr);
      it_age_lineages_for_one_chr++;
    }
  }
   
}

void
AgeFrequencyDistribution(const Data& data, CollapsedMatrix<float>& age_lineages, CollapsedMatrix<float>& age_frequency){

  const float log_Nsquared = log(data.N * data.N);
  for(int x = 0; x < data.N; x++){
    int epoche_index = 0;
    for(std::vector<float>::iterator it_age_frequency = age_frequency.rowbegin(x); it_age_frequency != age_frequency.rowend(x);){
      *it_age_frequency = 0.0;
      int num_lineages  = data.N;
      float x_fraction  = x/(float)data.N;
      //calculate E[A(t) * (A(t) - 1) * (1-x)^(A(t) - 2)]
      for(std::vector<float>::iterator it_age_lineages = age_lineages.rowbegin(epoche_index); it_age_lineages != age_lineages.rowend(epoche_index); it_age_lineages++){
        *it_age_frequency += *it_age_lineages * exp(log( num_lineages * (num_lineages-1.0) ) + (num_lineages - 2.0) * log( 1.0 - x_fraction ) - log_Nsquared); 
        num_lineages--;
      } 
      epoche_index++;
      it_age_frequency++;
    }

    //calculate normalizing constant
    float normalizing_constant = 0.0;
    for(std::vector<float>::iterator it_age_frequency = age_frequency.rowbegin(x); it_age_frequency != age_frequency.rowend(x); it_age_frequency++){
      normalizing_constant += *it_age_frequency;
    }
    assert(normalizing_constant > 0.0);

    //calculate cumulative distribution function
    std::vector<float>::iterator it_age_frequency_prev = age_frequency.rowbegin(x);
    *it_age_frequency_prev /= normalizing_constant;
    for(std::vector<float>::iterator it_age_frequency = std::next(age_frequency.rowbegin(x),1); it_age_frequency != age_frequency.rowend(x); it_age_frequency++){
      *it_age_frequency = *it_age_frequency_prev + *it_age_frequency/normalizing_constant;
      it_age_frequency_prev++;
    }

  }

}



int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options
  cxxopts::Options options("GetAgeFrequencyDistribution");
  options.add_options()
    ("help", "Print help")	
    ("S,num_haps", "Number of haplotypes", cxxopts::value<int>())
    ("s,chr_start", "First Chromosome", cxxopts::value<int>())
    ("e,chr_end", "Last Chromosome", cxxopts::value<int>())
    ("i,input", "Filename .mut file without file extension", cxxopts::value<std::string>())
    ("o,output", "Output file", cxxopts::value<std::string>());

  options.parse(argc, argv);
  bool help = false;
  if(!options.count("num_haps") || !options.count("chr_start") || !options.count("chr_end") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .mut file and creates .sele file containing info on age-frequency." << std::endl;
    exit(0);
  }  

  ///////
  Data data(options["num_haps"].as<int>(), 1);

  //create age_frequency matrix
  CollapsedMatrix<float> age_frequency;
  int num_epoches = 150;
  std::vector<float> epoches(num_epoches);
  for(int e = 0; e < num_epoches; e++){
    epoches[e] = std::exp(e/10.0);
  }

  int chr_start = options["chr_start"].as<int>();
  int chr_end   = options["chr_end"].as<int>();
  assert(chr_start <= chr_end);

  std::vector<std::string> filenames_bin, filenames_mut;
  for(int chr = chr_start; chr <= chr_end; chr++){
    filenames_mut.push_back(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
    filenames_bin.push_back(options["output"].as<std::string>() + "_age_lineages_chr" + std::to_string(chr) + ".bin");
  }

  //summarize age lineage distribution calculated previously
  CollapsedMatrix<float> age_lineages;
  age_lineages.resize(epoches.size(), data.N-1); //stores distribution over number of lineages at times given in epoches
  SummarizeAgeLineagesDistribution(filenames_bin, filenames_mut, age_lineages);

  age_frequency.resize(data.N, epoches.size());
  AgeFrequencyDistribution(data, age_lineages, age_frequency);

  std::ofstream os("age_frequency_distribution.txt");
  for(int i = 0; i < data.N; i++){
    for(int j = 0; j < epoches.size(); j++){
      os << age_frequency[i][j] << " ";
    }
    os << "\n";
  }

  //Dump result to file
  FILE* fp;
  std::string filename = options["output"].as<std::string>() + "_age_frequency.bin";
  fp = fopen(filename.c_str(), "wb");  
  age_frequency.DumpToFile(fp);
  fclose(fp);

  return 0;

}
