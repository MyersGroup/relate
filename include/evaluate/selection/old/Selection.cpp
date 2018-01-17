#include <iostream>
#include <utility>

#include "collapsed_matrix.hpp"
#include "arg.hpp"
#include "arg_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options
  cxxopts::Options options("Selection");
  options.add_options()
    ("help", "Print help")	
    ("S,num_haps", "Number of haplotypes", cxxopts::value<int>())
    ("L,num_SNPs", "Number of SNPs", cxxopts::value<int>())
    ("c,chr", "Chromosome index", cxxopts::value<int>())
    ("f,age_frequency", "filename of file containing age frequency distribution", cxxopts::value<std::string>())
    ("i,input", "Filename .mut file without file extension", cxxopts::value<std::string>())
    ("o,output", "Output file", cxxopts::value<std::string>());

  options.parse(argc, argv);
  bool help = false;
  if(!options.count("num_haps") || !options.count("num_SNPs") || !options.count("chr") || !options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .mut file and creates .sele file containing info on age-frequency." << std::endl;
    exit(0);
  }  


  ///////
  Data data(options["num_haps"].as<int>(), options["num_SNPs"].as<int>());

  //read mutations file
  Mutations mutations(data);
  std::cerr << options["input"].as<std::string>() + "_chr" + std::to_string(options["chr"].as<int>()) + ".mut" << std::endl;
  mutations.ReadLongFormatFromFile(options["input"].as<std::string>() + "_chr" + std::to_string(options["chr"].as<int>()) + ".mut");

  //read population labels
  std::ifstream is_sample;
  is_sample.open("population_labels.txt");
  if(is_sample.fail()){
    std::cerr << "Error while reading excluded samples." << std::endl;
    exit(1);
  } 
  
  std::vector<std::string> sample_label, unique_sample_label;
  std::map<std::string, int> population_id;
  std::vector<int> group_sizes;
  std::string line, read_label, read_size;
  int i, count_groups = 0;
  while(getline(is_sample, line)){

    i = 0;
    //while(line[i] != ' ') i++;
    //i++;
    read_label.clear();
    while(line[i] != ' '){
      read_label += line[i];
      i++;
    }
    i++;
    while(line[i] != ' ') i++;
    i++;
    read_size.clear();
    while(i < (int) line.size()){
      read_size += line[i];
      i++;
    }
    i++;

    sample_label.push_back(read_label); 
    if(population_id.find(read_label) == population_id.end() ){
      unique_sample_label.push_back(read_label);
      population_id[read_label]    = count_groups;
      group_sizes.push_back(std::stoi(read_size));
      count_groups++;
    }
    
  }

  ////////////////////

  //read age_frequency matrix
  CollapsedMatrix<float> age_frequency;
  int num_epoches = 150;
  std::vector<float> epoches(num_epoches);
  for(int e = 0; e < num_epoches; e++){
    epoches[e] = std::exp(e/10.0);
  }

  FILE* fp;
  fp = fopen(options["age_frequency"].as<std::string>().c_str(), "rb");  
  age_frequency.ReadFromFile(fp);
  fclose(fp);

 
  /* 
  for(int i = 0; i < age_frequency.size(); i++){
    for(int j = 0; j < age_frequency.subVectorSize(0) - 1; j++){
      std::cout << age_frequency[i][j+1] - age_frequency[i][j] << " ";
    }
    std::cout << std::endl;
  }
  */

  ///////////////////////////
  //create .sele file
  std::ofstream os(options["output"].as<std::string>() + "_chr" + std::to_string(options["chr"].as<int>()) + ".sele");

  os << "snp pos rs_id freq mean_startage all ";
  for(std::vector<std::string>::iterator it_unique_sample = unique_sample_label.begin(); it_unique_sample != unique_sample_label.end(); it_unique_sample++){
    //os << *it_unique_sample << " ";
  }
  os << "is_flipped\n";

  int snp = 0, epoche_index;
  std::vector<float> freq(population_id.size());
  int freq_all;
  float epoche_interpol, age;

  for(std::vector<SNPInfo>::iterator it_info = mutations.info.begin(); it_info != mutations.info.end(); it_info++){

    //if mapping
    if((*it_info).branch.size() == 1 && (*it_info).age_begin > 0.0){
      
      //calculate frequency of mutation
      std::fill(freq.begin(), freq.end(), 0);
      freq_all = 0.0;
      for(int k = 0; k < (int) (*it_info).freq.size(); k++){
        freq[population_id[sample_label[k]]] += (*it_info).freq[k];
        freq_all                             += (*it_info).freq[k];
      }
     
      //age = (*it_info).age_begin;
      //age = ((*it_info).age_begin + (*it_info).age_end)/2.0;
      age = (*it_info).age_end;
      if(age < 1.0) age = 1.0;

      //output freq, age
      os << (*it_info).snp_id << " " << (*it_info).pos << " " << (*it_info).rs_id << " " << freq_all << " " << age << " ";
   
      assert(freq_all <= data.N);
      epoche_index    = 10.0*log(age);
      assert(epoche_index >= 0);
      if(epoches[epoche_index] - age > 0.0) std::cerr << epoches[epoche_index] << " " << age << std::endl;
      //assert(epoches[epoche_index] - age <= 0.0);
      epoche_interpol = (age - epoches[epoche_index])/(epoches[epoche_index+1] - epoches[epoche_index]);

      
      //get p-value for this snp 
      if(freq_all == data.N){
        os << "1.0";
      }else{     
        os << age_frequency[freq_all][epoche_index] + epoche_interpol * (age_frequency[freq_all][epoche_index+1] - age_frequency[freq_all][epoche_index]);
      }
      os << " ";

      //get p-value when only considering subpopulation
      /*
      for(int k = 0; k < (int) freq.size(); k++){
        freq_all     = (freq[k] * data.N/((float) group_sizes[k])); //group size of group that I am looking at
        assert(freq_all <= data.N);
        if(!(*it_info).rs_id.compare("rs4988235")){
          std::cerr << freq[k] << " " << freq_all << std::endl;
        }
        if(freq_all == data.N){
          os << "1.0";
        }else{
          os << age_frequency[freq_all][epoche_index] + epoche_interpol * (age_frequency[freq_all][epoche_index+1] - age_frequency[freq_all][epoche_index]);
        }
        os << " ";
      } 
      */
      os <<  (*it_info).flipped << "\n";

    }
    snp++;

  }

  os.close();


  return 0;

}
