#include <iostream>
#include <utility>

#include "collapsed_matrix.hpp"
#include "arg.hpp"
#include "arg_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

void
AgeLineagesDistribution(const Data& data, const std::string& filename, CollapsedMatrix<float>& age_lineages, const std::vector<float>& epoches, const int first_pos, const int last_pos){

  age_lineages.resize(epoches.size(), data.N-1); //stores distribution over number of lineages at times given in epoches

  //go through trees, at each time in epoches, check how many lineages are remaining. Add to entry in age_lineages
  //Using age_lineages, calculate age_frequency

  Arg arg;
  float factor = 0.0;
  const float denominator = last_pos - first_pos;  
  std::vector<float> coordinates(2*data.N - 1);

  arg.ReadArgLong(filename);
  CorrTrees::iterator it_arg = arg.seq.begin();

  for(; it_arg != arg.seq.end(); it_arg++){

    if(std::next(it_arg,1) == arg.seq.end()){
      factor = (last_pos - (*it_arg).pos)/denominator; 
    }else{
      factor = ((*std::next(it_arg,1)).pos - (*it_arg).pos)/denominator; 
    }

    //find number of lineages at times in epoche:
    //determine coordinates, sort coordinates, count
    (*it_arg).tree.GetCoordinates(coordinates);
    std::sort(coordinates.begin(), coordinates.end());

    std::vector<float>::const_iterator it_epoches = epoches.begin();
    std::vector<float>::iterator it_coords  = std::next(coordinates.begin(), data.N);
    int epoche_index = 0;
    int coords_index = 0;

    for(; it_epoches != epoches.end();){ 
      while(*it_epoches >= *it_coords){
        if(it_coords != coordinates.end()){
          coords_index++;
          it_coords++;
        }else{
          break;
        }
      }

      if(*it_epoches < *it_coords){
        assert(coords_index < data.N-1);
        age_lineages[epoche_index][coords_index] += factor;
        epoche_index++;
        it_epoches++; 
      }else{
        break;
      }
    }

  }

}

int main(int argc, char* argv[]){

  //////////////////////////////////
  //Program options
  cxxopts::Options options("GetAgeLineagesDistribution");
  options.add_options()
    ("help", "Print help")	
    ("S,num_haps", "Number of haplotypes", cxxopts::value<int>())
    ("L,num_SNPs", "Number of SNPs", cxxopts::value<int>())
    ("c,chr", "Chromosome index", cxxopts::value<int>())
    ("i,input", "Filename .mut file without file extension", cxxopts::value<std::string>())
    ("o,output", "Output file", cxxopts::value<std::string>());

  options.parse(argc, argv);
  bool help = false;
  if(!options.count("num_haps") || !options.count("num_SNPs")  || !options.count("chr") || !options.count("input") || !options.count("output")){
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
  mutations.ReadLongFormatFromFile(options["input"].as<std::string>() + "_chr" + std::to_string(options["chr"].as<int>()) + ".mut");

  // parse pos
  int first_pos, last_pos;
  first_pos = mutations.info[0].pos;
  last_pos  = mutations.info[data.L-1].pos; 


  //create age_lineages matrix
  int num_epoches = 150;
  std::vector<float> epoches(num_epoches);
  for(int e = 0; e < num_epoches; e++){
    epoches[e] = std::exp(e/10.0);
  }

  CollapsedMatrix<float> age_lineages;
  std::string filename = options["input"].as<std::string>() + "_chr" + std::to_string(options["chr"].as<int>()) + ".arg";
  AgeLineagesDistribution(data, filename, age_lineages, epoches, first_pos, last_pos);

  //Dump result to file
  FILE* fp;
  filename = options["output"].as<std::string>() + "_age_lineages" + "_chr" + std::to_string(options["chr"].as<int>()) + ".bin";
  fp = fopen(filename.c_str(), "wb");  
  age_lineages.DumpToFile(fp);
  fclose(fp);

  return 0;

}
