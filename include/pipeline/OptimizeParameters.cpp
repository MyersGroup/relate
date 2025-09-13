/*! \file OptimizeParameters.cpp 
 *  \brief Optimize parameters for painting
 *
 *  Input: data
 *  Output: optimal parameters for painting (theta and recombination_factor)
 */

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <random>


#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "fast_painting.hpp"
#include "anc_builder.hpp"

#include "MakeChunks.cpp"
#include "Paint.cpp"
#include "Clean.cpp"

int OptimizeParameters(cxxopts::Options& options){

  bool help = false;
  if(!options.count("haps") || !options.count("sample") || !options.count("map") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: haps, sample, map, output. Optional: dist." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use to make smaller chunks from the data." << std::endl;
    exit(0);
  }

  std::cerr << "############" << std::endl;
  std::cerr << "Optimizing Parameters..." << std::endl;

  int N, L;
  double memory_size;
  int start_chunk, end_chunk;
  std::string file_out = options["output"].as<std::string>() + "/";

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Using:" << std::endl;
  std::cerr << "  " << options["haps"].as<std::string>() << std::endl;    
  std::cerr << "  " << options["sample"].as<std::string>() << std::endl;
  std::cerr << "  " << options["map"].as<std::string>() << std::endl;
  if(!options.count("coal")){
    std::cerr << "with mu = " << options["mutation_rate"].as<float>() << " and 2Ne = " << options["effectiveN"].as<float>() << "." << std::endl;
  }else{
    std::cerr << "with mu = " << options["mutation_rate"].as<float>() << " and coal = " << options["coal"].as<std::string>() << "." << std::endl;
  }

  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

  MakeChunks(options);

  FILE* fp = fopen((file_out + "parameters.bin").c_str(), "r");
  assert(fp != NULL);
  fread(&N, sizeof(int), 1, fp);
  fread(&L, sizeof(int), 1, fp);
  fread(&end_chunk, sizeof(int), 1, fp);
  fread(&memory_size, sizeof(double), 1, fp);
  fclose(fp);

  end_chunk--;
  start_chunk = 0;

  Data data(N,L);
  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Read " << data.N << " haplotypes with " << data.L << " SNPs per haplotype." << std::endl;
  std::cerr << "Expected minimum memory usage: " << memory_size << "Gb." << std::endl;
  std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

  std::vector<float> theta = {1e-3, 2.5e-3, 5e-3, 7.5e-3, 1e-2, 2.5e-2, 5e-2, 7.5e-2, 1e-1};
  std::vector<float> rec_factor = {0.0001,0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 10};

  float val;
  std::string line;
  if(options.count("input")){

    std::ifstream is(options["input"].as<std::string>());

    theta.clear();
    rec_factor.clear();

    getline(is, line);
    std::istringstream itheta;
    itheta.str(line);
    while(itheta >> val){
      if(val >= 1.0 || val <= 0){
        std::cerr << "Error: theta value has to be in (0,1)" << std::endl;
        exit(1);
      }
      theta.push_back(val);
    }

    getline(is, line);
    std::istringstream irec;
    irec.str(line);
    while(irec >> val){
      if(val <= 0){
        std::cerr << "Error: rho value has to be positive" << std::endl;
        exit(1);
      }
      rec_factor.push_back(val);
    }

    is.close();

  }

  std::vector<std::vector<int>> num_notmapping(theta.size());
  for(int i = 0; i < num_notmapping.size(); i++){
    num_notmapping[i].resize(rec_factor.size());
    std::fill(num_notmapping[i].begin(), num_notmapping[i].end(), 0);
  }

  for(int c = start_chunk; c <= end_chunk; c++){

    filesys f;
    f.MakeDir((file_out + "chunk_" + std::to_string(c) + "/").c_str());
    f.MakeDir((file_out + "chunk_" + std::to_string(c) + "/paint/").c_str());

    std::cerr << "---------------------------------------------------------" << std::endl;
    std::cerr << "Starting chunk " << c << " of " << end_chunk << "." << std::endl;
    std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

    std::string file_out = options["output"].as<std::string>() + "/";

    int N, L, num_windows;
    std::vector<int> window_boundaries;
    FILE* fp = fopen((file_out + "parameters_c" + std::to_string(c) + ".bin").c_str(), "r");
    assert(fp != NULL);
    fread(&N, sizeof(int), 1, fp);
    fread(&L, sizeof(int), 1, fp);
    fread(&num_windows, sizeof(int), 1, fp);
    window_boundaries.resize(num_windows);
    fread(&window_boundaries[0], sizeof(int), num_windows, fp);
    fclose(fp);
    num_windows--;

    bool use_transitions = true;

    Data data((file_out + "chunk_" + std::to_string(c) + ".hap").c_str(), (file_out + "chunk_" + std::to_string(c) + ".bp").c_str(), (file_out + "chunk_" + std::to_string(c) + ".dist").c_str(), (file_out + "chunk_" + std::to_string(c) + ".r").c_str(), (file_out + "chunk_" + std::to_string(c) + ".rpos").c_str(), (file_out + "chunk_" + std::to_string(c) + ".state").c_str()); //struct data is defined in data.hpp
    data.name = (file_out + "chunk_" + std::to_string(c) + "/paint/relate");
    const std::string dirname = file_out + "chunk_" + std::to_string(c) + "/";

    int k = 5;
    std::vector<int> random_sections;
    std::vector<int> pool(static_cast<std::size_t>(num_windows));
    std::iota(pool.begin(), pool.end(), 0);

    std::mt19937 gen(std::time(0) + getpid());

    if (pool.size() <= k){
      random_sections = pool;
    }else{
      random_sections.reserve(k); 
      std::shuffle(pool.begin(), pool.end(), gen);
      random_sections = std::vector<int>(pool.begin(), pool.begin() + static_cast<std::ptrdiff_t>(k));
      std::sort(random_sections.begin(), random_sections.end());
    }

    //I could change selected SNPs to 0 in data.sequence
    //and record their status in a separate vector.
    //This can be passed to OptimizeParameters and then remapped (including full BuildTopology machinery)
    
    std::vector<std::vector<char>> sequence; 
    std::vector<int> snps;
    int snp;
    int num_test_snps = 10000;
    for(int s = 0; s < random_sections.size(); s++){
      int section = random_sections[s];
      int section_startpos = window_boundaries[section];
      int section_endpos   = window_boundaries[section+1]-1;

      std::vector<int> random_snps;
      std::vector<int> pool(static_cast<std::size_t>(section_endpos-section_startpos));
      std::iota(pool.begin(), pool.end(), section_startpos);

      num_test_snps = std::min(10000, (int)(0.1*(section_endpos-section_startpos)));

      std::shuffle(pool.begin(), pool.end(), gen);

      for(int index = 0; index < pool.size(); index++){
        int AF = 0;
        snp = pool[index];
        for(int n = 0; n < data.N; n++){
          if(data.sequence[snp][n] == '1') AF++;
          if(AF > 1){
            random_snps.push_back(snp);
            break;
          }
        }
        if(random_snps.size() == num_test_snps) break;
      }

      std::sort(random_snps.begin(), random_snps.end());
      for(int i = 0; i < random_snps.size(); i++){
        snp = random_snps[i];
        snps.push_back(snp);
        sequence.emplace_back(data.sequence[snp], data.sequence[snp] + data.N);
        for(int j = 0; j < data.N; j++){
          data.sequence[snp][j] = '1';
        }
      }

    }

    Leaves sequences_carrying_mutation;
    sequences_carrying_mutation.member.resize(data.N);

    std::vector<double> rec_rate = data.r;
    for(int theta_index = 0; theta_index < (int) theta.size(); theta_index++){
      for(int rec_index = 0; rec_index < (int) rec_factor.size(); rec_index++){

        float mean_rec = 0.0;
        data.theta     = theta[theta_index];
        data.ntheta    = 1.0 - data.theta;

        for(int l = 0; l < (int)data.r.size(); l++){
          data.r[l] = rec_rate[l] * rec_factor[rec_index];
          mean_rec  += data.r[l];
        }

        //Paint stepping stones
        //Paint(options,c);
        char filename[1024];
        std::vector<FILE*> pfiles(num_windows);
        for(int w = 0; w < num_windows; w++){
          snprintf(filename, sizeof(char) * 1024, "%s_%i.bin", data.name.c_str(), w);
          pfiles[w] = fopen(filename, "wb");
          assert(pfiles[w] != NULL);
        }

        for(int hap = 0; hap < data.N; hap++){
          FastPainting painter(data);
          painter.PaintSteppingStones(data, window_boundaries, pfiles, hap);
        }

        for(int w = 0; w < num_windows; w++){  
          fclose(pfiles[w]);
        }

        int snpindex = 0;
        int section;
        for(int i = 0; i < random_sections.size(); i++){

          section = random_sections[i];

          AncesTree anc;
          AncesTreeBuilder ancbuilder(data);

          int section_startpos = window_boundaries[section];
          int section_endpos   = window_boundaries[section+1]-1;
          if(section_endpos >= data.L) section_endpos = data.L-1;

          int seed = c + section + std::time(0) + getpid();

          //I want to remove some SNPs, paint and build, then remap
          //num_notmapping[theta_index][rec_index] += ancbuilder.OptimizeParameters(section, section_startpos, section_endpos, data, seed);

          ///////////////////
          ancbuilder.BuildTopology(section, section_startpos, section_endpos, data, anc, rand(), true, 10000);
          //map mutations to anc.

          //std::cerr << "first: " << snpindex << " " << snps[snpindex] << " " << (*anc.seq.begin()).pos << " " << section << " " << section_startpos << " " << section_endpos << std::endl;
          while(snps[snpindex] < (*anc.seq.begin()).pos){
            //std::cerr << "skipping" << std::endl;
            snpindex++;
          }
          //if(snpindex < snps.size()){
           // std::cerr << "second: " << snpindex << " " << snps[snpindex] << " " << (*anc.seq.begin()).pos << " " << section << " " << section_startpos << " " << section_endpos << std::endl;
          //}
          for(CorrTrees::iterator it_seq = anc.seq.begin(); it_seq != std::prev(anc.seq.end(),1); it_seq++){

            while((*it_seq).pos <= snps[snpindex] && (*std::next(it_seq,1)).pos > snps[snpindex]){

              sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
              for(int n = 0; n < data.N; n++){
                if(sequence[snpindex][n] == '1'){
                  sequences_carrying_mutation.member[n] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
                  sequences_carrying_mutation.num_leaves++;
                }else{
                  sequences_carrying_mutation.member[n] = 0;
                }
              }

              if(ancbuilder.IsSNPMapping((*it_seq).tree, sequences_carrying_mutation, snpindex) != 1){
                num_notmapping[theta_index][rec_index]++;
              }
              snpindex++;

              if(snpindex == snps.size()) break;

            }
            if(snpindex == snps.size()) break;

          }

        }

        std::cerr << theta[theta_index] << " " << rec_factor[rec_index] << " " << num_notmapping[theta_index][rec_index] << std::endl;
      }

    }

  }

  Clean(options);

  std::ofstream os(options["output"].as<std::string>() + ".opt");
  for(int theta_index = 0; theta_index < (int) theta.size(); theta_index++){
    for(int rec_index = 0; rec_index < (int) rec_factor.size(); rec_index++){
      os << theta[theta_index] << " " << rec_factor[rec_index] << " " << num_notmapping[theta_index][rec_index] << std::endl;
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

  return 0;
}
