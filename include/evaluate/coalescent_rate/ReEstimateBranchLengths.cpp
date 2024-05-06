//Inference of branch lengths using EM and MCMC.

#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include "gzstream.hpp"
#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "branch_length_estimator.hpp"
#include "anc_builder.hpp"

void ShowProgress(int progress){

  std::cerr << "[" << progress << "%]\r";
  std::cerr.flush();
  /*
     int p = 0;
     std::cerr << "[";
     for(; p < progress; p++){
     std::cerr << "*";
     }
     for(; p < 100; p++){
     std::cerr << " ";
     }
     std::cerr << "]" << progress << "%\r";
     std::cerr.flush();
     */

}

int ReEstimateBranchLengths(cxxopts::Options& options){

  int seed;
  if(!options.count("seed")){
    seed = std::time(0) + getpid();
  }else{
    seed = options["seed"].as<int>();
    srand(seed);
    std::string name = options["input"].as<std::string>();
    int tmp = 0;
    for(int i = 0; i < name.size(); i++){
      if(std::isdigit(name[i])) tmp += name[i]-48; 
    }
    for(int i = 0; i < tmp; i++){
      seed = rand();
    }
  }
  srand(seed);

  int Ne = 3e4;
  double mutation_rate = options["mutation_rate"].as<float>();
  std::string line;
  double tmp;

  //parse data
  int N;
  igzstream is_N(options["input"].as<std::string>() + ".anc");
  if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
  if(is_N.fail()){
    std::cerr << "Error while opening .anc file." << std::endl;
    exit(1);
  } 
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  //make this more efficient
  int L = 0;
  if(options.count("dist")){
    igzstream is_L(options["dist"].as<std::string>());
    if(is_L.fail()){
      std::cerr << "Error while opening .dist file." << std::endl;
      exit(1);
    } 
    while(std::getline(is_L, line)){
      ++L;
    }
    L--;
    is_L.close();
  }else{
    igzstream is_L(options["input"].as<std::string>() + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz");
    if(is_L.fail()){
      std::cerr << "Error while opening .mut file." << std::endl;
      exit(1);
    } 
    while(std::getline(is_L, line)){
      ++L;
    }
    L--;
    is_L.close();
  }

  Data data(N, L, Ne, mutation_rate);
  Mutations mut(data);
  mut.Read(options["input"].as<std::string>() + ".mut");

  data.dist.resize(L);
  if(options.count("dist")){
    igzstream is_dist(options["dist"].as<std::string>());
    if(is_dist.fail()){
      std::cerr << "Error while opening " << options["dist"].as<std::string>() << std::endl;
      exit(1);
    }
    getline(is_dist, line); 
    int dtmp, snp = 0;
    while(std::getline(is_dist, line)){
      sscanf(line.c_str(), "%d %d", &dtmp, &data.dist[snp]);
      snp++;
    }
    is_dist.close();
  }else{
    std::vector<int>::iterator it_pos = data.dist.begin();
    for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
      *it_pos = (*it_mut).dist;
      it_pos++;
    }
  }


  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Reinferring branch lengths for " << options["input"].as<std::string>() << " ..." << std::endl;

  Sample sample;
  if(options.count("poplabels")){
    sample.Read(options["poplabels"].as<std::string>());
  }
  std::vector<std::vector<std::vector<double>>> coal_rate_pair;

  // read epochs and population size 
  igzstream is(options["coal"].as<std::string>()); 
  if(is.fail()){
    is.open(options["coal"].as<std::string>() + ".gz");
    if(is.fail()){ 
      std::cerr << "Error while opening " << options["coal"].as<std::string>() << "(.gz)." << std::endl;
      exit(1);
    }
  } 

  std::vector<double> epoch, coalescent_rate;
  getline(is, line);
  std::vector<std::string> groups;
  std::istringstream is_group(line);
  std::string name;
  while(is_group >> name){
    //is_group >> name;
    groups.push_back(name);
  }
  getline(is, line);
  std::istringstream is_epoch(line);
  while(is_epoch >> tmp){
    //is_epoch >> tmp;
    epoch.push_back(tmp/data.Ne);
  }

  if(options.count("poplabels")){

    if(groups.size() != sample.groups.size()){
      std::cerr << "Coal file doesn't contain all groups vs all groups rates" << std::endl;
      exit(1);
    }

    //convert from coal file to poplabels file
    std::vector<int> convert(groups.size(), -1);
    for(int g = 0; g < groups.size(); g++){
      for(int g2 = 0; g2 < sample.groups.size(); g2++){
        if(groups[g] == sample.groups[g2]){
          convert[g] = g2;
          break;
        }
      }
      if(convert[g] == -1){
        std::cerr << "Groups in coal file don't match poplabels file" << std::endl;
        exit(1);
      }
    }

    coal_rate_pair.resize(epoch.size());
    for(int e = 0; e < epoch.size(); e++){
      coal_rate_pair[e].resize(sample.groups.size());
      for(int g = 0; g < sample.groups.size(); g++){
        coal_rate_pair[e][g].resize(sample.groups.size());
      }
    }

    for(int g1 = 0; g1 < sample.groups.size(); g1++){
      for(int g2 = 0; g2 < sample.groups.size(); g2++){
        if(!getline(is,line)){
          std::cerr << "Coal file doesn't contain all groups vs all groups rates" << std::endl;
          exit(1);
        }

        std::istringstream is_pop_size(line);
        is_pop_size >> tmp;
        assert(tmp == g1);
        is_pop_size >> tmp;
        assert(tmp == g2);
        int ep = 0;
        while(is_pop_size){
          is_pop_size >> tmp;
          //if(ep == 0) std::cerr << convert[g1] << " " << convert[g2] << " " << tmp << " " << data.Ne << " " << tmp * data.Ne << std::endl;
          if(tmp == 0.0){
            coal_rate_pair[ep][convert[g1]][convert[g2]] = 5e-10 * data.Ne;
          }else{
            coal_rate_pair[ep][convert[g1]][convert[g2]] = tmp * data.Ne;
          }
          ep++;
          if(ep == epoch.size()) break;
        }

        if(0){
        for(int g1 = 0; g1 < sample.groups.size(); g1++){
          for(int g2 = 0; g2 < sample.groups.size(); g2++){
            for(int i = epoch.size()-1; i > 0; i--){
              if(coal_rate_pair[i-1][g1][g2] == 0){
                if(coal_rate_pair[i][g1][g2] > 0.0){
                  coal_rate_pair[i-1][g1][g2] = coal_rate_pair[i][g1][g2];
                }else{
                  coal_rate_pair[i-1][g1][g2] = 1.0;
                }
              } 
            }
          }
        }
        }

      }
    }

    /*
    for(int g1 = 0; g1 < sample.groups.size(); g1++){
      for(int g2 = 0; g2 < sample.groups.size(); g2++){
        std::cerr << coal_rate_pair[0][g1][g2] << " ";
      }
      std::cerr << std::endl;
    }
    */

  }else{

    getline(is, line);
    std::istringstream is_pop_size(line);
    is_pop_size >> tmp >> tmp;
    while(is_pop_size){
      is_pop_size >> tmp;
      //tmp = 1.0/data.Ne; 
      if(tmp == 0.0 && coalescent_rate.size() > 0){
        if(*std::prev(coalescent_rate.end(),1) > 0.0){
          coalescent_rate.push_back(*std::prev(coalescent_rate.end(),1));
        }
        //coalescent_rate.push_back(1);
      }else{
        coalescent_rate.push_back(tmp * data.Ne);
      }
    }

    for(int i = (int)coalescent_rate.size()-1; i > 0; i--){
      if(coalescent_rate[i-1] == 0){
        if(coalescent_rate[i] > 0.0){
          coalescent_rate[i-1] = coalescent_rate[i];
        }else{
          coalescent_rate[i-1] = 1.0;
        }
      } 
    }

    //std::cerr << coalescent_rate[0] << std::endl;

  }
  is.close();


  if(1){
    if(options.count("mrate")){
      //multiply by mutation rate
      is.open(options["mrate"].as<std::string>());
      double mepoch, mrate;
      int e = 0;
      while(getline(is, line)){
        sscanf(line.c_str(), "%lf %lf", &mepoch, &mrate);
        assert(mepoch/data.Ne == epoch[e]);
        if(mrate > 0){     
          //double diff = (data.mu - mrate)/data.mu;
          //if(diff > 1) diff = 1;
          //if(diff < -1) diff = 1;
          //coalescent_rate[e] *= exp(log(10)*diff);
          coalescent_rate[e] *= data.mu/mrate;
        }
        e++;
      }
    }
  }


  ///////////////////////////////////////// TMRCA Inference /////////////////////////
  //Infer Branchlengths

  AncesTree anc;
  anc.Read(options["input"].as<std::string>() + ".anc");

  int num_trees = anc.seq.size();
  int progress_interval = (int)(num_trees/100.0) + 1;
  int count_trees = 0, progress = 0, progress_step = 1;
  if(num_trees < 100){
    progress_step = 100/num_trees;
  }

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  CorrTrees::iterator it_seq   = anc.seq.begin();
  //Infer branch lengths
  if(anc.sample_ages.size() == 0){
    EstimateBranchLengthsWithSampleAge bl(data);
    for(; it_seq != anc.seq.end(); it_seq++){
      if(count_trees % progress_interval == 0){
        progress += progress_step;
        ShowProgress(progress); 
      }
      count_trees++; 
      if(options.count("poplabels")){
        bl.MCMCCoalRatesForRelate(data, (*it_seq).tree, sample.group_of_haplotype, epoch, coal_rate_pair, rand()); //this is estimating times     
      }else{
        bl.MCMCVariablePopulationSizeForRelate(data, (*it_seq).tree, epoch, coalescent_rate, rand()); //this is estimating times
      }
    }
  }else{
    EstimateBranchLengthsWithSampleAge bl(data, anc.sample_ages);
    for(; it_seq != anc.seq.end(); it_seq++){
      if(count_trees % progress_interval == 0){
        progress += progress_step;
        ShowProgress(progress); 
      }
      count_trees++; 
      if(options.count("poplabels")){
        bl.MCMCCoalRatesForRelate(data, (*it_seq).tree, sample.group_of_haplotype, epoch, coal_rate_pair, rand()); //this is estimating times     
      }else{
        bl.MCMCVariablePopulationSizeForRelate(data, (*it_seq).tree, epoch, coalescent_rate, rand()); //this is estimating times
      }
    }
  }

  ShowProgress(100);
  std::cerr << std::endl;
  //Dump to file
  anc.Dump(options["output"].as<std::string>() + ".anc");

  ////////////////////////// Update mutation file

  CorrTrees::iterator it_anc = anc.seq.begin();
  std::vector<float> coordinates(2*data.N-1);
  int num_tree = mut.info[0].tree;
  int root = 2*data.N-2;

  (*it_anc).tree.GetCoordinates(coordinates);

  std::vector<SNPInfo>::iterator it_mut = mut.info.begin();
  for(; it_mut != mut.info.end(); it_mut++){
    //need the tree such that snp_of_next_tree > (*it_mut).snp_id
    //and snp_of_current_tree <= (*it_mut).snp_id
    if((*it_mut).tree > num_tree){
      while((*it_mut).tree > num_tree){
        it_anc++;
        if(it_anc == anc.seq.end()){
          it_anc--;
          break;
        }
        num_tree++;
      }
      (*it_anc).tree.GetCoordinates(coordinates);
    }
    if((*it_mut).tree != num_tree) std::cerr << (*it_mut).tree << " " << num_tree << std::endl;
    if((*it_mut).branch.size() == 1){
      int branch = *(*it_mut).branch.begin();
      if(branch != root){
        (*it_mut).age_begin = coordinates[branch];
        (*it_mut).age_end   = coordinates[(*(*it_anc).tree.nodes[branch].parent).label]; 
      }else{
        (*it_mut).age_begin = coordinates[branch];
        (*it_mut).age_end   = coordinates[branch];
      }
    }
  }
  mut.Dump(options["output"].as<std::string>() + ".mut"); 

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

////////////////////////////

int SampleBranchLengths(cxxopts::Options& options){

  int seed;
  if(!options.count("seed")){
    seed = std::time(0) + getpid();
  }else{
    seed = options["seed"].as<int>();
    srand(seed);
    std::string name = options["input"].as<std::string>();
    int tmp = 0;
    for(int i = 0; i < name.size(); i++){
      if(std::isdigit(name[i])) tmp += name[i]-48; 
    }
    for(int i = 0; i < tmp; i++){
      seed = rand();
    }
  }
  srand(seed);

  int Ne = 2e4;
  double mutation_rate = options["mutation_rate"].as<float>();
  std::string line;
  double tmp;

  //parse data
  int N;
  igzstream is_N(options["input"].as<std::string>() + ".anc");
  if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
  if(is_N.fail()){
    std::cerr << "Error while opening .anc file." << std::endl;
    exit(1);
  } 
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  //make this more efficient
  int L = 0;
  if(options.count("dist")){
    igzstream is_L(options["dist"].as<std::string>());
    if(is_L.fail()){
      std::cerr << "Error while opening .dist file." << std::endl;
      exit(1);
    } 
    while(std::getline(is_L, line)){
      ++L;
    }
    L--;
    is_L.close();
  }else{
    igzstream is_L(options["input"].as<std::string>() + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz");
    if(is_L.fail()){
      std::cerr << "Error while opening .mut file." << std::endl;
      exit(1);
    } 
    while(std::getline(is_L, line)){
      ++L;
    }
    L--;
    is_L.close();
  }

  Data data(N, L, Ne, mutation_rate);
  int root = 2*N-2;

  Mutations mut(data);
  mut.Read(options["input"].as<std::string>() + ".mut");
  data.dist.resize(L);
  std::vector<int> bp(L);
  if(options.count("dist")){
    igzstream is_dist(options["dist"].as<std::string>());
    if(is_dist.fail()){
      std::cerr << "Error while opening " << options["dist"].as<std::string>() << std::endl;
      exit(1);
    }
    getline(is_dist, line); 
    int dtmp, snp = 0;
    while(std::getline(is_dist, line)){
      sscanf(line.c_str(), "%d %d", &bp[snp], &data.dist[snp]);
      snp++;
    }
    is_dist.close();
  }else{
    std::vector<int>::iterator it_pos = data.dist.begin();
    std::vector<int>::iterator it_bp  = bp.begin();
    for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
      *it_pos = (*it_mut).dist;
      *it_bp  = (*it_mut).pos;
      it_bp++;
      it_pos++;
    }
  }


  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Sampling branch lengths for " << options["input"].as<std::string>() << " ..." << std::endl;

  Sample sample;
  if(options.count("poplabels")){
    sample.Read(options["poplabels"].as<std::string>());
  }
  std::vector<std::vector<std::vector<double>>> coal_rate_pair;

  // read epochs and population size 
  igzstream is(options["coal"].as<std::string>()); 
  if(is.fail()){
    is.open(options["coal"].as<std::string>() + ".gz");
    if(is.fail()){ 
      std::cerr << "Error while opening " << options["coal"].as<std::string>() << "(.gz)." << std::endl;
      exit(1);
    }
  } 

  std::vector<double> epoch, coalescent_rate;
  getline(is, line);
  std::vector<std::string> groups;
  std::istringstream is_group(line);
  std::string name;
  while(is_group >> name){
    //is_group >> name;
    groups.push_back(name);
  }
  getline(is, line);
  std::istringstream is_epoch(line);
  while(is_epoch >> tmp){
    //is_epoch >> tmp;
    epoch.push_back(tmp/data.Ne);
  }

  if(options.count("poplabels")){

    if(groups.size() != sample.groups.size()){
      std::cerr << "Coal file doesn't contain all groups vs all groups rates" << std::endl;
      exit(1);
    }

    //convert from coal file to poplabels file
    std::vector<int> convert(groups.size(), -1);
    for(int g = 0; g < groups.size(); g++){
      for(int g2 = 0; g2 < sample.groups.size(); g2++){
        if(groups[g] == sample.groups[g2]){
          convert[g] = g2;
          break;
        }
      }
      if(convert[g] == -1){
        std::cerr << "Groups in coal file don't match poplabels file" << std::endl;
        exit(1);
      }
    }

    coal_rate_pair.resize(epoch.size());
    for(int e = 0; e < epoch.size(); e++){
      coal_rate_pair[e].resize(sample.groups.size());
      for(int g = 0; g < sample.groups.size(); g++){
        coal_rate_pair[e][g].resize(sample.groups.size());
      }
    }

    for(int g1 = 0; g1 < sample.groups.size(); g1++){
      for(int g2 = 0; g2 < sample.groups.size(); g2++){
        if(!getline(is,line)){
          std::cerr << "Coal file doesn't contain all groups vs all groups rates" << std::endl;
          exit(1);
        }

        std::istringstream is_pop_size(line);
        is_pop_size >> tmp;
        assert(tmp == g1);
        is_pop_size >> tmp;
        assert(tmp == g2);
        int ep = 0;
        while(is_pop_size){
          is_pop_size >> tmp;
          //if(ep == 0) std::cerr << convert[g1] << " " << convert[g2] << " " << tmp << " " << data.Ne << " " << tmp * data.Ne << std::endl;
          if(tmp == 0.0){
            coal_rate_pair[ep][convert[g1]][convert[g2]] = 5e-10 * data.Ne;
          }else{
            coal_rate_pair[ep][convert[g1]][convert[g2]] = tmp * data.Ne;
          }
          ep++;
          if(ep == epoch.size()) break;
        }

        if(0){
          for(int g1 = 0; g1 < sample.groups.size(); g1++){
            for(int g2 = 0; g2 < sample.groups.size(); g2++){
              for(int i = epoch.size()-1; i > 0; i--){
                if(coal_rate_pair[i-1][g1][g2] == 0){
                  if(coal_rate_pair[i][g1][g2] > 0.0){
                    coal_rate_pair[i-1][g1][g2] = coal_rate_pair[i][g1][g2];
                  }else{
                    coal_rate_pair[i-1][g1][g2] = 1.0;
                  }
                } 
              }
            }
          }
        }

      }
    }

  }else{

    getline(is, line);
    std::istringstream is_pop_size(line);
    is_pop_size >> tmp >> tmp;
    while(is_pop_size){
      is_pop_size >> tmp;
      //tmp = 1.0/data.Ne; 
      if(tmp == 0.0 && coalescent_rate.size() > 0){
        if(*std::prev(coalescent_rate.end(),1) > 0.0){
          coalescent_rate.push_back(*std::prev(coalescent_rate.end(),1));
        }
        //coalescent_rate.push_back(1);
      }else{
        coalescent_rate.push_back(tmp * data.Ne);
      }
    }

    for(int i = (int)coalescent_rate.size()-1; i > 0; i--){
      if(coalescent_rate[i-1] == 0){
        if(coalescent_rate[i] > 0.0){
          coalescent_rate[i-1] = coalescent_rate[i];
        }else{
          coalescent_rate[i-1] = 1.0;
        }
      } 
    }

    //std::cerr << coalescent_rate[0] << std::endl;

  }
  is.close();

  if(1){
    if(options.count("mrate")){
      //multiply by mutation rate
      is.open(options["mrate"].as<std::string>());
      double mepoch, mrate;
      int e = 0;
      while(getline(is, line)){
        sscanf(line.c_str(), "%lf %lf", &mepoch, &mrate);
        assert(mepoch/data.Ne == epoch[e]);
        if(mrate > 0){     
          //double diff = (data.mu - mrate)/data.mu;
          //if(diff > 1) diff = 1;
          //if(diff < -1) diff = 1;
          //coalescent_rate[e] *= exp(log(10)*diff);
          coalescent_rate[e] *= data.mu/mrate;
        }
        e++;
      }
    }
  }

  ///////////////////////////////////////// TMRCA Inference /////////////////////////
  //Infer Branchlengths

  AncesTree anc;
  anc.Read(options["input"].as<std::string>() + ".anc");

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  int num_trees = anc.seq.size();
  int progress_interval = (int)(num_trees/100.0) + 1;
  int count_trees = 0, progress = 0, progress_step = 1;
  if(num_trees < 100){
    progress_step = 100/num_trees;
  }

  //need to make these three variables to arguments
  int num_proposals = 1000*std::max(data.N/10.0, 10.0);
  if(options.count("num_proposals")){
    num_proposals = options["num_proposals"].as<int>();
  }
  int num_samples   = options["num_samples"].as<int>();
  std::string chrid = "chr";

  if(num_samples < 1){
    std::cerr << "Error: num_samples value < 1" << std::endl;
    exit(1);
  }
  if(num_proposals < 0){
    std::cerr << "Error: num_proposals value < 0" << std::endl;
    exit(1);
  }

  //////////output files

  std::string format = "a";
  if(options.count("format")){
    format = options["format"].as<std::string>();
    if(format != "a" && format != "n"){
      std::cerr << "Error: output format doesn't exist." << std::endl;
      exit(1);
    }
  }

  std::vector<Node>::iterator n_it;
  std::vector<std::vector<float>> branch_lengths(2*data.N-1);
  for(std::vector<std::vector<float>>::iterator it_branch_lengths = branch_lengths.begin(); it_branch_lengths != branch_lengths.end(); it_branch_lengths++){
    (*it_branch_lengths).resize(num_samples);
    std::fill((*it_branch_lengths).begin(), (*it_branch_lengths).end(), 0.0);
  }

  std::ofstream os, os_sites;
  std::string filename;
  if(format == "n"){
    filename = options["output"].as<std::string>() + ".newick";
    os.open(filename);
    os << "#chrom\tchromStart\tchromEnd\tMCMC_sample\ttree" << std::endl;
    os_sites.open(options["output"].as<std::string>() + ".sites");

    os_sites << "NAMES\t";
    for(int i = 0; i < data.N; i++){
      os_sites << i << "\t";
    }
    os_sites << "\n";
    if(mut.info.size() > 0){
      os_sites << "REGION\t" << chrid << "\t" << mut.info[0].pos << "\t" << mut.info[mut.info.size()-1].pos + 1 << "\n";
    }
  }else{
    filename = options["output"].as<std::string>() + ".anc";
    os.open(filename);
    os << "NUM_HAPLOTYPES " << ((*anc.seq.begin()).tree.nodes.size() + 1)/2 << " ";
    for(std::vector<double>::iterator it_sample_ages = anc.sample_ages.begin(); it_sample_ages != anc.sample_ages.end(); it_sample_ages++){
      os << *it_sample_ages << " ";
    }
    os << "\n";
    os << "NUM_TREES " << anc.seq.size() << "\n";
    if(num_samples > 1) os << "NUM_SAMPLES_PER_TREE " << num_samples << "\n";
  }

  std::vector<Leaves> leaves;
  CorrTrees::iterator it_seq   = anc.seq.begin();
  std::vector<SNPInfo>::iterator it_mut = mut.info.begin();

  //Infer branch lengths

  if(anc.sample_ages.size() == 0){
    EstimateBranchLengthsWithSampleAge bl(data);
    for(; it_seq != anc.seq.end(); it_seq++){

      if(count_trees % progress_interval == 0){
        progress += progress_step;
        ShowProgress(progress); 
      }

      for(std::vector<Node>::iterator it_node = (*it_seq).tree.nodes.begin(); it_node != (*it_seq).tree.nodes.end(); it_node++){
        (*it_node).branch_length /= (double) data.Ne;
      }

      int count = 0;
      if(count < num_samples){

        if(options.count("poplabels")){
          bl.MCMCCoalRatesSample(data, (*it_seq).tree, sample.group_of_haplotype, epoch, coal_rate_pair, num_proposals, 1, rand()); //this is estimating times     
        }else{
          bl.MCMCVariablePopulationSizeSample(data, (*it_seq).tree, epoch, coalescent_rate, num_proposals, 1, rand()); //this is estimating times
        }

        if(format == "n"){
          if(it_seq != std::prev(anc.seq.end(),1)){
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << bp[(*std::next(it_seq,1)).pos] << "\t" << count << "\t";
          }else{
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << (*std::prev(mut.info.end(),1)).pos + 1 << "\t" << count << "\t";
          } 
          (*it_seq).tree.WriteNewick(os, (double) data.Ne);
        }else{
          for(n_it = (*it_seq).tree.nodes.begin(); n_it != (*it_seq).tree.nodes.end(); n_it++){
            branch_lengths[(*n_it).label][count] = (*n_it).branch_length;
          }	
        }
      }
      count++;
      for(;count < num_samples; count++){
        if(options.count("poplabels")){
          bl.MCMCCoalRatesSample(data, (*it_seq).tree, sample.group_of_haplotype, epoch, coal_rate_pair, num_proposals, 0, rand()); //this is estimating times     
        }else{
          bl.MCMCVariablePopulationSizeSample(data, (*it_seq).tree, epoch, coalescent_rate, num_proposals, 0, rand()); //this is estimating times
        }

        if(format == "n"){
          if(it_seq != std::prev(anc.seq.end(),1)){
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << bp[(*std::next(it_seq,1)).pos] << "\t" << count << "\t";
          }else{
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << (*std::prev(mut.info.end(),1)).pos + 1 << "\t" << count << "\t";
          } 
          (*it_seq).tree.WriteNewick(os, (double) data.Ne);
        }else{
          for(n_it = (*it_seq).tree.nodes.begin(); n_it != (*it_seq).tree.nodes.end(); n_it++){
            branch_lengths[(*n_it).label][count] = (*n_it).branch_length;
          }	
        }
      }

      if(format == "n"){
        (*it_seq).tree.FindAllLeaves(leaves);

        if(it_mut != mut.info.end()){
          while((*it_mut).tree == count_trees){

            if((*it_mut).branch.size() == 1 && (*it_mut).flipped == false){

              //.sites file
              //get ancestral and derived allele
              char ancestral = (*it_mut).mutation_type[0];
              char derived   = (*it_mut).mutation_type[2]; 
              //get list of descendants and output string

              std::sort(leaves[*(*it_mut).branch.begin()].member.begin(), leaves[*(*it_mut).branch.begin()].member.end());

              std::vector<int>::iterator it_member = leaves[*(*it_mut).branch.begin()].member.begin();
              os_sites << (*it_mut).pos << "\t";
              for(int node = 0; node < data.N; node++){
                if(it_member == leaves[*(*it_mut).branch.begin()].member.end()){
                  os_sites << ancestral;           
                }else if(node == *it_member){
                  os_sites << derived;
                  it_member++;
                }else{
                  os_sites << ancestral;
                }
              }
              os_sites << "\n";

            }

            it_mut++;
            if(it_mut == mut.info.end()) break;

          }
        }
      }else{

        int parent;
        n_it = (*it_seq).tree.nodes.begin();
        os << (*it_seq).pos << ": ";
        for(; n_it != (*it_seq).tree.nodes.end(); n_it++){
          if((*n_it).parent == NULL){
            parent = -1;
          }else{
            parent = (*(*n_it).parent).label;
          }
          os << parent << ":(";
          for(std::vector<float>::iterator it_branch_length = branch_lengths[(*n_it).label].begin(); it_branch_length != branch_lengths[(*n_it).label].end(); it_branch_length++){
            os << std::fixed << std::setprecision(5) << (*it_branch_length * data.Ne) << " ";
          }
          os << std::setprecision(2) << (*n_it).num_events << " " << (*n_it).SNP_begin << " " << (*n_it).SNP_end << ") ";
        }

        os << "\n";

      }

      count_trees++; 

    }
  }else{
    //has sample ages
    EstimateBranchLengthsWithSampleAge bl(data, anc.sample_ages);
    for(; it_seq != anc.seq.end(); it_seq++){

      if(count_trees % progress_interval == 0){
        progress += progress_step;
        ShowProgress(progress); 
      }

      for(std::vector<Node>::iterator it_node = (*it_seq).tree.nodes.begin(); it_node != (*it_seq).tree.nodes.end(); it_node++){
        (*it_node).branch_length /= (double) data.Ne;
      }

      int count = 0;
      if(count < num_samples){
        if(options.count("poplabels")){
          bl.MCMCCoalRatesSample(data, (*it_seq).tree, sample.group_of_haplotype, epoch, coal_rate_pair, num_proposals, 1, rand()); //this is estimating times     
        }else{
          bl.MCMCVariablePopulationSizeSample(data, (*it_seq).tree, epoch, coalescent_rate, num_proposals, 1, rand()); //this is estimating times
        }

        if(format == "n"){
          if(it_seq != std::prev(anc.seq.end(),1)){
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << bp[(*std::next(it_seq,1)).pos] << "\t" << count << "\t";
          }else{
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << (*std::prev(mut.info.end(),1)).pos + 1 << "\t" << count << "\t";
          } 
          (*it_seq).tree.WriteNewick(os, (double) data.Ne);
        }else{
          for(n_it = (*it_seq).tree.nodes.begin(); n_it != (*it_seq).tree.nodes.end(); n_it++){
            branch_lengths[(*n_it).label][count] = (*n_it).branch_length;
          }	
        }
      }
      count++;
      for(;count < num_samples; count++){
        if(options.count("poplabels")){
          bl.MCMCCoalRatesSample(data, (*it_seq).tree, sample.group_of_haplotype, epoch, coal_rate_pair, num_proposals, 0, rand()); //this is estimating times     
        }else{
          bl.MCMCVariablePopulationSizeSample(data, (*it_seq).tree, epoch, coalescent_rate, num_proposals, 0, rand()); //this is estimating times
        }

        if(format == "n"){
          if(it_seq != std::prev(anc.seq.end(),1)){
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << bp[(*std::next(it_seq,1)).pos] << "\t" << count << "\t";
          }else{
            os << chrid << "\t" << bp[(*it_seq).pos] << "\t" << (*std::prev(mut.info.end(),1)).pos + 1 << "\t" << count << "\t";
          } 
          (*it_seq).tree.WriteNewick(os, (double) data.Ne);
        }else{
          for(n_it = (*it_seq).tree.nodes.begin(); n_it != (*it_seq).tree.nodes.end(); n_it++){
            branch_lengths[(*n_it).label][count] = (*n_it).branch_length;
          }	
        }
      }

      if(format == "n"){
        (*it_seq).tree.FindAllLeaves(leaves);

        if(it_mut != mut.info.end()){
          while((*it_mut).tree == count_trees){

            if((*it_mut).branch.size() == 1 && (*it_mut).flipped == false){

              //.sites file
              //get ancestral and derived allele
              char ancestral = (*it_mut).mutation_type[0];
              char derived   = (*it_mut).mutation_type[2]; 
              //get list of descendants and output string

              std::sort(leaves[*(*it_mut).branch.begin()].member.begin(), leaves[*(*it_mut).branch.begin()].member.end());

              std::vector<int>::iterator it_member = leaves[*(*it_mut).branch.begin()].member.begin();
              os_sites << (*it_mut).pos << "\t";
              for(int node = 0; node < data.N; node++){
                if(it_member == leaves[*(*it_mut).branch.begin()].member.end()){
                  os_sites << ancestral;           
                }else if(node == *it_member){
                  os_sites << derived;
                  it_member++;
                }else{
                  os_sites << ancestral;
                }
              }
              os_sites << "\n";

            }

            it_mut++;
            if(it_mut == mut.info.end()) break;

          }
        }
      }else{

        int parent;
        n_it = (*it_seq).tree.nodes.begin();
        os << (*it_seq).pos << ": ";
        for(; n_it != (*it_seq).tree.nodes.end(); n_it++){
          if((*n_it).parent == NULL){
            parent = -1;
          }else{
            parent = (*(*n_it).parent).label;
          }
          os << parent << ":(";
          for(std::vector<float>::iterator it_branch_length = branch_lengths[(*n_it).label].begin(); it_branch_length != branch_lengths[(*n_it).label].end(); it_branch_length++){
            os << std::fixed << std::setprecision(5) << (*it_branch_length * data.Ne) << " ";
          }
          os << std::setprecision(2) << (*n_it).num_events << " " << (*n_it).SNP_begin << " " << (*n_it).SNP_end << ") ";
        }

        os << "\n";

      }

      count_trees++; 

    }
  }
  os.close();

  ShowProgress(100);
  std::cerr << std::endl;

  os_sites.close(); 

  if(format == "a"){

    for(std::vector<double>::iterator it_sample_ages = anc.sample_ages.begin(); it_sample_ages != anc.sample_ages.end(); it_sample_ages++){
      *it_sample_ages /= data.Ne;
    }
    Mutations mut(data);
    mut.Read(options["input"].as<std::string>() + ".mut");

    CorrTrees::iterator it_anc = anc.seq.begin();
    std::vector<float> coordinates(2*data.N-1);
    int num_tree = mut.info[0].tree;

    (*it_anc).tree.GetCoordinates(coordinates);

    std::vector<SNPInfo>::iterator it_mut = mut.info.begin();
    for(; it_mut != mut.info.end(); it_mut++){
      //need the tree such that snp_of_next_tree > (*it_mut).snp_id
      //and snp_of_current_tree <= (*it_mut).snp_id
      if((*it_mut).tree > num_tree){
        while((*it_mut).tree > num_tree){
          it_anc++;
          if(it_anc == anc.seq.end()){
            it_anc--;
            break;
          }
          num_tree++;
        }
        (*it_anc).tree.GetCoordinates(coordinates);
      }
      if((*it_mut).tree != num_tree){
        std::cerr << (*it_mut).tree << " " << num_tree << std::endl;
        exit(1);
      }
      if((*it_mut).branch.size() == 1){
        int branch = *(*it_mut).branch.begin();
        if(branch != root){
          (*it_mut).age_begin = data.Ne*coordinates[branch];
          (*it_mut).age_end   = data.Ne*coordinates[(*(*it_anc).tree.nodes[branch].parent).label]; 
        }else{
          (*it_mut).age_begin = data.Ne*coordinates[branch];
          (*it_mut).age_end   = data.Ne*coordinates[branch]; 
        }
      }
    }
    mut.Dump(options["output"].as<std::string>() + ".mut"); 
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


/////////////////////////////////
float           
GetCoords(int node, Tree& tree, int branch, float Ne, char m, std::vector<float>::iterator& it_dertimes, std::vector<float>::iterator& it_anctimes){

  float coordinate = 0.0;
  if(tree.nodes[node].child_left != NULL){

    int child_left  = (*tree.nodes[node].child_left).label;
    int child_right = (*tree.nodes[node].child_right).label;

    if(child_left == branch || m == 'd'){
      coordinate = GetCoords(child_left, tree, branch, Ne, 'd', it_dertimes, it_anctimes);
    }else{
      coordinate = GetCoords(child_left, tree, branch, Ne, 'a', it_dertimes, it_anctimes);
    }

    if(child_right == branch || m == 'd'){
      coordinate = GetCoords(child_right, tree, branch, Ne, 'd', it_dertimes, it_anctimes);
    }else{
      coordinate = GetCoords(child_right, tree, branch, Ne, 'a', it_dertimes, it_anctimes);
    }
    coordinate += tree.nodes[child_right].branch_length;

    if(child_left != branch && child_right != branch){
      if(m == 'a'){
        *it_anctimes = Ne*coordinate;
        it_anctimes++;
      }else{
        *it_dertimes = Ne*coordinate;
        it_dertimes++;
      }
    }

  }else{

    if(tree.sample_ages != NULL){
      assert((*tree.sample_ages).size() > 0);
      coordinate = (*tree.sample_ages)[node];
    }

  }

  return coordinate;

}

int SampleBranchLengthsBinary(cxxopts::Options& options){

  int seed;
  if(!options.count("seed")){
    seed = std::time(0) + getpid();
  }else{
    seed = options["seed"].as<int>();
    srand(seed);
    std::string name = options["input"].as<std::string>();
    int tmp = 0;
    for(int i = 0; i < name.size(); i++){
      if(std::isdigit(name[i])) tmp += name[i] - 48; 
    }
    for(int i = 0; i < tmp; i++){
      seed = rand();
    }
  }
  srand(seed);

  int Ne = 3e4;
  double mutation_rate = options["mutation_rate"].as<float>();
  std::string line;
  double tmp;

  //parse data
  int N;
  igzstream is_N(options["input"].as<std::string>() + ".anc");
  if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
  if(is_N.fail()){
    std::cerr << "Error while opening .anc file." << std::endl;
    exit(1);
  } 
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  //make this more efficient
  int L = 0;
  if(options.count("dist")){
    igzstream is_L(options["dist"].as<std::string>());
    if(is_L.fail()){
      std::cerr << "Error while opening .dist file." << std::endl;
      exit(1);
    } 
    while(std::getline(is_L, line)){
      ++L;
    }
    L--;
    is_L.close();
  }else{
    igzstream is_L(options["input"].as<std::string>() + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz");
    if(is_L.fail()){
      std::cerr << "Error while opening .mut file." << std::endl;
      exit(1);
    } 
    while(std::getline(is_L, line)){
      ++L;
    }
    L--;
    is_L.close();
  }

  Data data(N, L, Ne, mutation_rate);

  Mutations mut(data);
  mut.Read(options["input"].as<std::string>() + ".mut");
  int num_mapping_SNPs = 0;
  for(Muts::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
    if((*it_mut).branch.size() <= 1) num_mapping_SNPs++;
  }
  if(num_mapping_SNPs == 0){
    std::cerr << "Error: No SNPs are mapping to tree" << std::endl;
    exit(1);
  }

  data.dist.resize(L);
  std::vector<int> bp(L);
  if(options.count("dist")){
    igzstream is_dist(options["dist"].as<std::string>());
    if(is_dist.fail()){
      std::cerr << "Error while opening " << options["dist"].as<std::string>() << std::endl;
      exit(1);
    }
    getline(is_dist, line); 
    int dtmp, snp = 0;
    while(std::getline(is_dist, line)){
      sscanf(line.c_str(), "%d %d", &bp[snp], &data.dist[snp]);
      snp++;
    }
    is_dist.close();
  }else{
    std::vector<int>::iterator it_pos = data.dist.begin();
    std::vector<int>::iterator it_bp  = bp.begin();
    for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
      *it_pos = (*it_mut).dist;
      *it_bp  = (*it_mut).pos;
      it_bp++;
      it_pos++;
    }
  }

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Sampling branch lengths for " << options["input"].as<std::string>() << " ..." << std::endl;

  // read epochs and population size 
  igzstream is(options["coal"].as<std::string>()); 
  if(is.fail()){
    is.open(options["coal"].as<std::string>() + ".gz");
    if(is.fail()){ 
      std::cerr << "Error while opening " << options["coal"].as<std::string>() << "(.gz)." << std::endl;
      exit(1);
    }
  } 

  std::vector<double> epoch, coalescent_rate;
  getline(is, line);
  getline(is, line);
  std::istringstream is_epoch(line);
  while(is_epoch){
    is_epoch >> tmp;
    epoch.push_back(tmp/data.Ne);
  }
  getline(is, line);
  is.close();


	getline(is, line);
	std::istringstream is_pop_size(line);
	is_pop_size >> tmp >> tmp;
	while(is_pop_size){
		is_pop_size >> tmp;
		//tmp = 1.0/data.Ne; 
		if(tmp == 0.0 && coalescent_rate.size() > 0){
			if(*std::prev(coalescent_rate.end(),1) > 0.0){
				coalescent_rate.push_back(*std::prev(coalescent_rate.end(),1));
			}
			//coalescent_rate.push_back(1);
		}else{
			coalescent_rate.push_back(tmp * data.Ne);
		}
	}

	for(int i = (int)coalescent_rate.size()-1; i > 0; i--){
		if(coalescent_rate[i-1] == 0){
			if(coalescent_rate[i] > 0.0){
				coalescent_rate[i-1] = coalescent_rate[i];
			}else{
				coalescent_rate[i-1] = 1.0;
			}
		} 
	}

  if(0){
    if(options.count("mrate")){
      //multiply by mutation rate
      is.open(options["mrate"].as<std::string>());
      double mepoch, mrate;
      int e = 0;
      while(getline(is, line)){
        sscanf(line.c_str(), "%lf %lf", &mepoch, &mrate);
        assert(mepoch/data.Ne == epoch[e]);
        if(mrate > 0){     
          //double diff = (data.mu - mrate)/data.mu;
          //if(diff > 1) diff = 1;
          //if(diff < -1) diff = 1;
          //coalescent_rate[e] *= exp(log(10)*diff);
          coalescent_rate[e] *= data.mu/mrate;
        }
        e++;
      }
    }
  }

  ///////////////////////////////////////// TMRCA Inference /////////////////////////
  //Infer Branchlengths

  MarginalTree mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
  Muts::iterator it_mut; //iterator for mut file
  float num_bases_tree_persists = 0.0;
  int root = 2*data.N - 2;

  AncMutIterators ancmut(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  int num_trees = ancmut.NumTrees();

  //////////////////////////////////////////// Read Tree ///////////////////////////////////

  int progress_interval = (int)(num_mapping_SNPs/100.0) + 1;
  int count_trees = 0, count_snps = 0, progress = 0, progress_step = 1;
  if(num_mapping_SNPs < 100){
    progress_step = 100/num_mapping_SNPs;
  }

  //need to make these three variables to arguments
  int num_proposals = 1000*std::max(data.N/10.0, 10.0);
  if(options.count("num_proposals")){
    num_proposals = options["num_proposals"].as<int>();
  }
  int num_samples   = options["num_samples"].as<int>();
  std::string chrid = "chr";

  //prepare files for output
  std::string filename = options["output"].as<std::string>() + ".timeb";
  FILE* fp = fopen(filename.c_str(), "wb");
  //write number of trees and number of proposals per SNP
  //std::cerr << num_mapping_SNPs << " " << sizeof(int) << std::endl;
  //std::cerr << num_samples << " " << sizeof(int) << std::endl;
  fwrite(&num_mapping_SNPs, sizeof(int), 1, fp);
  fwrite(&num_samples, sizeof(int), 1, fp);

  char anc_allele, der_allele;
  if(ancmut.sample_ages.size() == 0){

    //Infer branch lengths
    EstimateBranchLengthsWithSampleAge bl(data);
    std::vector<Leaves> leaves;
    bool first_snp = false;
    std::vector<Tree> sampled_trees(num_samples);
    //iterate through whole file
    while(num_bases_tree_persists >= 0.0){

      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      first_snp = true;

      if(it_mut != ancmut.mut_end()){
        while((*it_mut).tree == count_trees){

          if((*it_mut).branch.size() <= 1){

            if(count_snps % progress_interval == 0){
              progress += progress_step;
              ShowProgress(progress); 
            }

            if(first_snp){
              for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
                (*it_node).branch_length /= (double) data.Ne;
              }
              mtr.tree.FindAllLeaves(leaves);

              int i = 0;
              if(i < num_samples){
                sampled_trees[i] = mtr.tree;
                bl.MCMCVariablePopulationSizeSample(data, sampled_trees[i], epoch, coalescent_rate, num_proposals, 1, rand()); //this is estimating times
              }
              i++;
              for(; i < num_samples; i++){
                sampled_trees[i] = sampled_trees[i-1];
                bl.MCMCVariablePopulationSizeSample(data, sampled_trees[i], epoch, coalescent_rate, num_proposals, 0, rand()); //this is estimating times
              }
              first_snp = false;
            }

            //store matrix with
            //num_snps x anc_times x num_samples
            //num_snps x der_times x num_samples

            int branch, DAF;
            std::vector<float> anctimes, dertimes;
            if((*it_mut).branch.size() == 1){
              branch = *(*it_mut).branch.begin();
              DAF    = leaves[branch].num_leaves;
              anctimes.resize(num_samples*std::max(0,(data.N-DAF-1)));
              dertimes.resize(num_samples*std::max(0,(DAF-1)));
              std::vector<float>::iterator it_anctimes = anctimes.begin(), it_dertimes = dertimes.begin();

              int count = 0;
              if(count < num_samples){
                //store anc and dertimes
                std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
                if(branch != root){
                  GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'a', it_dertimes, it_anctimes);
                }else{
                  GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'd', it_dertimes, it_anctimes);
                }
                assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
                assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
                std::sort(it_anctimes_s, it_anctimes);
                std::sort(it_dertimes_s, it_dertimes);
              }
              count++;
              for(;count < num_samples; count++){
                //store anc and dertimes
                std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
                if(branch != root){
                  GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'a', it_dertimes, it_anctimes);
                }else{
                  GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'd', it_dertimes, it_anctimes);
                }
                assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
                assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
                std::sort(it_anctimes_s, it_anctimes);
                std::sort(it_dertimes_s, it_dertimes);
              }
            }else{
              DAF    = 0;
              anctimes.resize(num_samples*std::max(0,(data.N-1)));
              dertimes.resize(0);
              std::vector<float>::iterator it_anctimes = anctimes.begin(), it_dertimes = dertimes.begin();
              int count = 0;
              if(count < num_samples){
                //store anc and dertimes
                std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
                GetCoords(2*data.N-2, sampled_trees[count], root, data.Ne, 'a', it_dertimes, it_anctimes);
                assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
                assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
                std::sort(it_anctimes_s, it_anctimes);
                std::sort(it_dertimes_s, it_dertimes);
              }
              count++;
              for(;count < num_samples; count++){
                //store anc and dertimes
                std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
                GetCoords(2*data.N-2, sampled_trees[count], root, data.Ne, 'a', it_dertimes, it_anctimes);
                assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
                assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
                std::sort(it_anctimes_s, it_anctimes);
                std::sort(it_dertimes_s, it_dertimes);
              }

            }

            //WriteBinary(anctimes, dertimes, fp);
            //BP, DAF, N, 
            //dump
            int msize = (*it_mut).mutation_type.size();
            if(msize >= 1){
              anc_allele = (*it_mut).mutation_type[0];
              der_allele = 'N';
              int i = 1;
              while((*it_mut).mutation_type[i] != '/' && i < msize){
                i++;
                if(i == msize) break;
              }
              i++;
              if(i < msize) der_allele = (*it_mut).mutation_type[i];
            }

            int BP = (*it_mut).pos;
            fwrite(&BP, sizeof(int), 1, fp);
            fwrite(&anc_allele, sizeof(char), 1, fp);
            fwrite(&der_allele, sizeof(char), 1, fp);
            fwrite(&DAF, sizeof(int), 1, fp);
            fwrite(&data.N, sizeof(int), 1, fp);
            fwrite(&anctimes[0], sizeof(float), num_samples*std::max(0,(data.N-DAF-1)), fp);
            fwrite(&dertimes[0], sizeof(float), num_samples*std::max(0,(DAF-1)), fp);

            count_snps++;
          }

          it_mut++;
          if(it_mut == ancmut.mut_end()) break;
        }
      }

      if(it_mut == ancmut.mut_end()) break;
      count_trees++; 

    }
  }else{

    //Infer branch lengths
    EstimateBranchLengthsWithSampleAge bl(data, ancmut.sample_ages);
    std::vector<Leaves> leaves;
    bool first_snp = false;
    std::vector<Tree> sampled_trees(num_samples);
    //iterate through whole file
    while(num_bases_tree_persists >= 0.0){

      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      first_snp = true;

      if(it_mut != ancmut.mut_end()){
        while((*it_mut).tree == count_trees){

          if((*it_mut).branch.size() <= 1){

            if(count_snps % progress_interval == 0){
              progress += progress_step;
              ShowProgress(progress); 
            }

            if(first_snp){
              for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != mtr.tree.nodes.end(); it_node++){
                (*it_node).branch_length /= (double) data.Ne;
              }
              mtr.tree.FindAllLeaves(leaves);

              int i = 0;
              if(i < num_samples){
								//bl.MCMCVariablePopulationSizeSample(data, sampled_trees[i], epoch, coalescent_rate, num_proposals, 1, rand()); //this is estimating times
								bl.MCMCVariablePopulationSizeSample(data, mtr.tree, epoch, coalescent_rate, num_proposals, 1, rand()); //this is estimating times
								sampled_trees[i] = mtr.tree;
              }
              i++;
              for(; i < num_samples; i++){
								//bl.MCMCVariablePopulationSizeSample(data, sampled_trees[i], epoch, coalescent_rate, num_proposals, 0, rand()); //this is estimating times
								bl.MCMCVariablePopulationSizeSample(data, mtr.tree, epoch, coalescent_rate, num_proposals, 0, rand());
								sampled_trees[i] = mtr.tree;
              }
              first_snp = false;
            }

						//store matrix with
						//num_snps x anc_times x num_samples
						//num_snps x der_times x num_samples
						int branch, DAF;
						std::vector<float> anctimes, dertimes;
						if((*it_mut).branch.size() == 1){
							branch = *(*it_mut).branch.begin();
							DAF    = leaves[branch].num_leaves;
							anctimes.resize(num_samples*std::max(0,(data.N-DAF-1)));
							dertimes.resize(num_samples*std::max(0,(DAF-1)));
							std::vector<float>::iterator it_anctimes = anctimes.begin(), it_dertimes = dertimes.begin();

							int count = 0;
							if(count < num_samples){
								//store anc and dertimes
								std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
								if(branch != root){
									GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'a', it_dertimes, it_anctimes);
								}else{
									GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'd', it_dertimes, it_anctimes);
								}
								assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
								assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
								std::sort(it_anctimes_s, it_anctimes);
								std::sort(it_dertimes_s, it_dertimes);
							}
							count++;
							for(;count < num_samples; count++){
								//store anc and dertimes
								std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
								if(branch != root){
									GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'a', it_dertimes, it_anctimes);
								}else{
									GetCoords(2*data.N-2, sampled_trees[count], branch, data.Ne, 'd', it_dertimes, it_anctimes);
								}
								assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
								assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
								std::sort(it_anctimes_s, it_anctimes);
								std::sort(it_dertimes_s, it_dertimes);
							}
						}else{
							DAF    = 0;
							anctimes.resize(num_samples*std::max(0,(data.N-1)));
							dertimes.resize(0);
							std::vector<float>::iterator it_anctimes = anctimes.begin(), it_dertimes = dertimes.begin();
							int count = 0;
							if(count < num_samples){
								//store anc and dertimes
								std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
								GetCoords(2*data.N-2, sampled_trees[count], root, data.Ne, 'a', it_dertimes, it_anctimes);
								assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
								assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
								std::sort(it_anctimes_s, it_anctimes);
								std::sort(it_dertimes_s, it_dertimes);
							}
							count++;
							for(;count < num_samples; count++){
								//store anc and dertimes
								std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
								GetCoords(2*data.N-2, sampled_trees[count], root, data.Ne, 'a', it_dertimes, it_anctimes);
								assert(std::next(it_anctimes_s, std::max(0,data.N-DAF-1)) == it_anctimes);
								assert(std::next(it_dertimes_s, std::max(0,DAF-1)) == it_dertimes);
								std::sort(it_anctimes_s, it_anctimes);
								std::sort(it_dertimes_s, it_dertimes);
							}

						}

						//WriteBinary(anctimes, dertimes, fp);
						//BP, DAF, N, 
						//dump
						int msize = (*it_mut).mutation_type.size();
						if(msize >= 1){
							anc_allele = (*it_mut).mutation_type[0];
							der_allele = 'N';
							int i = 1;
							while((*it_mut).mutation_type[i] != '/' && i < msize){
								i++;
								if(i == msize) break;
							}
							i++;
							if(i < msize) der_allele = (*it_mut).mutation_type[i];
						}

						int BP = (*it_mut).pos;
						fwrite(&BP, sizeof(int), 1, fp);
						fwrite(&anc_allele, sizeof(char), 1, fp);
						fwrite(&der_allele, sizeof(char), 1, fp);
						fwrite(&DAF, sizeof(int), 1, fp);
						fwrite(&data.N, sizeof(int), 1, fp);
						fwrite(&anctimes[0], sizeof(float), num_samples*std::max(0,(data.N-DAF-1)), fp);
						fwrite(&dertimes[0], sizeof(float), num_samples*std::max(0,(DAF-1)), fp);

						count_snps++;
					}

					it_mut++;
					if(it_mut == ancmut.mut_end()) break;
				}
			}

			if(it_mut == ancmut.mut_end()) break;
			count_trees++; 

		}


	}


	ShowProgress(100);
	std::cerr << std::endl;

	fclose(fp);

	//Resource Usage
	/////////////////////////////////////////////

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
