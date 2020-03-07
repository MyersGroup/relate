#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#include "gzstream.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

//////// functions for estimating pairwise coalescence rate /////////

float
GetCoalescentRate(Node n, float factor, std::vector<float>& epoch, std::vector<CollapsedMatrix<float>>& coalescence_rate_data, std::vector<int>& leaves){

  //go through tree and fill in the matrix
  if(n.child_left != NULL){

    std::vector<int> leaves_child_left, leaves_child_right;

    Node child_left  = *n.child_left;
    Node child_right = *n.child_right;

    int num_children_left, num_children_right;

    float coalescent_time = GetCoalescentRate(child_left, factor, epoch, coalescence_rate_data, leaves_child_left) + child_left.branch_length;
    GetCoalescentRate(child_right, factor, epoch, coalescence_rate_data, leaves_child_right);

    leaves.resize(leaves_child_left.size() + leaves_child_right.size());
    std::vector<int>::iterator it_leaves = leaves.begin();

    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      *it_leaves = *it_leaves_child_left;
      it_leaves++;
    }
    for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){
      *it_leaves = *it_leaves_child_right;
      it_leaves++;
    }

    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){

        assert(*it_leaves_child_left != *it_leaves_child_right);

        if(*it_leaves_child_left < *it_leaves_child_right){
          for(int e = 0; e < (int) epoch.size() - 1; e++){
            if(coalescent_time < epoch[e+1]){
              assert(coalescent_time >= epoch[e]);
              coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
              coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - epoch[e]);
              break;
            }else{
              assert(coalescent_time >= epoch[e+1]);
              coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - epoch[e]);
            }
          }
        }else{
          for(int e = 0; e < (int) epoch.size() - 1; e++){
            if(coalescent_time < epoch[e+1]){
              assert(coalescent_time >= epoch[e]);
              coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
              coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - epoch[e]);
              break;
            }else{
              assert(coalescent_time >= epoch[e+1]);
              coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - epoch[e]);
            }
          } 
        }

      }
    }

    return coalescent_time;

  }else{

    leaves.push_back(n.label);
    return 0.0;

  }


}

float
GetCoalescentRate(Node n, float factor, std::vector<float>& epoch, std::vector<double>& sample_ages, std::vector<CollapsedMatrix<float>>& coalescence_rate_data, std::vector<int>& leaves){

  //go through tree and fill in the matrix
  if(n.child_left != NULL){

    std::vector<int> leaves_child_left, leaves_child_right;

    Node child_left  = *n.child_left;
    Node child_right = *n.child_right;

    int num_children_left, num_children_right;

    float coalescent_time = GetCoalescentRate(child_left, factor, epoch, sample_ages, coalescence_rate_data, leaves_child_left) + child_left.branch_length;
    GetCoalescentRate(child_right, factor, epoch, sample_ages, coalescence_rate_data, leaves_child_right);

    leaves.resize(leaves_child_left.size() + leaves_child_right.size());
    std::vector<int>::iterator it_leaves = leaves.begin();

    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      *it_leaves = *it_leaves_child_left;
      it_leaves++;
    }
    for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){
      *it_leaves = *it_leaves_child_right;
      it_leaves++;
    }

    double max_sample_age = 0.0;
    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){

        assert(*it_leaves_child_left != *it_leaves_child_right);

        max_sample_age = sample_ages[*it_leaves_child_left];
        if(max_sample_age < sample_ages[*it_leaves_child_right]){
          max_sample_age = sample_ages[*it_leaves_child_right];
        }

        if(max_sample_age == 0.0){
          if(*it_leaves_child_left < *it_leaves_child_right){
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(coalescent_time < epoch[e+1]){
                assert(coalescent_time >= epoch[e]);
                coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
                coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - epoch[e]);
                break;
              }else{
                assert(coalescent_time >= epoch[e+1]);
                coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - epoch[e]);
              }
            }
          }else{
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(coalescent_time < epoch[e+1]){
                assert(coalescent_time >= epoch[e]);
                coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
                coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - epoch[e]);
                break;
              }else{
                assert(coalescent_time >= epoch[e+1]);
                coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - epoch[e]);
              }
            } 
          }
        }else{        
          if(*it_leaves_child_left < *it_leaves_child_right){
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(max_sample_age < epoch[e+1]){
                if(max_sample_age >= epoch[e]){                
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - max_sample_age);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - max_sample_age);
                  }
                }else{
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor;
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (coalescent_time - epoch[e]);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor * (epoch[e+1] - epoch[e]);
                  }
                }
              }
            }
          }else{
            for(int e = 0; e < (int) epoch.size() - 2; e++){
              if(max_sample_age < epoch[e+1]){
                if(max_sample_age >= epoch[e]){                
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - max_sample_age);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - max_sample_age);
                  }
                }else{
                  if(coalescent_time < epoch[e+1]){
                    assert(coalescent_time >= epoch[e]);
                    coalescence_rate_data[e][*it_leaves_child_right][*it_leaves_child_left] += factor;
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (coalescent_time - epoch[e]);
                    break;
                  }else{
                    assert(coalescent_time >= epoch[e+1]);
                    coalescence_rate_data[e][*it_leaves_child_left][*it_leaves_child_right] += factor * (epoch[e+1] - epoch[e]);
                  }
                }
              }
            }
          }
        }

      }
    }

    return coalescent_time;

  }else{

    leaves.push_back(n.label);
    return sample_ages[n.label];

  }


}

int CoalescentRateForSection(cxxopts::Options& options, int chr = -1){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: years_per_gen, num_bins." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .anc file and calculates pairwise coalescence rate. Output is bin file. Use SummarizeCoalesecntRate to obtain coalescence rates." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  if(chr == -1){
    std::cerr << "Calculating coalescence rate for " << options["input"].as<std::string>() << " ..." << std::endl;
  }else{
    std::cerr << "Calculating coalescence rate for " << options["input"].as<std::string>() << "_chr" << chr << " ..." << std::endl;
  }

  ////////////////////////
  //read in anc file

  AncMutIterators ancmut;
  if(chr == -1){
    ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  }else{
    ancmut.OpenFiles(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc", options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
  }   
  
  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0; 
  Data data(N,L);

  //Read mut
  Mutations mut;

  if(chr == -1){
    mut.Read(options["input"].as<std::string>() + ".mut");
  }else{
    mut.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
  }

  //calculate epoch times

  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epoch(num_epochs);
  std::vector<CollapsedMatrix<float>> coalescence_rate_data(num_epochs);
  epoch[0] = 0.0;
  coalescence_rate_data[0].resize(data.N, data.N);
  std::fill(coalescence_rate_data[0].vbegin(), coalescence_rate_data[0].vend(), 0.0);
  epoch[1] = 1e3/years_per_gen;
  coalescence_rate_data[1].resize(data.N, data.N);
  std::fill(coalescence_rate_data[1].vbegin(), coalescence_rate_data[1].vend(), 0.0);

  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){

    epoch[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    coalescence_rate_data[e].resize(data.N, data.N);
    std::fill(coalescence_rate_data[e].vbegin(), coalescence_rate_data[e].vend(), 0.0);

  }
  epoch[num_epochs-1] = 1e8/years_per_gen;
  coalescence_rate_data[num_epochs-1].resize(data.N, data.N);
  std::fill(coalescence_rate_data[num_epochs-1].vbegin(), coalescence_rate_data[num_epochs-1].vend(), 0.0);


  ////////////////////////////////
  //Pairwise coalescence rate
  //In each tree, find the coalescent time. Then update count_per_epoch and coalescent_time_in_epoch. 

  if(ancmut.sample_ages.size() > 0){
    float factor = 0.0;
    num_bases_tree_persists = 0.0;
    while(num_bases_tree_persists >= 0.0){
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      std::vector<int> leaves;
      factor = 1.0;
      GetCoalescentRate(*std::prev(mtr.tree.nodes.end(),1), factor, epoch, ancmut.sample_ages, coalescence_rate_data, leaves);
    }  
  }else{
    float factor = 0.0;
    num_bases_tree_persists = 0.0;
    while(num_bases_tree_persists >= 0.0){
      num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
      std::vector<int> leaves;
      factor = 1.0;
      GetCoalescentRate(*std::prev(mtr.tree.nodes.end(),1), factor, epoch, coalescence_rate_data, leaves);
    }  
  }

  //output as bin
  FILE* fp;
  if(chr == -1){
    fp = fopen((options["output"].as<std::string>() + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + std::to_string(chr) + ".bin" ).c_str(), "wb");  
  }

  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epoch[0], sizeof(float), epoch.size(), fp);
  for(std::vector<CollapsedMatrix<float>>::iterator it_coalescence_rate_data = coalescence_rate_data.begin(); it_coalescence_rate_data != coalescence_rate_data.end();){
    (*it_coalescence_rate_data).DumpToFile(fp);
    it_coalescence_rate_data++;
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

  return 0;

}

//TODO: the functions below need to be debugged and adapted to ancient samples
//////// functions for estimating directional migration /////////

int 
MatrixIndex(int i, int j, int k, int N){

  //need to make this function as efficient as possible

  int f = 0;

  int tmp, tmp2;
  if(i < k && k < j){

    //case i < k < j
    f = 2;

    tmp = j;
    j = k;
    k = tmp;

  }else if(k < i && i < j){

    //case k < i < j
    f = 2;

    tmp  = j;
    tmp2 = k;
    j = i;
    i = tmp2;
    k = tmp;

  }else if(k < j && j < i){

    //case k < j < i
    f = 1;

    tmp  = k;
    k = i;
    i = tmp;

  }else if(j < i && i < k){

    //case j < i < k

    tmp  = j;
    j = i;
    i = tmp;

  }else if(j < k && k < i){

    //case j < k < i
    f = 1;

    tmp  = j;
    tmp2 = k;
    k = i;
    i = tmp;
    j = tmp2;

  }

  assert(i >= 0);
  assert(i < j);
  assert(j < k);
  assert(k < N);

  int res = 0;

  if(i > 0){
    res  = N* ( (N-3.0)/2.0 - (i-1.0)/2.0 - 1.0 ) + 3.0/4.0*(i-1.0) + 1.0;
    res *= i;
  }
  if(j > 0){
    res += j * (N - (j+1.0)/2.0 - 1.0);
  }
  res   += (i+1.0)*(i+2.0)/2.0 + k - 1.0 - N;

  if(i-1 >= 0){
    tmp = 0;
    for(int n = 0; n < i; n++){
      tmp += n*n;
    }
    res += 0.5 * tmp;
  }

  res *= 3;
  res += f;

  return(res);

}

float
GetDirCoalescentRate(Node n, int N, std::vector<float>& epoch, std::vector<CollapsedMatrix<float>>& coalescence_rate_data, CollapsedMatrix<float>& ctime, std::vector<int>& leaves){

  //go through tree and fill in the matrix
  if(n.child_left != NULL){

    std::vector<int> leaves_child_left, leaves_child_right;

    Node child_left  = *n.child_left;
    Node child_right = *n.child_right;

    int num_children_left, num_children_right;

    float coalescent_time = GetDirCoalescentRate(child_left, N, epoch, coalescence_rate_data, ctime, leaves_child_left) + child_left.branch_length;
    GetDirCoalescentRate(child_right, N, epoch, coalescence_rate_data, ctime, leaves_child_right);

    leaves.resize(leaves_child_left.size() + leaves_child_right.size());
    std::vector<int>::iterator it_leaves = leaves.begin();

    for(std::vector<int>::iterator it_leaves_child_left = leaves_child_left.begin(); it_leaves_child_left != leaves_child_left.end(); it_leaves_child_left++){
      *it_leaves = *it_leaves_child_left;
      it_leaves++;
    }
    for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){
      *it_leaves = *it_leaves_child_right;
      it_leaves++;
    }

    int m;
    for(std::vector<int>::iterator it_leaves_child_left1 = leaves_child_left.begin(); it_leaves_child_left1 != leaves_child_left.end(); it_leaves_child_left1++){
      for(std::vector<int>::iterator it_leaves_child_left2 = std::next(it_leaves_child_left1, 1); it_leaves_child_left2 != leaves_child_left.end(); it_leaves_child_left2++){
        for(std::vector<int>::iterator it_leaves_child_right = leaves_child_right.begin(); it_leaves_child_right != leaves_child_right.end(); it_leaves_child_right++){

          //(*it_leaves_child_left1, *it_leaves_child_left2), *it_leaves_child_right)
          //I need coalescence time of the two children to the left and the current coalescence time.
          m = MatrixIndex(*it_leaves_child_left1, *it_leaves_child_left2, *it_leaves_child_right, N);
          int e = 0;
          for(; e < (int) epoch.size() - 1; e++){
            if(ctime[*it_leaves_child_left1][*it_leaves_child_left2] < epoch[e]){
              break;
            }
          }

          if(coalescent_time < epoch[e+1]){
            coalescence_rate_data[e][0][m] += 1.0;
            coalescence_rate_data[e][1][m] += (coalescent_time - ctime[*it_leaves_child_left1][*it_leaves_child_left2]);
          }else{
            coalescence_rate_data[e][1][m] += (epoch[e+1] - ctime[*it_leaves_child_left1][*it_leaves_child_left2]);
            e++;
            for(; e < (int) epoch.size() - 1; e++){
              if(coalescent_time < epoch[e+1]){
                assert(coalescent_time >= epoch[e]);
                coalescence_rate_data[e][0][m] += 1.0;
                coalescence_rate_data[e][1][m] += (coalescent_time - epoch[e]);
                break;
              }else{
                assert(coalescent_time >= epoch[e+1]);
                coalescence_rate_data[e][1][m] += (epoch[e+1] - epoch[e]);
              }
            } 
          }


        }
      }
      for(std::vector<int>::iterator it_leaves_child_right1 = leaves_child_right.begin(); it_leaves_child_right1 != leaves_child_right.end(); it_leaves_child_right1++){
        for(std::vector<int>::iterator it_leaves_child_right2 = std::next(it_leaves_child_right1, 1); it_leaves_child_right2 != leaves_child_right.end(); it_leaves_child_right2++){

          //(*it_leaves_child_right1, *it_leaves_child_right2), *it_leaves_child_left1)
          //I need coalescence time of the two children to the right and the current coalescence time.
          m = MatrixIndex(*it_leaves_child_right1, *it_leaves_child_right2, *it_leaves_child_left1, N);
          int e = 0;
          for(; e < (int) epoch.size() - 1; e++){
            if(ctime[*it_leaves_child_right1][*it_leaves_child_right2] < epoch[e]){
              break;
            }
          }

          if(coalescent_time < epoch[e+1]){
            coalescence_rate_data[e][0][m] += 1.0;
            coalescence_rate_data[e][1][m] += (coalescent_time - ctime[*it_leaves_child_right1][*it_leaves_child_right2]);
          }else{
            coalescence_rate_data[e][1][m] += (epoch[e+1] - ctime[*it_leaves_child_right1][*it_leaves_child_right2]);
            e++;
            for(; e < (int) epoch.size() - 1; e++){
              if(coalescent_time < epoch[e+1]){
                assert(coalescent_time >= epoch[e]);
                coalescence_rate_data[e][0][m] += 1.0;
                coalescence_rate_data[e][1][m] += (coalescent_time - epoch[e]);
                break;
              }else{
                assert(coalescent_time >= epoch[e+1]);
                coalescence_rate_data[e][1][m] += (epoch[e+1] - epoch[e]);
              }
            } 
          }


        }

        ctime[*it_leaves_child_left1][*it_leaves_child_right1] = coalescent_time;
        ctime[*it_leaves_child_right1][*it_leaves_child_left1]  = coalescent_time;

      }
    }



    return coalescent_time;

  }else{

    leaves.push_back(n.label);
    return 0.0;

  }


}

int CoalescentRateDir(cxxopts::Options& options, int chr = -1){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: years_per_gen, num_bins." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .anc file and calculates coalescence rate between three haplotypes. Output is bin file. Use SummarizeCoalesecntRate to obtain coalescence rates." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  if(chr == -1){
    std::cerr << "Calculating coalescence rate for " << options["input"].as<std::string>() << " ..." << std::endl;
  }else{
    std::cerr << "Calculating coalescence rate for " << options["input"].as<std::string>() << "_chr" << chr << " ..." << std::endl;
  }

  ////////////////////////
  //read in anc file

  int N;
  igzstream is_N;
  if(chr == -1){
    is_N.open(options["input"].as<std::string>() + ".anc");
    if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
  }else{
    is_N.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc");
    if(is_N.fail()) is_N.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc.gz");
  }
  if(is_N.fail()){
    std::cerr << "Error while opening .anc file." << std::endl;
    exit(1);
  } 
  is_N.ignore(256, ' ');
  is_N >> N;
  is_N.close();

  int L = 0;
  igzstream is_L;
  if(chr == -1){
    is_L.open(options["input"].as<std::string>() + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz");
  }else{
    is_L.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
    if(is_L.fail()) is_L.open(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut.gz");
  }
  if(is_L.fail()){
    std::cerr << "Error while opening .mut file." << std::endl;
    exit(1);
  } 

  std::string unused;
  std::getline(is_L, unused); 
  while ( std::getline(is_L, unused) ){
    ++L;
  }
  is_L.close();
  Data data(N,L);
  int dim = data.N * (data.N-1.0) * (data.N-2.0)/2.0;

  //Read anc
  AncesTree anc;
  Mutations mut;

  if(chr == -1){
    anc.Read(options["input"].as<std::string>() + ".anc");
    mut.Read(options["input"].as<std::string>() + ".mut");
  }else{
    anc.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".anc");
    mut.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
  }

  //calculate epoch times

  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epoch(num_epochs);
  std::vector<CollapsedMatrix<float>> coalescence_rate_data(num_epochs);
  epoch[0] = 0.0;
  coalescence_rate_data[0].resize(2, dim);
  std::fill(coalescence_rate_data[0].vbegin(), coalescence_rate_data[0].vend(), 0.0);
  epoch[1] = 1e3/years_per_gen;
  coalescence_rate_data[1].resize(2, dim);
  std::fill(coalescence_rate_data[1].vbegin(), coalescence_rate_data[1].vend(), 0.0);

  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){

    epoch[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
    coalescence_rate_data[e].resize(2, dim);
    std::fill(coalescence_rate_data[e].vbegin(), coalescence_rate_data[e].vend(), 0.0);

  }
  epoch[num_epochs-1] = 1e8/years_per_gen;
  coalescence_rate_data[num_epochs-1].resize(2, dim);
  std::fill(coalescence_rate_data[num_epochs-1].vbegin(), coalescence_rate_data[num_epochs-1].vend(), 0.0);


  ////////////////////////////////
  //Pairwise coalescence rate
  //In each tree, find the coalescent time. Then update count_per_epoch and coalescent_time_in_epoch. 

  int prog_step = anc.seq.size()/100.0+1;
  int percent = 0, tree_count = 0;

  CollapsedMatrix<float> ctime;
  ctime.resize(data.N, data.N);
  CorrTrees::iterator it_anc = anc.seq.begin();
  for(; it_anc != std::prev(anc.seq.end(),1); it_anc++){


    if(tree_count % prog_step == 0){
      std::cerr << "[" << percent << "%]\r";
      std::cerr.flush();
      percent++;
    }

    std::vector<int> leaves;
    std::fill(ctime.vbegin(), ctime.vend(), 0); //not needed - check when debugging
    GetDirCoalescentRate(*std::prev((*it_anc).tree.nodes.end(),1), data.N, epoch, coalescence_rate_data, ctime, leaves);

    tree_count++;

  }


  //output as bin
  FILE* fp;
  if(chr == -1){
    fp = fopen((options["output"].as<std::string>() + ".bin" ).c_str(), "wb");  
  }else{
    fp = fopen((options["output"].as<std::string>() + "_chr" + std::to_string(chr) + ".bin" ).c_str(), "wb");  
  }

  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epoch[0], sizeof(float), epoch.size(), fp);
  for(std::vector<CollapsedMatrix<float>>::iterator it_coalescence_rate_data = coalescence_rate_data.begin(); it_coalescence_rate_data != coalescence_rate_data.end();){
    (*it_coalescence_rate_data).DumpToFile(fp);
    it_coalescence_rate_data++;
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

////////////////////////////////////////////////////////////////
//Conditional coalescence rate

//Calculate the coalescence rate between labelled branch and any other branch in the genealogy.
//For each labelled branch, I need the branchID and start_time in the target tree

void
GetConditionalCoalescenceRate(std::vector<float>& start_times, std::vector<float>& end_times, std::vector<int>& flag, std::vector<int>& num_repeat, std::vector<int>& ind1, std::vector<int>& ind2, std::vector<int>& ind3, std::vector<float>& epoch, std::vector<std::vector<std::vector<float>>>& numer, std::vector<std::vector<std::vector<float>>>& denom){

  assert(num_repeat.size() == start_times.size());
  assert(num_repeat.size() == end_times.size());
  assert(ind1.size() == flag.size());
  assert(ind2.size() == flag.size());
  assert(ind3.size() == flag.size());
  int total_num_repeat = 0;
  for(std::vector<int>::iterator it_num_repeat = num_repeat.begin(); it_num_repeat != num_repeat.end(); it_num_repeat++){
    total_num_repeat += *it_num_repeat;
  }
  assert(flag.size() == total_num_repeat);

  //assuming numer and denom are already initialsed

  std::vector<float>::iterator it_start_times = start_times.begin(), it_end_times = end_times.begin();
  std::vector<int>::iterator it_ind1 = ind1.begin(), it_ind2 = ind2.begin(), it_ind3 = ind3.begin(), it_flag = flag.begin();
  std::vector<int>::iterator it_num_repeat = num_repeat.begin();

  for(; it_start_times != start_times.end();){

    if(*it_start_times <= *it_end_times){

      std::vector<float>::iterator it_epoch = epoch.begin();

      int e = 0;
      while(*it_epoch <= *it_start_times){
        it_epoch++;
        e++;
        if(it_epoch == epoch.end()) break;
      }
      //it_epoch--;
      e--;

      //assert(*it_epoch <= *it_start_times);
      //assert(*it_ind1 < numer.size());
      //assert(*it_ind2 < numer[*it_ind1].size());
      //assert(*it_ind1 < denom.size());
      //assert(*it_ind2 < denom[*it_ind1].size());

      //////
      float lower_time = *it_start_times;
      while(lower_time < *it_end_times){

        if(it_epoch == epoch.end()){
          //assert((*it_ind3) * epoch.size()+e-1 < numer[*it_ind1][*it_ind2].size());
          //assert((*it_ind3) * epoch.size()+e-1 < denom[*it_ind1][*it_ind2].size());

          for(int i = 0; i < *it_num_repeat; i++){
            denom[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] += (*it_end_times - lower_time);
            assert(denom[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] >= 0.0);
            if(*it_flag == 10 || *it_flag == 11){
              numer[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] += 1.0;
            }
          }

          //done
          break;
        }

        //assert((*it_ind3) * epoch.size()+e < numer[*it_ind1][*it_ind2].size());
        //assert((*it_ind3) * epoch.size()+e < denom[*it_ind1][*it_ind2].size());
        if(*it_epoch < *it_end_times){
          for(int i = 0; i < *it_num_repeat; i++){
            denom[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] += (*it_epoch - lower_time);
            assert(denom[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] >= 0.0);
          }
          lower_time = *it_epoch;
        }else{
          for(int i = 0; i < *it_num_repeat; i++){
            denom[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] += (*it_end_times - lower_time);
            assert(denom[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] >= 0.0);
            if(*std::next(it_flag,i) == 10 || *std::next(it_flag,i) == 11){
              numer[*std::next(it_ind1,i)][*std::next(it_ind2,i)][(*std::next(it_ind3,i))*epoch.size()+e] += 1.0;
            }
          }
          //done
          break;
        }
        e++;
        it_epoch++;

      }

    }

    it_ind1 = std::next(it_ind1, *it_num_repeat);
    it_ind2 = std::next(it_ind2, *it_num_repeat);
    it_ind3 = std::next(it_ind3, *it_num_repeat);
    it_flag = std::next(it_flag, *it_num_repeat);

    it_start_times++;
    it_end_times++;
    it_num_repeat++;

  }

}

void
GetBranchesForIntrogression(Tree& tree, std::vector<int>& flag, std::vector<float>& start_times, std::vector<float>& end_times, std::vector<float>& sorted_coords){

  std::vector<float> coords;
  flag.clear();
  start_times.clear();
  end_times.clear();
  sorted_coords.clear();

  int N = ((int) (tree.nodes.size() + 1)/2.0);
  tree.GetCoordinates(coords);

  for(int i = 0; i < tree.nodes.size(); i++){
    if(coords[i] <= 5000/28.0 && coords[(*tree.nodes[i].parent).label] >= 5000/28.0){
      end_times.push_back(coords[(*tree.nodes[i].parent).label]);
      flag.push_back(10);
    }
  }
  start_times.resize(flag.size());
  std::fill(start_times.begin(), start_times.end(), 5000/28.0);

  //for(int i = 0; i < tree.nodes.size()-1; i++){
  //  //branches.push_back(i);
  //  end_times.push_back(coords[(*tree.nodes[i].parent).label]);
  //  start_times.push_back(coords[i]);
  //  flag.push_back(10);
  //}

  sort(coords.begin(), coords.end());
  sorted_coords.resize(N);
  std::vector<float>::iterator it_sorted_coords = sorted_coords.begin();
  for(std::vector<float>::iterator it_coords = std::next(coords.begin(), N-1); it_coords != coords.end(); it_coords++){
    *it_sorted_coords = *it_coords;
    it_sorted_coords++;
  }

}

void
Populate(Tree& subtree, int child, float etime, int i_sample, std::vector<Leaves>& leaves, std::vector<float>& coords, std::vector<int>& convert_back, std::vector<float>& start_times, std::vector<float>& end_times, std::vector<int>& num_repeat, std::vector<int>& ind1, std::vector<int>& ind2, std::vector<int>& ind3, std::vector<int>& flag){

  if(subtree.nodes[child].child_left != NULL){

    start_times.push_back(coords[child]);
    end_times.push_back(etime);
    int child1 = (*subtree.nodes[child].child_left).label;
    int child2 = (*subtree.nodes[child].child_right).label;
    num_repeat.push_back(2*leaves[child1].num_leaves*leaves[child2].num_leaves);

    for(std::vector<int>::iterator it_member1 = leaves[child1].member.begin(); it_member1 != leaves[child1].member.end(); it_member1++){
      for(std::vector<int>::iterator it_member2 = leaves[child2].member.begin(); it_member2 != leaves[child2].member.end(); it_member2++){
        ind1.push_back(convert_back[*it_member1]);
        ind2.push_back(convert_back[*it_member2]);
        ind3.push_back(i_sample); 
        flag.push_back(10);

        ind1.push_back(convert_back[*it_member2]);
        ind2.push_back(convert_back[*it_member1]);
        ind3.push_back(i_sample); 
        flag.push_back(10);
      }
    }

    Populate(subtree, child1, etime, i_sample, leaves, coords, convert_back, start_times, end_times, num_repeat, ind1, ind2, ind3, flag);
    Populate(subtree, child2, etime, i_sample, leaves, coords, convert_back, start_times, end_times, num_repeat, ind1, ind2, ind3, flag);

  }

}

void
GetBranchesForMigration(Tree& tree, Tree& subtree, Sample& sample, std::vector<int>& flag, std::vector<float>& start_times, std::vector<float>& end_times, std::vector<int>& num_repeat, std::vector<int>& ind1, std::vector<int>& ind2, std::vector<int>& ind3, const int pop1, const int pop2){

  //int pop1 = 0, pop2 = 1; //need to make these to arguments
  flag.clear();
  start_times.clear();
  end_times.clear();
  num_repeat.clear();
  //vectors saying which triplet this is
  ind1.clear();
  ind2.clear();
  ind3.clear();

  Sample sample_orig = sample;
  sample.AssignPopOfInterest(sample.groups[pop2]);
  //get subtree
  std::vector<int> convert_index, number_in_subpop;
  //convert_index: maps old branch id to new branch id (equals -1 if it doesn't exists)
  //number_in_subpop: number of descendants in subpop (takes old branch ids)
  tree.GetSubTree(sample, subtree, convert_index, number_in_subpop); 
  int N_sub = (subtree.nodes.size() + 1.0)/2.0, N = (tree.nodes.size() + 1.0)/2.0;    

  std::vector<int> convert_back(N_sub);
  int i = 0;
  for(std::vector<int>::iterator it_convert = convert_index.begin(); it_convert != std::next(convert_index.begin(), N); it_convert++){
    if(*it_convert != -1){
      assert(*it_convert < N_sub);
      convert_back[*it_convert] = i;
    }
    i++;
  }

  std::vector<Leaves> leaves;
  subtree.FindAllLeaves(leaves);

  std::vector<float> coords;
  subtree.GetCoordinates(coords);
  std::vector<float> coords_full; //need for start time
  tree.GetCoordinates(coords_full);

  int root = 2*N_sub-2;
  int node, parent, child;
  for(int i_sample = 0; i_sample < N; i_sample++){

    if(sample_orig.group_of_haplotype[i_sample] == pop1){

      ///////////////////
      //get branches for ((pop1, pop2), pop2)

      node = (*tree.nodes[i_sample].parent).label;
      if(pop1 != pop2){
        while(number_in_subpop[node] == 0){
          node = (*tree.nodes[node].parent).label;
          if(node == tree.nodes.size()) break;
        }
      }else{
        while(number_in_subpop[node] == 1){
          node = (*tree.nodes[node].parent).label;
          if(node == tree.nodes.size()) break;
        }
      }

      //node is the node in the full tree from where we want to measure the coalescence rate.
      if(node < tree.nodes.size() - 1){

        std::vector<float> stimes;
        std::vector<int> ntimes, ind2_index;

        //include node
        int num_in_subpop = number_in_subpop[node];
        if(pop1 == pop2) num_in_subpop--;

        if(pop1 == pop2){
          //need to go one up and then down
          child = (*tree.nodes[node].child_left).label;
          assert(convert_index[child] >= 0);
          if(leaves[convert_index[child]].member[0] != convert_index[i_sample]){
            ind2_index.push_back(convert_index[child]);
            assert(ind2_index[0] >= 0);
            assert(ind2_index[0] != convert_index[i_sample]);
          }else{
            child = (*tree.nodes[node].child_right).label;
            ind2_index.push_back(convert_index[child]);
            assert(leaves[convert_index[child]].member[0] != convert_index[i_sample]);
            assert(ind2_index[0] >= 0);
            assert(ind2_index[0] != convert_index[i_sample]);
          }
        }else{
          ind2_index.push_back(convert_index[node]); //ind2
          assert(ind2_index[0] >= 0);
        }

        assert(num_in_subpop > 0);
        float stime = coords_full[node];
        node        = convert_index[node];
        assert(node >= 0);
        assert(node < leaves.size());
        if(pop1 == pop2) assert(num_in_subpop+1 == leaves[node].num_leaves);
        if(pop1 != pop2) assert(num_in_subpop == leaves[node].num_leaves);
        if(node != root){
          assert(node >= 0);
          //assume we want to look at all triplets
          //therefore, we want to add each branch as often as it would appear (ntimes x num_leaves)
          stimes.push_back(stime);
          ntimes.push_back(num_in_subpop); 
          while(subtree.nodes[node].parent != NULL){

            assert(node >= 0);
            assert(node < leaves.size());
            int num_leaves = leaves[node].num_leaves;
            parent = (*subtree.nodes[node].parent).label;               
            assert(parent >= 0);
            assert(parent < leaves.size());
            num_leaves = leaves[parent].num_leaves - num_leaves;   

            child = (*subtree.nodes[parent].child_left).label;
            if(child != node){
              ind2_index.push_back(child);
            }else{
              child = (*subtree.nodes[parent].child_right).label;
              ind2_index.push_back(child);
            }
            assert(child >= 0);
            assert(child < leaves.size());
            assert(leaves[child].num_leaves == num_leaves);

            //ntimes.size(): lower coalescence point
            //ntimes[k]: number of descendants of lower coalescence point, minus already accounted for
            //num_leaves: number of tips on the other side of parent
            for(int k = 0; k < ntimes.size(); k++){  
              assert(k < stimes.size());
              assert(ind2_index[k] < leaves.size());
              assert(ind2_index[k] >= 0);
              assert(leaves[ind2_index[k]].member.size() == ntimes[k]);
              start_times.push_back(stimes[k]);
              end_times.push_back(coords[parent]);
              num_repeat.push_back(ntimes[k] * num_leaves);
              for(int l = 0; l < ntimes[k]; l++){
                for(int m = 0; m < num_leaves; m++){
                  ind1.push_back(i_sample);
                  ind2.push_back(convert_back[leaves[ind2_index[k]].member[l]]); //need to convert back to full tree
                  ind3.push_back(convert_back[leaves[child].member[m]]);  //need to convert back to full tree
                  flag.push_back(10);
                }
              }              
            }  
            ntimes.push_back(num_leaves);
            stimes.push_back(coords[parent]);

            node = parent; 
          }

        }

      }

      if(pop1 != pop2){

        ///////////////////
        //get branches for ((pop2, pop2), pop1)
        //strategy: i_sample is now ind3, need to get pairs of pop2
        node = i_sample;
        int prev_node = -1;
        while(number_in_subpop[node] <= 1){
          prev_node = node;
          node = (*tree.nodes[node].parent).label;
          if(node == tree.nodes.size()-1) break;
        }
        //now there are at least two samples from pop2 under node
        child = (*tree.nodes[node].child_left).label;
        if(child == prev_node){
          child = (*tree.nodes[node].child_right).label;
        }
        child = convert_index[child];
        //child is root of subtree underneath node
        Populate(subtree, child, coords_full[node], i_sample, leaves, coords, convert_back, start_times, end_times, num_repeat, ind1, ind2, ind3, flag);

        node = convert_index[node];
        while(node < root){
          prev_node = node;
          node = (*subtree.nodes[node].parent).label;
          child = (*subtree.nodes[node].child_left).label;
          if(child == prev_node){
            child = (*subtree.nodes[node].child_right).label;
          }
          Populate(subtree, child, coords[node], i_sample, leaves, coords, convert_back, start_times, end_times, num_repeat, ind1, ind2, ind3, flag);
        
        }

      }//if pop1 != pop2

    }//end of if to check if i_sample is in pop1
  }// end for i_sample

}

int 
ConditionalCoalescenceRateForSection(cxxopts::Options& options, int chr = -1){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output. Optional: years_per_gen, num_bins, poplabels." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Reads .anc file and calculates pairwise coalescence rate. Output is bin file. Use SummarizeCoalesecntRate to obtain coalescence rates." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  if(chr == -1){
    std::cerr << "Calculating conditional coalescence rate for " << options["input"].as<std::string>() << " ..." << std::endl;
  }else{
    std::cerr << "Calculating conditional coalescence rate for " << options["input"].as<std::string>() << "_chr" << chr << " ..." << std::endl;
  }

  ////////////////////////
  //Parse anc/mut

  MarginalTree mtr;
  Muts::iterator it_mut;
  double num_bases_of_tree = 0.0; //Starts halfway between previous SNP and SNP where tree starts and goes to halfway of SNP where tree ends 

  Sample sample;
  std::string label;
  if(options.count("poplabels") > 0){
    sample.Read(options["poplabels"].as<std::string>());
  }
  int num_pops = sample.groups.size(); //adjust this eventually
  if(num_pops != 2){
    std::cerr << "Error: This function assumes two populations in the poplabels file" << std::endl;
    exit(1);
  }
  assert(num_pops == 2);
  int N = sample.group_of_haplotype.size();

  ///////////////////////////////////////////////////////////
  //calculate epoch times

  float years_per_gen = 28.0;
  if(options.count("years_per_gen")){
    years_per_gen = options["years_per_gen"].as<float>();
  }

  int num_epochs = 30;
  if(options.count("num_bins") > 0){
    num_epochs = options["num_bins"].as<int>();
  }
  num_epochs++;
  std::vector<float> epoch(num_epochs);
  epoch[0] = 0.0;
  epoch[1] = 1e3/years_per_gen;
  float log_10 = std::log(10);
  for(int e = 2; e < num_epochs-1; e++){
    epoch[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
  }
  epoch[num_epochs-1] = 1e8/years_per_gen;

  ////////////////////////////////
  //Triplet coalescence rate
  //In each tree, find the coalescent time. Then update count_per_epoch and coalescent_time_in_epoch. 

  std::vector<std::vector<std::vector<float>>> numer, denom;
  numer.resize(N);
  denom.resize(N);
  for(int i = 0; i < numer.size(); i++){
    numer[i].resize(N);
    denom[i].resize(N);
    for(int j = 0; j < numer[i].size(); j++){
      numer[i][j].resize(epoch.size()*N);
      std::fill(numer[i][j].begin(), numer[i][j].end(), 0.0);

      denom[i][j].resize(epoch.size()*N);
      std::fill(denom[i][j].begin(), denom[i][j].end(), 0.0);
    }
  }


  //////////////////////
  //populate numer and denom

  for(int pop1 = 0; pop1 < num_pops; pop1++){
    for(int pop2 = 0; pop2 < num_pops; pop2++){

      std::cerr << pop1 << " " << pop2 << std::endl;
      AncMutIterators ancmut(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
      num_bases_of_tree = 0.0;

      Tree subtree;
      std::vector<float> start_times, end_times;
      std::vector<int> num_repeat,flag, ind1, ind2, ind3;

      while(num_bases_of_tree >= 0.0){
        num_bases_of_tree = ancmut.NextTree(mtr, it_mut);

        ////////////////////////////////////
        //find branches satisfying condition

        //populate branches and start_times
        //GetBranchesForIntrogression(mtr.tree, flag, start_times, end_times, sorted_coords);
        GetBranchesForMigration(mtr.tree, subtree, sample, flag, start_times, end_times, num_repeat, ind1, ind2, ind3, pop1, pop2);

        //////////////////////////////////////
        //flag: 01 - invalid external, 00 - invalid internal, 10 - valid internal, 11 - valid external
        GetConditionalCoalescenceRate(start_times, end_times, flag, num_repeat, ind1, ind2, ind3, epoch, numer, denom);

      }

    }
  }

  //TODO: for across chromosome rates, need to code up a summarize and move the below into a finalize function
  ///////////
  //now calculate group-wise coalescence rates

  std::vector<std::vector<float>> group_rate(6);
  std::vector<std::vector<int>> group_count(6);
  for(std::vector<std::vector<float>>::iterator it_group_rate = group_rate.begin(); it_group_rate != group_rate.end(); it_group_rate++){
    (*it_group_rate).resize(epoch.size());
    std::fill((*it_group_rate).begin(), (*it_group_rate).end(), 0.0);
  }
  for(std::vector<std::vector<int>>::iterator it_group_count = group_count.begin(); it_group_count != group_count.end(); it_group_count++){
    (*it_group_count).resize(epoch.size());
    std::fill((*it_group_count).begin(), (*it_group_count).end(), 0.0);
  }

  for(int i = 0; i < numer.size(); i++){
    for(int j = 0; j < numer[i].size(); j++){
      for(int k = 0; k < numer[i][j].size(); k++){
        if(denom[i][j][k] != 0.0){
          int l = (k - k % epoch.size())/epoch.size();
          int entry = sample.group_of_haplotype[i] * 2 + sample.group_of_haplotype[j];
          if(entry == 0){
            if(sample.group_of_haplotype[l] == 1) entry = 4;
          }
          if(entry == 3){
            if(sample.group_of_haplotype[l] == 0) entry = 5;
          }
          group_rate[entry][k % epoch.size()] += numer[i][j][k]/denom[i][j][k];
          group_count[entry][k % epoch.size()]++; 
        }
      }
    }
  }

  std::ofstream os;
  for(int pop1 = 0; pop1 < num_pops; pop1++){
    for(int pop2 = 0; pop2 < num_pops; pop2++){

      int entry = pop1 * 2 + pop2;

      os.open(options["output"].as<std::string>() + std::to_string(pop1) + std::to_string(pop2) + ".coal");
      os << "group" << pop1 << pop2 << std::endl;
      for(int e = 0; e < epoch.size(); e++){
        os << epoch[e] << " ";
      }
      os << std::endl;
      os << "0 0 ";
      for(int e = 0; e < epoch.size(); e++){
        os << group_rate[entry][e]/group_count[entry][e] << " ";
      }
      os << std::endl;
      os.close();

    }
  }

  int pop1 = 0, pop2 = 1;

  os.open(options["output"].as<std::string>() + std::to_string(pop1) + std::to_string(pop1) + std::to_string(pop2) + ".coal");
  os << "group" << pop1 << pop1 << pop1 << std::endl;
  for(int e = 0; e < epoch.size(); e++){
    os << epoch[e] << " ";
  }
  os << std::endl;
  os << "0 0 ";
  for(int e = 0; e < epoch.size(); e++){
    os << group_rate[4][e]/group_count[4][e] << " ";
  }
  os << std::endl;
  os.close();

  os.open(options["output"].as<std::string>() + std::to_string(pop2) + std::to_string(pop2) + std::to_string(pop1) + ".coal");
  os << "group" << pop2 << pop2 << pop2 << std::endl;
  for(int e = 0; e < epoch.size(); e++){
    os << epoch[e] << " ";
  }
  os << std::endl;
  os << "0 0 ";
  for(int e = 0; e < epoch.size(); e++){
    os << group_rate[5][e]/group_count[5][e] << " ";
  }
  os << std::endl;
  os.close();

  /*
  //output as bin
  FILE* fp;
  if(chr == -1){
  fp = fopen((options["output"].as<std::string>() + ".bin" ).c_str(), "wb");  
  }else{
  fp = fopen((options["output"].as<std::string>() + "_chr" + std::to_string(chr) + ".bin" ).c_str(), "wb");  
  }

  fwrite(&num_epochs, sizeof(int), 1, fp);
  fwrite(&epoch[0], sizeof(float), epoch.size(), fp);
  for(std::vector<CollapsedMatrix<float>>::iterator it_coalescence_rate_data = coalescence_rate_data.begin(); it_coalescence_rate_data != coalescence_rate_data.end();){
  (*it_coalescence_rate_data).DumpToFile(fp);
  it_coalescence_rate_data++;
  }

  fclose(fp);
  */

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

