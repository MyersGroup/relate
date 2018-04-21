#ifndef SAMPLE_HPP
#define SAMPLE_HPP

#include "gzstream.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>

struct Sample{

  std::vector<std::string> groups;
  std::vector<int> group_sizes;
  std::map<std::string, std::string> region_of_group;

  std::vector<int> group_of_interest;
  int group_of_interest_size;

  std::vector<int> group_of_haplotype;
  
  void Read(const std::string& filename); //read from .sample file
  std::string AssignPopOfInterest(const std::string& pop); //comma separated string with pop of interest

};



#endif //SAMPLE_HPP
