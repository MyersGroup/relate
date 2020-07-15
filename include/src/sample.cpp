#include "sample.hpp"

void
Sample::Read(const std::string& filename){

	bool diploid = true;

  std::string line, read, ploidy;
  bool exists;

  //Read all possible labels
  igzstream is(filename);
  getline(is, line);
  while(getline(is, line)){

    int i = 0;
    while(line[i] != ' ' && line[i] != '\t') i++;
    i++;
    read.clear();
    while(line[i] != ' ' && line[i] != '\t'){
      read += line[i];
      i++;
    }
    i++;
		ploidy.clear();
		while(line[i] != ' ' && line[i] != '\t'){
			i++;
		}
		i++;
		ploidy.clear();
		while(line[i] != ' ' && line[i] != '\t'){
			ploidy += line[i];
			i++;
			if(i == line.size()) break;
		}
		i++;
   
		if(ploidy != "NA"){
      if(ploidy == "1"){
				diploid = false;
			}else{
        if(!diploid){
          std::cerr << "Error: Detected both haploid and diploid samples." << std::endl;
					exit(1);
				}
			}
		}
    exists = false;
    for(std::vector<std::string>::iterator it_groups = groups.begin(); it_groups != groups.end(); it_groups++){
      if(!read.compare(*it_groups)){
        exists = true;
        break;
      }
    }
    if(!exists){
      groups.push_back(read);

      read.clear();
      while(line[i] != ' ' && line[i] != '\t'){
        read += line[i];
        i++;
      }

      region_of_group[*std::prev(groups.end(),1)] = read;
    }

  }
  is.close();
  std::sort(groups.begin(), groups.end());

  //read group_of_haplotype
  is.open(filename);
  if(is.fail()){
    std::cerr << "Error while opening file." << std::endl;
    exit(1);
  }
  getline(is, line);
  int ind;
  while(getline(is, line)){

    int i = 0;
    while(line[i] != ' ' && line[i] != '\t') i++;
    i++;
    read.clear();
    while(line[i] != ' ' && line[i] != '\t'){
      read += line[i];
      i++;
    }

    ind = 0;
    for(std::vector<std::string>::iterator it_groups = groups.begin(); it_groups != groups.end(); it_groups++){
      if(!read.compare(*it_groups)){
        group_of_haplotype.push_back(ind);
        if(diploid) group_of_haplotype.push_back(ind);
        break;
      }
      ind++;
    }

  }
  is.close();

  //calculate group sizes
  group_sizes.resize(groups.size());
  std::fill(group_sizes.begin(), group_sizes.end(), 0);
  for(std::vector<int>::iterator it_group_of_haplotype = group_of_haplotype.begin(); it_group_of_haplotype != group_of_haplotype.end(); it_group_of_haplotype++){
    group_sizes[*it_group_of_haplotype]++;
  }

}


std::string
Sample::AssignPopOfInterest(const std::string& s_pops){

  //read group_of_interest
  std::string label;
  group_of_interest.clear();

  if(s_pops != "All"){

    std::string pop;
    int i = 0;
    while(i < s_pops.size()){
      pop.clear();
      while(s_pops[i] != ','){
        pop += s_pops[i];
        i++;
        if(i == s_pops.size()) break;
      } 
      i++;

      bool exists = false;
      int index = 0;
      for(std::vector<std::string>::iterator it_groups = groups.begin(); it_groups != groups.end(); it_groups++){
        if(!(*it_groups).compare(pop)){
          exists = true;
          break;
        }
        index++;
      }
      if(!exists){
        std::cerr << "Group label does not exist." << std::endl;
        exit(1);
      }

      exists = false;
      for(std::vector<int>::iterator it_index = group_of_interest.begin(); it_index != group_of_interest.end(); it_index++){
        if(index == *it_index){
          exists = true;
          break;
        }
      }

      if(!exists){
        group_of_interest.push_back(index);
        label += pop;
      }
    }
    std::sort(group_of_interest.begin(), group_of_interest.end());
    
  }else{
    group_of_interest.resize(groups.size());
    for(int i = 0; i < (int) group_of_interest.size(); i++){
      group_of_interest[i] = i;
    }
    label = "All";
  }
  std::sort(group_of_interest.begin(), group_of_interest.end());

  group_of_interest_size = 0;
  for(std::vector<int>::iterator it_group_of_interest = group_of_interest.begin(); it_group_of_interest != group_of_interest.end(); it_group_of_interest++){
    group_of_interest_size += group_sizes[*it_group_of_interest];
  }

  return label;

}

