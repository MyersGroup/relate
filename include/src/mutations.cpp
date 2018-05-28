#include "mutations.hpp"


///////////////////////////////////////////////////

Mutations::Mutations(Data& data){
  N       = data.N;
  L       = data.L;

  info.resize(L);
  num_flips = 0; 
  num_notmappingmutations = 0; 
}

void
Mutations::Init(Data& data){
  N       = data.N;
  L       = data.L;

  info.resize(L);
  num_flips = 0; 
  num_notmappingmutations = 0; 
}

//////////////////// Public Members //////////

void
Mutations::GetAge(AncesTree& anc){

  CorrTrees::iterator it_seq = anc.seq.begin();
  int count_tree = 0;
  for(std::vector<SNPInfo>::iterator it = info.begin(); it != info.end(); it++){

    if((*it).tree > count_tree){
      count_tree++;
      it_seq++;
      assert((*it).tree == count_tree);
    }

    //find the branch on the tree and get age
    if((*it).branch.size() == 1){
      Node n = (*it_seq).tree.nodes[*(*it).branch.begin()];
      //traverse this node down to the bottom
      (*it).age_begin = 0.0;
      (*it).age_end = n.branch_length;
      while(n.child_left != NULL){
        n = *n.child_left;
        (*it).age_begin += n.branch_length;
      }
      (*it).age_end += (*it).age_begin;
    }
  }

}

void
Mutations::ReadShortFormat(const std::vector<std::string>& filenames){

  info.resize(L);
  std::string line;
  int snp = 0;

  int add_tree_index  = 0;
  num_flips           = 0;
  num_notmappingmutations = 0;

  for(int i = 0; i < (int)filenames.size(); i++){

    igzstream is(filenames[i]);
    if(is.fail()) is.open(filenames[i] + ".gz");
    if(is.fail()){
      std::cerr << "Error while reading " << filenames[i] << "(.gz)."  << std::endl;
      exit(1);
    }

    std::getline(is, line);

    while(std::getline(is, line)){

      //structure:
      //tree index;branches;flipped;age_begin;age_end;
      int i = 0;
      std::string inread;
      while(line[i] != ';'){
        inread += line[i];
        i++;
      }
      info[snp].tree = std::stoi(inread) + add_tree_index;
      inread.clear();

      do{
        i++;
        while(line[i] != ' ' && line[i] != ';'){
          inread += line[i];
          i++;
        }
        if(inread.size() > 0){ 
          info[snp].branch.push_back(std::stoi(inread));
          inread.clear();
        }
      }while(line[i] != ';');

      i++;
      while(line[i] != ';') i++;
      i++;
      while(line[i] != ';'){
        inread += line[i];
        i++;
      } 
      info[snp].flipped = std::stoi(inread);
      inread.clear();
      i++;

      while(line[i] != ';'){
        inread += line[i];
        i++;
      } 
      info[snp].age_begin = std::stof(inread);
      inread.clear();
      i++;

      while(line[i] != ';'){
        inread += line[i];
        i++;
      } 
      info[snp].age_end = std::stof(inread);
      inread.clear();

      if(info[snp].flipped) num_flips++;
      if(info[snp].branch.size() > 1) num_notmappingmutations++;
      snp++;

    }

    is.close();
    add_tree_index = info[snp-1].tree + 1;
  }

}

void
Mutations::Read(igzstream& is){

  info.clear();
  info.resize(L);
  std::string line;
  int snp = 0;

  num_flips           = 0;
  num_notmappingmutations = 0;

  std::string inread;
  
  //file structure:
  //snp;pos_of_snp;rs-id;tree_index;branch_indices;is_mapping;is_flipped;(age_begin,age_end);ancestral_allele/alternative_allele;downstream_allele;upstream_allele;ACB;ASW;BEB;CDX;CEU;CHB;CHS;CLM;ESN;FIN;GBR;GIH;GWD;IBS;ITU;JPT;KHV;LWK;MSL;MXL;PEL;PJL;PUR;STU;TSI;YRI

   int i = 0;
   while(getline(is, line)){ 

    i = 0;

    //snp-id
    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].snp_id = std::stoi(inread);

    //pos_of_snp
    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].pos = std::stoi(inread);

    //dist_to_next_snp
    
    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].dist = std::stoi(inread);
    
    //rs-id
    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    i++;
    info[snp].rs_id = inread;

    inread.clear();
    while(line[i] != ';'){
      inread += line[i];
      i++;
    }
    info[snp].tree = std::stoi(inread);
    inread.clear();

    do{
      i++;
      while(line[i] != ' ' && line[i] != ';'){
        inread += line[i];
        i++;
      }
      if(inread.size() > 0){ 
        info[snp].branch.push_back(std::stoi(inread));
        inread.clear();
      }
    }while(line[i] != ';');

    i++;
    while(line[i] != ';') i++;
    i++;
    while(line[i] != ';'){
      inread += line[i];
      i++;
    } 
    info[snp].flipped = std::stoi(inread);
    inread.clear();
    i++;

    while(line[i] != ';'){
      inread += line[i];
      i++;
    } 
    info[snp].age_begin = std::stof(inread);
    inread.clear();
    i++;

    while(line[i] != ';'){
      inread += line[i];
      i++;
    } 
    info[snp].age_end = std::stof(inread);
    inread.clear();
    i++;

    if(i < line.size()){

      //allele type
      while(line[i] != ';' && i < line.size()){
        inread += line[i];
        i++;
      } 
      info[snp].mutation_type = inread;
      i++;
      if(i < line.size()){
        inread.clear();
        while(line[i] != ';'){
          inread += line[i];
          i++;
        } 
        info[snp].upstream_base = inread;
        i++;
        inread.clear();
        while(line[i] != ';'){
          inread += line[i];
          i++;
        } 
        info[snp].downstream_base = inread;
        i++;
        inread.clear();

        while(i < line.size()){ 
          inread.clear();
          while(line[i] != ';' && i < line.size()){
            inread += line[i];
            i++;
          } 
          i++;
          info[snp].freq.push_back(std::stoi(inread));
        }
      }

    }

    if(info[snp].flipped) num_flips++;
    if(info[snp].branch.size() > 1) num_notmappingmutations++;
    snp++;

  }

  info.resize(snp);
  is.close();

}

void
Mutations::Read(const std::string& filename){

  std::string line;
  igzstream is(filename);
  if(is.fail()) is.open(filename + ".gz");
  if(is.fail()){
    std::cerr << "Error while reading " << filename << "(.gz)." << std::endl;
    exit(1);
  }
  std::getline(is, header); 
  
  std::string unused;
  L = 0;
  while( std::getline(is, unused) ){
    ++L;
  }
  is.close();

  is.open(filename);
  if(is.fail()) is.open(filename + ".gz");
  std::getline(is, line);
  Read(is);

}

void 
Mutations::Dump(const std::string& filename){

  std::ofstream os;
  os.open(filename);

  if(os.fail()){
    std::cerr << "Error while writing to " << filename << "." << std::endl; 
  }else{

    if(header.size() > 0){
      os << header;
    }else{
      os << "snp;pos_of_snp;dist;rs-id;tree_index;branch_indices;is_not_mapping;is_flipped;age_begin;age_end;ancestral_allele/alternative_allele;upstream_allele;downstream_allele;"; 
    }
    os << "\n";

    for(std::vector<SNPInfo>::iterator it = info.begin(); it != info.end(); it++){
      os << (*it).snp_id << ";" << (*it).pos << ";" << (*it).dist << ";" << (*it).rs_id << ";" << (*it).tree << ";";
      std::deque<int>::iterator it_branch = (*it).branch.begin();
      if((*it).branch.size() > 0){
        os << *it_branch;
        it_branch++;
      }
      for(; it_branch != (*it).branch.end(); it_branch++){
        os << " " << *it_branch;
      }
      if((*it).branch.size() > 1){
        os << ";1;"; 
      }else{
        os << ";0;";
      }
      os << (*it).flipped << ";" << (*it).age_begin << ";" << (*it).age_end << ";";
      os << (*it).mutation_type << ";";

      if((*it).freq.size() > 0){
        os << (*it).upstream_base << ";" << (*it).downstream_base << ";";
        for(std::vector<int>::iterator it_freq = (*it).freq.begin(); it_freq != (*it).freq.end(); it_freq++){
          os << *it_freq << ";";
        }

      }

      os << "\n";

    }

    os.close();

  }

}


void 
Mutations::DumpShortFormat(const std::string& filename){

  std::ofstream os;
  os.open(filename);

  if(os.fail()){
    std::cerr << "Error while writing to " << filename << "." << std::endl; 
  }else{

    os << "tree_index;branch_index;is_mapping;is_flipped;age_of_mutation" << "\n";

    for(std::vector<SNPInfo>::iterator it = info.begin(); it != info.end(); it++){
      os << (*it).tree << ";";
      std::deque<int>::iterator it_branch = (*it).branch.begin();
      if((*it).branch.size() > 0){
        os << *it_branch;
        it_branch++;
      }
      for(; it_branch != (*it).branch.end(); it_branch++){
        os << " " << *it_branch;
      }
      if((*it).branch.size() > 1){
        os << ";1;"; 
      }else{
        os << ";0;";
      }
      os << (*it).flipped << ";" << (*it).age_begin << ";" << (*it).age_end << ";" << "\n";
    }

    os.close();

  }

}

void 
Mutations::DumpShortFormat(const std::string& filename, const int section_startpos, const int section_endpos){

  std::ofstream os;
  os.open(filename);

  if(os.fail()){
    std::cerr << "Error while writing to " << filename << "." << std::endl; 
  }else{

    os << "tree_index;branch_index;is_mapping;is_flipped;age_of_mutation" << "\n";

    for(std::vector<SNPInfo>::iterator it = std::next(info.begin(),section_startpos); it != std::next(info.begin(), section_endpos+1); it++){
      os << (*it).tree << ";";
      std::deque<int>::iterator it_branch = (*it).branch.begin();
      if((*it).branch.size() > 0){
        os << *it_branch;
        it_branch++;
      }
      for(; it_branch != (*it).branch.end(); it_branch++){
        os << " " << *it_branch;
      }
      if((*it).branch.size() > 1){
        os << ";1;"; 
      }else{
        os << ";0;";
      }
      os << (*it).flipped << ";" << (*it).age_begin << ";" << (*it).age_end << ";" << "\n";
    }

    os.close();

  }

}


