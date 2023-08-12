//Build topology for trees in regions of 500 SNPs.
// Input: .anc and .mut file
// Output: .anc and .mut file 

//Have it as a final step.
//Have another function that works after FindEquivalentBranches.

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"


void
Relabel(Tree& tree, int N){

  bool relabel = true;

  //currently not most efficient
  Node dummy, dummy2;
  int parent;
  while(relabel){
    relabel = false;
    for(int i = N; i < tree.nodes.size()-1; i++){

      if( i > (*tree.nodes[i].parent).label ){

        dummy = tree.nodes[i];
        parent = (*tree.nodes[i].parent).label;

        tree.nodes[i] = tree.nodes[parent];
        tree.nodes[parent] = dummy;
        tree.nodes[i].label = i;
        tree.nodes[parent].label = parent;

        if((*tree.nodes[i].child_left).label == i) tree.nodes[i].child_left = &tree.nodes[parent];
        if((*tree.nodes[i].child_right).label == i) tree.nodes[i].child_right = &tree.nodes[parent];
        if((*tree.nodes[parent].parent).label == parent) tree.nodes[parent].parent = &tree.nodes[i];

        (*tree.nodes[i].child_left).parent = &tree.nodes[i]; 
        (*tree.nodes[i].child_right).parent = &tree.nodes[i];
        if( (*(*tree.nodes[i].parent).child_left).label == parent ) (*tree.nodes[i].parent).child_left = &tree.nodes[i];
        if( (*(*tree.nodes[i].parent).child_right).label == parent ) (*tree.nodes[i].parent).child_right = &tree.nodes[i];

        (*tree.nodes[parent].child_left).parent = &tree.nodes[parent]; 
        (*tree.nodes[parent].child_right).parent = &tree.nodes[parent]; 
        if( (*(*tree.nodes[parent].parent).child_left).label == i ) (*tree.nodes[parent].parent).child_left = &tree.nodes[parent];
        if( (*(*tree.nodes[parent].parent).child_right).label == i ) (*tree.nodes[parent].parent).child_right = &tree.nodes[parent];
        relabel = true;

        assert((*(*tree.nodes[parent].child_left).parent).label == parent);
        assert((*(*tree.nodes[parent].child_right).parent).label == parent);
        assert((*(*tree.nodes[i].child_left).parent).label == i);
        assert((*(*tree.nodes[i].child_right).parent).label == i);
        assert((*(*tree.nodes[i].parent).child_left).label == i || (*(*tree.nodes[i].parent).child_right).label == i);
        //std::cerr << i << " " << parent << " " << (*tree.nodes[parent].parent).label << " " << (*(*tree.nodes[parent].parent).child_left).label << " " << (*(*tree.nodes[parent].parent).child_right).label << std::endl;
        assert((*(*tree.nodes[parent].parent).child_left).label == parent || (*(*tree.nodes[parent].parent).child_right).label == parent);


      }

    }
  }

}


int
Map(std::vector<Leaves>& desc, int k, CollapsedMatrix<char>& seq, int DAF, std::vector<int>& nodes, int num_desc, int thr){

  //check if seq[k] maps onto any branch in desc

  std::vector<char>::iterator it_seq = seq.rowbegin(k);

  int matching = 0, non_matching = 0; 
  if(thr <= 1 || DAF < 4){
    //has to map exactly
    if(DAF == num_desc){
      for(std::vector<int>::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++){
        for(std::vector<int>::iterator it_mem = desc[*it_n].member.begin(); it_mem != desc[*it_n].member.end(); it_mem++){
          if(*std::next(it_seq,*it_mem) != '1'){
            return(thr);
          }else{
            matching++;
          }
        }
      }
    }else{
      return(thr);
    }
    return(0);
  }else{
    for(std::vector<int>::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++){
      for(std::vector<int>::iterator it_mem = desc[*it_n].member.begin(); it_mem != desc[*it_n].member.end(); it_mem++){
        if(*std::next(it_seq,*it_mem) != '1'){
          non_matching++;
          if(non_matching >= thr) return(thr);
        }else{
          matching++;
        }
      }
    }
    assert(non_matching < thr);
    if(DAF == matching && non_matching == 0) return(0);

    if(DAF - matching + non_matching >= thr) return(thr);
    if(matching <= 0.7*DAF) return(thr); //! matching/DAF > 0.7;
    if(non_matching >= 0.3*(seq.size() - DAF)) return(thr); //!non_matching/(N-DAF) < 0.3;
  }

  //std::cerr << "debug3" << std::endl;
  if(matching > 0.7*num_desc && (seq.size()-DAF-non_matching) > 0.7 * (seq.size() - num_desc)){
    return(DAF - matching + non_matching);
  }

  //int correct_carriers = matching;
  //int incorrect_carriers = DAF - matching;
  //int correct_noncarriers = seq.size()-DAF-non_matching;
  //int incorrect_noncarriers = non_matching;
  //assert(incorrect_carriers >= 0);
  //assert(correct_noncarriers >= 0);

  //if(matching/(matching + non_matching) > 0.7 && correct_noncarriers/(incorrect_carriers + correct_noncarriers) > 0.7){
  //	return(incorrect_carriers + incorrect_noncarriers);
  //}

  return(thr);

}


int
Map(std::vector<Leaves>& desc, std::vector<char>& seq, int DAF, std::vector<int>& nodes, int thr){

  int num_leaves = 0;
  for(std::vector<int>::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++){
    num_leaves += desc[(*it_n)].num_leaves;
  }


  int diff = DAF - num_leaves;
  if(diff < 0) diff *= -1;
  if( diff < thr){
    int matching = 0, non_matching = 0; 
    if(thr <= 1 || DAF < 4){
      int ret = thr;
      if(DAF < 4) thr = 1;
      for(std::vector<int>::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++){
        for(std::vector<int>::iterator it_mem = desc[*it_n].member.begin(); it_mem != desc[*it_n].member.end(); it_mem++){
          if(seq[*it_mem] != '1'){
            return(ret);
          }else{
            matching++;
          }
        }
      }
      if(DAF == matching) return(0);
      return(ret);
    }else{
      for(std::vector<int>::iterator it_n = nodes.begin(); it_n != nodes.end(); it_n++){
        for(std::vector<int>::iterator it_mem = desc[*it_n].member.begin(); it_mem != desc[*it_n].member.end(); it_mem++){
          if(seq[*it_mem] != '1'){
            non_matching++;
            if(non_matching >= thr) return(thr);
          }else{
            matching++;
          }
        }
      }
      assert(non_matching < thr);
      if(DAF == matching && non_matching == 0) return(0);

      if(DAF - matching + non_matching >= thr) return(thr);
      if(matching <= 0.7*DAF) return(thr); //! matching/DAF > 0.7;
      if(non_matching >= 0.3*(seq.size() - DAF)) return(thr); //!non_matching/(N-DAF) < 0.3;
    }

    //std::cerr << "debug3" << std::endl;
    if(matching > 0.7*num_leaves && (seq.size()-DAF-non_matching) > 0.7 * (seq.size() - num_leaves)){
      return(DAF - matching + non_matching);
    }
    //int correct_carriers = matching;
    //int incorrect_carriers = DAF - matching;
    //int correct_noncarriers = seq.size()-DAF-non_matching;
    //int incorrect_noncarriers = non_matching;
    //assert(incorrect_carriers >= 0);
    //assert(correct_noncarriers >= 0);

    //if(matching/(matching + non_matching) > 0.7 && correct_noncarriers/(incorrect_carriers + correct_noncarriers) > 0.7){
    //	return(incorrect_carriers + incorrect_noncarriers);
    //}

  }

  return(thr);

}

bool
CheckBranch(Data& data, std::vector<Leaves>& desc, std::vector<int>& DAF, int thr, int node1, int node2, int node3, int& closest_event12, int& closest_event13, int& closest_event23, int dist, int k, bool approx = true){

  int score12, score13, score23;
  std::vector<int> nodes(2);
  bool mapped;
  int num_desc;
  int threshold = 1e6;

  mapped = false;
  nodes[0] = node1;
  nodes[1] = node2;

  num_desc = desc[node1].num_leaves + desc[node2].num_leaves;
  if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
    score12 = Map(desc, k, data.sequence, DAF[k], nodes, num_desc, thr);
  }else{
    score12 = thr;
  }
  if(score12 == 0){
    mapped = true;
    if( dist < closest_event12 ) closest_event12 = dist;
  }else{
    nodes[1] = node3;

    num_desc = desc[node1].num_leaves + desc[node3].num_leaves;
    if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
      score13 = Map(desc, k, data.sequence, DAF[k], nodes, num_desc, thr);
    }else{
      score13 = thr;
    }
    if(score13 == 0){
      mapped = true;
      if( dist < closest_event13 ) closest_event13 = dist;
    }else{
      nodes[0] = node2;
      nodes[1] = node3;
      num_desc = desc[node2].num_leaves + desc[node3].num_leaves;
      if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
        score23 = Map(desc, k, data.sequence, DAF[k], nodes, num_desc, thr);
      }else{
        score23 = thr;
      }
      if(score23 == 0){
        mapped = true;
        if( dist < closest_event23 ) closest_event23 = dist;
      }
    }
  }

  if(approx && !mapped && thr > 1){

    if(closest_event12 > threshold && closest_event13 > threshold && closest_event23 > threshold){

      dist += threshold;
      if(score12 < thr || score13 < thr || score23 < thr){

        int min = thr;
        nodes = {node1};
        num_desc = desc[node1].num_leaves;
        if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
          min  = std::min(min, Map(desc, k, data.sequence, DAF[k], nodes, num_desc, thr));
        }
        nodes[0] = node2;
        num_desc = desc[node2].num_leaves;
        if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
          min  = std::min(min, Map(desc, k, data.sequence, DAF[k], nodes, num_desc, thr));
        }
        nodes[0] = node3;
        num_desc = desc[node3].num_leaves;
        if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
          min  = std::min(min, Map(desc, k, data.sequence, DAF[k], nodes, num_desc, thr));
        }
        nodes = {node1,node2,node3};
        num_desc = desc[node1].num_leaves + desc[node2].num_leaves + desc[node3].num_leaves;
        if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
          min  = std::min(min, Map(desc, k, data.sequence, DAF[k], nodes, num_desc, thr));
        }
        nodes.resize(2);

        //map to the branch that fits best.
        if(score12 < score13 && score12 < score23 && score12 < min){
          if( dist < closest_event12 ){
            closest_event12 = dist;
            mapped = true;
          }
        }else if(score13 < score12 && score13 < score23 && score13 < min){
          if( dist < closest_event13 ){
            closest_event13 = dist;
            mapped = true;
          }
        }else if(score23 < score12 && score23 < score13 && score23 < min){
          if( dist < closest_event23 ){
            closest_event23 = dist;
            mapped = true;
          }
        }

      }

    }

  }

  return(mapped);

}

void
PostProcess(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if(!options.count("input") || !options.count("haps") || !options.count("sample") || !options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, haps, sample, output, Optional: map, randomise, transversion." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Post process anc/mut files." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Post-processing " << options["input"].as<std::string>() << "..." << std::endl;

  bool randomise = false;
  if(options.count("randomise")){
    randomise = true;
  }

  bool use_transitions = true;
  if(options.count("transversion")) use_transitions = false;


  std::mt19937 rng;
  int seed = 1;
  if(options.count("seed")){
    seed = options["seed"].as<int>();
  }
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0.0,1.0);

  AncesTree anc; //anc.seq is a list of MarginalTrees (type CorrTrees)
  Mutations mut; //mut.info is a vector of SNPInfo (type Muts)

  anc.Read(options["input"].as<std::string>() + ".anc");
  mut.Read(options["input"].as<std::string>() + ".mut");

  int N_total = 2*anc.N-1;
  int N = anc.N;
  int root = 2*anc.N-2;

  std::vector<float> rdist(mut.info.size());
  float threshold = 10e6, threshold2 = 1e6;
  if(options.count("map")){

    map map(options["map"].as<std::string>().c_str());

    int i = 0;
    float r = 0;

    threshold = 10, threshold2 = 1;
    Muts::iterator it_mut = mut.info.begin();
    for(std::vector<float>::iterator it_rdist = rdist.begin(); it_rdist != rdist.end(); it_rdist++){

      if(i < map.bp.size()){
        if(i == 0 && map.bp[i] > (*it_mut).pos){
          r = map.gen_pos[i]/map.bp[i]*(*it_mut).pos;
        }else{
          while(map.bp[i] < (*it_mut).pos){
            i++;
            if(i == map.bp.size()){
              break;
            }
          }
          assert(i > 0);
          assert((*it_mut).pos <= map.bp[i]);
          assert((*it_mut).pos >= map.bp[i-1]);
          if(i < map.bp.size()){
            r = (map.gen_pos[i] - map.gen_pos[i-1])/(map.bp[i] - map.bp[i-1])*((*it_mut).pos - map.bp[i-1]) + map.gen_pos[i-1];
          }else{
            r = (map.gen_pos[i-1] - map.gen_pos[i-2])/(map.bp[i-1] - map.bp[i-2])*((*it_mut).pos - map.bp[i-2]) + map.gen_pos[i-1];
          }
        }
      }else{
        r = (map.gen_pos[i-1] - map.gen_pos[i-2])/(map.bp[i-1] - map.bp[i-2])*((*it_mut).pos - map.bp[i-2]) + map.gen_pos[i-1];
      }

      *it_rdist = r;
      //std::cerr << (*it_mut).pos << " " << *it_rdist << " " << *it_rdist/(*it_mut).pos*1e6 << std::endl;
      it_mut++;
    }

  }else{

    std::cerr << "Using physical distance." << std::endl;
    Muts::iterator it_mut = mut.info.begin();
    for(std::vector<float>::iterator it_rdist = rdist.begin(); it_rdist != rdist.end(); it_rdist++){
      *it_rdist = (*it_mut).pos;
      it_mut++;
    }

  }

  haps m_haps(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
  int buf = 5000;
  std::vector<std::vector<char>> seq(buf);
  std::vector<int> bps(buf), index(buf), DAF(buf,0);
  int L = m_haps.GetL();

  if(L != mut.info.size()){
    std::cerr << "Error: Haps file is likely not the one used to infer tree" << std::endl;
    exit(1);
  }

  int snp = 0, snp2 = 0;
  while(snp2 < buf && snp < L){
    seq[snp2].resize(N);
    m_haps.ReadSNP(seq[snp2], bps[snp2]);
    index[snp2] = snp;
    for(std::vector<char>::iterator it_seq = seq[snp2].begin(); it_seq != seq[snp2].end(); it_seq++){
      if(*it_seq == '1') DAF[snp2]++;
    }
    if(DAF[snp2] > 1) snp2++;
    //snp2++;
    snp++;
  }

  std::vector<Leaves> desc;
  std::vector<float> coords;

  float val;
  std::vector<int> nodes(2);
  std::vector<int> remaining;

  int thr = (int) (0.03 * N) + 1;
  //thr = 1.0;

  CorrTrees::iterator it_seq;
  Muts::iterator it_mut = mut.info.begin(); //iterator for mut file
  int tree_index = 0;
  int bp_init = rdist[mut.info.size()-1];
  float tree_bp, dist;
  int score12, score13, score23, tmp;
  int num_desc;
  bool mapped;

  for(it_seq = anc.seq.begin(); it_seq != anc.seq.end(); it_seq++){

    if(snp < L && it_mut != mut.info.end()){

      if( (*it_seq).pos - index[snp2 % buf] > buf/2){
        tmp = index[snp2 % buf];
        while( (*it_seq).pos - tmp > buf/2 ){
          if(snp % (int)(L/100.0) == 0){
            std::cerr << "[" << (int)(snp/((int)(L/100.0))) << "%]" << "\r";
          }
          m_haps.ReadSNP(seq[snp2 % buf], bps[snp2 % buf]);
          index[snp2 % buf] = snp;
          DAF[snp2 % buf] = 0;
          for(std::vector<char>::iterator it_seq = seq[snp2 % buf].begin(); it_seq != seq[snp2 % buf].end(); it_seq++){
            if(*it_seq == '1') DAF[snp2 % buf]++;
          }
          if(DAF[snp2 % buf] > 1){
            snp2++;
            tmp = index[snp2 % buf];
          }
          snp++;
          it_mut++;
          if(snp == L) break;
          if(it_mut == mut.info.end()) break;
        }
      }
      //assert( index[snp2 % buf] <= (*it_seq).pos && index[(snp2-1) % buf] >= (*it_seq).pos );

    }

    (*it_seq).tree.FindAllLeaves(desc);
    (*it_seq).tree.GetCoordinates(coords);
    tree_bp = rdist[(*it_seq).pos];
    for(int iter = 0; iter < 5; iter++){

      bool is_updated = false;

      if(0){
        std::vector<int> nodes;
        if(iter == 0){
          for(int i = N; i < root; i++){
            nodes.push_back(i);
          }
        }else{
          for(int i = root-1; i >= N; i--){
            nodes.push_back(i);
          }
        }
      }

      for(int i = root-1; i >= N; i--){

        int node1 = (*(*it_seq).tree.nodes[i].child_left).label;
        int node2 = (*(*it_seq).tree.nodes[i].child_right).label;

        int parent = (*(*it_seq).tree.nodes[i].parent).label;
        int node3 = (*(*it_seq).tree.nodes[parent].child_left).label;
        if(node3 == i){
          node3 = (*(*it_seq).tree.nodes[parent].child_right).label; 
        }

        assert(node1 != i);
        assert(node2 != i);

        if( (*it_seq).tree.nodes[i].num_events < 1.0 ){

          //check if there is evidence for node1+node2, node1 + node3, node2 + node3
          int num_events12 = 0, num_events13 = 0, num_events23 = 0;
          int closest_event12 = bp_init, closest_event13 = bp_init, closest_event23 = bp_init;

          for(int k = 0; k < buf; k++){

            //std::cerr << (*it_seq).pos << " " << index[k] << " " << rdist[index[k]] << " " << tree_bp << std::endl;
            dist = rdist[index[k]] - tree_bp;
            if(dist < 0) dist *= -1.0;

            if( dist < threshold){

              //  report.num_incorrect_carriers + report.num_incorrect_noncarriers
              //  thr = (int) (0.03 * N) + 1;
              //
              //  necessary conditions:
              //  (((float) report.num_incorrect_carriers)/total_carriers < 0.3);
              //  (((float) report.num_incorrect_noncarriers)/total_noncarriers < 0.3);
              //  (((float) report.num_correct_carriers)/(report.num_correct_carriers + report.num_incorrect_noncarriers) > 0.7
              //  (((float) report.num_correct_noncarriers)/(report.num_incorrect_carriers + report.num_correct_noncarriers) > 0.7)
              //
              //  total_carriers = DAF
              //  total_non_carriers = N-DAF
              //
              //  data:1, tree:1 (matching) - correct_carriers
              //  data:1, tree:0 (DAF - matching) - incorrect_carriers
              //  data:0, tree:1 (non_matching) - incorrect_noncarriers
              //  data:0, tree:0 (N-DAF-non_matching) - correct_noncarriers
              //

              mapped = false;
              nodes[0] = node1;
              nodes[1] = node2;
              num_desc = desc[node1].num_leaves + desc[node2].num_leaves;
              if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
                score12 = Map(desc, seq[k], DAF[k], nodes, thr);
              }else{
                score12 = thr;
              }
              if(score12 == 0){
                num_events12++;
                mapped = true;
                if( dist < closest_event12 ) closest_event12 = dist;
              }else{
                nodes[1] = node3;
                num_desc = desc[node1].num_leaves + desc[node3].num_leaves;
                if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
                  score13 = Map(desc, seq[k], DAF[k], nodes, thr);
                }else{
                  score13 = thr;
                }
                if(score13 == 0){
                  num_events13++;
                  mapped = true;
                  if( dist < closest_event13 ) closest_event13 = dist;
                }else{
                  nodes[0] = node2;
                  nodes[1] = node3;
                  num_desc = desc[node2].num_leaves + desc[node3].num_leaves;
                  if( DAF[k] - num_desc < thr && num_desc - DAF[k] < thr ){
                    score23 = Map(desc, seq[k], DAF[k], nodes, thr);
                  }else{
                    score23 = thr;
                  }
                  if(score23 == 0){
                    num_events23++;
                    mapped = true;
                    if( dist < closest_event23 ) closest_event23 = dist;
                  }
                }
              }

              if(!mapped && thr > 1){

                if(closest_event12 > threshold && closest_event13 > threshold && closest_event23 > threshold){

                  dist += threshold;
                  if(score12 < thr || score13 < thr || score23 < thr){

                    int min = thr;
                    nodes = {node1};
                    min  = std::min(min, Map(desc, seq[k], DAF[k], nodes, thr));
                    nodes[0] = node2;
                    min  = std::min(min, Map(desc, seq[k], DAF[k], nodes, thr));
                    nodes[0] = node3;
                    min  = std::min(min, Map(desc, seq[k], DAF[k], nodes, thr));
                    nodes = {node1,node2,node3};
                    min  = std::min(min, Map(desc, seq[k], DAF[k], nodes, thr));
                    nodes.resize(2);

                    //map to the branch that fits best.
                    if(score12 < score13 && score12 < score23 && score12 < min){
                      num_events12++;
                      if( dist < closest_event12 ) closest_event12 = dist;
                    }else if(score13 < score12 && score13 < score23 && score13 < min){
                      num_events13++;
                      if( dist < closest_event13 ) closest_event13 = dist;
                    }else if(score23 < score12 && score23 < score13 && score23 < min){
                      num_events23++;
                      if( dist < closest_event23 ) closest_event23 = dist;
                    }

                  }

                }

              }

            }
          }

          if( (closest_event13 < closest_event12 && closest_event13 <= closest_event23) || (closest_event13 <= closest_event12 && closest_event13 < closest_event23) ){

            is_updated = true;

            (*it_seq).tree.nodes[i].child_left       = &(*it_seq).tree.nodes[node1];
            (*it_seq).tree.nodes[i].child_right      = &(*it_seq).tree.nodes[node3];
            (*it_seq).tree.nodes[node1].parent       = &(*it_seq).tree.nodes[i];
            (*it_seq).tree.nodes[node3].parent       = &(*it_seq).tree.nodes[i];

            (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[i];
            (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[node2];
            (*it_seq).tree.nodes[i].parent           = &(*it_seq).tree.nodes[parent];
            (*it_seq).tree.nodes[node2].parent       = &(*it_seq).tree.nodes[parent];

            if(coords[node3] >= coords[i]) coords[i] = (coords[parent] + coords[node3])/2.0;

            (*it_seq).tree.nodes[i].num_events    = 1.0;

            (*it_seq).tree.nodes[node1].branch_length = coords[i] - coords[node1];
            (*it_seq).tree.nodes[node3].branch_length = coords[i] - coords[node3];
            (*it_seq).tree.nodes[node2].branch_length = coords[parent] - coords[node2];
            (*it_seq).tree.nodes[i].branch_length     = coords[parent] - coords[i];

            //update desc[i]
            desc[i].member = desc[node1].member;
            desc[i].member.insert(desc[i].member.end(),desc[node3].member.begin(), desc[node3].member.end());
            desc[i].num_leaves = desc[node1].num_leaves + desc[node3].num_leaves;
            std::sort(desc[i].member.begin(), desc[i].member.end());

          }else if( (closest_event23 < closest_event12 && closest_event23 <= closest_event13) || (closest_event23 <= closest_event12 && closest_event23 < closest_event13) ){

            is_updated = true;

            (*it_seq).tree.nodes[i].child_left       = &(*it_seq).tree.nodes[node2];
            (*it_seq).tree.nodes[i].child_right      = &(*it_seq).tree.nodes[node3];
            (*it_seq).tree.nodes[node2].parent       = &(*it_seq).tree.nodes[i];
            (*it_seq).tree.nodes[node3].parent       = &(*it_seq).tree.nodes[i];

            (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[i];
            (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[node1];
            (*it_seq).tree.nodes[i].parent           = &(*it_seq).tree.nodes[parent];
            (*it_seq).tree.nodes[node1].parent       = &(*it_seq).tree.nodes[parent];

            if(coords[node3] >= coords[i]) coords[i] = (coords[parent] + coords[node3])/2.0;

            (*it_seq).tree.nodes[i].num_events    = 1.0;

            (*it_seq).tree.nodes[node2].branch_length = coords[i] - coords[node2];
            (*it_seq).tree.nodes[node3].branch_length = coords[i] - coords[node3];
            (*it_seq).tree.nodes[node1].branch_length = coords[parent] - coords[node1];
            (*it_seq).tree.nodes[i].branch_length     = coords[parent] - coords[i];

            desc[i].member = desc[node2].member;
            desc[i].member.insert(desc[i].member.end(),desc[node3].member.begin(), desc[node3].member.end());
            desc[i].num_leaves = desc[node2].num_leaves + desc[node3].num_leaves;
            std::sort(desc[i].member.begin(), desc[i].member.end());

          }else if( (closest_event12 < closest_event23 && closest_event12 <= closest_event13) || (closest_event12 <= closest_event23 && closest_event12 < closest_event13) ){
            (*it_seq).tree.nodes[i].num_events    = 1.0;
          }
        }
      }
      if(!is_updated){
        break;
      }
    }

    if(randomise){

      for(int i = root-1; i >= N; i--){

        int node1 = i;
        int parent = (*(*it_seq).tree.nodes[i].parent).label;
        int node2 = (*(*it_seq).tree.nodes[parent].child_left).label;
        if(node2 == i){
          node2 = (*(*it_seq).tree.nodes[parent].child_right).label; 
        }

        if( (*it_seq).tree.nodes[node1].num_events < 1.0 ){
          if( (*it_seq).tree.nodes[node2].num_events < 1.0 || (*it_seq).tree.nodes[parent].num_events < 1.0 ){

            int child1 = (*(*it_seq).tree.nodes[node1].child_left).label;
            int child2 = (*(*it_seq).tree.nodes[node1].child_right).label;
            int child3, child4;

            remaining = {child1, child2, node2, -1};

            bool shuffle_four = false;

            if( (*it_seq).tree.nodes[node2].num_events < 1.0 && (*it_seq).tree.nodes[node2].child_left != NULL){
              child3 = (*(*it_seq).tree.nodes[node2].child_left).label;
              child4 = (*(*it_seq).tree.nodes[node2].child_right).label;
              remaining = {child1, child2, child3, child4};
              shuffle_four = true;
            }else if( (*it_seq).tree.nodes[parent].num_events < 1.0 && parent != root ){
              if(0){
                child3 = node2;
                node2  = parent;
                parent = (*(*it_seq).tree.nodes[node2].parent).label;
                child4 = (*(*it_seq).tree.nodes[parent].child_left).label;
                if(child4 == node2){
                  child4 = (*(*it_seq).tree.nodes[parent].child_right).label;
                }
                remaining = {child1, child2, child3, child4};
                shuffle_four = true;
              }
            }

            if(shuffle_four){

              if(shuffle_four){

                if(coords[child1] >= coords[node1]) coords[node1] = (coords[parent] + coords[child1])/2.0;
                if(coords[child2] >= coords[node1]) coords[node1] = (coords[parent] + coords[child2])/2.0;
                if(coords[child3] >= coords[node1]) coords[node1] = (coords[parent] + coords[child3])/2.0;
                if(coords[child4] >= coords[node1]) coords[node1] = (coords[parent] + coords[child4])/2.0;
                if(coords[child1] >= coords[node2]) coords[node2] = (coords[parent] + coords[child1])/2.0;
                if(coords[child2] >= coords[node2]) coords[node2] = (coords[parent] + coords[child2])/2.0;
                if(coords[child3] >= coords[node2]) coords[node2] = (coords[parent] + coords[child3])/2.0;
                if(coords[child4] >= coords[node2]) coords[node2] = (coords[parent] + coords[child4])/2.0;

                //make sure that node2 < node1 and coords[node2] < coords[node1]
                if(node2 > node1){
                  int foo = node1;
                  node1 = node2;
                  node2 = foo;
                }
                if(coords[node2] > coords[node1]){
                  float foo = coords[node1];
                  coords[node1] = coords[node2];
                  coords[node2] = foo;
                }

                val = dist_unif(rng);
                if(val < 1.0/6.0){
                  nodes[0] = child1;
                  nodes[1] = child2;
                  remaining[0] = node2;
                  remaining[1] = remaining[3];
                  remaining[3] = -1;
                }else if(val < 2.0/6.0){
                  nodes[0] = child1;
                  nodes[1] = child3;
                  remaining[0] = node2;
                  remaining[2] = remaining[3];
                  remaining[3] = -1;
                }else if(val < 3.0/6.0){
                  nodes[0] = child1;
                  nodes[1] = child4;
                  remaining[0] = node2;
                  remaining[3] = -1;
                }else if(val < 4.0/6.0){
                  nodes[0] = child2;
                  nodes[1] = child3;
                  remaining[1] = node2;
                  remaining[2] = remaining[3];
                  remaining[3] = -1;
                }else if(val < 5.0/6.0){
                  nodes[0] = child2;
                  nodes[1] = child4;
                  remaining[1] = node2;
                  remaining[3] = -1;
                }else{ 
                  nodes[0] = child3;
                  nodes[1] = child4;
                  remaining[2] = node2;
                  remaining[3] = -1;
                }

                (*it_seq).tree.nodes[node2].child_left       = &(*it_seq).tree.nodes[nodes[0]];
                (*it_seq).tree.nodes[node2].child_right      = &(*it_seq).tree.nodes[nodes[1]];
                (*it_seq).tree.nodes[nodes[0]].parent        = &(*it_seq).tree.nodes[node2];
                (*it_seq).tree.nodes[nodes[1]].parent        = &(*it_seq).tree.nodes[node2];
                (*it_seq).tree.nodes[nodes[0]].branch_length = coords[node2] - coords[nodes[0]];
                (*it_seq).tree.nodes[nodes[1]].branch_length = coords[node2] - coords[nodes[1]];
                assert((*it_seq).tree.nodes[nodes[0]].branch_length >= 0.0);
                assert((*it_seq).tree.nodes[nodes[1]].branch_length >= 0.0);
                //(*it_seq).tree.nodes[node2].num_events = 1.0;
              }else{
                if(coords[node2] >= coords[node1]) coords[node1] = (coords[parent] + coords[node2])/2.0;
              }

              val = dist_unif(rng);
              if(val < 1.0/3.0){
                nodes[0]     = remaining[0];
                nodes[1]     = remaining[1];
                remaining[0] = node1;
                remaining[1] = remaining[2];
                remaining[2] = -1;
              }else if(val < 2.0/3.0){
                nodes[0]     = remaining[0];
                nodes[1]     = remaining[2];
                remaining[0] = node1;
                remaining[2] = -1;
              }else{ 
                nodes[0]     = remaining[1];
                nodes[1]     = remaining[2];
                remaining[1] = node1;
                remaining[2] = -1;
              }

              if(coords[nodes[0]] >= coords[node1]) coords[node1] = (coords[parent] + coords[nodes[0]])/2.0;
              if(coords[nodes[1]] >= coords[node1]) coords[node1] = (coords[parent] + coords[nodes[1]])/2.0;

              (*it_seq).tree.nodes[node1].child_left       = &(*it_seq).tree.nodes[nodes[0]];
              (*it_seq).tree.nodes[node1].child_right      = &(*it_seq).tree.nodes[nodes[1]];
              (*it_seq).tree.nodes[nodes[0]].parent        = &(*it_seq).tree.nodes[node1];
              (*it_seq).tree.nodes[nodes[1]].parent        = &(*it_seq).tree.nodes[node1];
              (*it_seq).tree.nodes[nodes[0]].branch_length = coords[node1] - coords[nodes[0]];
              (*it_seq).tree.nodes[nodes[1]].branch_length = coords[node1] - coords[nodes[1]];
              assert((*it_seq).tree.nodes[nodes[0]].branch_length >= 0.0);
              assert((*it_seq).tree.nodes[nodes[1]].branch_length >= 0.0);
              //(*it_seq).tree.nodes[node1].num_events = 1.0;

              (*it_seq).tree.nodes[parent].child_left          = &(*it_seq).tree.nodes[remaining[0]];
              (*it_seq).tree.nodes[parent].child_right         = &(*it_seq).tree.nodes[remaining[1]];
              (*it_seq).tree.nodes[remaining[0]].parent        = &(*it_seq).tree.nodes[parent];
              (*it_seq).tree.nodes[remaining[1]].parent        = &(*it_seq).tree.nodes[parent];
              (*it_seq).tree.nodes[remaining[0]].branch_length = coords[parent] - coords[remaining[0]];
              (*it_seq).tree.nodes[remaining[1]].branch_length = coords[parent] - coords[remaining[1]];
              assert((*it_seq).tree.nodes[remaining[0]].branch_length >= 0.0);
              assert((*it_seq).tree.nodes[remaining[1]].branch_length >= 0.0);

            }

          }

        }
      }
    }

    Relabel((*it_seq).tree, N);

    for(int i = 0; i < N_total; i++){
      (*it_seq).tree.nodes[i].SNP_begin  = (*it_seq).pos;
      if(std::next(it_seq,1) != anc.seq.end()){
        (*it_seq).tree.nodes[i].SNP_end  = (*std::next(it_seq,1)).pos;
      }else{
        (*it_seq).tree.nodes[i].SNP_end  = L-1;
      }
      (*it_seq).tree.nodes[i].num_events = 0;
    }

    tree_index++;

  }

  //Map mutations
  haps mhaps(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str()); 
  Data data(mhaps.GetN(), mhaps.GetL());

  AncesTreeBuilder ab(data);

  Leaves sequences_carrying_mutation;
  sequences_carrying_mutation.member.resize(data.N);

  int bp;
  std::vector<char> sequence(data.N);

  it_seq = anc.seq.begin();
  (*it_seq).tree.GetCoordinates(coords);
  tree_index = 0;
  for(snp = 0; snp < data.L; snp++){

    mhaps.ReadSNP(sequence, bp);
    assert(bp == mut.info[snp].pos);

    //map mutation onto tree
    std::fill(sequences_carrying_mutation.member.begin(), sequences_carrying_mutation.member.end(), 0);
    sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
    for(int i = 0; i < data.N; i++){
      if(sequence[i] == '1'){
        sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
        sequences_carrying_mutation.num_leaves++;
      }else{
        sequences_carrying_mutation.member[i] = 0;
      }
    }

    if(mut.info[snp].tree > tree_index){
      it_seq++;
      (*it_seq).tree.GetCoordinates(coords);
      tree_index = mut.info[snp].tree;
    }

    mut.info[snp].branch.clear();
    if(sequences_carrying_mutation.num_leaves == data.N){
      (*it_seq).tree.nodes[root].num_events += 1.0;
      mut.info[snp].tree      = tree_index;
      mut.info[snp].branch.resize(1);
      mut.info[snp].branch[0] = root;
      mut.info[snp].age_begin = coords[root];
      mut.info[snp].age_end   = coords[root];
    }else{

      if(use_transitions){
        ab.IsSNPMapping((*it_seq).tree, sequences_carrying_mutation, snp);
      }else{
        if( (strcmp(mhaps.ancestral,"C") == 0 && strcmp(mhaps.alternative,"T") == 0) || (strcmp(mhaps.ancestral,"T") == 0 && strcmp(mhaps.alternative,"C") == 0) || 
            (strcmp(mhaps.ancestral,"G") == 0 && strcmp(mhaps.alternative,"A") == 0) || (strcmp(mhaps.ancestral,"A") == 0 && strcmp(mhaps.alternative,"G") == 0) ){
          ab.IsSNPMapping((*it_seq).tree, sequences_carrying_mutation, snp, 0);
        }else{
          ab.IsSNPMapping((*it_seq).tree, sequences_carrying_mutation, snp, 1);
        }
      }

      mut.info[snp].tree      = tree_index;
      mut.info[snp].branch    = ab.mutations.info[snp].branch;
      if(mut.info[snp].branch.size() == 1){	
        int branch = mut.info[snp].branch[0];
        if(branch < root){
          mut.info[snp].age_begin = coords[branch];
          mut.info[snp].age_end   = coords[(*(*it_seq).tree.nodes[branch].parent).label];
        }else{
          mut.info[snp].age_begin = coords[branch];
          mut.info[snp].age_end   = coords[branch];
        }
      }else{
        mut.info[snp].age_begin = 0.0;
        mut.info[snp].age_end   = 0.0;
      }
    }
    mut.info[snp].flipped = ab.mutations.info[snp].flipped;
  }

  //Associate equivalent branches (defined in anc.hpp) 
  anc.AssociateEquivalentBranches();

  //anc for dumping these files:
  anc.Dump(options["output"].as<std::string>() + ".anc");
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

}

int PostProcess(cxxopts::Options& options, int chunk_index){

  bool help = false;
  if(!options.count("output")){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Use after FindEquivalentBranches to improve tree topologies. Output written as bin file." << std::endl;
    exit(0);
  }

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Post-processing topology " << options["input"].as<std::string>() << "..." << std::endl;

  bool randomise = false;
  if(options.count("randomise")){
    randomise = true;
  }

  bool use_transitions = true;
  if(options.count("transversion")) use_transitions = false;

  std::mt19937 rng;
  int seed = 1;
  if(options.count("seed")){
    seed = options["seed"].as<int>();
  }
  rng.seed(seed);
  std::uniform_real_distribution<double> dist_unif(0.0,1.0);

  std::string file_out = options["output"].as<std::string>() + "/";

  int N, L, num_windows;
  std::vector<int> window_boundaries;
  FILE* fp = fopen((file_out + "parameters_c" + std::to_string(chunk_index) + ".bin").c_str(), "r");
  assert(fp != NULL);
  fread(&N, sizeof(int), 1, fp);
  fread(&L, sizeof(int), 1, fp);
  fread(&num_windows, sizeof(int), 1, fp);
  window_boundaries.resize(num_windows);
  fread(&window_boundaries[0], sizeof(int), num_windows, fp);
  fclose(fp);
  num_windows--;

  Data data((file_out + "chunk_" + std::to_string(chunk_index) + ".hap").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".bp").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".dist").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".r").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".rpos").c_str(), (file_out + "chunk_" + std::to_string(chunk_index) + ".state").c_str()); //struct data is defined in data.hpp
  data.name = (file_out + "chunk_" + std::to_string(chunk_index) + "/paint/relate");
  const std::string dirname = file_out + "chunk_" + std::to_string(chunk_index) + "/";

  std::vector<int> DAF(data.L, 0);
  for(int snp = 0; snp < data.L; snp++){
    for(int i = 0; i < data.N; i++){
      if(data.sequence[snp][i] == '1'){
        DAF[snp]++;
      }
    }
  }

  AncesTree anc;
  Mutations mut(data);
  AncesTreeBuilder ancbuilder(data);
  ancbuilder.PreCalcPotentialBranches(); // precalculating the number of decendants a branch needs to be equivalent (narrowing search space)
  Leaves sequences_carrying_mutation;
  sequences_carrying_mutation.member.resize(data.N);

  /////
  //First anc
  std::string filename;

  std::vector<Leaves> desc;

  float val;
  std::vector<int> nodes(2);
  std::vector<int> remaining;

  int thr = (int) (0.03 * N) + 1; //add 1 so I can use <

  std::vector<std::string> mut_filenames(num_windows);
  for(int i = 0; i < num_windows; i++){
    mut_filenames[i] = dirname + options["output"].as<std::string>() + "_" + std::to_string(i) + ".mut";
  }
  //mut.ReadShortFormat(mut_filenames);

	int count = 0;
  for(int anc_index = 0; anc_index < num_windows; anc_index++){

    int section_startpos = window_boundaries[anc_index];
    int section_endpos;
		if(anc_index < num_windows-1){
			section_endpos = window_boundaries[anc_index+1]-1;
		}else{
			section_endpos = data.L-1;
		}
		if(section_endpos >= data.L) section_endpos = data.L-1;

    filename = dirname + options["output"].as<std::string>() + "_" + std::to_string(anc_index) + ".anc";
    anc.ReadBin(filename);

		mut.ReadShortFormat(mut_filenames[anc_index], section_startpos, section_endpos);

    //do the post processing here
    int N_total = 2*anc.N-1;
    int N = anc.N;
    int root = 2*anc.N-2;
    float threshold = 1e6;

    CorrTrees::iterator it_seq;
    int tree_index = 0;
    int bp_init = data.bp_pos[data.L-1];
    float tree_bp, dist;
    int score12, score13, score23, tmp;
    bool mapped;
    int num_desc;

    int snp_start = 0, snp_end = 0;

    std::vector<int> equivalent_branches;
    std::vector<int> is_mapped(data.L, 0);
    for(it_seq = anc.seq.begin(); it_seq != anc.seq.end(); it_seq++){

      tree_bp = data.bp_pos[(*it_seq).pos];
      //std::cerr << (*it_seq).pos << " " << data.L << " " << ((float)(*it_seq).pos)/data.L << std::endl;
      if(snp_start < data.L){
        while(data.bp_pos[snp_start] < tree_bp - threshold){
          snp_start++;
          if(snp_start == data.L) break;
        }
      }
      if(snp_end < data.L){
        while(data.bp_pos[snp_end] < tree_bp + threshold){
          snp_end++;
          if(snp_end == data.L) break;
        }
      }
      (*it_seq).tree.FindAllLeaves(desc);

      if(1){
        if(it_seq != anc.seq.begin()){
          //check which branches are not associated with previous tree and only check those
          ancbuilder.BranchAssociation((*std::prev(it_seq,1)).tree, (*it_seq).tree, equivalent_branches); //O(N^2)  
          for(int i = 0; i < (*it_seq).tree.nodes.size(); i++){
            if(equivalent_branches[i] != -1){
              //if((*std::prev(it_seq,1)).tree.nodes[equivalent_branches[i]].num_events >= 1.0){
              (*it_seq).tree.nodes[i].num_events = 1.0;
              //}
            }
          }
        }
      }

      std::fill(is_mapped.begin(), is_mapped.end(), 0);
      int iter;
      for(iter = 0; iter < 5; iter++){

        if(iter > 0){
          if(randomise){

            for(int i = root-1; i >= N; i--){

              int node1 = i;
              int parent = (*(*it_seq).tree.nodes[i].parent).label;
              int node2 = (*(*it_seq).tree.nodes[parent].child_left).label;
              if(node2 == i){
                node2 = (*(*it_seq).tree.nodes[parent].child_right).label; 
              }

              if( (*it_seq).tree.nodes[node1].num_events < 1.0 ){
                if( (*it_seq).tree.nodes[node2].num_events < 1.0 || (*it_seq).tree.nodes[parent].num_events < 1.0 ){

                  int child1 = (*(*it_seq).tree.nodes[node1].child_left).label;
                  int child2 = (*(*it_seq).tree.nodes[node1].child_right).label;
                  int child3, child4;

                  remaining = {child1, child2, node2, -1};

                  bool shuffle_four = false;
                  if( (*it_seq).tree.nodes[node2].num_events < 1.0 && (*it_seq).tree.nodes[node2].child_left != NULL){
                    child3 = (*(*it_seq).tree.nodes[node2].child_left).label;
                    child4 = (*(*it_seq).tree.nodes[node2].child_right).label;
                    remaining = {child1, child2, child3, child4};
                    shuffle_four = true;
                  }else if( (*it_seq).tree.nodes[parent].num_events < 1.0 && parent != root ){
                    if(0){
                      child3 = node2;
                      node2  = parent;
                      parent = (*(*it_seq).tree.nodes[node2].parent).label;
                      child4 = (*(*it_seq).tree.nodes[parent].child_left).label;
                      if(child4 == node2){
                        child4 = (*(*it_seq).tree.nodes[parent].child_right).label;
                      }
                      remaining = {child1, child2, child3, child4};
                      shuffle_four = true;
                    }
                  }

                  if(1){

                    if(shuffle_four){

                      //make sure that node2 < node1 and coords[node2] < coords[node1]
                      if(node2 > node1){
                        int foo = node1;
                        node1 = node2;
                        node2 = foo;
                      }

                      val = dist_unif(rng);
                      if(val < 1.0/6.0){
                        nodes[0] = child1;
                        nodes[1] = child2;
                        remaining[0] = node2;
                        remaining[1] = remaining[3];
                        remaining[3] = -1;
                      }else if(val < 2.0/6.0){
                        nodes[0] = child1;
                        nodes[1] = child3;
                        remaining[0] = node2;
                        remaining[2] = remaining[3];
                        remaining[3] = -1;
                      }else if(val < 3.0/6.0){
                        nodes[0] = child1;
                        nodes[1] = child4;
                        remaining[0] = node2;
                        remaining[3] = -1;
                      }else if(val < 4.0/6.0){
                        nodes[0] = child2;
                        nodes[1] = child3;
                        remaining[1] = node2;
                        remaining[2] = remaining[3];
                        remaining[3] = -1;
                      }else if(val < 5.0/6.0){
                        nodes[0] = child2;
                        nodes[1] = child4;
                        remaining[1] = node2;
                        remaining[3] = -1;
                      }else{ 
                        nodes[0] = child3;
                        nodes[1] = child4;
                        remaining[2] = node2;
                        remaining[3] = -1;
                      }

                      (*it_seq).tree.nodes[node2].child_left       = &(*it_seq).tree.nodes[nodes[0]];
                      (*it_seq).tree.nodes[node2].child_right      = &(*it_seq).tree.nodes[nodes[1]];
                      (*it_seq).tree.nodes[nodes[0]].parent        = &(*it_seq).tree.nodes[node2];
                      (*it_seq).tree.nodes[nodes[1]].parent        = &(*it_seq).tree.nodes[node2];
                      //(*it_seq).tree.nodes[node2].num_events = 1.0;
                    }

                    val = dist_unif(rng);
                    if(val < 1.0/3.0){
                      nodes[0]     = remaining[0];
                      nodes[1]     = remaining[1];
                      remaining[0] = node1;
                      remaining[1] = remaining[2];
                      remaining[2] = -1;
                    }else if(val < 2.0/3.0){
                      nodes[0]     = remaining[0];
                      nodes[1]     = remaining[2];
                      remaining[0] = node1;
                      remaining[2] = -1;
                    }else{ 
                      nodes[0]     = remaining[1];
                      nodes[1]     = remaining[2];
                      remaining[1] = node1;
                      remaining[2] = -1;
                    }

                    (*it_seq).tree.nodes[node1].child_left       = &(*it_seq).tree.nodes[nodes[0]];
                    (*it_seq).tree.nodes[node1].child_right      = &(*it_seq).tree.nodes[nodes[1]];
                    (*it_seq).tree.nodes[nodes[0]].parent        = &(*it_seq).tree.nodes[node1];
                    (*it_seq).tree.nodes[nodes[1]].parent        = &(*it_seq).tree.nodes[node1];
                    //(*it_seq).tree.nodes[node1].num_events = 1.0;

                    (*it_seq).tree.nodes[parent].child_left          = &(*it_seq).tree.nodes[remaining[0]];
                    (*it_seq).tree.nodes[parent].child_right         = &(*it_seq).tree.nodes[remaining[1]];
                    (*it_seq).tree.nodes[remaining[0]].parent        = &(*it_seq).tree.nodes[parent];
                    (*it_seq).tree.nodes[remaining[1]].parent        = &(*it_seq).tree.nodes[parent];

                  }

                }

              }
            }
          }
        }

        bool is_updated = false;
        for(int i = root-1; i >= N; i--){

          int node1 = (*(*it_seq).tree.nodes[i].child_left).label;
          int node2 = (*(*it_seq).tree.nodes[i].child_right).label;

          int parent = (*(*it_seq).tree.nodes[i].parent).label;
          int node3 = (*(*it_seq).tree.nodes[parent].child_left).label;
          if(node3 == i){
            node3 = (*(*it_seq).tree.nodes[parent].child_right).label; 
          }

          assert(node1 != i);
          assert(node2 != i);

          if( (*it_seq).tree.nodes[i].num_events < 1.0 ){

            //check if there is evidence for node1+node2, node1 + node3, node2 + node3
            int closest_event12 = bp_init, closest_event13 = bp_init, closest_event23 = bp_init;

            //need data.seq, data.bp_pos, DAF
            assert((*it_seq).pos >= snp_start);
            assert((*it_seq).pos <= snp_end);

            if(0){
              for(int k = snp_start; k < snp_end; k++){
                dist = data.bp_pos[k] - tree_bp;
                if(dist < 0) dist *= -1.0;
                if( (desc[node1].num_leaves + desc[node2].num_leaves + desc[node3].num_leaves > DAF[k] - thr) &&
                    (desc[node1].num_leaves - DAF[k] < thr || desc[node2].num_leaves - DAF[k] < thr || desc[node3].num_leaves - DAF[k] < thr)){
                  CheckBranch(data, desc, DAF, thr, node1, node2, node3, closest_event12, closest_event13, closest_event23, dist, k);
                }
              }
            }else{

              mapped = false;
              int k = (*it_seq).pos;
              if(is_mapped[k] == 0 && DAF[k] > 1){
                if( (desc[node1].num_leaves + desc[node2].num_leaves + desc[node3].num_leaves > DAF[k] - thr) &&
                    (desc[node1].num_leaves - DAF[k] < thr || desc[node2].num_leaves - DAF[k] < thr || desc[node3].num_leaves - DAF[k] < thr)){
                  dist = data.bp_pos[k] - tree_bp;
                  if(dist < 0) dist *= -1.0;
                  //std::cerr << "1: " << k << " " << data.L << " " << data.sequence.size() << std::endl;
                  mapped = CheckBranch(data, desc, DAF, thr, node1, node2, node3, closest_event12, closest_event13, closest_event23, dist, k);
                }
              }

              if(!mapped){
                for(int l = 1; l < std::max((*it_seq).pos - snp_start, snp_end - (*it_seq).pos); l++){
                  //std::cerr << l << std::endl;
                  k = (*it_seq).pos - l;
                  if(k > 0){
                    if(is_mapped[k] == 0 && DAF[k] > 1){
                      if( (desc[node1].num_leaves + desc[node2].num_leaves + desc[node3].num_leaves > DAF[k] - thr) &&
                          (desc[node1].num_leaves - DAF[k] < thr || desc[node2].num_leaves - DAF[k] < thr || desc[node3].num_leaves - DAF[k] < thr)){
                        dist = data.bp_pos[k] - tree_bp;
                        if(dist < 0) dist *= -1.0;
                        //std::cerr << "2: " << k << " " << data.L << " " << data.sequence.size() << std::endl;
                        mapped = CheckBranch(data, desc, DAF, thr, node1, node2, node3, closest_event12, closest_event13, closest_event23, dist, k);
                      }
                    }
                  }
                  if(mapped){
                    if(k >= 0) is_mapped[k] = 1;
                    break;
                  }
                  k = (*it_seq).pos + l;
                  if(k < data.L){
                    if(is_mapped[k] == 0 && DAF[k] > 1){
                      if( (desc[node1].num_leaves + desc[node2].num_leaves + desc[node3].num_leaves > DAF[k] - thr) &&
                          (desc[node1].num_leaves - DAF[k] < thr || desc[node2].num_leaves - DAF[k] < thr || desc[node3].num_leaves - DAF[k] < thr)){
                        dist = data.bp_pos[k] - tree_bp;
                        if(dist < 0) dist *= -1.0;
                        //std::cerr << "3: " << k << " " << data.L << " " << data.sequence.size() << std::endl;
                        mapped = CheckBranch(data, desc, DAF, thr, node1, node2, node3, closest_event12, closest_event13, closest_event23, dist, k);
                      }
                    }
                  }
                  if(mapped){
                    if(k < data.L) is_mapped[k] = 1;
                    break;
                  }
                }
              }else{
                is_mapped[k] = 1;
              }
            }

            //SNP is out of question if its DAF is 
            if( (closest_event13 < closest_event12 && closest_event13 <= closest_event23) || (closest_event13 <= closest_event12 && closest_event13 < closest_event23) ){

              is_updated = true;

              (*it_seq).tree.nodes[i].child_left       = &(*it_seq).tree.nodes[node1];
              (*it_seq).tree.nodes[i].child_right      = &(*it_seq).tree.nodes[node3];
              (*it_seq).tree.nodes[node1].parent       = &(*it_seq).tree.nodes[i];
              (*it_seq).tree.nodes[node3].parent       = &(*it_seq).tree.nodes[i];

              (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[i];
              (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[node2];
              (*it_seq).tree.nodes[i].parent           = &(*it_seq).tree.nodes[parent];
              (*it_seq).tree.nodes[node2].parent       = &(*it_seq).tree.nodes[parent];

              (*it_seq).tree.nodes[i].num_events    = 1.0;

              //update desc[i]
              desc[i].member = desc[node1].member;
              desc[i].member.insert(desc[i].member.end(),desc[node3].member.begin(), desc[node3].member.end());
              desc[i].num_leaves = desc[node1].num_leaves + desc[node3].num_leaves;
              std::sort(desc[i].member.begin(), desc[i].member.end());

            }else if( (closest_event23 < closest_event12 && closest_event23 <= closest_event13) || (closest_event23 <= closest_event12 && closest_event23 < closest_event13) ){

              is_updated = true;

              (*it_seq).tree.nodes[i].child_left       = &(*it_seq).tree.nodes[node2];
              (*it_seq).tree.nodes[i].child_right      = &(*it_seq).tree.nodes[node3];
              (*it_seq).tree.nodes[node2].parent       = &(*it_seq).tree.nodes[i];
              (*it_seq).tree.nodes[node3].parent       = &(*it_seq).tree.nodes[i];

              (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[i];
              (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[node1];
              (*it_seq).tree.nodes[i].parent           = &(*it_seq).tree.nodes[parent];
              (*it_seq).tree.nodes[node1].parent       = &(*it_seq).tree.nodes[parent];

              (*it_seq).tree.nodes[i].num_events    = 1.0;

              desc[i].member = desc[node2].member;
              desc[i].member.insert(desc[i].member.end(),desc[node3].member.begin(), desc[node3].member.end());
              desc[i].num_leaves = desc[node2].num_leaves + desc[node3].num_leaves;
              std::sort(desc[i].member.begin(), desc[i].member.end());

            }else if( (closest_event12 < closest_event23 && closest_event12 <= closest_event13) || (closest_event12 <= closest_event23 && closest_event12 < closest_event13) ){
              (*it_seq).tree.nodes[i].num_events    = 1.0;
            }
          }
        }
        if(randomise){
          if(iter > 0 && !is_updated) break;
        }else{
          if(!is_updated) break;
        }
      }

      //std::cerr << iter << std::endl;

      Relabel((*it_seq).tree, N);

      for(int i = 0; i < N_total; i++){
        (*it_seq).tree.nodes[i].SNP_begin  = (*it_seq).pos;
        if(std::next(it_seq,1) != anc.seq.end()){
          (*it_seq).tree.nodes[i].SNP_end  = (*std::next(it_seq,1)).pos;
        }else{
          (*it_seq).tree.nodes[i].SNP_end  = L-1;
        }
        (*it_seq).tree.nodes[i].num_events = 0;
      }

      for(int snp = (*it_seq).pos; snp < (*std::next(it_seq,1)).pos; snp++){
        //map mutation onto tree
        std::fill(sequences_carrying_mutation.member.begin(), sequences_carrying_mutation.member.end(), 0);
        sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
        for(int i = 0; i < data.N; i++){
          if(data.sequence[snp][i] == '1'){
            sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
            sequences_carrying_mutation.num_leaves++;
          }else{
            sequences_carrying_mutation.member[i] = 0;
          }
        }

        mut.info[snp].branch.clear();
        ancbuilder.mutations.info[snp].flipped = false;
        if(sequences_carrying_mutation.num_leaves == data.N){
          (*it_seq).tree.nodes[root].num_events += 1.0;
          mut.info[snp].branch.resize(1);
          mut.info[snp].branch[0] = root;
        }else{
          if(use_transitions){
            ancbuilder.IsSNPMapping((*it_seq).tree, sequences_carrying_mutation, snp);
          }else{
            ancbuilder.IsSNPMapping((*it_seq).tree, sequences_carrying_mutation, snp, data.state[snp]);
          }
          mut.info[snp].branch  = ancbuilder.mutations.info[snp].branch;
          mut.info[snp].flipped = ancbuilder.mutations.info[snp].flipped;

        }
      }

			count++;
    }

    filename = dirname + options["output"].as<std::string>() + "_" + std::to_string(anc_index) + ".anc";
    anc.DumpBin(filename);
    filename = dirname + options["output"].as<std::string>() + "_" + std::to_string(anc_index) + ".mut";
    mut.DumpShortFormat(filename, section_startpos, section_endpos);
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


