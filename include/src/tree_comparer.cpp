#include "tree_comparer.hpp"

float 
DistanceUsingPearsonCorrelation(Tree& tr1, Tree& tr2){

  //number of leaves
  int N = (tr1.nodes.size() + 1)/2;
  std::vector<Leaves> tr1_leaves;
  std::vector<Leaves> tr2_leaves;

  Correlation cor(N);

  tr1.FindAllLeaves(tr1_leaves);
  tr2.FindAllLeaves(tr2_leaves);

  float correlation = 0.0;

  //iterate through all internal nodes in tr1
  for(int i = N; i < (int) tr1.nodes.size(); i++){

    float correlation_max = 0.0;

    if(tr1.nodes[i].parent != NULL){
      for(int j = N; j < (int) tr2.nodes.size(); j++){

        if(tr2.nodes[j].parent != NULL){
          float correlation_tmp = cor.Pearson(tr1_leaves[i], tr2_leaves[j]);
          if(correlation_tmp > correlation_max){
            correlation_max = correlation_tmp;
          }
          assert(correlation_max >= 0.0);
          if(correlation_max == 1.0) break;
        }
      }
      correlation += correlation_max * correlation_max;
    }
  }

  assert(correlation/((float) N - 2.0) <= 1.0);
  return correlation/((float) N - 2.0);

}

float
PartitionMetric(Tree& tr1, Tree& tr2, float threshold){

  //number of clades defined by one but not the other

  //number of leaves
  int N = (tr1.nodes.size() + 1)/2;
  std::vector<Leaves> tr1_leaves;
  std::vector<Leaves> tr2_leaves;

  Correlation cor(N);

  tr1.FindAllLeaves(tr1_leaves);
  tr2.FindAllLeaves(tr2_leaves);

  float correlation_tmp;
  int distance = 0;

  //iterate through all internal nodes in tr1
  for(int i = N; i < (int) tr1.nodes.size(); i++){

    if(tr1.nodes[i].parent != NULL){ //not interested in root

      correlation_tmp = 0.0;
      for(int j = N; j < (int) tr2.nodes.size(); j++){
        if(tr2.nodes[j].parent != NULL){ //not interested in root
          correlation_tmp = cor.Pearson(tr1_leaves[i], tr2_leaves[j]);
          if(correlation_tmp >= threshold) break;
        }
      }
      if(correlation_tmp < threshold){
        //std::cerr << correlation_tmp << " " << threshold << std::endl;
        distance++;
      }
    }

  }


  //iterate through all internal nodes in tr1
  for(int i = N; i < (int) tr2.nodes.size(); i++){

    if(tr2.nodes[i].parent != NULL){ //not interested in root

      correlation_tmp = 0.0;
      for(int j = N; j < (int) tr1.nodes.size(); j++){
        if(tr1.nodes[j].parent != NULL){ //not interested in root
          correlation_tmp = cor.Pearson(tr2_leaves[i], tr1_leaves[j]);
          if(correlation_tmp >= threshold) break;
        }
      }
      if(correlation_tmp < threshold) distance++;
    }

  }

  return distance/(2.0 * N - 4.0);

}

float
BranchScoreMetric(Tree& tr1, Tree& tr2, float threshold){

  //number of leaves
  int N = (tr1.nodes.size() + 1)/2;
  std::vector<Leaves> tr1_leaves;
  std::vector<Leaves> tr2_leaves;

  Data data(N,1);
  Correlation cor(N);

  tr1.FindAllLeaves(tr1_leaves);
  tr2.FindAllLeaves(tr2_leaves);

  float correlation_tmp;
  float distance = 0.0;
  int equiv_branch = 0;
  //iterate through all internal nodes in tr1
  for(int i = N; i < (int) tr1.nodes.size(); i++){

    if(tr1.nodes[i].parent != NULL){ //not interested in root

      correlation_tmp = 0.0;
      for(int j = N; j < (int) tr2.nodes.size(); j++){
        if(tr2.nodes[j].parent != NULL){ //not interested in root
          correlation_tmp = cor.Pearson(tr1_leaves[i], tr2_leaves[j]);
          if(correlation_tmp >= threshold){
            equiv_branch = j;
            break;
          }
        }
      }
      if(correlation_tmp < threshold){
        distance += tr1.nodes[i].branch_length/data.Ne * tr1.nodes[i].branch_length/data.Ne;
      }else{
        float distance_tmp = (tr1.nodes[i].branch_length/data.Ne - tr2.nodes[equiv_branch].branch_length/data.Ne);
        distance += distance_tmp * distance_tmp;
      }
    }

  }


  //iterate through all internal nodes in tr1
  for(int i = N; i < (int) tr2.nodes.size(); i++){

    if(tr2.nodes[i].parent != NULL){ //not interested in root

      correlation_tmp = 0.0;
      for(int j = N; j < (int) tr1.nodes.size(); j++){
        if(tr1.nodes[j].parent != NULL){ //not interested in root
          correlation_tmp = cor.Pearson(tr2_leaves[i], tr1_leaves[j]);
          if(correlation_tmp >= threshold) break;
        }
      }
      if(correlation_tmp < threshold){
        distance += tr2.nodes[i].branch_length/data.Ne * tr2.nodes[i].branch_length/data.Ne;
      }else{
        float distance_tmp = (tr2.nodes[i].branch_length/data.Ne - tr1.nodes[equiv_branch].branch_length/data.Ne);
        distance += distance_tmp * distance_tmp;
      }
    }

  }

  return distance/(2.0 * N - 4.0);

}

float
TimeWhileKAncestorsDistance(Tree& tr1, Tree& tr2){

  //number of leaves
  int N = (tr1.nodes.size() + 1)/2;

  std::vector<double> coordinates_tr1(tr1.nodes.size());
  std::vector<double> coordinates_tr2(tr2.nodes.size());

  Data data(N,1);
  InferBranchLengths bl(data);

  int root = (int) tr1.nodes.size() - 1;
  if(tr1.nodes[root].parent != NULL){
    for(int i = 0; i < (int) tr1.nodes.size(); i++){
      if(tr1.nodes[i].parent == NULL){
        root = i; //root
        break;
      }
    }
  }
  bl.GetCoordinates(tr1.nodes[root], coordinates_tr1);
  root = (int) tr2.nodes.size() - 1;
  if(tr2.nodes[root].parent != NULL){
    for(int i = 0; i < (int) tr2.nodes.size(); i++){
      if(tr2.nodes[i].parent == NULL){
        root = i; //root
        break;
      }
    }
  }
  bl.GetCoordinates(tr2.nodes[root], coordinates_tr2);

  std::sort(coordinates_tr1.begin(), coordinates_tr1.end());
  std::sort(coordinates_tr2.begin(), coordinates_tr2.end());

  float distance = 0.0;
  float tmp_distance;
  float num_lineages;
  for(int i = N; i < tr1.nodes.size(); i++){
    num_lineages = 2*N-i;
    assert(coordinates_tr1[i] >= coordinates_tr1[i-1]);
    assert(coordinates_tr2[i] >= coordinates_tr2[i-1]);
    float tmp_distance = ((coordinates_tr1[i] - coordinates_tr1[i-1]) - (coordinates_tr2[i] - coordinates_tr2[i-1])) * (num_lineages * (num_lineages-1.0))/2.0;
    distance += tmp_distance * tmp_distance;
  }

  //std::cerr << coordinates_tr1[tr1.nodes.size() -1] << " " << coordinates_tr2[tr2.nodes.size() - 1] << std::endl;

  return sqrt(distance)/(N-1.0);

}



float
GetTotalBranchLength(Tree& tr){

  float total_branch_length = 0.0;

  for(int i = 0; i < (int) tr.nodes.size(); i++){
    if(tr.nodes[i].parent != NULL){
      total_branch_length += tr.nodes[i].branch_length;
    }
  }

  return total_branch_length;

}

float
GetTMRCA(Tree& tr){

  int root = (int) tr.nodes.size() - 1;
  for(int i = 0; i < (int) tr.nodes.size(); i++){
    if(tr.nodes[i].parent == NULL){
      root = i; //root
      break;
    }
  }

  float total_tree_height = 0.0;
  Node parent = tr.nodes[root];
  while(parent.child_left != NULL){
    parent = *parent.child_left;
    total_tree_height += parent.branch_length;
  }

  return total_tree_height;

}

void
PairwiseTMRCA(Tree& tr, std::vector<float>& pairwiseTMRCA){

  //initialize vector of length N choose 2
  //start at root and record pairwise TMRCAs in this vector.
  int N = (tr.nodes.size() + 1)/2;
  //std::vector<float> pairwiseTMRCA(N*N);
  pairwiseTMRCA.resize(N*N);

  int root = (int) tr.nodes.size() - 1;
  for(int i = 0; i < (int) tr.nodes.size(); i++){
    if(tr.nodes[i].parent == NULL){
      root = i; //root
      break;
    }
  }

  std::vector<Leaves> leaves;
  tr.FindAllLeaves(leaves);

  //std::cerr << pairwiseTMRCA.size() << std::endl;
  GetPairwiseTMRCA(tr.nodes[root], pairwiseTMRCA, leaves);

  /*
     for(int i = 0; i < (int)pairwiseTMRCA.size(); i++){
     std::cerr << (int) (i/N) << " " << i % N << " " << pairwiseTMRCA[i] << std::endl;
     }
     */

}

float
GetPairwiseTMRCA(Node& n, std::vector<float>& pairwiseTMRCA, std::vector<Leaves>& leaves){

  float height = 0.0;
  int N = std::sqrt(pairwiseTMRCA.size());

  if(n.child_left != NULL){

    Node child1 = *n.child_left;
    Node child2 = *n.child_right;

    height = GetPairwiseTMRCA(child1, pairwiseTMRCA, leaves) + child1.branch_length;
    GetPairwiseTMRCA(child2, pairwiseTMRCA, leaves);

    for(std::vector<int>::iterator it = leaves[child1.label].member.begin(); it != leaves[child1.label].member.end(); it++){
      for(std::vector<int>::iterator jt = leaves[child2.label].member.begin(); jt != leaves[child2.label].member.end(); jt++){
        pairwiseTMRCA[(*it)*N+(*jt)] = height;
        pairwiseTMRCA[(*jt)*N+(*it)] = height;
      }
    }

  }

  return height;

}

