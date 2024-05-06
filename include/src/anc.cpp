#include "anc.hpp"

//////////////////////////////////

void
Tree::GetMsPrime(igzstream& is, int num_nodes){

  nodes.resize(2*num_nodes-1);

  std::string line;
  int count_lines = 0;
  while(count_lines < 2*num_nodes-1){
    getline(is, line); 

    //std::cerr << count_lines << " " << 2*num_nodes << ": " << line << std::endl;

    float node, child_node_left, child_node_right, branch_length_left, branch_length_right; 
    std::stringstream lineStream(line);

    lineStream >> node; //if leaf then only one entry
    if(lineStream.peek() != decltype(lineStream)::traits_type::eof()){ //if there is more in lineStream
      lineStream >> child_node_left >> child_node_right >> branch_length_left >> branch_length_right;
      nodes[node].child_left         = &nodes[child_node_left];
      nodes[node].child_right        = &nodes[child_node_right];
      nodes[child_node_left].parent  = &nodes[node];
      nodes[child_node_right].parent = &nodes[node];
      nodes[child_node_left].branch_length  = branch_length_left;
      nodes[child_node_right].branch_length = branch_length_right;
    }
    nodes[node].label = node;

    count_lines++;

  }

}

void
Tree::ReadTree(const char* line, int N){

  //clock_t tbegin = clock();
  int i = 0;
  nodes.clear();
  nodes.resize(2*N-1);

  std::vector<Node>::iterator it_node = nodes.begin(); 
  Node* p_parent;
  int num_node = 0, parent;

  while(num_node < 2*N-1){

    sscanf(&line[i], "%d:(%lf %f %d %d)", &parent, &(*it_node).branch_length, &(*it_node).num_events, &(*it_node).SNP_begin, &(*it_node).SNP_end); 

    while(line[i] != ')' && line[i] != '\n') i++;
    i += 2;

    if(parent != -1){
      p_parent   = &nodes[parent];
      (*it_node).parent         = p_parent;
      (*it_node).label          = num_node;
      if((*p_parent).child_left == NULL){
        (*p_parent).child_left  = &(*it_node);
      }else{
        (*p_parent).child_right = &(*it_node);       
      }
    }else{
      (*it_node).parent         = NULL;
      (*it_node).label          = num_node;
    }       

    num_node++;
    it_node++;

  }

  //clock_t tend = clock();
  //double elapsed_secs = double(tend - tbegin) / CLOCKS_PER_SEC;
  //std::cerr << elapsed_secs << std::endl;

}

void
Tree::ReadTreeBin(FILE* pfile, int N){

  //clock_t tbegin = clock();
  nodes.clear();
  nodes.resize(2*N-1);

  std::vector<Node>::iterator it_node = nodes.begin(); 
  Node* p_parent;
  int num_node = 0, parent;

  while(num_node < 2*N-1){

    fread(&parent, sizeof(int), 1, pfile);

    if(parent != -1){
      p_parent   = &nodes[parent];
      (*it_node).parent         = p_parent;
      (*it_node).label          = num_node;
      if((*p_parent).child_left == NULL){
        (*p_parent).child_left  = &(*it_node);
      }else{
        (*p_parent).child_right = &(*it_node);       
      }
    }else{
      (*it_node).parent         = NULL;
      (*it_node).label          = num_node;
    }       

    fread(&(*it_node).branch_length, sizeof(double), 1, pfile);
    fread(&(*it_node).num_events, sizeof(float), 1, pfile);
    fread(&(*it_node).SNP_begin, sizeof(int), 1, pfile);
    fread(&(*it_node).SNP_end, sizeof(int), 1, pfile);

    num_node++;
    it_node++;

  }

  //clock_t tend = clock();
  //double elapsed_secs = double(tend - tbegin) / CLOCKS_PER_SEC;
  //std::cerr << elapsed_secs << std::endl;

}


void
Tree::WriteNewick(const std::string& filename_newick, double factor, const bool add) const{
 
  std::ofstream os;
  if(!add){ 
    os.open(filename_newick);
  }else{
    os.open(filename_newick, std::ofstream::app);
  }

	WriteNewick(os, factor);

  os.close();

}

void
Tree::WriteNewick(std::ofstream& os, double factor) const{

	//coordinates.clear();
	//maybe not the most efficient convertion algorithm but it works
	int root = (int) nodes.size() - 1;
	for(int i = 0; i < (int) nodes.size(); i++){
		if(nodes[i].parent == NULL){
			root = i; //root
			break;
		}
	}

	std::list<Node> todo_nodes;
	float l1 = ((*nodes[root].child_left).branch_length) * factor;
	float l2 = ((*nodes[root].child_right).branch_length) * factor;

	//if(l1 < 0.0) std::cerr << root << ": " << coordinates[root].tree << ", " << children[root][0] << ": " << coordinates[children[root][0]].tree << std::endl;    
	//if(l2 < 0.0) std::cerr << root << ": " << coordinates[root].tree << ", " << children[root][1] << ": " << coordinates[children[root][1]].tree << std::endl;

	std::string newick = "(" + std::to_string((*nodes[root].child_left).label) + ":" + std::to_string(l1)  + "," + std::to_string((*nodes[root].child_right).label) + ":" + std::to_string(l2) + ")";
	todo_nodes.push_back(*nodes[root].child_left);
	todo_nodes.push_back(*nodes[root].child_right);

	while(todo_nodes.size() > 0){

		std::list<Node>::iterator node = todo_nodes.begin();

		if((*node).child_left == NULL){
			todo_nodes.erase(node);
		}else{
			int index = 0;
			std::string node_index;
			for(; index < (int) newick.size(); index++){
				while( !isdigit(newick[index]) ){
					//std::cerr << newick[index] << " " << isdigit(newick[index]) << std::endl;
					index++;
				}
				node_index.clear();
				while( newick[index] != ',' && newick[index] != ')' && newick[index] != ':'){
					node_index = node_index + newick[index];
					index++;
					//std::cerr << node_index << " " << index << " " << newick[index] << std::endl;
				}
				if((float)(*node).label == stof(node_index) && newick[index] == ':') break;
			}
			assert(index < (int) newick.size());

			float l1 = ((*(*node).child_left).branch_length) * factor;
			float l2 = ((*(*node).child_right).branch_length) * factor;

			std::string new_brackets = "(" + std::to_string((*(*node).child_left).label) + ":" + std::to_string(l1) + "," + std::to_string((*(*node).child_right).label) + ":" + std::to_string(l2) + ")";
			newick.replace(index-node_index.size(), node_index.size(), new_brackets);
			todo_nodes.push_back(*(*node).child_left);
			todo_nodes.push_back(*(*node).child_right);
			todo_nodes.erase(node);

		}

		//std::cerr << "newick: " << newick << " " << todo_nodes.size() << std::endl;

	}
	newick += ";";

	os << newick << std::endl;

}

void
Tree::WriteNHX(const std::string& filename_nhx, std::vector<std::string>& property, const bool add) const{

  if(property.size() != nodes.size()){
    std::cerr << "Property vector has wrong size." << std::endl;
    exit(1);
  }

  //maybe not the most efficient convertion algorithm but it works
  int root = (int) nodes.size() - 1;
  for(int i = 0; i < (int) nodes.size(); i++){
    if(nodes[i].parent == NULL){
      root = i; //root
      break;
    }
  }

  std::ofstream os_new;
  if(!add){ 
    os_new.open(filename_nhx);
  }else{
    os_new.open(filename_nhx, std::ofstream::app);
  }

  std::list<Node> todo_nodes; //stores node indices in newick string
  float l1 = ((*nodes[root].child_left).branch_length);
  float l2 = ((*nodes[root].child_right).branch_length);

  std::string newick = "(" + std::to_string((*nodes[root].child_left).label) + ":" + std::to_string(l1) + "[&&NHX:S=" + property[(*nodes[root].child_left).label] + "]," + std::to_string((*nodes[root].child_right).label) + ":" + std::to_string(l2) + "[&&NHX:S=" + property[(*nodes[root].child_right).label] + "])";
  todo_nodes.push_back(*nodes[root].child_left);
  todo_nodes.push_back(*nodes[root].child_right);

  while(todo_nodes.size() > 0){

    std::list<Node>::iterator node = todo_nodes.begin();

    if((*node).child_left == NULL){
      todo_nodes.erase(node); //erase if it is a tip
    }else{
      int index = 0;
      std::string node_index;
      for(; index < (int) newick.size(); index++){
        while( !isdigit(newick[index]) ){
          index++;
        }
        node_index.clear();
        while( newick[index] != ',' && newick[index] != ')' && newick[index] != ':'){
          node_index = node_index + newick[index];
          index++;
        }
        if((float)(*node).label == stof(node_index) && newick[index] == ':') break;
      }
      assert(index < (int) newick.size());

      float l1 = ((*(*node).child_left).branch_length);
      float l2 = ((*(*node).child_right).branch_length);

      std::string new_brackets = "(" + std::to_string((*(*node).child_left).label) + ":" + std::to_string(l1) + "[&&NHX:S=" + property[(*(*node).child_left).label]  + "]," + std::to_string((*(*node).child_right).label) + ":" + std::to_string(l2) + "[&&NHX:S=" + property[(*(*node).child_right).label] + "])";
      newick.replace(index-node_index.size(), node_index.size(), new_brackets); //replace node index by branching event
      todo_nodes.push_back(*(*node).child_left);
      todo_nodes.push_back(*(*node).child_right);
      todo_nodes.erase(node);

    }

  }
  newick += "[&&NHX:S=" + property[root] + "];";

  os_new << newick << std::endl;

  os_new.close();

}

void
Tree::WriteOrientedTree(const std::string& filename, const bool add){

  //file format:
  //value at position u is parent of node u and value after colon is branch length to parent
  //5:0.5 6:0.2 etc 

  std::ofstream os_tree;
  if(!add){ 
    os_tree.open(filename);
  }else{
    os_tree.open(filename, std::ofstream::app);
  }

  std::string oriented_tree;
  //oriented_tree.reserve() TODO: argument is number of char

  std::vector<Node>::iterator n_it = nodes.begin();
  int parent;
  for(; n_it != nodes.end(); n_it++){
    if((*n_it).parent == NULL){ 
      parent = -1;
    }else{
      parent = (*(*n_it).parent).label;
    }
    oriented_tree += std::to_string(parent) + ":" + std::to_string((*n_it).branch_length) + " ";
  }

  os_tree << oriented_tree << std::endl;

  os_tree.close();

}

//for debugging
void
Tree::PrintTree(){

  for(int n = 0; n < (int) nodes.size(); n++){
    if(nodes[n].child_left != NULL) std::cerr << nodes[n].label << " " << (*(nodes[n].child_left)).label << " " << (*(nodes[n].child_right)).label << " ";
    std::cerr << nodes[n].branch_length << std::endl;
  }

}

//Align two trees (can be inefficient as its only used for plotting small trees)
void
Tree::DetermineOrderOfLeaves(Node root, std::vector<int>& leaves){

  //std::cerr << "det order" << std::endl;

  leaves.clear();

  std::list<Node> todo_nodes;
  todo_nodes.push_back(root);
  while(todo_nodes.size() > 0){
    std::list<Node>::iterator node = todo_nodes.begin();
    if((*node).child_left == NULL){
      leaves.push_back((*node).label);
      todo_nodes.erase(node);
    }else{
      Node child_left = *(*node).child_left, child_right = *(*node).child_right;
      todo_nodes.erase(node);
      todo_nodes.push_front(child_left);
      todo_nodes.push_front(child_right);
    }
  }

  /*
     for(std::vector<Node>::iterator it = leaves.begin(); it != leaves.end(); it++){
     std::cerr << (*it).label << " ";
     }
     std::cerr << std::endl;
     */

}

int
Tree::CountCrossings(std::vector<int>& leaves1, std::vector<int>& leaves2){

  int num_crossings = 0;

  std::vector<int> dispos(leaves1.size());
  for(int i = 0; i < (int)leaves1.size(); i++){
    //search for same value in leaves2 (this can be done more efficiently using binary search or hashmaps)
    int j = 0;
    for(; j < (int)leaves2.size(); j++){
      if(leaves1[i] == leaves2[j]) break;
    }
    dispos[i] = j - i;
  }

  for(int i = 0; i < (int)leaves1.size(); i++){
    for(int j = i; j < (int)leaves1.size(); j++){
      if(j + dispos[j] < i + dispos[i]) num_crossings += i + dispos[i] - j - dispos[j];
    }
  }

  return num_crossings;

}

void
Tree::AlignTrees(Tree& ref_tree){

  //this algorithm also changes reference_tree

  int N_total = nodes.size();
  int N = (N_total + 1)/2;
  int root = nodes[N_total - 1].label;
  int ref_root = ref_tree.nodes[N_total - 1].label;
  if(nodes[root].parent != NULL){
    for(int i = N; i < N_total; i++){
      if(nodes[i].parent == NULL){
        root = nodes[i].label;
        break;
      }
    }
  }
  if(ref_tree.nodes[ref_root].parent != NULL){
    for(int i = N; i < N_total; i++){
      if(ref_tree.nodes[i].parent == NULL){
        ref_root = ref_tree.nodes[i].label;
        break;
      }
    }
  }

  //tree_leaves and reference_tree_leaves will contain the order of leaves in which they are printed
  std::vector<int> tree_leaves, ref_tree_leaves;

  //determine initial order of leaves
  DetermineOrderOfLeaves(nodes[root], tree_leaves);
  ref_tree.DetermineOrderOfLeaves(ref_tree.nodes[ref_root], ref_tree_leaves);  

  //count number of crossings  
  int num_crossings = CountCrossings(tree_leaves, ref_tree_leaves); 
  int num_crossings_after_swap;
  for(int iterate = 0; iterate < 1; iterate++){
    //swap and update number of crossings
    for(int i = N; i < N_total; i++){
      //swap
      int child1 = (*nodes[i].child_left).label;
      int child2 = (*nodes[i].child_right).label;

      nodes[i].child_left  = &nodes[child2];
      nodes[i].child_right = &nodes[child1];
      DetermineOrderOfLeaves(nodes[root], tree_leaves);
      num_crossings_after_swap = CountCrossings(tree_leaves, ref_tree_leaves);

      if(num_crossings_after_swap > num_crossings){
        nodes[i].child_left  = &nodes[child1];
        nodes[i].child_right = &nodes[child2];        
      }else{
        num_crossings = num_crossings_after_swap;
      }
    }
  }

  std::cerr << "number of crossings in alignment " << num_crossings << std::endl;

}

void 
Tree::FindAllLeaves(std::vector<Leaves>& leaves) const{

  int N_total = nodes.size();
  int N = (N_total + 1)/2;
  leaves.resize(N_total);

  Node root = nodes[N_total - 1];
  if(root.parent != NULL){
    for(int i = N; i < N_total; i++){
      if(nodes[i].parent == NULL){
        root = nodes[i];
        break;
      }
    }
  }

  FindLeaves(root, leaves);

}

void 
Tree::FindLeaves(Node& node, std::vector<Leaves>& leaves) const{

  if(node.child_left != NULL){

    Node child1 = *node.child_left;
    Node child2 = *node.child_right;

    FindLeaves(child1, leaves);
    FindLeaves(child2, leaves);

    leaves[node.label].member.resize( leaves[child1.label].member.size() + leaves[child2.label].member.size() );

    std::vector<int>::iterator it_member = leaves[node.label].member.begin();
    std::vector<int>::iterator it_child1_member = leaves[child1.label].member.begin();
    std::vector<int>::iterator it_child2_member = leaves[child2.label].member.begin();

    const std::vector<int>::iterator it_child1_member_end = leaves[child1.label].member.end();
    const std::vector<int>::iterator it_child2_member_end = leaves[child2.label].member.end();

    for(; it_member != leaves[node.label].member.end();){

      if(it_child1_member != it_child1_member_end && it_child2_member != it_child2_member_end){
        if(*it_child1_member < *it_child2_member){
          *it_member = *it_child1_member;
          it_child1_member++;
          it_member++;
        }else{
          *it_member = *it_child2_member;
          it_child2_member++;
          it_member++;
        }
      }else if(it_child1_member != it_child1_member_end){
        *it_member = *it_child1_member;
        it_child1_member++;
        it_member++;
      }else{
        *it_member = *it_child2_member;
        it_child2_member++;
        it_member++;
      }

    }

    leaves[node.label].num_leaves = leaves[child1.label].num_leaves + leaves[child2.label].num_leaves;

  }else{
    leaves[node.label].member.resize(1);    //resizing bitset and filling it with FALSE
    leaves[node.label].member[0]  = node.label; //setting position node.label to TRUE
    leaves[node.label].num_leaves = 1;
  }

}

void
Tree::TraverseTreeToGetCoordinates(Node& n, std::vector<float>& coordinates){

  if(n.child_left != NULL){

    TraverseTreeToGetCoordinates(*n.child_left, coordinates);
    TraverseTreeToGetCoordinates(*n.child_right, coordinates);
    coordinates[n.label] = std::max(coordinates[(*n.child_right).label] + (*n.child_right).branch_length, coordinates[(*n.child_left).label] + (*n.child_left).branch_length);  

  }else{
    coordinates[n.label] = 0.0;
  }

}

void
Tree::TraverseTreeToGetCoordinates_sample_age(Node& n, std::vector<float>& coordinates){

  if(n.child_left != NULL){

    TraverseTreeToGetCoordinates_sample_age(*n.child_left, coordinates);
    TraverseTreeToGetCoordinates_sample_age(*n.child_right, coordinates);
    coordinates[n.label] = std::max(coordinates[(*n.child_right).label] + (*n.child_right).branch_length, coordinates[(*n.child_left).label] + (*n.child_left).branch_length);  

  }else{
    coordinates[n.label] = (*sample_ages)[n.label];
  }

}

void 
Tree::GetCoordinates(std::vector<float>& coordinates){
  coordinates.resize(nodes.size());
  if(sample_ages == NULL){
    TraverseTreeToGetCoordinates(nodes[nodes.size() - 1], coordinates);
  }else if((*sample_ages).size() > 0){
    assert(2*(*sample_ages).size() - 1 == coordinates.size());
    TraverseTreeToGetCoordinates_sample_age(nodes[nodes.size() - 1], coordinates);
  }else{
    TraverseTreeToGetCoordinates(nodes[nodes.size() - 1], coordinates);
  }
}

void 
Tree::GetCoordinates(int node, std::vector<float>& coordinates){
  coordinates.resize(nodes.size());
  if(sample_ages == NULL){
    TraverseTreeToGetCoordinates(nodes[node], coordinates);
  }else if((*sample_ages).size() > 0){
    assert(2*(*sample_ages).size() - 1 == coordinates.size());
    TraverseTreeToGetCoordinates_sample_age(nodes[node], coordinates);
  }else{ 
    TraverseTreeToGetCoordinates(nodes[node], coordinates);
  }
}


/////////////////////////////////

void
Tree::GetNumberOfLeavesInSubpop(const Node& n, std::vector<int>& subpop, std::vector<int>& number_in_subpop) const{

  if(n.child_left != NULL){

    GetNumberOfLeavesInSubpop((*n.child_left), subpop, number_in_subpop);
    GetNumberOfLeavesInSubpop((*n.child_right), subpop, number_in_subpop);
    number_in_subpop[n.label] = number_in_subpop[(*n.child_left).label] + number_in_subpop[(*n.child_right).label];

  }else{

    for(std::vector<int>::iterator it_subpop = subpop.begin(); it_subpop != subpop.end(); it_subpop++){
      if(*it_subpop == n.label){
        number_in_subpop[*it_subpop] = 1;
        break;
      }
    }

  }

}

void
Tree::GetSubTree(Sample& sample, Tree& subtree) const{

  std::vector<int> subpop;

  int hap = 0;
  for(std::vector<int>::iterator it_group_of_haplotype = sample.group_of_haplotype.begin(); it_group_of_haplotype != sample.group_of_haplotype.end(); it_group_of_haplotype++){
    bool exists = false;
    for(std::vector<int>::iterator it_group_of_interest = sample.group_of_interest.begin(); it_group_of_interest != sample.group_of_interest.end(); it_group_of_interest++){
      if(*it_group_of_haplotype == *it_group_of_interest){
        exists = true;
        break;
      }
    }
    if(exists){
      subpop.push_back(hap);
    }
    hap++;
  }

  std::vector<int> convert_index, number_in_subpop;
  GetSubTree(subpop, subtree, convert_index, number_in_subpop);

}

void
Tree::GetSubTree(Sample& sample, Tree& subtree,  std::vector<int>& convert_index, std::vector<int>& number_in_subpop) const{

  std::vector<int> subpop;

  int hap = 0;
  for(std::vector<int>::iterator it_group_of_haplotype = sample.group_of_haplotype.begin(); it_group_of_haplotype != sample.group_of_haplotype.end(); it_group_of_haplotype++){
    bool exists = false;
    for(std::vector<int>::iterator it_group_of_interest = sample.group_of_interest.begin(); it_group_of_interest != sample.group_of_interest.end(); it_group_of_interest++){
      if(*it_group_of_haplotype == *it_group_of_interest){
        exists = true;
        break;
      }
    }
    if(exists){
      subpop.push_back(hap);
    }
    hap++;
  }

  GetSubTree(subpop, subtree, convert_index, number_in_subpop);

}

void
Tree::GetSubTree(std::vector<int>& subpop, Tree& subtree, std::vector<int>& convert_index, std::vector<int>& number_in_subpop) const{

  convert_index.resize(nodes.size(), -1);
  std::fill(convert_index.begin(), convert_index.end(), -1);

  //for each node, calculate the number of leaves in subpop below it
  number_in_subpop.resize(nodes.size(), 0.0);
  std::fill(number_in_subpop.begin(), number_in_subpop.end(), 0.0);
  GetNumberOfLeavesInSubpop(*std::prev(nodes.end(),1), subpop, number_in_subpop);

  if(subpop.size() >= (int)(nodes.size() + 1.0)/2.0){

    subtree = *this;
    for(int i = 0; i < (int)nodes.size(); i++){
      convert_index[i] = i;
    }

  }else{

    int node = 0;
    subtree.nodes.resize(2*subpop.size() - 1);

    for(node = 0; node < subpop.size(); node++){
      subtree.nodes[node] = nodes[subpop[node]];
      subtree.nodes[node].label = node;
      convert_index[subpop[node]] = node;
    }

    for(int i = (int)(nodes.size() + 1.0)/2.0; i < (int) nodes.size(); i++){

      int child_left  = (*nodes[i].child_left).label;
      int child_right = (*nodes[i].child_right).label;

			//std::cerr << i << " " << child_left << " " << child_right << " " << number_in_subpop[child_left] << " " << number_in_subpop[child_right] << std::endl;

      if(number_in_subpop[child_left] > 0 && number_in_subpop[child_right] > 0){

				//both children exist in subtree
        assert(convert_index[child_left] != -1);
        assert(convert_index[child_right] != -1);
        subtree.nodes[node] = nodes[i];
        //connect to children.
        subtree.nodes[node].label       = node;
        subtree.nodes[node].child_left  = &subtree.nodes[convert_index[child_left]];
        subtree.nodes[node].child_right = &subtree.nodes[convert_index[child_right]];
        subtree.nodes[convert_index[child_left]].parent  = &subtree.nodes[node];
        subtree.nodes[convert_index[child_right]].parent = &subtree.nodes[node];

        convert_index[i]    = node;
        assert(number_in_subpop[i] > 0);
        node++;

      }else if(number_in_subpop[child_left] > 0){

				//only child left exists in subtree, so i inherits child_left's label
        assert(convert_index[child_left] != -1);
        convert_index[i]                               = convert_index[child_left];
        subtree.nodes[convert_index[i]].branch_length += nodes[i].branch_length;
        subtree.nodes[convert_index[i]].num_events    += nodes[i].num_events;
        assert(number_in_subpop[i] > 0);

      }else if(number_in_subpop[child_right] > 0){

				//only child right exists in subtree, so i inherits child_left's label
        assert(convert_index[child_right] != -1);
        convert_index[i]                               = convert_index[child_right];
        subtree.nodes[convert_index[i]].branch_length += nodes[i].branch_length;
        subtree.nodes[convert_index[i]].num_events    += nodes[i].num_events;
        assert(number_in_subpop[i] > 0);

      }else{

        assert(convert_index[i] == -1);
        assert(number_in_subpop[i] == 0);

      }

    }

    assert(node == 2*subpop.size()-1);
    subtree.nodes[node-1].parent = NULL;
  }

}


/////////////////////////////////

void
MarginalTree::Read(const std::string& line, int N){

  int i = 0;
  while(line[i] != ':'){
    i++;
  }
  i += 2;

  sscanf(line.c_str(), "%d: ", &pos);
  tree.ReadTree(&(line.c_str())[i], N);

}

void
MarginalTree::Read(const std::string& line, int N, std::vector<double>& sample_ages){

  int i = 0;
  while(line[i] != ':'){
    i++;
  }
  i += 2;

  sscanf(line.c_str(), "%d: ", &pos);
  tree.ReadTree(&(line.c_str())[i], N);
  tree.sample_ages = &sample_ages;

}


void
MarginalTree::Dump(std::ofstream& os){

  int parent;
  std::vector<Node>::iterator n_it;
  //std::string line;
  std::stringstream stream;

  n_it = tree.nodes.begin();
  stream << pos << ": ";
  for(; n_it != tree.nodes.end(); n_it++){
    if((*n_it).parent == NULL){
      parent = -1;
    }else{
      parent = (*(*n_it).parent).label;
    }
    stream << parent << ":(" << std::setprecision(5) << (*n_it).branch_length << " "  << std::setprecision(2) << (*n_it).num_events << " " << (*n_it).SNP_begin << " " << (*n_it).SNP_end << ")";
  }
  stream << "\n";

  os << stream.str();

}

void
MarginalTree::Dump(FILE *pfile){

  int parent;
  std::vector<Node>::iterator n_it;

  n_it = tree.nodes.begin();
  fprintf(pfile, "%d: ", pos);
  for(; n_it != tree.nodes.end(); n_it++){
    if((*n_it).parent == NULL){
      parent = -1;
    }else{
      parent = (*(*n_it).parent).label;
    }
    fprintf(pfile, "%d:(%.5f %.3f %d %d) ", parent, (*n_it).branch_length, (*n_it).num_events, (*n_it).SNP_begin, (*n_it).SNP_end);
  }

  fprintf(pfile, "\n");

}



///////////////////////////////////

float 
Correlation::Pearson(const Leaves& set1, const Leaves& set2){

	if(set1.num_leaves == N || set2.num_leaves == N){
		if(set1.num_leaves == set2.num_leaves) return 1;
		return 0;
	}

	float prod = 0.0; 
	std::vector<int>::const_iterator it_set1_member = set1.member.begin();
	std::vector<int>::const_iterator it_set2_member = set2.member.begin();

	const std::vector<int>::const_iterator it_set1_member_end = set1.member.end();
	const std::vector<int>::const_iterator it_set2_member_end = set2.member.end();

	while(it_set1_member != it_set1_member_end && it_set2_member != it_set2_member_end){
		if(*it_set1_member == *it_set2_member){
			prod += 1.0;
			it_set1_member++;
			it_set2_member++;
		}else if(*it_set1_member < *it_set2_member){
			it_set1_member++;
		}else{
			it_set2_member++;
		}
	}

	if(prod == set1.num_leaves && prod == set2.num_leaves) return 1.0;
	float r = prod - set1.num_leaves * (((float) set2.num_leaves)/N_float);
	if(r <= 0.0) return 0.0;

	r /= sqrt( ((((float)set1.num_leaves)/N_float) * (N_float - set1.num_leaves)) * ((((float)set2.num_leaves)/N_float) * (N_float - set2.num_leaves)) );
	assert(!std::isnan(r));

	return r;

}

////////////////////////////////


void 
AncesTree::Read(igzstream& is){

  std::string line;
  seq.clear();
  seq.emplace_back();
  CorrTrees::iterator it_seq = seq.begin();

  int num_tree = 0;
  int i;
  while(num_tree < L){

    getline(is, line);

    i = 0;
    while(line[i] != ':'){
      i++;
    }
    i += 2;
    sscanf(line.c_str(), "%d: ", &(*it_seq).pos);
    (*it_seq).tree.ReadTree(&(line.c_str())[i], N);
    (*it_seq).tree.sample_ages = &sample_ages;

    seq.emplace_back();
    it_seq++;
    num_tree++;
  }

  seq.pop_back();

}

void
AncesTree::Read(const std::string& filename){

  double start_time = time(NULL);
  clock_t begin = clock();

  igzstream is(filename);
  if(is.fail()) is.open(filename + ".gz");
  if(is.fail()){ 
    std::cerr << "Error while opening file " << filename << "(.gz)." << std::endl;
    exit(1);
  }

  std::istringstream is_header;

  char sdummy[30];
  std::string line, tmp;
  //read num_haplotypes
  getline(is, line);
  //sscanf(line.c_str(), "%s %d", sdummy, &N);
  is_header.str(line);
  is_header >> tmp;
  is_header >> N;

  sample_ages.resize(N);
  std::vector<double>::iterator it_sample_ages = sample_ages.begin();
  int i = 0;
  while(is_header >> *it_sample_ages){
    it_sample_ages++;
    i++;
    if(it_sample_ages == sample_ages.end()) break;
  }
  if(i != N) sample_ages.clear();
  getline(is, line);
  sscanf(line.c_str(), "%s %d", sdummy, &L);

  Read(is);
  
  clock_t end = clock();
  double end_time = time(NULL);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //std::cerr << "anc read in: " << elapsed_secs << " CPU secs and " << end_time - start_time << " real secs." << std::endl;

}

void 
AncesTree::ReadBin(FILE* pfile){

  bool has_sample_ages;
  fread(&has_sample_ages, sizeof(bool), 1, pfile);
  fread(&N, sizeof(unsigned int), 1, pfile);
  if(has_sample_ages){
    sample_ages.resize(N);
    fread(&sample_ages[0], sizeof(double), N, pfile);
  }
  fread(&L, sizeof(unsigned int), 1, pfile);

  seq.clear();
  seq.emplace_back();
  CorrTrees::iterator it_seq = seq.begin();

  int num_tree = 0;
  while(num_tree < L){
    fread(&(*it_seq).pos, sizeof(int), 1, pfile);
    (*it_seq).tree.ReadTreeBin(pfile, N);

    seq.emplace_back();
    it_seq++;
    num_tree++;
  }

  seq.pop_back();

}

void
AncesTree::ReadBin(const std::string& filename){

  double start_time = time(NULL);
  clock_t begin = clock();

  FILE* pfile = fopen(filename.c_str(), "rb");
  assert(pfile);
  ReadBin(pfile);
  fclose(pfile);

  clock_t end = clock();
  double end_time = time(NULL);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //std::cerr << "anc read in: " << elapsed_secs << " CPU secs and " << end_time - start_time << " real secs." << std::endl;

}



void
AncesTree::Dump(FILE *pfile){

  double start_time = time(NULL);
  clock_t begin = clock();

  int parent;
  std::vector<Node>::iterator n_it;
  for(CorrTrees::iterator it_seq = seq.begin(); it_seq != seq.end(); it_seq++){

    n_it = (*it_seq).tree.nodes.begin();
    fprintf(pfile, "%d: ", (*it_seq).pos);
    for(; n_it != (*it_seq).tree.nodes.end(); n_it++){
      if((*n_it).parent == NULL){
        parent = -1;
      }else{
        parent = (*(*n_it).parent).label;
      }
      fprintf(pfile, "%d:(%.5f %.3f %d %d) ", parent, (*n_it).branch_length, (*n_it).num_events, (*n_it).SNP_begin, (*n_it).SNP_end);
    }

    fprintf(pfile, "\n");

  }

  clock_t end = clock();
  double end_time = time(NULL);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //std::cerr << "anc dumped in: " << elapsed_secs  << " CPU secs and " << end_time - start_time << " real secs." << std::endl;

}

void
AncesTree::Dump(std::ofstream& os){

  int parent;
  std::vector<Node>::iterator n_it;
  //std::string line;
  std::stringstream stream;
  for(CorrTrees::iterator it_seq = seq.begin(); it_seq != seq.end(); it_seq++){

    n_it = (*it_seq).tree.nodes.begin();
    stream << (*it_seq).pos << ": ";
    for(; n_it != (*it_seq).tree.nodes.end(); n_it++){
      if((*n_it).parent == NULL){
        parent = -1;
      }else{
        parent = (*(*n_it).parent).label;
      }
      stream << parent << ":(" << std::setprecision(5) << (*n_it).branch_length << " "  << std::setprecision(2) << (*n_it).num_events << " " << (*n_it).SNP_begin << " " << (*n_it).SNP_end << ")";
    }

    stream << "\n";
  }

  os << stream.str();

}

void
AncesTree::Dump(const std::string& filename){

  FILE *pfile = std::fopen(filename.c_str(), "w");

  if(pfile == NULL){

    std::cerr << "Error while writing to " << filename << "." << std::endl;

  }else{

    fprintf(pfile, "NUM_HAPLOTYPES %ld ", ((*seq.begin()).tree.nodes.size() + 1)/2);
    for(std::vector<double>::iterator it_sample_ages = sample_ages.begin(); it_sample_ages != sample_ages.end(); it_sample_ages++){
      fprintf(pfile, "%f ", *it_sample_ages);
    }
    fprintf(pfile, "\n");
    fprintf(pfile, "NUM_TREES %ld\n", seq.size());

    Dump(pfile);

    fclose(pfile);

  }

}

void
AncesTree::DumpStream(const std::string& filename){

  std::ofstream os(filename);

  if(os.fail()){

    std::cerr << "Error while opening file " << filename << "." << std::endl;
    exit(1);

  }else{

    os << "NUM_HAPLOTYPES " << ((*seq.begin()).tree.nodes.size() + 1)/2 << " ";
    for(std::vector<double>::iterator it_sample_ages = sample_ages.begin(); it_sample_ages != sample_ages.end(); it_sample_ages++){
      os << *it_sample_ages << " ";
    }
    os << "\n";
    os << "NUM_TREES " << seq.size() << "\n";

    Dump(os);

    os.close();

  }

}


void
AncesTree::DumpBin(FILE *pfile){

  double start_time = time(NULL);
  clock_t begin = clock();

  int parent;
  std::vector<Node>::iterator n_it;
  for(CorrTrees::iterator it_seq = seq.begin(); it_seq != seq.end(); it_seq++){

    n_it = (*it_seq).tree.nodes.begin();
    fwrite(&(*it_seq).pos, sizeof(int), 1,pfile);
    for(; n_it != (*it_seq).tree.nodes.end(); n_it++){
      if((*n_it).parent == NULL){
        parent = -1;
      }else{
        parent = (*(*n_it).parent).label;
      }
      fwrite(&parent, sizeof(int), 1, pfile);
      fwrite(&(*n_it).branch_length, sizeof(double), 1, pfile);
      fwrite(&(*n_it).num_events, sizeof(float), 1, pfile);
      fwrite(&(*n_it).SNP_begin, sizeof(int), 1, pfile);
      fwrite(&(*n_it).SNP_end, sizeof(int), 1, pfile);
    }

  }

  clock_t end = clock();
  double end_time = time(NULL);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  //std::cerr << "anc dumped in: " << elapsed_secs  << " CPU secs and " << end_time - start_time << " real secs." << std::endl;

}

void
AncesTree::DumpBin(const std::string& filename){

  FILE *pfile = std::fopen(filename.c_str(), "wb");
  assert(pfile);

  if(pfile == NULL){

    std::cerr << "Error while writing to " << filename << "." << std::endl;

  }else{

    unsigned int N = ((*seq.begin()).tree.nodes.size() + 1)/2;
    unsigned int num_trees = seq.size();

    bool has_sample_ages = (sample_ages.size() > 0);
    fwrite(&has_sample_ages, sizeof(bool), 1, pfile);
    fwrite(&N, sizeof(unsigned int), 1, pfile);
    if(has_sample_ages){
      fwrite(&sample_ages[0], sizeof(double), N, pfile);
    }
    fwrite(&num_trees, sizeof(unsigned int), 1, pfile);

    DumpBin(pfile);

    fclose(pfile);

  }


}


////////////////////////////////////

void 
AncesTree::ReadMsPrime(const std::string& filename_el){

  igzstream is(filename_el);

  if(is.fail()){
    std::cerr << "Error while opening file " << filename_el << "." << std::endl;
    exit(1);
  } 


  //////
  //#comment
  //N L
  //pos1
  //Tree1
  //pos2
  //Tree2
  //etc

  std::string line;
  getline(is, line); //skip header
  getline(is, line); //read N and L
  std::stringstream posline(line);
  int num_nodes, num_snp;
  posline >> num_nodes >> num_snp;

  CorrTrees::iterator it_seq;
  for(int snp = 0; snp < num_snp; snp++){
    seq.emplace_back();
    it_seq = std::prev(seq.end(),1);
    //is >> (*it_seq).pos; 
    getline(is, line);
    (*it_seq).pos = std::stoi(line); 
    //std::cerr << "SNP:" << snp << ",POS: " << (*it_seq).pos << std::endl;
    (*it_seq).tree.GetMsPrime(is, num_nodes);
  }

  is.close();

}

void
AncesTree::ReadArgweaverSMC(const std::string& filename){

  int N       = 0;
  int N_total = 0;

  igzstream is(filename);
  std::string line;

  std::string sdummy, newick;
  int idummy;
  int node_index;
  std::vector<int> convert_index;

  //get number of nodes
  std::getline(is, line);
  std::stringstream readN(line);
  readN >> sdummy;
  while(readN >> idummy){
    convert_index.push_back(idummy-1);
    N++;
  }
 
  N_total = 2*N-1;
  convert_index.resize(N_total);
  for(int n = N; n < N_total;n++){
    convert_index[n] = n;
  }
  //std::cerr << N << std::endl;

  //std::getline(is, line);
  int num_tree = 0;
  seq.emplace_back();
  CorrTrees::iterator it_seq = seq.begin();
  while(std::getline(is, line)){ 

    assert(std::getline(is, line));

    std::stringstream tree_stream(line); 
    (*it_seq).tree.nodes.resize(N_total);

    tree_stream >> sdummy;
    tree_stream >> (*it_seq).pos;
    tree_stream >> idummy;

    //std::cerr << (*it_seq).pos << std::endl; 
    //(*it_seq).pos -= 1.0;
    tree_stream >> newick; 

    //need to convert newick into my tree data structure
    int i = 0;
    while(newick.size() > 0){ 
      int startpos, endpos;
      std::string index_child1, index_child2, index_parent;
      std::string branch_length1, branch_length2;

      //std::cerr << newick << std::endl << std::endl;
      while(newick[i] == '(') i++;
      startpos = i-1;
      while(newick[i] != ':'){
        index_child1 += newick[i];
        i++;
      }
      i++;
      while(newick[i] != '['){
        branch_length1 += newick[i];
        i++;
      }
      while(newick[i] != ',') i++;
      i++;
      if(newick[i] != '('){
        while(newick[i] != ':'){
          index_child2 += newick[i];
          i++;
        } 
        i++;
        while(newick[i] != '['){
          branch_length2 += newick[i];
          i++;
        }
        while(newick[i] != ')') i++;
        i++;
        endpos = i;
        while( newick[i] != ':' && newick[i] != '['){
          index_parent += newick[i];
          i++;
        }
        //std::cerr << newick[startpos] << " " << newick[endpos-1] << " " << newick[i] << std::endl; 
        //std::cerr << index_child1 << " " << branch_length1 << " " << index_child2 << " " << branch_length2 << " " << index_parent << std::endl;

        int child_left = convert_index[stoi(index_child1)], child_right = convert_index[stoi(index_child2)], parent = convert_index[stoi(index_parent)];
        //std::cerr << child_left << " " << bl1 << " " << child_right << " " << bl2 << " " << parent << std::endl << std::endl;

        (*it_seq).tree.nodes[child_left].label  = child_left;
        (*it_seq).tree.nodes[child_right].label = child_right;
        (*it_seq).tree.nodes[parent].label      = parent;

        (*it_seq).tree.nodes[child_left].parent  = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[child_right].parent = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[child_left];
        (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[child_right];

        (*it_seq).tree.nodes[child_left].branch_length = stof(branch_length1);
        (*it_seq).tree.nodes[child_right].branch_length = stof(branch_length2);

        if(newick[i] == '['){ 
          break;
        }
        newick.replace(startpos, endpos - startpos, "");  
        i = 0;
      }

    }

    //make sure that node root is root
    //find root
    int root = N_total-1;
    if((*it_seq).tree.nodes[root].parent != NULL){

      int real_root = 0;
      for(std::vector<Node>::iterator it_node = (*it_seq).tree.nodes.begin(); it_node != (*it_seq).tree.nodes.end(); it_node++){
        if((*it_node).parent == NULL){
          break;
        }
        real_root++;
      }

      Node n1 = (*it_seq).tree.nodes[real_root];
      Node n2 = (*it_seq).tree.nodes[root];

      bool left = false;
      if((*(*n2.parent).child_left).label == root){
        left = true;
      }else{ 
        assert((*(*n2.parent).child_right).label == root);
      }

      if((*n2.parent).label == real_root){

        (*it_seq).tree.nodes[real_root]       = n2;
        (*it_seq).tree.nodes[real_root].label = real_root; 
        (*it_seq).tree.nodes[root]       = n1;
        (*it_seq).tree.nodes[root].label = root;

        if(left){
          (*it_seq).tree.nodes[root].child_left  = &(*it_seq).tree.nodes[real_root];
        }else{ 
          (*it_seq).tree.nodes[root].child_right = &(*it_seq).tree.nodes[real_root];
        }
        (*(*it_seq).tree.nodes[real_root].child_left).parent  = &(*it_seq).tree.nodes[real_root];
        (*(*it_seq).tree.nodes[real_root].child_right).parent = &(*it_seq).tree.nodes[real_root];
        (*(*it_seq).tree.nodes[root].child_left).parent       = &(*it_seq).tree.nodes[root];
        (*(*it_seq).tree.nodes[root].child_right).parent      = &(*it_seq).tree.nodes[root];

      }else{

        (*it_seq).tree.nodes[real_root]       = n2;
        (*it_seq).tree.nodes[real_root].label = real_root; 
        if(left){
          (*(*it_seq).tree.nodes[real_root].parent).child_left  = &(*it_seq).tree.nodes[real_root];
        }else{ 
          (*(*it_seq).tree.nodes[real_root].parent).child_right = &(*it_seq).tree.nodes[real_root];
        }
        (*(*it_seq).tree.nodes[real_root].child_left).parent  = &(*it_seq).tree.nodes[real_root];
        (*(*it_seq).tree.nodes[real_root].child_right).parent = &(*it_seq).tree.nodes[real_root];

        (*it_seq).tree.nodes[root]       = n1;
        (*it_seq).tree.nodes[root].label = root;
        (*(*it_seq).tree.nodes[root].child_left).parent  = &(*it_seq).tree.nodes[root];
        (*(*it_seq).tree.nodes[root].child_right).parent = &(*it_seq).tree.nodes[root];
   
      }
      
      /* 
      std::cerr << (*it_seq).tree.nodes[root].label << " " << ((*it_seq).tree.nodes[root].parent == NULL) << " " << (*(*(*it_seq).tree.nodes[root].child_left).parent).label << " " << (*(*(*it_seq).tree.nodes[root].child_right).parent).label << std::endl;
  
      std::cerr << (*it_seq).tree.nodes[real_root].label << " " << ((*it_seq).tree.nodes[real_root].parent == NULL) << " " << (*(*(*it_seq).tree.nodes[real_root].child_left).parent).label << " " << (*(*(*it_seq).tree.nodes[real_root].child_right).parent).label << std::endl;

      std::cerr << (*(*(*it_seq).tree.nodes[real_root].parent).child_left).label << " " << (*(*(*it_seq).tree.nodes[real_root].parent).child_right).label << std::endl;

      std::cerr << "------------" << std::endl;
      if(num_tree == 789) exit(1);     
      */
      

    }




    num_tree++;
    seq.emplace_back();
    it_seq++;
    //std::getline(is, line); //skipping every second line
  }

  seq.pop_back();

  //std::cerr << seq.size() << std::endl;
}

void
AncesTree::ReadRent(const std::string& filename, float Ne){

  igzstream is(filename);
  std::string line;

  std::string newick; 

  int N = -1;
  int N_total;
  int i;

  std::vector<float> coordinates;
  int num_tree = 1;
  seq.resize(1);
  CorrTrees::iterator it_seq = seq.begin();
  while(std::getline(is, line)){ 

    i = 0;
    if(N == -1){
      N = 0;
      while(i < line.size()){
        if(line[i] == ',') N++;
        i++;
      }
      N += 1;
      N_total = 2*N-1;

    }

    //std::cerr << line << " " << N_total << std::endl;
    //exit(1);

    std::stringstream tree_stream(line); 
    (*it_seq).tree.nodes.resize(N_total);
    tree_stream >> (*it_seq).pos;
    tree_stream >> newick; 

    //need to convert newick into my tree data structure
    i = 0;
    int node = N+1;
    int count_bracket = 0, count_comma = 0;
    while(node <= N_total){ 
      int startpos, endpos;
      std::string index_child1, index_child2, index_parent;
      std::string branch_length1, branch_length2;

      while(newick[i] == '(') i++;
      startpos = i;

      while(newick[i] != ':'){
        index_child1 += newick[i];
        i++;
      }
      i++;
      while(newick[i] != ',' && i < newick.size()){
        branch_length1 += newick[i];
        i++;
      }
      i++;
      if(newick[i] != '(' && i < newick.size()){
        while(newick[i] != ':'){
          index_child2 += newick[i];
          i++;
        } 
        i++;
        while(newick[i] != ')' && i < newick.size()){
          branch_length2 += newick[i];
          i++;
        }
        i++;
        endpos = i;

        //std::cerr << index_child1 << " " << index_child2 << " " << branch_length1 << " " << branch_length2 << std::endl;
        int child_left = stoi(index_child1)-1, child_right = stoi(index_child2)-1, parent = node-1;

        (*it_seq).tree.nodes[child_left].label  = child_left;
        (*it_seq).tree.nodes[child_right].label = child_right;
        (*it_seq).tree.nodes[parent].label      = parent;

        (*it_seq).tree.nodes[child_left].parent  = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[child_right].parent = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[child_left];
        (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[child_right];

        (*it_seq).tree.nodes[child_left].branch_length  = stof(branch_length1) * Ne;
        (*it_seq).tree.nodes[child_right].branch_length = stof(branch_length2) * Ne;

        //std::cerr << newick << std::endl;
        //std::cerr << newick.substr(startpos-1, endpos-startpos+1) << std::endl;
        newick.replace(startpos-1, endpos - startpos + 1, std::to_string(node)); 
        count_bracket = 0; 
        count_comma = 0;
        i = 0;
        for(; i < newick.size(); i++){
          if(newick[i] == '(') count_bracket++;
          if(newick[i] == ',') count_comma++;
        }
        //std::cerr << count_comma << " " << count_bracket << std::endl;
        //assert(count_comma == count_bracket);
        if(count_comma != count_bracket){
          break;
        }
        //std::cerr << newick << std::endl; 
        i = 0;
        node++;
      }
    }

    bool everyone_has_parent = true;
    for(int i = 0; i < N_total - 1; i++){
      if((*it_seq).tree.nodes[i].parent == NULL){
        everyone_has_parent = false;
        break;
      } 
    }
    //std::cerr << everyone_has_parent << std::endl;
    if(node != N_total + 1 || count_comma != count_bracket || !everyone_has_parent){
      seq.pop_back();
      num_tree--;
      it_seq--;
    }

    /*
       coordinates.resize(N_total); 
       (*it_seq).tree.GetCoordinates(coordinates);
       std::cerr << coordinates[N_total-1] << " " << coordinates[N_total-1]/Ne << std::endl; 
       */

    num_tree++;
    seq.emplace_back();
    it_seq = std::prev(seq.end(),1);
    //it_seq++;
  }

  seq.pop_back();

  //std::cerr << num_tree << " " << (*std::prev(seq.end(),1)).pos << std::endl;
}

void
AncesTree::ReadNewick(const std::string& filename, float Ne){

  igzstream is(filename);
  std::string line;

  std::string newick; 

  int N = -1;
  int N_total;
  int i;

  std::vector<float> coordinates;
  int num_tree = 1;
  seq.resize(1);
  CorrTrees::iterator it_seq = seq.begin();
  while(std::getline(is, line)){ 

    i = 0;
    if(N == -1){
      N = 0;
      while(i < line.size()){
        if(line[i] == ',') N++;
        i++;
      }
      N += 1;
      N_total = 2*N-1;

    }

    //exit(1);

    std::stringstream tree_stream(line); 
    (*it_seq).tree.nodes.resize(N_total);
    tree_stream >> (*it_seq).pos;
    tree_stream >> newick; 

    //need to convert newick into my tree data structure
    i = 0;
    int node = N;
    int count_bracket = 0, count_comma = 0;
    while(node < N_total){

			//std::cerr << newick << std::endl;

      int startpos, endpos;
      std::string index_child1, index_child2, index_parent;
      std::string branch_length1, branch_length2;

      while(newick[i] == '(') i++;
      startpos = i;

      while(newick[i] != ':'){
        index_child1 += newick[i];
        i++;
      }
      i++;
      while(newick[i] != ',' && i < newick.size()){
        branch_length1 += newick[i];
        i++;
      }
      i++;
      if(newick[i] != '(' && i < newick.size()){
        while(newick[i] != ':'){
          index_child2 += newick[i];
          i++;
        } 
        i++;
        while(newick[i] != ')' && i < newick.size()){
          branch_length2 += newick[i];
          i++;
        }
        i++;
        endpos = i;

        //std::cerr << index_child1 << " " << index_child2 << " " << branch_length1 << " " << branch_length2 << std::endl;
        int child_left = stoi(index_child1), child_right = stoi(index_child2), parent = node;

        (*it_seq).tree.nodes[child_left].label  = child_left;
        (*it_seq).tree.nodes[child_right].label = child_right;
        (*it_seq).tree.nodes[parent].label      = parent;

        (*it_seq).tree.nodes[child_left].parent  = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[child_right].parent = &(*it_seq).tree.nodes[parent];
        (*it_seq).tree.nodes[parent].child_left  = &(*it_seq).tree.nodes[child_left];
        (*it_seq).tree.nodes[parent].child_right = &(*it_seq).tree.nodes[child_right];

        (*it_seq).tree.nodes[child_left].branch_length  = stof(branch_length1) * Ne;
        (*it_seq).tree.nodes[child_right].branch_length = stof(branch_length2) * Ne;

        //std::cerr << newick << std::endl;
        //std::cerr << newick.substr(startpos-1, endpos-startpos+1) << std::endl;
        newick.replace(startpos-1, endpos - startpos + 1, std::to_string(node)); 
        count_bracket = 0; 
        count_comma = 0;
        i = 0;
        for(; i < newick.size(); i++){
          if(newick[i] == '(') count_bracket++;
          if(newick[i] == ',') count_comma++;
        }
        //assert(count_comma == count_bracket);
        if(count_comma != count_bracket){
          break;
        }
        //std::cerr << newick << std::endl; 
        i = 0;
        node++;
      }
    }

    bool everyone_has_parent = true;
    for(int i = 0; i < N_total - 1; i++){
      if((*it_seq).tree.nodes[i].parent == NULL){
        everyone_has_parent = false;
        break;
      } 
    }
    //std::cerr << everyone_has_parent << std::endl;
    if(node != N_total || count_comma != count_bracket || !everyone_has_parent){
      seq.pop_back();
      num_tree--;
      it_seq--;
    }

    /*
       coordinates.resize(N_total); 
       (*it_seq).tree.GetCoordinates(coordinates);
       std::cerr << coordinates[N_total-1] << " " << coordinates[N_total-1]/Ne << std::endl; 
       */

    num_tree++;
    seq.emplace_back();
    it_seq = std::prev(seq.end(),1);
    //it_seq++;
  }

  seq.pop_back();

}



//Functions for assiciating equivalent branches across adjacent trees
void
AncesTree::BranchAssociation(const Tree& ref_tree, const Tree& tree, std::vector<int>& equivalent_branches, std::vector<std::vector<int>>& potential_branches, int N, int N_total, float threshold_brancheq){

	////
	equivalent_branches.resize(N_total);
	std::fill(equivalent_branches.begin(), equivalent_branches.end(), -1);
	std::vector<int> equivalent_branches_ref(N_total, -1);

	Correlation cor(N);

	//1. calculate leaves
	std::vector<Leaves> tr_leaves;
	std::vector<Leaves> rtr_leaves;

	tree.FindAllLeaves(tr_leaves);
	ref_tree.FindAllLeaves(rtr_leaves);

	//2. Sort nodes according to number of leaves and store in std::vector of size(N_total)

	std::vector<int> sorted_branches(N_total);
	std::size_t n(0);
	std::generate(std::begin(sorted_branches), std::end(sorted_branches), [&]{ return n++; });
	std::sort(std::begin(sorted_branches), std::end(sorted_branches), [&](int i1, int i2) { return rtr_leaves[i1].num_leaves < rtr_leaves[i2].num_leaves; } );

	std::vector<int> index_sorted_branches(N,0); //this vector is for accessing sorted_branches
	for(std::vector<Leaves>::iterator it = rtr_leaves.begin(); it != std::prev(rtr_leaves.end(),1); it++){
		index_sorted_branches[(*it).num_leaves]++;
	}
	int cum = 0;
	for(std::vector<int>::iterator it = index_sorted_branches.begin(); it != index_sorted_branches.end(); it++){
		*it += cum;
		cum = *it;
	}

	//3. For each branch, search for exact equivalent branch in ref_tree
	//   record which branches are unpaired

	std::vector<int> unpaired_branches;

	Node parent, ref_parent;

	//treat leaves separately (this is faster)
	for(int i = 0; i < N; i++){

		if(equivalent_branches[i] == -1){

			parent = *tree.nodes[i].parent;
			ref_parent = *ref_tree.nodes[i].parent;

			if((*parent.child_left).label == i){

				int child_right = (*parent.child_right).label;
				if(child_right < N){
					if( child_right == (*ref_parent.child_right).label ){ 
						equivalent_branches[i]     = i;
						equivalent_branches_ref[i] = i;
						equivalent_branches[child_right]     = child_right;
						equivalent_branches_ref[child_right] = child_right;
					}else if( child_right == (*ref_parent.child_left).label ){ 
						equivalent_branches[i]     = i;
						equivalent_branches_ref[i] = i;
						equivalent_branches[child_right]     = child_right;
						equivalent_branches_ref[child_right] = child_right;
					} 
				}else{
					if(cor.Pearson(tr_leaves[parent.label], rtr_leaves[ref_parent.label]) >= threshold_brancheq){
						equivalent_branches[i]     = i;
						equivalent_branches_ref[i] = i;
					}
				}

			}else{

				int child_left = (*parent.child_left).label;
				if(child_left < N){
					if( child_left == (*ref_parent.child_left).label ){ 
						equivalent_branches[i]     = i;
						equivalent_branches_ref[i] = i;
						equivalent_branches[child_left]     = child_left;
						equivalent_branches_ref[child_left] = child_left;
					}else if( child_left == (*ref_parent.child_right).label ){ 
						equivalent_branches[i]     = i;
						equivalent_branches_ref[i] = i;
						equivalent_branches[child_left]     = child_left;
						equivalent_branches_ref[child_left] = child_left;
					} 
				}else{
					if(cor.Pearson(tr_leaves[parent.label], rtr_leaves[ref_parent.label]) >= threshold_brancheq){
						equivalent_branches[i]     = i;
						equivalent_branches_ref[i] = i;
					}
				}

			}

		}

	}

	for(int i = N; i < N_total-1; i++){ 

		if(cor.Pearson(tr_leaves[i], rtr_leaves[i]) >= 0.9999 && cor.Pearson(tr_leaves[(*tree.nodes[i].parent).label], rtr_leaves[(*ref_tree.nodes[i].parent).label]) >= 0.9999){     
			//branches i and i are equivalent
			equivalent_branches[i] = i;
			equivalent_branches_ref[i] = i;
		}

		if(equivalent_branches[i] == -1){
			int num_leaves = tr_leaves[i].num_leaves;
			for(std::vector<int>::iterator it = std::next(sorted_branches.begin(),index_sorted_branches[num_leaves-1]); it != std::next(sorted_branches.begin(), index_sorted_branches[num_leaves]); it++ ){
				if(cor.Pearson(tr_leaves[i], rtr_leaves[*it]) >= 0.9999 && cor.Pearson(tr_leaves[(*tree.nodes[i].parent).label], rtr_leaves[(*ref_tree.nodes[*it].parent).label]) >= 0.9999){     
					//branches i and *it are equivalent
					equivalent_branches[i]       = *it;
					equivalent_branches_ref[*it] = i;

					break;  
				}
			}
		}

		if(equivalent_branches[i] == -1){
			unpaired_branches.push_back(i);
		}

	}

	//4. For the nodes that are not paired up, find all possible pairs above threshold 
	std::vector<EquivalentNode> possible_pairs;

	for(std::vector<int>::iterator it_unpaired = unpaired_branches.begin(); it_unpaired != unpaired_branches.end(); it_unpaired++){

		int num_leaves = tr_leaves[*it_unpaired].num_leaves - 1;
		for(std::vector<int>::iterator it_k = potential_branches[num_leaves].begin(); it_k != potential_branches[num_leaves].end(); it_k++){

			for(std::vector<int>::iterator it = std::next(sorted_branches.begin(),index_sorted_branches[*it_k-1]); it != std::next(sorted_branches.begin(), index_sorted_branches[*it_k]); it++ ){
				if(equivalent_branches_ref[*it] == -1){
					float score = cor.Pearson(tr_leaves[*it_unpaired], rtr_leaves[*it]);
					if(score >= threshold_brancheq && cor.Pearson(tr_leaves[(*tree.nodes[*it_unpaired].parent).label], rtr_leaves[(*ref_tree.nodes[*it].parent).label]) >= threshold_brancheq){     
						//branches i and *it are equivalent
						possible_pairs.push_back(EquivalentNode(*it_unpaired, *it, score));
					}
				}
			}

		}

	}


	//5. Sort possible_pairs by score
	std::sort(std::begin(possible_pairs), std::end(possible_pairs), std::greater<EquivalentNode>());

	//6. Pair up approximate matches
	for(std::vector<EquivalentNode>::iterator it = possible_pairs.begin(); it != possible_pairs.end(); it++){
		if(equivalent_branches[(*it).node1] == -1 && equivalent_branches_ref[(*it).node2] == -1){
			equivalent_branches[(*it).node1]     = (*it).node2;
			equivalent_branches_ref[(*it).node2] = (*it).node1;
		}
	}

}

//Carry over information across equivalent branches
void
AncesTree::AssociateEquivalentBranches(){

	CorrTrees::iterator it_seq_prev;
	CorrTrees::iterator it_seq; 

	int N_total = (*seq.begin()).tree.nodes.size();
	int N       = (N_total+1.0)/2.0;

	//Associate branches
	//Pre calculate how many descendants a branch needs to be equivalent
	float threshold_brancheq = 0.95;
	//float threshold_brancheq = 1.0;
	std::vector<std::vector<int>> potential_branches;
	//the number of leaves a branch needs to be equivalent
	potential_branches.resize(N);
	float threshold_inv = 1/(threshold_brancheq * threshold_brancheq);
	float N_float = N;
	for(int i = 1; i <= N; i++){
		potential_branches[i-1].push_back(i);
		//for branches with i+1 leaves, list the number of leaves a potential equivalent branch needs
		for(int j = i+1; j <= N; j++){
			if(threshold_inv >= j/(N_float-j) * ((N_float-i)/i) ){
				potential_branches[i-1].push_back(j);
				potential_branches[j-1].push_back(i);
			}
		}
	}

	/////
	// Find equivalent branches

	it_seq_prev = seq.begin();
	it_seq      = std::next(it_seq_prev,1); 

	std::vector<std::vector<int>> equivalent_branches;
	std::vector<std::vector<int>>::iterator it_equivalent_branches;
	std::vector<std::vector<int>>::reverse_iterator rit_equivalent_branches;

	for(; it_seq != seq.end();){
		equivalent_branches.emplace_back();
		it_equivalent_branches = std::prev(equivalent_branches.end(),1);
		BranchAssociation((*it_seq_prev).tree, (*it_seq).tree, *it_equivalent_branches, potential_branches, N, N_total, threshold_brancheq); //O(N^2) 
		it_seq++;
		it_seq_prev++;
	}  

	///////////////////////////////////////////
	//Now carry over information on branches, starting from the first tree.

	it_equivalent_branches = equivalent_branches.begin();
	std::vector<Node>::iterator it_nodes;

	it_seq_prev = seq.begin();
	it_seq      = std::next(it_seq_prev,1); 

	for(; it_seq != seq.end();){
		it_nodes = (*it_seq).tree.nodes.begin();
		for(std::vector<int>::iterator it = (*it_equivalent_branches).begin(); it != (*it_equivalent_branches).end(); it++){
			if(*it != -1){
				(*it_nodes).num_events += (*it_seq_prev).tree.nodes[*it].num_events;
				(*it_nodes).SNP_begin   = (*it_seq_prev).tree.nodes[*it].SNP_begin;
			}
			it_nodes++;
		}
		it_equivalent_branches++;

		it_seq++;
		it_seq_prev++;
	}
	assert(it_equivalent_branches == equivalent_branches.end());

	///////////////////////////////////////////
	//Now go from the last tree to the first

	rit_equivalent_branches = equivalent_branches.rbegin();

	CorrTrees::reverse_iterator rit_seq_next;
	CorrTrees::reverse_iterator rit_seq; 
	rit_seq_next = seq.rbegin();
	rit_seq      = std::next(rit_seq_next,1); 
	for(; rit_seq != seq.rend();){
		it_nodes = (*rit_seq_next).tree.nodes.begin();
		for(std::vector<int>::iterator it = (*rit_equivalent_branches).begin(); it != (*rit_equivalent_branches).end(); it++){
			if(*it != -1){
				(*rit_seq).tree.nodes[*it].num_events = (*it_nodes).num_events;
				(*rit_seq).tree.nodes[*it].SNP_end    = (*it_nodes).SNP_end;
			}
			it_nodes++;
		}
		rit_equivalent_branches++;

		rit_seq++;
		rit_seq_next++;
	}

	assert(rit_equivalent_branches == equivalent_branches.rend());

}





