#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"
#include <ctime>

float           
GetCoords(int node, Tree& tree, int branch, char m, std::vector<float>::iterator& it_dertimes, std::vector<float>::iterator& it_anctimes){

	float coordinate = 0.0;
	if(tree.nodes[node].child_left != NULL){

		int child_left  = (*tree.nodes[node].child_left).label;
		int child_right = (*tree.nodes[node].child_right).label;

		if(child_left == branch || m == 'd'){
			coordinate = GetCoords(child_left, tree, branch, 'd', it_dertimes, it_anctimes);
		}else{
			coordinate = GetCoords(child_left, tree, branch, 'a', it_dertimes, it_anctimes);
		}

		if(child_right == branch || m == 'd'){
			coordinate = GetCoords(child_right, tree, branch, 'd', it_dertimes, it_anctimes);
		}else{
			coordinate = GetCoords(child_right, tree, branch, 'a', it_dertimes, it_anctimes);
		}
		coordinate += tree.nodes[child_right].branch_length;

		if(child_left != branch && child_right != branch){
			if(m == 'a'){
				*it_anctimes = coordinate;
				it_anctimes++;
			}else{
				*it_dertimes = coordinate;
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

int
ReadNewick(std::string line, int& bp_start, int& bp_end, Tree& tree){

	std::string newick, dummy; 

	int N = -1;
	int N_total;
	int i;

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

	std::stringstream tree_stream(line); 
	tree.nodes.resize(N_total);
	tree_stream >> dummy;
	tree_stream >> bp_start;
	tree_stream >> bp_end;
	tree_stream >> dummy;
	tree_stream >> newick; 

	//need to convert newick into my tree data structure
	i = 0;
	int node = N;
	int count_bracket = 0, count_comma = 0;
	while(node < N_total){ 
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

			int child_left = stoi(index_child1), child_right = stoi(index_child2), parent = node;

			tree.nodes[child_left].label  = child_left;
			tree.nodes[child_right].label = child_right;
			tree.nodes[parent].label      = parent;

			tree.nodes[child_left].parent  = &tree.nodes[parent];
			tree.nodes[child_right].parent = &tree.nodes[parent];
			tree.nodes[parent].child_left  = &tree.nodes[child_left];
			tree.nodes[parent].child_right = &tree.nodes[child_right];

			tree.nodes[child_left].branch_length  = stof(branch_length1);
			tree.nodes[child_right].branch_length = stof(branch_length2);

			newick.replace(startpos-1, endpos - startpos + 1, std::to_string(node)); 
			count_bracket = 0; 
			count_comma = 0;
			i = 0;
			for(; i < newick.size(); i++){
				if(newick[i] == '(') count_bracket++;
				if(newick[i] == ',') count_comma++;
			}
			if(count_comma != count_bracket){
				break;
			}
			i = 0;
			node++;
		}
	}

	bool everyone_has_parent = true;
	for(int i = 0; i < N_total - 1; i++){
		if(tree.nodes[i].parent == NULL){
			everyone_has_parent = false;
			break;
		} 
	}
	if(node != N_total || count_comma != count_bracket || !everyone_has_parent){
		std::cerr << "Failed to read tree at bp " << bp_start << std::endl;
		return 1;
	}

	return 0;

}

void 
ConvertNewickToTimeb(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("input") || !options.count("anc_genome") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, anc_genome, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Converts newick/sites to timeb." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Convert newick/sites to timeb for file" << options["input"].as<std::string>() << ".newick/sites..." << std::endl;


	fasta anc_genome;
	anc_genome.Read(options["anc_genome"].as<std::string>());

	std::ifstream is_newick(options["input"].as<std::string>() + ".newick");
	if(is_newick.fail()){
		is_newick.open(options["input"].as<std::string>() + ".newick.gz");
	}
	if(is_newick.fail()){
		std::cerr << "Failed to open file " << options["input"].as<std::string>() + ".newick" << std::endl;
		exit(1);
	}

	std::ifstream is_sites(options["input"].as<std::string>() + ".sites");
	if(is_sites.fail()){
		is_sites.open(options["input"].as<std::string>() + ".sites.gz");
	}
	if(is_sites.fail()){
		std::cerr << "Failed to open file " << options["input"].as<std::string>() + ".sites" << std::endl;
		exit(1);
	}
	std::string line_newick, line_sites, dummy;

  //////////
	//read entire sites file
	//header
	assert(getline(is_sites, line_sites));
	std::stringstream is_line(line_sites);
  is_line >> dummy;
	int N = 0;
	while(is_line >> dummy) N++;
	assert(getline(is_sites, line_sites));
  
	int i = 0;
	int L = 1000;
	std::vector<int> pos(L);
	std::vector<char> anc_allele(L), der_allele(L);
	std::vector<Leaves> hap(L);
	while(getline(is_sites, line_sites)){

		if(pos.size() - i < 100){
      pos.resize(pos.size() + L);
			hap.resize(hap.size() + L);
			anc_allele.resize(anc_allele.size() + L);
			der_allele.resize(der_allele.size() + L);
		}

		std::string line_hap;
		std::stringstream isites(line_sites); 
		isites >> pos[i];
		isites >> line_hap;
	
		char a0 = anc_genome.seq[pos[i]-1];
		hap[i].member.resize(N);
		anc_allele[i] = a0;
		for(int j = 0; j < N; j++){
      if(line_hap[j] == a0){
        hap[i].member[j] = 0;
			}else{
				der_allele[i] = line_hap[j];
        hap[i].member[j] = 1;
				hap[i].num_leaves++;
			}
		}
    i++;

	}
	pos.resize(i);
	hap.resize(i);
	anc_allele.resize(i);
	der_allele.resize(i);
	L = i;

	//////////

	//file format: chr bp_start bp_end sample newick_string
	AncesTree anc;
	std::vector<int> num_samples_per_tree, tree_start_index, tree_pos;
	int bp = -1;
	i = -1;
	int j = 0;
	assert(getline(is_newick, line_newick)); //header
	while(getline(is_newick, line_newick)){
		Tree tree;
		int bp_start;
		int bp_end;
		int ret = ReadNewick(line_newick, bp_start, bp_end, tree);
		if(ret == 0){
			anc.seq.emplace_back(MarginalTree(bp_start, tree));
			if(bp_start > bp){
				num_samples_per_tree.push_back(1);
				tree_start_index.push_back(j);
				tree_pos.push_back(bp_start);
				bp = bp_start;
				i++;
			}else if(bp_start == bp){
				num_samples_per_tree[i]++; 
			}else{
				std::cerr << "Trees are not sorted by bp in newick file" << std::endl;
				exit(1);
			}
			j++;
		}
	}

	int num_samples = num_samples_per_tree[0];
  for(i = 0; i < num_samples_per_tree.size(); i++){
    if(num_samples_per_tree[i] != num_samples){
      std::cerr << "Error: num_samples per tree is not the same." << std::endl;
			exit(1);
		}
	}

	Data data(N, L);
	AncesTreeBuilder ancbuilder(data);
	CorrTrees::iterator it_seq = anc.seq.begin();

	//use anc and hap to generate timeb file
	std::string filename = options["output"].as<std::string>() + ".timeb";
	FILE* fp = fopen(filename.c_str(), "wb");

	int num_mapping_SNPs = pos.size();
	fwrite(&num_mapping_SNPs, sizeof(int), 1, fp);
	fwrite(&num_samples, sizeof(int), 1, fp);

	for(i = 0; i < pos.size(); i++){

		//find tree corresponding to pos[i]
		j = 0;
		if(j < tree_pos.size()){
			while(tree_pos[j] <= pos[i]){
				j++;
				if(j == tree_pos.size()) break;
			}
		}
		j--;

		int DAF = hap[i].num_leaves;
		std::vector<float> anctimes(num_samples*(data.N-DAF-1), 0.0);
		std::vector<float> dertimes(num_samples*(DAF-1), 0.0);
		std::vector<float>::iterator it_anctimes = anctimes.begin(), it_dertimes = dertimes.begin();

		it_seq = std::next(anc.seq.begin(), tree_start_index[j]);
		for(int k = 0; k < num_samples; k++){
			assert(ancbuilder.IsSNPMapping((*it_seq).tree, hap[i], 0) == 1);
			int branch = *ancbuilder.mutations.info[0].branch.begin();

			std::vector<float>::iterator it_anctimes_s = it_anctimes, it_dertimes_s = it_dertimes;
			GetCoords(2*data.N-2, (*it_seq).tree, branch, 'a', it_dertimes, it_anctimes);
			assert(std::next(it_anctimes_s, data.N-DAF-1) == it_anctimes);
			assert(std::next(it_dertimes_s, DAF-1) == it_dertimes);
			std::sort(it_anctimes_s, it_anctimes);
			std::sort(it_dertimes_s, it_dertimes);

			it_seq++;
		}

		//WriteBinary(anctimes, dertimes, fp);
		//BP, DAF, N, 
		//dump
		fwrite(&pos[i], sizeof(int), 1, fp);
		fwrite(&anc_allele[i], sizeof(char), 1, fp);
		fwrite(&der_allele[i], sizeof(char), 1, fp);
		fwrite(&DAF, sizeof(int), 1, fp);
		fwrite(&data.N, sizeof(int), 1, fp);
		fwrite(&anctimes[0], sizeof(float), num_samples*(data.N-DAF-1), fp);
		fwrite(&dertimes[0], sizeof(float), num_samples*(DAF-1), fp);

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

}
