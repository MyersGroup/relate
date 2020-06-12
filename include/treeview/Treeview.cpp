#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"
#include <ctime>


void TraverseTree(Node& node, Tree& tree, double& x_coordinate, std::vector<float>& coordinates, std::vector<Leaves>& leaves, std::ofstream& os, int& counter){

	int N = (tree.nodes.size() + 1)/2;
	if(node.child_left != NULL){

		Node child_left, child_right;
		child_left = *node.child_left;
		child_right = *node.child_right;
		if(0){
			if(leaves[(*node.child_left).label].num_leaves > leaves[(*node.child_right).label].num_leaves){
				child_left = *node.child_left;
				child_right = *node.child_right;
			}else{
				child_left = *node.child_right;
				child_right = *node.child_left;
			}
		}

		int num_grandchildren_left = 0, num_grandchildren_right = 0;
		if(child_left.child_left != NULL){
			num_grandchildren_left  = leaves[(*child_left.child_left).label].num_leaves;
		}
		if(child_left.child_right != NULL){
			num_grandchildren_right = leaves[(*child_left.child_right).label].num_leaves;
		}

		double x_coord_child_left, x_coord_child_right;

		TraverseTree(child_left, tree, x_coord_child_left, coordinates, leaves, os, counter);
		TraverseTree(child_right, tree, x_coord_child_right, coordinates, leaves, os, counter);

		x_coordinate = (x_coord_child_left + x_coord_child_right)/2;

		os << x_coord_child_left << " " << x_coordinate << " " << coordinates[node.label] << " " << coordinates[node.label] << " " << child_left.label << " " << "h" << "\n"; //horizontal line
		if(child_left.label < N){
			os << x_coord_child_left << " " << x_coord_child_left << " " << coordinates[child_left.label] << " " << coordinates[node.label] << " " << child_left.label << " " << "t" << "\n"; //vertical line
		}else{
			os << x_coord_child_left << " " << x_coord_child_left << " " << coordinates[child_left.label] << " " << coordinates[node.label] << " " << child_left.label << " " << "v" << "\n"; //vertical line
		}

		int num_events = child_left.num_events;
		for(int i = 0; i < num_events; i++){
			double coord = coordinates[child_left.label] + child_left.branch_length/(num_events + 1.0) * (i+1.0);
			os << x_coord_child_left << " " << x_coord_child_left << " " << coord  << " " << coord << " " << child_left.label << " " << "m" << "\n"; //mutation
		}

		os << x_coord_child_right << " " << x_coordinate << " " << coordinates[node.label] << " " << coordinates[node.label] << " " << child_right.label << " " << "h" << "\n"; //horizontal line
		if(child_right.label < N){
			os << x_coord_child_right << " " << x_coord_child_right << " " << coordinates[child_right.label] << " " << coordinates[node.label] << " " << child_right.label << " " << "t" << "\n"; //vertical line   
		}else{
			os << x_coord_child_right << " " << x_coord_child_right << " " << coordinates[child_right.label] << " " << coordinates[node.label] << " " << child_right.label << " " << "v" << "\n"; //vertical line
		}

		num_events = child_right.num_events;
		for(int i = 0; i < num_events; i++){
			double coord = coordinates[child_right.label] + child_right.branch_length/(num_events + 1.0) * (i+1.0);
			os << x_coord_child_right << " " << x_coord_child_right << " " << coord  << " " << coord << " " << child_right.label << " " << "m" << "\n"; //mutation
		}


	}else{
		counter++;
		x_coordinate = counter;
	}


}

void ExtractPlotCoordinates(Tree& tree, std::ofstream& os){

	//calculate number of leaves below each coalescent event
	//calculate y-coordinates of each coalescent event
	//
	//start at root, get number of leaves of children and number of leaves of grandchildren
	//calculate x coordinate from this.
	//add a segment to output.
	//to each segment, output branch id

	std::vector<float> coordinates;
	tree.GetCoordinates(coordinates);
	std::vector<Leaves> leaves;
	tree.FindAllLeaves(leaves);

	int root = tree.nodes.size()-1;

	os << "x_begin x_end y_begin y_end branchID seg_type\n";
	int counter = 0;
	double x_coordinate = 0.0;
	TraverseTree(tree.nodes[root], tree, x_coordinate, coordinates, leaves, os, counter);

	os << x_coordinate << " " << x_coordinate << " " << coordinates[root] << " " << coordinates[root] << " " << root << " " << "v" << "\n"; //vertical line

}

void
TreeView(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("snp_of_interest") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, snp_of_interest, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Outputs tree of interest as .newick file." << std::endl;
		exit(0);
	}  

	int snp_of_interest = options["snp_of_interest"].as<int>();
	int index_of_snp_of_interest;

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Get tree at BP " << snp_of_interest << "..." << std::endl;

	//////////////////////////////////
	//Parse Data
	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N = ancmut.NumTips();
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0;
	int i; 
	std::string line, line2, read;

	//////////////////////////////////////////// Read Tree ///////////////////////////////////

	Mutations mut;
	mut.Read(options["mut"].as<std::string>());
	index_of_snp_of_interest = 0;
	for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
		if((*it_mut).pos >= snp_of_interest) break;
		index_of_snp_of_interest++;
	}
	if(index_of_snp_of_interest == mut.info.size()) index_of_snp_of_interest--;
	int tree_index_of_interest = mut.info[index_of_snp_of_interest].tree;

	//getline(is_anc, line2);
	int count_trees = 0;
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	while(num_bases_tree_persists >= 0.0){

		if(count_trees == tree_index_of_interest){
			//output plot coordinates
			std::ofstream os(options["output"].as<std::string>() + ".plotcoords");
			ExtractPlotCoordinates(mtr.tree, os);
			os.close();

			break; 
		}
		count_trees++;
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

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

}

////////////////////////////////////

void
TraverseTreeSample(Node& n, std::vector<double>& sample_ages, std::vector<std::vector<double>>& coords, std::vector<std::vector<double>>& ages){

	if(n.child_left != NULL){

		Node child_left = (*n.child_left), child_right = (*n.child_right);
		TraverseTreeSample(child_left, sample_ages, coords, ages);
		TraverseTreeSample(child_right, sample_ages, coords, ages);

		for(int i = 0; i < ages[n.label].size(); i++){
			coords[n.label][i] = coords[child_left.label][i] + ages[child_left.label][i];
		}

	}else{

		if(sample_ages.size() > 0){
			for(int i = 0; i < ages[n.label].size(); i++){
				coords[n.label][i] = sample_ages[n.label];
			}
		}else{
			for(int i = 0; i < ages[n.label].size(); i++){
				coords[n.label][i] = 0.0;
			}
		}

	}

}

/*
void
TreeViewSample(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("snp_of_interest") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, snp_of_interest, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Outputs tree of interest as .newick file." << std::endl;
		exit(0);
	}  

	int snp_of_interest = options["snp_of_interest"].as<int>();
	int index_of_snp_of_interest;

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Get tree at BP " << snp_of_interest << "..." << std::endl;

	//////////////////////////////////
	//Parse Data

	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N_ref = ancmut.NumTips();
	MarginalTree mtr_ref;
	Muts::iterator it_mut_ref;
	float num_bases_tree_persists = 0.0;


	int N, num_trees, num_samples;
	MarginalTree mtr;
	Muts::iterator it_mut;

	std::string filename_anc = options["input"].as<std::string>() + ".anc";
	igzstream is(filename_anc);
	if(is.fail()) is.open(filename_anc + ".gz");
	if(is.fail()){
		std::cerr << "Failed to open file " << filename_anc << "(.gz)" << std::endl;
		exit(1);
	}

	std::istringstream is_header;
	std::string line, tmp, line2, read;
	//read num_haplotypes
	getline(is, line);
	is_header.str(line);
	is_header >> tmp;
	is_header >> N;

	//read sample ages
	std::vector<double> sample_ages(N);
	std::vector<double>::iterator it_sample_ages = sample_ages.begin();
	int i = 0;
	while(is_header >> *it_sample_ages){
		it_sample_ages++;
		i++;
		if(it_sample_ages == sample_ages.end()) break;
	}
	if(i != N) sample_ages.clear();

	//read num trees
	getline(is, line);
	is_header.str(line);
	is_header.clear();
	is_header >> tmp;
	is_header >> num_trees;

	//read num trees
	getline(is, line);
	is_header.str(line);
	is_header.clear();
	is_header >> tmp;
	is_header >> num_samples;

	assert(N == N_ref);
	assert(num_trees == ancmut.NumTrees());

	//////////////////////////////////////////// Read Tree ///////////////////////////////////

	Mutations mut;
	mut.Read(options["mut"].as<std::string>());
	index_of_snp_of_interest = 0;
	for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
		if((*it_mut).pos >= snp_of_interest) break;
		index_of_snp_of_interest++;
	}
	int tree_index_of_interest = mut.info[index_of_snp_of_interest].tree;

	std::vector<float> coords;

	//getline(is_anc, line2);
	int count_trees = 0;
	num_bases_tree_persists = ancmut.NextTree(mtr_ref, it_mut_ref);
	while(getline(is, line)){

		if(count_trees == tree_index_of_interest){
			//output plot coordinates

			//output format
			//branchID age

			std::vector<std::vector<double>> ages(2*N-1);
			std::vector<std::vector<double>> coords(2*N-1);
			std::vector<std::vector<double>>::iterator it_ages = ages.begin();
			std::vector<double>::iterator it2_ages;
			for(; it_ages != ages.end(); it_ages++){
				(*it_ages).resize(num_samples);
			}
			for(it_ages = coords.begin(); it_ages != coords.end(); it_ages++){
				(*it_ages).resize(num_samples);
			}

			i = 0;
			for(it_ages = ages.begin(); it_ages != std::prev(ages.end(),1); it_ages++){
				while(line[i] != '(' && i != line.size()-1){
					i++;
				}
				if(i == line.size()-1) break;
				i++;
				tmp.clear();
				it2_ages = (*it_ages).begin();
				for(int k = 0; k < num_samples; k++){
					while(line[i] != ' ' && line[i] != '\t'){
						tmp += line[i];
						assert(line[i] != '(');
						assert(line[i] != ')');
						i++;
					}
					//std::cerr << tmp << std::endl;
					*it2_ages = stof(tmp);
					tmp.clear();
					it2_ages++;
					i++;
				}
				if(i == line.size()-1) break;
			}

			TraverseTreeSample(mtr_ref.tree.nodes[2*N-2], sample_ages, coords, ages);

			std::ofstream os(options["output"].as<std::string>() + ".plotcoords");
			int i = 0;
			os << "branchID age\n";
			for(it_ages = coords.begin(); it_ages != coords.end(); it_ages++){
				for(it2_ages = (*it_ages).begin(); it2_ages != (*it_ages).end(); it2_ages++){
					os << i << " " << *it2_ages << "\n"; 
				}
				i++;
			}
			os.close();

			break; 
		}

		num_bases_tree_persists = ancmut.NextTree(mtr_ref, it_mut_ref);
		count_trees++;

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

}
*/


void
TreeViewSample(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("snp_of_interest") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, snp_of_interest, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Outputs tree of interest as .newick file." << std::endl;
		exit(0);
	}  

	int snp_of_interest = options["snp_of_interest"].as<int>();
	int index_of_snp_of_interest;

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Get tree at BP " << snp_of_interest << "..." << std::endl;

	//////////////////////////////////
	//Parse Data

	int N, num_trees, num_samples;
	MarginalTree mtr;
	Muts::iterator it_mut;

	std::string filename_anc = options["anc"].as<std::string>();
	igzstream is(filename_anc);
	if(is.fail()) is.open(filename_anc + ".gz");
	if(is.fail()){
		std::cerr << "Failed to open file " << filename_anc << "(.gz)" << std::endl;
		exit(1);
	}

	std::istringstream is_header;
	std::string line, tmp, line2, read;
	//read num_haplotypes
	getline(is, line);
	is_header.str(line);
	is_header >> tmp;
	is_header >> N;

	//read sample ages
	std::vector<double> sample_ages(N);
	std::vector<double>::iterator it_sample_ages = sample_ages.begin();
	int i = 0;
	while(is_header >> *it_sample_ages){
		it_sample_ages++;
		i++;
		if(it_sample_ages == sample_ages.end()) break;
	}
	if(i != N) sample_ages.clear();

	//read num trees
	getline(is, line);
	is_header.str(line);
	is_header.clear();
	is_header >> tmp;
	is_header >> num_trees;

	//read num trees
	getline(is, line);
	is_header.str(line);
	is_header.clear();
	is_header >> tmp;
  if(tmp != "NUM_SAMPLES_PER_TREE"){
    std::cerr << "Error: need anc/mut with at least two sampled branch lengths." << std::endl;
    exit(1);
  }
	is_header >> num_samples;

	//////////////////////////////////////////// Read Tree ///////////////////////////////////

	Mutations mut;
	mut.Read(options["mut"].as<std::string>());
	if(mut.info.size() == 0){
    std::cerr << "Error: anc/mut needs to span at least one mutation" << std::endl;
		exit(1); 
	}
	index_of_snp_of_interest = 0;
	for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
		if((*it_mut).pos >= snp_of_interest) break;
		index_of_snp_of_interest++;
	}
	if(index_of_snp_of_interest == mut.info.size()) index_of_snp_of_interest--;
	int tree_index_of_interest = mut.info[index_of_snp_of_interest].tree;

	std::vector<float> coords;

	//getline(is_anc, line2);
	int count_trees = 0;
	while(getline(is, line)){

		if(count_trees == tree_index_of_interest){
			//output plot coordinates

			//output format
			//branchID age

			std::vector<std::vector<double>> ages(2*N-1);
			std::vector<std::vector<double>> coords(2*N-1);
			std::vector<std::vector<double>>::iterator it_ages = ages.begin();
			std::vector<double>::iterator it2_ages;
			for(; it_ages != ages.end(); it_ages++){
				(*it_ages).resize(num_samples);
			}
			for(it_ages = coords.begin(); it_ages != coords.end(); it_ages++){
				(*it_ages).resize(num_samples);
			}

			mtr.tree.sample_ages = &sample_ages;
			mtr.tree.nodes.resize(2*N-1);
			i = 0;
			tmp.clear();
			while(line[i] != ':' && i != line.size()-1){
				tmp += line[i];
				i++;
			}
			mtr.pos = stoi(tmp);
			i += 2;
			int node = 0;
			for(it_ages = ages.begin(); it_ages != std::prev(ages.end(),1); it_ages++){
	
        tmp.clear();
				while(line[i] != ':' && i != line.size()-1){
					tmp += line[i];
					i++;
				}
				int parent = stoi(tmp);
				if(parent != -1){
					mtr.tree.nodes[node].parent = &mtr.tree.nodes[parent];
					mtr.tree.nodes[node].label  = node;
					if(mtr.tree.nodes[parent].child_left == NULL){
						mtr.tree.nodes[parent].child_left = &mtr.tree.nodes[node];
					}else{
						mtr.tree.nodes[parent].child_right = &mtr.tree.nodes[node];
					}
				}else{
          mtr.tree.nodes[node].parent = NULL;
					mtr.tree.nodes[node].label  = node;
				}

				i += 2;
				if(i >= line.size()-1) break;
				tmp.clear();
				it2_ages = (*it_ages).begin();
				double mean_bl = 0.0;
				for(int k = 0; k < num_samples; k++){
					while(line[i] != ' ' && line[i] != '\t'){
						tmp += line[i];
						assert(line[i] != '(');
						assert(line[i] != ')');
						i++;
					}
					*it2_ages = stof(tmp);
					mean_bl += *it2_ages;
					tmp.clear();
					it2_ages++;
					i++;
				}
        mtr.tree.nodes[node].branch_length = mean_bl/num_samples;
				tmp.clear();
				while(line[i] != ' ' && line[i] != '\t'){
          tmp += line[i];
					i++;
					if(i == line.size()) break;
				}
				mtr.tree.nodes[node].num_events = stof(tmp);
				tmp.clear();
				i++;
				while(line[i] != ' ' && line[i] != '\t'){
					tmp += line[i];
					i++;
					if(i == line.size()) break;
				}
				mtr.tree.nodes[node].SNP_begin = stof(tmp);
				tmp.clear();
				i++;
				while(line[i] != ')'){
					tmp += line[i];
          i++;
					if(i == line.size()) break;
				}
				mtr.tree.nodes[node].SNP_end = stof(tmp);
				node++;
				i++;
				i++;
				if(i >= line.size()) break;
			}

			mtr.tree.nodes[2*N-2].parent = NULL;
			mtr.tree.nodes[2*N-2].label  = node;
			
		  AncesTree anc;
			anc.sample_ages = sample_ages;
			anc.seq.push_back(mtr);
			anc.Dump(options["output"].as<std::string>() + ".anc");

			Mutations mut_avg;
			mut_avg.info.resize(1);
			mut_avg.info[0] = mut.info[index_of_snp_of_interest];
			mut_avg.info[0].tree = 0;
			mut_avg.Dump(options["output"].as<std::string>() + ".mut");

			TraverseTreeSample(mtr.tree.nodes[2*N-2], sample_ages, coords, ages);

			std::ofstream os(options["output"].as<std::string>() + ".plotcoords");
			int i = 0;
			os << "branchID age\n";
			for(it_ages = coords.begin(); it_ages != coords.end(); it_ages++){
				for(it2_ages = (*it_ages).begin(); it2_ages != (*it_ages).end(); it2_ages++){
					os << i << " " << *it2_ages << "\n"; 
				}
				i++;
			}
			os.close();

			break; 
		}

		count_trees++;
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

}


//////////////////////////////////////////////

void
MutationsOnBranches(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("snp_of_interest") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, haps, sample, snp_of_interest, output. Optional: dist, mask." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Outputs tree of interest as .newick file." << std::endl;
		exit(0);
	}  

	int snp_of_interest = options["snp_of_interest"].as<int>();
	int index_of_snp_of_interest;

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Get list of mutations that map to tree at BP " << snp_of_interest << "..." << std::endl;

	//read in ancestor
	fasta mask;
	bool is_mask = false;
	if(options.count("mask")){
		mask.Read(options["mask"].as<std::string>());
		is_mask = true;
	}

	//////////////////////////////////
	//Parse Data
	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N = ancmut.NumTips();
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0;

	int i; 
	std::string line, line2, read;

	//////////////////////////////////////////// Read Tree ///////////////////////////////////

	Mutations mut;
	mut.Read(options["mut"].as<std::string>());
	index_of_snp_of_interest = 0;
	for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
		if((*it_mut).pos >= snp_of_interest) break;
		index_of_snp_of_interest++;
	}
	if(index_of_snp_of_interest == mut.info.size()) index_of_snp_of_interest--;
	int tree_index_of_interest = mut.info[index_of_snp_of_interest].tree;

	// read pos
	int L = 0;
	igzstream is_L;
	std::vector<int> pos;
	if(options.count("dist")){
		is_L.open(options["dist"].as<std::string>());
		if(is_L.fail()) is_L.open(options["dist"].as<std::string>() + ".gz");
		if(is_L.fail()){
			std::cerr << "Error while opening file " << options["dist"].as<std::string>() << std::endl;
			exit(1);
		} 
	}else{
		is_L.open(options["mut"].as<std::string>());
		if(is_L.fail()) is_L.open(options["mut"].as<std::string>() + ".gz");
		if(is_L.fail()){
			std::cerr << "Error while opening file " << options["mut"].as<std::string>() << std::endl;
			exit(1);
		} 
	}

	while(std::getline(is_L, line)){
		++L;
	}
	L--;
	is_L.close();

	pos.resize(L);
	if(options.count("dist")){
		igzstream is_dist(options["dist"].as<std::string>());
		if(is_dist.fail()) is_dist.open(options["dist"].as<std::string>() + ".gz");
		if(is_dist.fail()){
			std::cerr << "Error while opening file " << options["dist"].as<std::string>() << std::endl;
			exit(1);
		}
		getline(is_dist, line); 
		int dtmp, snp = 0;
		while(std::getline(is_dist, line)){
			sscanf(line.c_str(), "%d %d", &pos[snp], &dtmp);
			snp++;
		}
		is_dist.close();
	}else{
		int count = 0;
		std::vector<int>::iterator it_pos = pos.begin();
		for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
			*it_pos = (*it_mut).pos;
			it_pos++;
			count++;
		}
	}

	/////
	int count_trees = 0;
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	while(num_bases_tree_persists >= 0.0){

		if(count_trees == tree_index_of_interest){

			//output list of SNPs mapping to each branch

			//get min_snp and max_snp, read in all sites between these two.
			//then try to map each of then and record which branch they mapped to (use assertions)

			std::vector<std::vector<int>> mut_on_branches(mtr.tree.nodes.size());

			int min_snp = mtr.tree.nodes[0].SNP_begin, max_snp = mtr.tree.nodes[0].SNP_end;
			for(int i = 0; i < mtr.tree.nodes.size(); i++){
				if(min_snp > mtr.tree.nodes[i].SNP_begin) min_snp = mtr.tree.nodes[i].SNP_begin;
				if(max_snp < mtr.tree.nodes[i].SNP_end) max_snp = mtr.tree.nodes[i].SNP_end;
			}
			int min_bp = pos[min_snp], max_bp = pos[max_snp];

			haps m_hap(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str());
			Data data(m_hap.GetN(), m_hap.GetL());
			if(data.N != ancmut.NumTips()){
				std::cerr << "Haps file and anc/mut have different number of samples" << std::endl;
				exit(1);
			}
			Leaves sequences_carrying_mutation;
			sequences_carrying_mutation.member.resize(data.N);
			AncesTreeBuilder ancbuilder(data);

			std::vector<char> sequence(data.N);
			int bp, num_carriers;
			int snp = 0;
			Tree tr = mtr.tree;
			do{
				m_hap.ReadSNP(sequence, bp);
				snp++;
			}while(bp < min_bp && snp < data.L);
			/*
				 if(bp != min_bp){
				 std::cerr << bp << " " << min_bp << std::endl;
				 std::cerr << "SNP is not in haps file" << std::endl;
				 exit(1);
				 } 
				 */

			while(bp <= max_bp && snp < data.L){

				//map sequence to tree
				sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
				for(int i = 0; i < data.N; i++){
					if(sequence[i] == '1'){
						sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
						sequences_carrying_mutation.num_leaves++;
					}else{
						sequences_carrying_mutation.member[i] = 0;
					}
				}

				if(sequences_carrying_mutation.num_leaves > 0 && sequences_carrying_mutation.num_leaves < data.N){
					//if(bp == 236591955) std::cerr << "SNP exists" << std::endl; 
					if(ancbuilder.IsSNPMapping(tr, sequences_carrying_mutation, snp) == 1){

						if(is_mask){
							if(std::toupper(mask.seq[bp-1]) == 'P'){           
								int branch = ancbuilder.mutations.info[snp].branch[0];
								if(pos[mtr.tree.nodes[branch].SNP_begin] <= bp && pos[mtr.tree.nodes[branch].SNP_end] >= bp && mtr.tree.nodes[branch].num_events > 0){
									mut_on_branches[branch].push_back(bp);
								}
							}
						}else{ 
							int branch = ancbuilder.mutations.info[snp].branch[0];
							if(pos[mtr.tree.nodes[branch].SNP_begin] <= bp && pos[mtr.tree.nodes[branch].SNP_end] >= bp){
								mut_on_branches[branch].push_back(bp);
							}
						}
					}
				}

				m_hap.ReadSNP(sequence, bp);  

				snp++;
			}
			m_hap.CloseFile();

			std::ofstream os(options["output"].as<std::string>() + ".plotcoords.mut");

			os << "pos branchID\n";

			i = 0;
			for(std::vector<std::vector<int>>::iterator it_mut_on_branches = mut_on_branches.begin(); it_mut_on_branches != mut_on_branches.end(); it_mut_on_branches++){
				for(std::vector<int>::iterator it_branch = (*it_mut_on_branches).begin(); it_branch != (*it_mut_on_branches).end(); it_branch++){
					os << *it_branch << " " << i << "\n";
				}
				i++;
			} 

			os.close();

			break; 

		}
		count_trees++;
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

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

}

void TraverseTree(Node& node, Tree& tree, std::ofstream& os){

	os << node.label << "\n";
	if(node.child_left != NULL){

		TraverseTree((*node.child_left), tree, os);
		TraverseTree((*node.child_right), tree, os);

	}

}

void
BranchesBelowMutation(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("snp_of_interest") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, haps, sample, snp_of_interest, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Outputs tree of interest as .newick file." << std::endl;
		exit(0);
	}  

	int snp_of_interest    = options["snp_of_interest"].as<int>();
	int index_of_snp_of_interest;
	int branch_of_interest;

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Get list of branches below mutation at BP " << snp_of_interest << "..." << std::endl;


	//////////////////////////////////
	//Parse Data
	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N = ancmut.NumTips();
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0;

	int i; 
	std::string line, line2, read;

	//////////////////////////////////////////// Read Tree ///////////////////////////////////

	Mutations mut;
	mut.Read(options["mut"].as<std::string>());
	index_of_snp_of_interest = 0;
	for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
		if((*it_mut).pos >= snp_of_interest) break;
		index_of_snp_of_interest++;
	}
	if(index_of_snp_of_interest == mut.info.size()) index_of_snp_of_interest--;
	int tree_index_of_interest = mut.info[index_of_snp_of_interest].tree;

	if(mut.info[index_of_snp_of_interest].branch.size() == 1){
		branch_of_interest = mut.info[index_of_snp_of_interest].branch[0];
	}else{
		std::cerr << "SNP is not mapping to a unique branch." << std::endl;
		exit(1);
	}

	int count_trees = 0;
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	while(num_bases_tree_persists >= 0.0){

		if(count_trees == tree_index_of_interest){

			//output list of SNPs mapping to each branch
			std::ofstream os(options["output"].as<std::string>() + ".plotcoords.mut");
			os << "branchID\n";
			TraverseTree(mtr.tree.nodes[branch_of_interest], mtr.tree, os);
			os.close();

			break; 

		}

		count_trees++; 
		ancmut.NextTree(mtr, it_mut);

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

}


