#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"
#include <ctime>


void
GetTreeOfInterest(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || (!options.count("bp_of_interest") && (!options.count("first_bp") || !options.count("last_bp")) ) || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, output, bp_of_interest or (first_bp and last_bp). Optional: years_per_gen (Default: 28)" << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Outputs tree of interest as .newick file." << std::endl;
		exit(0);
	}  

	int first_bp, last_bp;
	if(options.count("bp_of_interest")){
		first_bp = options["bp_of_interest"].as<int>();
		last_bp  = first_bp;
	}else if(options.count("first_bp") && options.count("last_bp")){
		first_bp = options["first_bp"].as<int>();
		last_bp  = options["last_bp"].as<int>();
	}else{
		std::cerr << "Error: need either --bp_of_interest or both --first_bp and --last_bp" << std::endl;
		exit(1);
	}

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Get tree from " << options["anc"].as<std::string>() << " in region [" << first_bp << "," << last_bp  << "]..." << std::endl;

	double years_per_gen = 28;
	if(options.count("years_per_gen")) years_per_gen = options["years_per_gen"].as<float>();

	//////////////////////////////////
	//Parse Data
	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N = ancmut.NumTips();
	int L = ancmut.NumSnps();
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0;
	Data data(N,L);

	int i; 
	std::string line, read;

	//////////////////////////////////////////// Read Tree ///////////////////////////////////

	Mutations mut;
	mut.Read(options["mut"].as<std::string>());
	int index_of_first_bp = -1;
	for(it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
		index_of_first_bp++;
		if((*it_mut).pos >= first_bp) break;
	}
	if(index_of_first_bp == -1){
		std::cerr << "BP not covered by anc/mut" << std::endl;	
		exit(1);
	}
	int tree_index_start = mut.info[index_of_first_bp].tree;

	int index_of_last_bp = index_of_first_bp;
	if(last_bp > first_bp && it_mut != mut.info.end()){
		if((*it_mut).pos < last_bp){
			for(; it_mut != mut.info.end(); it_mut++){
				index_of_last_bp++;
				if((*it_mut).pos >= last_bp) break;
			}
			if(it_mut == mut.info.end()) index_of_last_bp = mut.info.size() - 1;
		}
	}
	int tree_index_end = mut.info[index_of_last_bp].tree;

	std::string filename = options["output"].as<std::string>() + ".newick";
	std::ofstream os(filename);
	std::ofstream os_pos(options["output"].as<std::string>() + ".pos");

	int count_trees = 0;
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	while(num_bases_tree_persists >= 0.0){

		if(count_trees >= tree_index_start && count_trees <= tree_index_end){
			//tree is in line
			os_pos << mut.info[mtr.pos].pos << "\n";
			mtr.tree.WriteNewick(os, years_per_gen);
		}
		if(count_trees == tree_index_end) break;

		count_trees++;
		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	}
	os.close();
	os_pos.close();

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


void
MapMutation(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("haps") || !options.count("sample") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, haps, sample, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Map mutations to tree." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Mapping mutations to " << options["anc"].as<std::string>() << "..." << std::endl;
  std::cerr << "Mutations mapping at already existing positions will be skipped." << std::endl;

	haps mhaps(options["haps"].as<std::string>().c_str(), options["sample"].as<std::string>().c_str()); 
	Data data(mhaps.GetN(), mhaps.GetL());

	int bp;
	std::vector<char> sequence(data.N);

	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N = ancmut.NumTips();
	int N_total = 2*N-1;
	int root = N_total - 1;
	MarginalTree mtr, mtr_prev;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0;

	Leaves sequences_carrying_mutation;
	sequences_carrying_mutation.member.resize(data.N);
	std::vector<float> coordinates(N_total, 0.0);  

	//read SNPs from haps/sample
	//if mutation at this position exists, throw warning and skip
	//otherwise, map to tree
	//output new mut file (but not anc file, or have this as option)

	AncesTreeBuilder ab(data);

	Mutations mut;
	mut.info.resize(ancmut.NumSnps() + data.L);

	num_bases_tree_persists = ancmut.FirstSNP(mtr, it_mut);
	mtr_prev = mtr;
	mtr_prev.tree.GetCoordinates(coordinates);	
	int snp = 0, num_not_mapping = 0, num_flipped = 0, snp_mut = 0, count_tree = 1;
	while(snp < data.L){

		mhaps.ReadSNP(sequence, bp);

		//fast forward to tree onto which I want to map this SNP, copying muts along the way
		if(num_bases_tree_persists >= 0){
			while(bp > (*it_mut).pos){
				mut.info[snp_mut] = (*it_mut);
				//mut.info[snp_mut].snp_id = snp_mut;
				snp_mut++;
				if(count_tree < (*it_mut).tree){
					mtr_prev = mtr;
					count_tree = (*it_mut).tree;
					mtr_prev.tree.GetCoordinates(coordinates);	
				}
				num_bases_tree_persists = ancmut.NextSNP(mtr, it_mut);
				if(num_bases_tree_persists < 0.0) break;
			}
		}

		//map new SNP to tree
		if(bp != (*it_mut).pos && bp != mut.info[snp_mut].pos){

			//map mutation onto tree
			sequences_carrying_mutation.num_leaves = 0; //this stores the number of nodes with a mutation at this snp.
			for(int i = 0; i < data.N; i++){
				if(sequence[i] == '1'){
					sequences_carrying_mutation.member[i] = 1; //member stores a sequence of 0 and 1, where 1 at position i means that i carries a mutation.
					sequences_carrying_mutation.num_leaves++;
				}else{
					sequences_carrying_mutation.member[i] = 0;
				}
			}

			if(sequences_carrying_mutation.num_leaves == data.N){

				mut.info[snp_mut].tree      = count_tree-1;
				mut.info[snp_mut].branch.resize(1);
				mut.info[snp_mut].branch[0] = root;
				mut.info[snp_mut].age_begin = coordinates[root];
				mut.info[snp_mut].age_end   = std::numeric_limits<float>::infinity();

			}else{
				if(ab.IsSNPMapping(mtr_prev.tree, sequences_carrying_mutation, snp) == 2){
					num_not_mapping++;
				}
				mut.info[snp_mut].tree      = count_tree-1;
				mut.info[snp_mut].branch    = ab.mutations.info[snp].branch;
				if(mut.info[snp_mut].branch.size() == 1){	
					int branch = mut.info[snp_mut].branch[0];
					if(branch < root){
						mut.info[snp_mut].age_begin = coordinates[branch];
						mut.info[snp_mut].age_end   = coordinates[(*mtr_prev.tree.nodes[branch].parent).label];
					}else{
						mut.info[snp_mut].age_begin = coordinates[branch];
						mut.info[snp_mut].age_end   = std::numeric_limits<float>::infinity();
					}
				}else{
					mut.info[snp_mut].age_begin = 0.0;
					mut.info[snp_mut].age_end   = 0.0;
				}
			}

			mut.info[snp_mut].flipped = ab.mutations.info[snp].flipped;
			if(mut.info[snp_mut].flipped) num_flipped++;

			mut.info[snp_mut].rs_id  = mhaps.rsid;
			//mut.info[snp_mut].snp_id = snp_mut;
			mut.info[snp_mut].snp_id = -1;
			mut.info[snp_mut].pos    = bp;
			mut.info[snp_mut].dist = 0.0;

			mut.info[snp_mut].mutation_type = mhaps.ancestral;
			mut.info[snp_mut].mutation_type += '/';
			mut.info[snp_mut].mutation_type += mhaps.alternative;

			snp_mut++;
		}

		snp++;
	}

	while(num_bases_tree_persists >= 0){
		mut.info[snp_mut] = (*it_mut);
		//mut.info[snp_mut].snp_id = snp_mut;
		snp_mut++;
		mtr_prev = mtr;
		num_bases_tree_persists = ancmut.NextSNP(mtr, it_mut);
	}

	mut.info.resize(snp_mut);
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

void
UnlinkTips(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, output" << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Unlink Tips." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Unlink tips of " << options["anc"].as<std::string>() << "..." << std::endl;

	bool use_transitions = true;
	if(options.count("transversion")) use_transitions = false;

	std::string line;

	std::ofstream os(options["output"].as<std::string>() + ".anc");

	igzstream is(options["anc"].as<std::string>());
	if(is.fail()){
		std::cerr << "Error opening .anc file" << std::endl;
		exit(1);
	}
	assert(getline(is,line));
	os << line << "\n";
	assert(getline(is,line));
	os << line << "\n";
	is.close();

	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N = ancmut.NumTips();
	int N_total = 2*N-1;
	int root = N_total - 1;
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0;


	igzstream is_in(options["input"].as<std::string>());
	if(is_in.fail()){
		std::cerr << "Error opening input file" << std::endl;
		exit(1);
	}
	std::vector<int> tips;
	while(getline(is_in,line)){
		int i = std::stoi(line);
		assert(i < 2*ancmut.NumTips()-1);
		tips.push_back(i);
	}
	std::sort(tips.begin(), tips.end());
	is_in.close();



	std::vector<Node>::iterator n_it;
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	int snp = 0;
	while(num_bases_tree_persists >= 0){
  
		int tree_index = (*it_mut).tree;
		float dist = 0.0;
		int SNP_begin = (*it_mut).snp_id;

		for(std::vector<int>::iterator it_tip = tips.begin(); it_tip != tips.end(); it_tip++){
      mtr.tree.nodes[*it_tip].num_events = 0.0;
			mtr.tree.nodes[*it_tip].SNP_begin = SNP_begin;
		}

		while((*it_mut).tree == tree_index){
		  dist += (*it_mut).dist;
			if((*it_mut).branch.size() == 1){
				if((*it_mut).branch[0] < ancmut.NumTips()){

					bool use = true;
					if(!use_transitions){
						if((*it_mut).mutation_type == "C/T" || (*it_mut).mutation_type == "T/C" || (*it_mut).mutation_type == "G/A" || (*it_mut).mutation_type == "A/G") use = false;
					}

					if(use){
						for(std::vector<int>::iterator it_tip = tips.begin(); it_tip != tips.end(); it_tip++){
							if((*it_mut).branch[0] == *it_tip){
								mtr.tree.nodes[(*it_mut).branch[0]].num_events += 1.0;
								break;
							}
						}
					}
				}
			}
			it_mut++;
		}
    int SNP_end   = (*it_mut).snp_id;
		for(std::vector<int>::iterator it_tip = tips.begin(); it_tip != tips.end(); it_tip++){
			mtr.tree.nodes[*it_tip].SNP_end = SNP_end;
		}

		int parent;
		n_it = mtr.tree.nodes.begin();
		os << mtr.pos << ": ";
		for(; n_it != mtr.tree.nodes.end(); n_it++){
			if((*n_it).parent == NULL){
				parent = -1;
			}else{
				parent = (*(*n_it).parent).label;
			}
			os << parent << ":(";
			os << std::fixed << std::setprecision(5) << (*n_it).branch_length << " ";
			os << std::setprecision(2) << (*n_it).num_events << " " << (*n_it).SNP_begin << " " << (*n_it).SNP_end << ") ";
		}
		os << "\n";

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
