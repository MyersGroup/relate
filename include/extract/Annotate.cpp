#include "filesystem.hpp"
#include "cxxopts.hpp"
#include <string>

void
GenerateSNPAnnotationsUsingTree(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if( !options.count("anc") || !options.count("mut") || !options.count("poplabels") || !options.count("output") ){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, poplabels, output. Optional: ancestor." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Generate additional annotation for SNPs." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Generating additional annotation for SNPs... " << std::endl;

	bool is_ancestor = false;
	if(options.count("ancestor")) is_ancestor = true;

	fasta ancestor;
	if(is_ancestor) ancestor.Read(options["ancestor"].as<std::string>());

	Sample sample;
	sample.Read(options["poplabels"].as<std::string>());

	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	MarginalTree mtr;
	Muts::iterator it_mut; //iterator for mut file
	float num_bases_SNP_persists = 0.0;

	int N = ancmut.NumTips();
	int num_trees = ancmut.NumTrees();
	int L = ancmut.NumSnps();
	Data data(N, L);
	assert(N == sample.group_of_haplotype.size());

	Mutations mut;
	bool is_mut = false;
	if(options.count("mut")){
		is_mut = true;
		mut.Read(options["mut"].as<std::string>());
	}

	char nucl;
	std::vector<int> carriers_by_pop(sample.groups.size());
	int percentage = 0;
	std::cerr << "[" << percentage << "%]\r";
	std::cerr.flush();

	int snp = 0;
	int bp = 0;
	int tree = 0;
	//get first SNP (Only necessary when using NextSNP, for NextTree I don't need to call this)
	num_bases_SNP_persists = ancmut.FirstSNP(mtr, it_mut);
	std::vector<Leaves> desc(2*N-1);
	mtr.tree.FindAllLeaves(desc);
	tree = (*it_mut).tree;

	//iterate through whole file
	while(num_bases_SNP_persists >= 0.0){

		if(tree < (*it_mut).tree){
			mtr.tree.FindAllLeaves(desc);
			tree = (*it_mut).tree;
		}
		bp = (*it_mut).pos;

		//std::cerr << (int)(data.L/100) << " " << snp % (int)(data.L/100) << std::endl;
		if((snp % (int)(data.L/100)) == 0){
			std::cerr << "[" << percentage << "%]\r";
			percentage++;
			std::cerr.flush();
		}

		//get nucleotide before and after in sequence
		if(is_ancestor){
			if(bp > 1){
				nucl = std::toupper(ancestor.seq[bp-2]);
				if(nucl == 'A' || nucl == 'C' || nucl == 'G' || nucl == 'T'){
					if(is_mut) mut.info[snp].upstream_base = nucl;
				}	
			}
			if(bp < ancestor.seq.size()){
				nucl = std::toupper(ancestor.seq[bp]);
				if(nucl == 'A' || nucl == 'C' || nucl == 'G' || nucl == 'T'){
					if(is_mut) mut.info[snp].downstream_base = nucl;
				}
			}
		}

		//calculate number of carriers in each population
		std::fill(carriers_by_pop.begin(), carriers_by_pop.end(), 0);
		if((*it_mut).branch.size() == 1){
			//std::cerr << (*it_mut).pos << " " << *(*it_mut).branch.begin() << " " << desc[*(*it_mut).branch.begin()].num_leaves << std::endl;
			for(std::vector<int>::iterator it_mem = desc[*(*it_mut).branch.begin()].member.begin(); it_mem != desc[*(*it_mut).branch.begin()].member.end(); it_mem++){
				carriers_by_pop[sample.group_of_haplotype[*it_mem]]++;
			}
		}

		if(is_mut){
			mut.info[snp].freq.resize(carriers_by_pop.size());
			int pop = 0;
			for(std::vector<int>::iterator it_carriers = carriers_by_pop.begin(); it_carriers != carriers_by_pop.end(); it_carriers++){
				mut.info[snp].freq[pop] = *it_carriers;
				pop++;
			} 
		}

		//it_mut points to a SNP, mtr stores the marginal tree corresponding to that SNP.
		num_bases_SNP_persists = ancmut.NextSNP(mtr, it_mut);
		snp++;
	}

	mut.header  = "snp;pos_of_snp;dist;rs-id;tree_index;branch_indices;is_not_mapping;is_flipped;age_begin;age_end;ancestral_allele/alternative_allele;"; 
	mut.header += "upstream_allele;downstream_allele;";
	for(std::vector<std::string>::iterator it_groups = sample.groups.begin(); it_groups != sample.groups.end(); it_groups++){
		mut.header += *it_groups + ";";
	}
	mut.Dump(options["output"].as<std::string>() + ".mut");
	std::cerr << "Output written to " << options["output"].as<std::string>() + ".mut" << std::endl;

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
PropagateMutations(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Example code for converting from tree sequence file format." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Propagate mutations " << options["anc"].as<std::string>() << " and " << options["mut"].as<std::string>() << "..." << std::endl;

	//Read anc/mut 
	AncesTree anc;
	anc.Read(options["anc"].as<std::string>());
	Mutations mut;
	mut.Read(options["mut"].as<std::string>());


	CorrTrees::iterator it_seq_prev;
	CorrTrees::iterator it_seq; 

	int N_total = (*anc.seq.begin()).tree.nodes.size();
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

	//for each tree, I want a
	std::vector<std::vector<std::vector<int>>> tree_mutations;
	std::vector<std::vector<std::vector<int>>>::iterator it_muts, it_muts_prev;

	tree_mutations.resize(anc.L);
	int snp = 0;
	int tree_index = mut.info[snp].tree;

	for(std::vector<std::vector<std::vector<int>>>::iterator it_muts = tree_mutations.begin(); it_muts != tree_mutations.end(); it_muts++){
		(*it_muts).resize(2*N-1);
	}

	for(int snp = 0; snp < mut.info.size(); snp++){
		tree_mutations[mut.info[snp].tree][*mut.info[snp].branch.begin()].push_back(snp);
	}

	it_seq_prev = anc.seq.begin();
	it_seq      = std::next(it_seq_prev,1); 

	std::vector<std::vector<int>> equivalent_branches;
	std::vector<std::vector<int>>::iterator it_equivalent_branches;
	std::vector<std::vector<int>>::reverse_iterator rit_equivalent_branches;

	for(; it_seq != anc.seq.end();){
		equivalent_branches.emplace_back();
		it_equivalent_branches = std::prev(equivalent_branches.end(),1);
		anc.BranchAssociation((*it_seq_prev).tree, (*it_seq).tree, *it_equivalent_branches, potential_branches, N, N_total, threshold_brancheq); //O(N^2) 
		it_seq++;
		it_seq_prev++;
	}  

	///////////////////////////////////////////
	//Now carry over information on branches, starting from the first tree.

	it_equivalent_branches = equivalent_branches.begin();
	std::vector<Node>::iterator it_nodes;

	it_seq_prev = anc.seq.begin();
	it_seq      = std::next(it_seq_prev,1); 
	it_muts_prev = tree_mutations.begin();
	it_muts      = std::next(it_muts_prev,1); 

	for(; it_seq != anc.seq.end();){
		it_nodes = (*it_seq).tree.nodes.begin();
		for(std::vector<int>::iterator it = (*it_equivalent_branches).begin(); it != (*it_equivalent_branches).end(); it++){
			if(*it != -1){
				(*it_nodes).num_events += (*it_seq_prev).tree.nodes[*it].num_events;
				(*it_nodes).SNP_begin   = (*it_seq_prev).tree.nodes[*it].SNP_begin;
				for(std::vector<int>::iterator it_snps = (*it_muts_prev)[*it].begin(); it_snps != (*it_muts_prev)[*it].end(); it_snps++){
					(*it_muts)[(*it_nodes).label].push_back(*it_snps);
				}
			}
			it_nodes++;
		}
		it_equivalent_branches++;

		it_seq++;
		it_seq_prev++;
		it_muts++;
		it_muts_prev++;
	}
	assert(it_equivalent_branches == equivalent_branches.end());

	///////////////////////////////////////////
	//Now go from the last tree to the first

	rit_equivalent_branches = equivalent_branches.rbegin();

	CorrTrees::reverse_iterator rit_seq_next;
	CorrTrees::reverse_iterator rit_seq; 
	std::vector<std::vector<std::vector<int>>>::reverse_iterator rit_muts, rit_muts_next;



	rit_seq_next = anc.seq.rbegin();
	rit_seq      = std::next(rit_seq_next,1); 
	rit_muts_next = tree_mutations.rbegin();
	rit_muts     = std::next(rit_muts_next,1);

	for(; rit_seq != anc.seq.rend();){
		it_nodes = (*rit_seq_next).tree.nodes.begin();
		for(std::vector<int>::iterator it = (*rit_equivalent_branches).begin(); it != (*rit_equivalent_branches).end(); it++){
			if(*it != -1){
				(*rit_seq).tree.nodes[*it].num_events = (*it_nodes).num_events;
				(*rit_seq).tree.nodes[*it].SNP_end    = (*it_nodes).SNP_end;
				(*rit_muts)[*it] = (*rit_muts_next)[(*it_nodes).label];	
				std::sort((*rit_muts)[*it].begin(), (*rit_muts)[*it].end());
			}
			it_nodes++;
		}
		rit_equivalent_branches++;

		rit_seq++;
		rit_seq_next++;
		rit_muts++;
		rit_muts_next++;
	}

	assert(rit_equivalent_branches == equivalent_branches.rend());


	std::ofstream os(options["output"].as<std::string>() + ".allmuts");
	os << "treeID branchID SNPID\n";
	int tree = 0;
	for(std::vector<std::vector<std::vector<int>>>::iterator it_muts = tree_mutations.begin(); it_muts != tree_mutations.end(); it_muts++){
		for(int b = 0; b < 2*N-1; b++){
			if((*it_muts)[b].size() > 0){
				for(std::vector<int>::iterator it_snps = (*it_muts)[b].begin(); it_snps != (*it_muts)[b].end(); it_snps++){
					os << tree << " " << b << " " << *it_snps << "\n";
				}
			}
		}
		tree++;
	}
	os.close();

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
PrintMutonBranches(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("anc") || !options.count("mut") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: anc, mut, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Example code for converting from tree sequence file format." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Output num mutations mapping to each branch " << options["anc"].as<std::string>() << " and " << options["mut"].as<std::string>() << "..." << std::endl;

	AncMutIterators ancmut(options["anc"].as<std::string>(), options["mut"].as<std::string>());
	int N = ancmut.NumTips();
	int num_trees = ancmut.NumTrees();
	int L = ancmut.NumSnps();
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0;

	Data data(N,L);
	data.dist.resize(L);
	std::string line;
	if(options.count("dist")){
		igzstream is_dist(options["dist"].as<std::string>());
		if(is_dist.fail()){
			std::cerr << "Error while opening " << options["dist"].as<std::string>() << std::endl;
			exit(1);
		}
		getline(is_dist, line); 
		int dtmp, snp = 0;
		while(std::getline(is_dist, line)){
			sscanf(line.c_str(), "%d %d", &dtmp, &data.dist[snp]);
			snp++;
		}
		is_dist.close();
	}else{
		std::vector<int>::iterator it_pos = data.dist.begin();
		for(std::vector<SNPInfo>::iterator it_mut = ancmut.mut.info.begin(); it_mut != ancmut.mut.info.end(); it_mut++){
			*it_pos = (*it_mut).dist;
			it_pos++;
		}
	}

	std::ofstream os(options["output"].as<std::string>() + ".allmuts");
  os << "treeID branchID pos_start pos_end dist branch_length num_muts\n";

	int count_trees = 0;
	int treeID, pos_start, pos_end, snp_begin, snp_end;
	float dist;
	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	treeID = (*it_mut).tree;
	while(num_bases_tree_persists >= 0.0){

		treeID = (*it_mut).tree;

		while((*it_mut).tree == treeID){
			pos_end = (*it_mut).pos;
			it_mut++;
      if(it_mut == ancmut.mut_end()) break;
		}
    pos_end = ((*it_mut).pos + pos_end)/2.0;

     
			for(std::vector<Node>::iterator it_n = mtr.tree.nodes.begin(); it_n != mtr.tree.nodes.end(); it_n++){

				int snp_begin = (*it_n).SNP_begin;
				int snp_end   = (*it_n).SNP_end;

				assert(snp_end < data.dist.size());
				dist          = 0.0;
				for(int snp = snp_begin; snp < snp_end; snp++){
					dist       += data.dist[snp];
				}

				if(snp_begin > 0){
					snp_begin--;
					pos_start   = (ancmut.mut.info[snp_begin].pos + ancmut.mut.info[snp_begin+1].pos)/2.0;
					dist       += 0.5 * data.dist[snp_begin];
				}else{
					pos_start   = ancmut.mut.info[snp_begin].pos;
				}
				if(snp_end < data.L-1){
					pos_end     = (ancmut.mut.info[snp_end].pos + ancmut.mut.info[snp_end+1].pos)/2.0;
					dist       += 0.5 * data.dist[snp_end];
				}else{
					pos_end     = ancmut.mut.info[snp_end].pos;
				}

				os << treeID << " " << (*it_n).label << " " << pos_start << " " << pos_end << " " << dist << " " << (*it_n).branch_length << " " << (int) (*it_n).num_events << "\n"; 
			}


		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	}
	os.close();








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


