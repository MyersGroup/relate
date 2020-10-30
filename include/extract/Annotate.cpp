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
		num_bases_SNP_persists = ancmut.NextTree(mtr, it_mut);
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




