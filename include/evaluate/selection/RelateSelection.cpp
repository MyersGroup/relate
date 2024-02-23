#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>

#include "gzstream.hpp"
#include "collapsed_matrix.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

// Function to sort the vector and create an index vector
void sortAndGetIndices(std::vector<float>& vec, std::vector<int>& index){

	// Create an index vector with initial values 0, 1, 2, ...
	std::vector<int> indices(vec.size());
	for (int i = 0; i < vec.size(); ++i) {
		indices[i] = i;
	}

	// Sort the index vector based on the values in the original vector
	std::sort(indices.begin(), indices.end(), [&vec](int a, int b) {
			if(vec[a] == vec[b]){
        return a < b;
			}else{
			  return vec[a] < vec[b];
			}
			});

	// Create the sorted vector using the index vector
	std::vector<float> sortedVec(vec.size());
	for(int i = 0; i < vec.size(); ++i) {
		sortedVec[i] = vec[indices[i]];
	}
	index = indices;
	vec = sortedVec;

}


void 
logFactorial(std::vector<float>& logF, int N){

	logF.resize(N+1); 
	std::vector<float>::iterator it_logF = logF.begin();
	*it_logF = 0;
	it_logF++;
	std::vector<float>::iterator it_logF_prev = logF.begin();
	for(int k = 1; k < N+1; k++){
		*it_logF = *it_logF_prev + std::log(k);
		it_logF++;
		it_logF_prev++;
	}

}

void
CopyCoordinates(int b, std::vector<float>& coordinates_mutation, const std::vector<float>& coordinates_unsrt, Tree& tr, int& DAF){

	//if( coordinates_unsrt[tr.nodes[b].label] != 0.0 ){
	coordinates_mutation[b] = coordinates_unsrt[b];
	if(tr.nodes[b].child_left != NULL){
		CopyCoordinates((*tr.nodes[b].child_left).label, coordinates_mutation, coordinates_unsrt, tr, DAF);
		CopyCoordinates((*tr.nodes[b].child_right).label, coordinates_mutation, coordinates_unsrt, tr, DAF);
	}
	//}
	if(tr.nodes[b].child_left == NULL) DAF++;

}

/*
	 float
	 log_pvalue_old(int k_now, int k_then, int z, std::vector<float>& logF){

	 if(k_then == -1) return 1;

	 assert(k_then >= 2);
	 assert(k_then < k_now);
	 assert(z < k_now);

	 float log_pvalue, tmp_pvalue;

	 if(0){
	 int z_tmp = 2;
	 log_pvalue  = logF[k_now - z_tmp - 1] - logF[k_now - z_tmp - k_then + 1] - logF[k_then - 2];
	 log_pvalue += log(z_tmp - 1);
	 log_pvalue -= logF[k_now - 1] - logF[k_now - k_then - 1] - logF[k_then];

	 for(int i = z_tmp + 1; i <= k_now - k_then + 1; i++){
	 tmp_pvalue  = logF[k_now - i - 1] - logF[k_now - i - k_then + 1] - logF[k_then - 2];
	 tmp_pvalue += log(i-1);
	 tmp_pvalue -= logF[k_now - 1] - logF[k_now - k_then - 1] - logF[k_then];

	 if(tmp_pvalue - log_pvalue <= 0.0){
	 log_pvalue  = log( 1.0 + exp(tmp_pvalue - log_pvalue) ) + log_pvalue;
	 }else{
	 log_pvalue  = log( 1.0 + exp(log_pvalue - tmp_pvalue) ) + tmp_pvalue;
	 }
	 }
	 assert(log_pvalue < 1e-2);
	 }

	 log_pvalue  = logF[k_now - z - 1] - logF[k_now - z - k_then + 1] - logF[k_then - 2];
	 log_pvalue += log(z-1);
	 log_pvalue -= logF[k_now - 1] - logF[k_now - k_then - 1] - logF[k_then];

	 for(int i = z+1; i <= k_now - k_then + 1; i++){
	 tmp_pvalue  = logF[k_now - i - 1] - logF[k_now - i - k_then + 1] - logF[k_then - 2];
	 tmp_pvalue += log(i-1);
	 tmp_pvalue -= logF[k_now - 1] - logF[k_now - k_then - 1] - logF[k_then];

	 if(tmp_pvalue - log_pvalue <= 0.0){
	 log_pvalue  = log( 1.0 + exp(tmp_pvalue - log_pvalue) ) + log_pvalue;
	 }else{
	 log_pvalue  = log( 1.0 + exp(log_pvalue - tmp_pvalue) ) + tmp_pvalue;
	 }
	 }

	 if(std::fabs(log_pvalue) < 1e-2) log_pvalue = 0.0;


	 if(log_pvalue > 0.0) std::cerr << log_pvalue << " " << k_now << " " << k_then << " " << z << std::endl;

	 assert(log_pvalue <= 0.0);

	 return log_pvalue/log(10);

	 }
	 */

float log_10 = std::log(10);

float 
log_pvalue(int k, float fk, int N, float fN, std::vector<float>& logF){

	float logp = 0.0;
	float px = 0.0;

	if(fk < 2) return 1.0;
	if(k == -1) return 1.0;

	assert(fN < N);
	assert(fk < k);
	assert(fN > 0);

	//Calculate P(fN | N, k, fk)
	px  = logF[N-fN-1] - logF[k-fk-1] - logF[N-k+fk-fN];
	px += logF[fN-1]   - logF[fk-1]   - logF[fN-fk];
	px -= logF[N-1]    - logF[k-1]    - logF[N-k];

	//std::cerr << px << std::endl;

	//Sum over fN <= f <= N-k+fk 
	logp = px;
	float x = fN - fk;
	int y = N - k;
	int c = N - 1;
	int var;
	while(x < N-k){
		var  = fk + x;
		px  += std::log( (y-x)/(x+1.0) * var/((float) (c - var)) );
		//std::cerr << px << " " << logp << " " << c << " " << var << " " << x << std::endl;
		//logp = std::log( 1.0 + exp( logp - px ) ) + px;
		logp = std::log( 1.0 + exp( px - logp ) ) + logp;
		x++;
	}

	//std::cerr << logp << std::endl;
	if(logp > 0.0) logp = 0.0;

	logp /= log_10;

	return(logp);

}

struct Qual{

	std::string id;
	int pos;
	float num_snps_on_tree;
	float frac_branches_with_mut;
	float frac_not_mapping;

};

void 
Selection(cxxopts::Options& options){

	//////////////////////////////////
	//Program options
	bool help = false;
	if(!options.count("input") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Calculating pvalue for selection using output of mode Frequency." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating evidence of selection for " + options["input"].as<std::string>() + ".\n";

	std::string line_freq, line_lin, read;

	//read in .freq and .lin files
	//obtain N, fN, k, fk for every SNP.
	//Calculate pvalue

	igzstream is_freq(options["input"].as<std::string>() + ".freq");
	if(is_freq.fail()) is_freq.open(options["input"].as<std::string>() + ".freq.gz");
	if(is_freq.fail()){
		std::cerr << "Error while opening file " << options["input"].as<std::string>() + ".freq(.gz)" << std::endl;
		exit(1);
	}
	igzstream is_lin(options["input"].as<std::string>() + ".lin");
	if(is_lin.fail()) is_lin.open(options["input"].as<std::string>() + ".freq.gz");
	if(is_lin.fail()){
		std::cerr << "Error while opening file " << options["input"].as<std::string>() + ".lin(.gz)" << std::endl;
		exit(1);
	}
	std::ofstream os(options["output"].as<std::string>() + ".sele");
	if(os.fail()){
		std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".sele" << std::endl;
		exit(1);
	}

	//skip headers
	getline(is_freq, line_freq);
	getline(is_lin, line_lin);

	os << line_lin << "\n";

	//precalculate logF, where logF[k] = log(k!)
	std::vector<float> logF;

	//read line by line
	int N, k;
	float fN, fk;
	std::vector<float> num_lin, num_freq;

	while(getline(is_freq, line_freq)){
		getline(is_lin, line_lin);

		//get fN, N, k, fk from the line
		std::stringstream s_freq(line_freq);
		std::stringstream s_lin(line_lin);

		s_freq >> read;
		os << read << " ";
		s_freq >> read;
		os << read << " ";
		s_lin >> read;
		s_lin >> read;

		int add_entries = 2;

		//read in k from s_lin and fk from s_freq
		if(logF.size() == 0){

			float foo;
			while(s_lin >> foo) num_lin.push_back(foo);
			num_freq.resize(num_lin.size()-add_entries);
			for(int i = 0; i < num_freq.size(); i++){
				s_freq >> num_freq[i];
			}
			N = (int)num_lin[num_lin.size() - add_entries - 1];
			logFactorial(logF, N);

		}else{

			for(int i = 0; i < num_lin.size(); i++){
				s_lin >> num_lin[i];
			}
			for(int i = 0; i < num_freq.size(); i++){
				s_freq >> num_freq[i];
			}

		}

		fN = num_freq[num_freq.size() - 1]; //frequency when N lineages are remaining

		if(fN <= 2){
			for(int i = 0; i < num_freq.size(); i++){
				os << "1 ";
			}
		}else{
			for(int i = 0; i < num_freq.size(); i++){
				os << log_pvalue(num_lin[i], num_freq[i], N, fN, logF) << " ";
			}
		}

		if(fN > 2){
			os << log_pvalue(num_lin[num_lin.size() - add_entries], (int)((fN+1.0)/2.0), N, fN, logF) << " ";
			os << log_pvalue(num_lin[num_lin.size() - add_entries+1], 2.0, N, fN, logF) << "\n";
		}else{
			os << "1 1\n";
		}

	}

	is_freq.close();
	is_lin.close();




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
Frequency(cxxopts::Options& options){

	//////////////////////////////////
	//Program options
	bool help = false;
	if(!options.count("input") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output. Optional: years_per_gen, first_snp, last_snp." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "..." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating frequency through time for " + options["input"].as<std::string>() + ".\n";

	std::string line, read;

	////////// PARSE DATA //////////

	AncMutIterators ancmut;
	ancmut.OpenFiles(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");

	int N = ancmut.NumTips();
	int L = ancmut.NumSnps();
	MarginalTree mtr;
	Muts::iterator it_mut;
	float num_bases_tree_persists = 0.0; 

	Data data(N,L);
	int N_total = 2*data.N-1;

	std::vector<float> logF;
	logFactorial(logF, data.N);

	int first_snp, last_snp;
	if(!options.count("first_snp")){
		first_snp = 0;
	}else{
		first_snp = options["first_snp"].as<int>();
	}
	if(!options.count("last_snp")){
		last_snp = data.L-1;
	}else{
		last_snp = options["last_snp"].as<int>();
	}

	///////// EPOCHES /////////
	float years_per_gen = 28.0;
	if(options.count("years_per_gen")){
		years_per_gen = options["years_per_gen"].as<float>();
	}

	int num_epochs;
	std::vector<float> epochs;
	float log_10 = std::log(10);
	if(options.count("bins")){

		double log_age = std::log(0);
		double age = 0; 

		double epoch_lower, epoch_upper, epoch_step;
		std::string str_epochs = options["bins"].as<std::string>();
		std::string tmp;
		int i = 0;
		tmp = "";
		while(str_epochs[i] != ','){
			tmp += str_epochs[i];
			i++;
			if(i == str_epochs.size()) break;
		}
		epoch_lower = std::stof(tmp);
		i++;
		if(i >= str_epochs.size()){
			std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
			exit(1);
		}
		tmp = "";
		while(str_epochs[i] != ','){
			tmp += str_epochs[i];
			i++;
			if(i == str_epochs.size()) break;
		}
		epoch_upper = std::stof(tmp);
		i++;
		if(i >= str_epochs.size()){
			std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
			exit(1);
		}
		tmp = "";
		while(str_epochs[i] != ','){
			tmp += str_epochs[i];
			i++;
			if(i == str_epochs.size()) break;
		}
		epoch_step = std::stof(tmp);

		int ep = 0;
		epochs.resize(1);
		epochs[ep] = 0.0;
		ep++; 
		double epoch_boundary = 0.0;
		if(log_age < epoch_lower && age != 0.0){
			epochs.push_back(age);
			ep++;
		}
		epoch_boundary = epoch_lower;
		while(epoch_boundary < epoch_upper){
			if(log_age < epoch_boundary){
				if(ep == 1 && age != 0.0) epochs.push_back(age);
				epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
				ep++;
			}
			epoch_boundary += epoch_step;
		}
		epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
		epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();

	}else{

		num_epochs = 31;
		epochs.resize(num_epochs);
		epochs[0] = 0.0;
		epochs[1] = 1e3/years_per_gen;
		for(int e = 2; e < num_epochs-1; e++){
			epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
		}
		epochs[num_epochs-1] = 1e8/years_per_gen;

	}


	//read mutations file
	Mutations mutations(data);
	mutations.Read(options["input"].as<std::string>() + ".mut");

	////////////////////
	//for a mutation, record how its frequency is changing

	//1. Get coordinates vector and sort. Use this to determine the number of lineages.
	//2. From branch on which mutation sits, record time of all coalescent events below.
	//3. Count.

	std::vector<int> index, index_mut;
	std::vector<float> coordinates_tree(N_total), coordinates_tree_unsrt(N_total), coordinates_mutation(N_total);
	int root = N_total-1;
	int i = 0;
	int snp_of_next_tree;

	std::ofstream os_freq(options["output"].as<std::string>() + ".freq");
	if(os_freq.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	}
	std::ofstream os_lin(options["output"].as<std::string>() + ".lin");
	if(os_lin.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	}

	int count_tree = 0;

	os_freq << "pos rs_id ";
	for(int ep = num_epochs-1; ep >= 0; ep--){
		os_freq << std::to_string(epochs[ep]) << " ";
	}
	os_freq << "TreeFreq DataFreq\n";

	os_lin << "pos rs_id ";
	for(int ep = num_epochs-1; ep >= 0; ep--){
		os_lin << std::to_string(epochs[ep]) << " ";
	}
	//os_lin << "when_mutation_has_freq2 when_mutation_has_freq3 when_mutation_has_freq4 when_mutation_has_freq5 when_mutation_has_freq6 when_mutation_has_freq7\n";
	os_lin << "when_DAF_is_half when_mutation_has_freq2\n";


	//read tree
	num_bases_tree_persists = ancmut.FirstSNP(mtr, it_mut);
	mtr.tree.GetCoordinates(coordinates_tree);
	coordinates_tree_unsrt = coordinates_tree;
	sortAndGetIndices(coordinates_tree, index);
	//std::sort(coordinates_tree.begin(), coordinates_tree.end());

	//std::vector<Leaves> leaves;
	//mtr.tree.FindAllLeaves(leaves);
	//std::cerr << leaves[*mutations.info[0].branch.begin()].num_leaves << std::endl;

	int current_tree = (*it_mut).tree;
	SNPInfo snp_info;
	int freq;
	float rec;
	for(int snp = first_snp; snp <= last_snp; snp++){

		snp_info = (*it_mut);

		int freq = 3;
		if(snp_info.freq.size() > 0){
			freq = 0;
			for(std::vector<int>::iterator it_freq = snp_info.freq.begin(); it_freq != snp_info.freq.end(); it_freq++){
				freq += *it_freq;
				if(freq > 2) break;
			}
		}

		if(snp_info.branch.size() == 1 && freq > 2 && !snp_info.flipped){

			if((*it_mut).tree != current_tree){
				current_tree = (*it_mut).tree;
				mtr.tree.GetCoordinates(coordinates_tree);
				coordinates_tree_unsrt = coordinates_tree; //unsorted coordinates of tree
				sortAndGetIndices(coordinates_tree, index);
				//std::sort(coordinates_tree.begin(), coordinates_tree.end()); //sorted coordinates of tree
			}
			assert(current_tree == snp_info.tree);

			if(snp_info.age_begin <= coordinates_tree[root]){

				int b = *snp_info.branch.begin();

				if(b != -1 && b != root){

					os_freq << snp_info.pos << " " << snp_info.rs_id << " ";
					os_lin  << snp_info.pos << " " << snp_info.rs_id << " ";

					int DAF = 0, DAF_half, num_lin_half = -1;
					std::fill(coordinates_mutation.begin(), coordinates_mutation.end(), -1.0);
					//get nodes below branch b
					//get coordinates from coordinates_tree_unsrt
					CopyCoordinates(b, coordinates_mutation, coordinates_tree_unsrt, mtr.tree, DAF); //get coordinates of branches below mutation
					DAF_half = (DAF+1)/2.0;
					coordinates_mutation[(*mtr.tree.nodes[b].parent).label] = coordinates_tree_unsrt[(*mtr.tree.nodes[b].parent).label]; //add to that list the coordinate of the parent of branch b

					//std::sort(coordinates_mutation.begin(), coordinates_mutation.end()); //sort coordinates_mutation
					sortAndGetIndices(coordinates_mutation, index_mut);

					std::vector<int> current_branches(data.N, -2);

					int num_carriers = 0, num_lineages = 1;
					//int k_when_mutation_appears = -1, k_when_mutation_has_freq2 = -1, k_when_mutation_has_freq3 = -1, k_when_mutation_has_freq4 = -1, k_when_mutation_has_freq5 = -1, k_when_mutation_has_freq6 = -1, k_when_mutation_has_freq7 = -1;
					int k_when_mutation_appears = -1, k_when_mutation_has_freq2 = -1, has_disappeared = -2, on_first_branch = -1;
					int n_mut = root;
					int n_tree = root;
					int ep = num_epochs-1; 

					//while epoch[ep] is larger than root of tree, number of lineages is 0
					while(coordinates_tree[n_tree] < epochs[ep]){
						os_freq << 0 << " ";
						os_lin  << 0 << " ";
						ep--;
					}

          //std::cerr << "start" << std::endl;

					do{

						if(num_carriers >= DAF_half && DAF_half > 1 && num_lin_half == -1){
							num_lin_half = num_lineages;
						}

						if(n_tree >= 0){
							//std::cerr << n_tree << " " << coordinates_tree.size() << " " << ep << " " << epochs.size() << std::endl;
							while(coordinates_tree[n_tree] <= epochs[ep]){ //n_tree+1 is above epoch[ep], n_tree is just below, but could span multiple epochs
								float num_muts = 0.0;
								if(k_when_mutation_appears != -1){ 

									if(has_disappeared == 1){
										os_freq << "0 ";
										os_lin  << num_lineages << " ";
									}else{
										os_freq << (num_carriers) << " ";
										os_lin  << num_lineages << " "; 
									}

								}else{
									os_freq << 0 << " ";
									os_lin  << num_lineages << " ";
								}

								ep--;
								if(ep == -1) break;
							}
						}

						assert(num_carriers >= 0);
						assert(num_lineages >= 0);
						int count = 0;
						for(int k = 0; k < num_carriers; k++){
							if(current_branches[k] == -1) count++;
						}
						assert(count == 0);
						count = 0;

						float coords = coordinates_tree[n_tree];
            if( coords != coordinates_mutation[n_mut] || (has_disappeared == 1) ){ //if n_tree is not next node of a branch onto which mutation falls, increase num_lineages
						
							//next coalescence does not involve carrier branches
							while(coords == coordinates_tree[n_tree]){
								if(index[n_tree] < data.N){
									assert(mtr.tree.nodes[index[n_tree]].child_left == NULL);
									num_lineages--;
									count++;
								}else{
									assert(mtr.tree.nodes[index[n_tree]].child_left != NULL);
									num_lineages++;
								}
								n_tree--;
								if(n_tree < 0) break;
							}

						}else{
						
							//next coalescence involves carrier branches
							assert(coordinates_tree[n_tree] == coordinates_mutation[n_mut]); //else, n_tree is the node of a branch onto which mutation falls

							while(coordinates_tree[n_tree] == coords){

                //if index don't match or coordinates_mutation == -1, then next event actually isn't a carrier (can happen if dates are identical among carriers and non-carriers)
								if(index[n_tree] != index_mut[n_mut] || coordinates_mutation[n_mut] == -1){
                  if(index[n_tree] < data.N){
										assert(mtr.tree.nodes[index[n_tree]].child_left == NULL);
                    num_lineages--;
										count++;
									}else{
										assert(mtr.tree.nodes[index[n_tree]].child_left != NULL);
                    num_lineages++;
									}
									n_tree--;
								}else{
                  //std::cerr << snp_info.pos << " " << coordinates_mutation[n_mut] << " " << coords << std::endl;
									assert(coordinates_mutation[n_mut] == coords);

									bool check_if_branch_is_found = false;
									if(k_when_mutation_appears == -1){
										num_lineages++;
										k_when_mutation_appears = num_lineages;

										assert(index[n_tree] == index_mut[n_mut]);
										assert((*mtr.tree.nodes[b].parent).label == index_mut[n_mut]);
										current_branches[0]     = mtr.tree.nodes[b].label; //current branches stores lineages that carry derived allele

										num_carriers = 1;
										on_first_branch = 1;
										has_disappeared = -1;

										check_if_branch_is_found = true;
									}else{
										on_first_branch = -1;
										for(int k = 0; k < num_carriers; k++){
											if(current_branches[k] >= 0){
												if(current_branches[k] == index_mut[n_mut]){
													check_if_branch_is_found = true;
													if(mtr.tree.nodes[index_mut[n_mut]].child_left == NULL){
														current_branches[k] = -1;
													}else{
														current_branches[k]              = (*mtr.tree.nodes[index_mut[n_mut]].child_left).label; //replace current_branch[k] by its children
														current_branches[num_carriers]   = (*mtr.tree.nodes[index_mut[n_mut]].child_right).label; 
														num_lineages++;
														num_carriers++;
													}
												}
											}
										}

									}

									assert(check_if_branch_is_found);
									n_tree--;
									n_mut--;	
								}
                
								if(n_tree < 0) break;
								if(n_mut < 0) break;
							}

							if(n_tree >= 0 && n_mut >= 0){
								assert(coordinates_mutation[n_mut] != coords);
							  assert(coordinates_tree[n_tree] != coords);
							}

						} 

						//num_carriers == 1 means (1+num_carriers) have the derived allele
						if(num_carriers >= 2){
							if(k_when_mutation_has_freq2 == -1){ //when this happens for the first time, record num_lineages at this point
								k_when_mutation_has_freq2 = num_lineages;
								//if(num_carriers > 1) k_when_mutation_has_freq2 -= num_carriers - 1; //this is needed because there might be many branches with exactly the same coordinates
								assert(k_when_mutation_has_freq2);
							}
						}

						for(int k = 0; k < num_carriers; k++){
							for(int l = num_carriers-1; l >= 0; l--){
								if(current_branches[l] != -1) break;
								if(current_branches[l] == -1){
									num_carriers--;
									num_lineages--;
								}
								if(num_carriers == 0) break;
							}
							if(k < num_carriers){
								if(current_branches[k] == -1){
									current_branches[k] = current_branches[num_carriers-1];
									num_carriers--;
									num_lineages--;
								}
							}
						}
						if(has_disappeared == -1 && num_carriers == 0){
							has_disappeared = 1;
						}

						if(0){
						if(n_tree >= 0){
							while(coordinates_tree[n_tree] <= epochs[ep]){ //n_tree+1 is above epoch[ep], n_tree is just below, but could span multiple epochs
								float num_muts = 0.0;
								if(k_when_mutation_appears != -1){ 

									if(has_disappeared == 1){
										os_freq << "0 ";
										os_lin  << num_lineages << " ";
									}else{
										os_freq << (num_carriers) << " ";
										os_lin  << num_lineages << " "; 
									}

								}else{
									os_freq << 0 << " ";
									os_lin  << num_lineages << " ";
								}

								ep--;
								if(ep == -1) break;
								}
							}
						}

						//std::cerr << "carriers: " << snp_info.pos << " " << num_carriers << " " << num_lineages << " " << b << " " << current_branches[0] << std::endl;
						assert(num_carriers <= num_lineages);
            assert(num_lineages >= 0);
						assert(num_carriers >= 0);

					}while(n_tree >= 0 && ep >= 0);

					os_freq << " " << num_carriers << " ";

					//assert(num_carriers == DAF);
					int carriers = 0;
					for(std::vector<int>::iterator it_freq = snp_info.freq.begin(); it_freq != snp_info.freq.end(); it_freq++){
						carriers += *it_freq;
					}
					os_freq << carriers << "\n";

					os_lin  << num_lin_half << " ";
					os_lin  << k_when_mutation_has_freq2 << "\n";
				}

			}

		}
		ancmut.NextSNP(mtr, it_mut);
	}

	os_freq.close();
	os_lin.close();

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
SDS(cxxopts::Options& options){

	//////////////////////////////////
	//Program options
	bool help = false;
	if(!options.count("input") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output. Optional: years_per_gen, first_snp, last_snp." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "..." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating frequency through time for " + options["input"].as<std::string>() + ".\n";

	std::string line, read;

	////////// PARSE DATA //////////

	int N;
	igzstream is_N;
	is_N.open(options["input"].as<std::string>() + ".anc");
	if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
	if(is_N.fail()){
		std::cerr << "Error while opening .anc file." << std::endl;
		exit(1);
	} 
	is_N.ignore(256, ' ');
	is_N >> N;
	is_N.close();

	int L = 0;
	igzstream is_L;
	is_L.open(options["input"].as<std::string>() + ".mut");
	if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz"); 
	if(is_L.fail()){
		std::cerr << "Error while opening .mut file." << std::endl;
		exit(1);
	} 

	std::string unused;
	std::getline(is_L, unused); 
	while ( std::getline(is_L, unused) ){
		++L;
	}

	Data data(N,L);
	int N_total = 2*data.N-1; 

	///////// EPOCHES /////////
	float years_per_gen = 28.0;
	if(options.count("years_per_gen")){
		years_per_gen = options["years_per_gen"].as<float>();
	}

	int num_epochs;
	std::vector<float> epochs;
	float log_10 = std::log(10);
	if(options.count("bins")){

		double log_age = std::log(0);
		double age = 0;

		double epoch_lower, epoch_upper, epoch_step;
		std::string str_epochs = options["bins"].as<std::string>();
		std::string tmp;
		int i = 0;
		tmp = "";
		while(str_epochs[i] != ','){
			tmp += str_epochs[i];
			i++;
			if(i == str_epochs.size()) break;
		}
		epoch_lower = std::stof(tmp);
		i++;
		if(i >= str_epochs.size()){
			std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
			exit(1);
		}
		tmp = "";
		while(str_epochs[i] != ','){
			tmp += str_epochs[i];
			i++;
			if(i == str_epochs.size()) break;
		}
		epoch_upper = std::stof(tmp);
		i++;
		if(i >= str_epochs.size()){
			std::cerr << "Error: epochs format is wrong. Specify x,y,stepsize." << std::endl;
			exit(1);
		}
		tmp = "";
		while(str_epochs[i] != ','){
			tmp += str_epochs[i];
			i++;
			if(i == str_epochs.size()) break;
		}
		epoch_step = std::stof(tmp);

		int ep = 0;
		epochs.resize(1);
		epochs[ep] = 0.0;
		ep++; 
		double epoch_boundary = 0.0;
		if(log_age < epoch_lower && age != 0.0){
			epochs.push_back(age);
			ep++;
		}
		epoch_boundary = epoch_lower;
		while(epoch_boundary < epoch_upper){
			if(log_age < epoch_boundary){
				if(ep == 1 && age != 0.0) epochs.push_back(age);
				epochs.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
				ep++;
			}
			epoch_boundary += epoch_step;
		}
		epochs.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
		epochs.push_back( std::max(1e8, 10.0*epochs[epochs.size()-1])/years_per_gen );
		num_epochs = epochs.size();

	}else{

		num_epochs = 31;
		epochs.resize(num_epochs);
		epochs[0] = 0.0;
		epochs[1] = 1e3/years_per_gen;
		for(int e = 2; e < num_epochs-1; e++){
			epochs[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs-3.0) ))/years_per_gen;
		}
		epochs[num_epochs-1] = 1e8/years_per_gen;

	}


	//read mutations file
	Mutations mutations(data);
	mutations.Read(options["input"].as<std::string>() + ".mut");

	int first_snp, last_snp;
	if(!options.count("first_snp")){
		first_snp = 0;
	}else{
		first_snp = options["first_snp"].as<int>();
	}
	if(!options.count("last_snp")){
		last_snp = data.L-1;
	}else{
		last_snp = options["last_snp"].as<int>();
	}


	////////////////////
	//for a mutation, record how its frequency is changing

	//1. Get coordinates vector and sort. Use this to determine the number of lineages.
	//2. From branch on which mutation sits, record time of all coalescent events below.
	//3. Count.

	MarginalTree mtr;
	int root = N_total-1;
	int i = 0;
	int snp_of_next_tree;

	igzstream is_anc(options["input"].as<std::string>() + ".anc");
	if(is_anc.fail()) is_anc.open(options["input"].as<std::string>() + ".anc.gz");
	if(is_anc.fail()){
		std::cerr << "Error while opening .anc file." << std::endl;
		exit(1);
	}
	std::ofstream os(options["output"].as<std::string>() + ".SDS");
	if(os.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	}

	int count_tree = 0;

	os << "pos rs_id ";
	os << "rSDS\n";

	getline(is_anc,line);
	getline(is_anc,line);
	getline(is_anc,line);

	//read tree
	std::vector<Leaves> leaves;
	mtr.Read(line, N);
	mtr.tree.FindAllLeaves(leaves);

	SNPInfo snp_info;
	int freq;
	float rec;
	for(int snp = first_snp; snp <= last_snp; snp++){

		snp_info = mutations.info[snp];

		int freq = 3;
		for(std::vector<int>::iterator it_freq = snp_info.freq.begin(); it_freq != snp_info.freq.end(); it_freq++){
			freq += *it_freq;
			if(freq > 2) break;
		}

		if(snp_info.branch.size() == 1 && freq > 2 && !snp_info.flipped){

			if(count_tree < snp_info.tree){
				while(count_tree < snp_info.tree){
					if(!getline(is_anc,line)){
						break; 
					};
					count_tree++;
				} 
				assert(count_tree == snp_info.tree);

				//read tree
				mtr.Read(line, N);
				mtr.tree.FindAllLeaves(leaves);
			}

			int b = *snp_info.branch.begin();

			if(b != -1 && b != root){

				os << snp_info.pos << " " << snp_info.rs_id << " ";

				std::sort(leaves[b].member.begin(), leaves[b].member.end());
				//calculate SDS
				double aSDS = 0.0, dSDS = 0.0;
				std::vector<int>::iterator it_member = leaves[b].member.begin();
				int k = 0;
				for(std::vector<Node>::iterator it_node = mtr.tree.nodes.begin(); it_node != std::next(mtr.tree.nodes.begin(), data.N); it_node++){
					
					if(it_member != leaves[b].member.end()){
						if((*it_node).label == *it_member){
							dSDS += (*it_node).branch_length; 
							it_member++;
							k++;
						}else{
							aSDS += (*it_node).branch_length;
						}
					}else{
						aSDS += (*it_node).branch_length;
					}
				}
				assert(it_member == leaves[b].member.end());
				os << log( (aSDS/dSDS) * leaves[b].num_leaves)/(data.N - leaves[b].num_leaves) << "\n";


			}

		}
	}

	is_anc.close();
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
FreqDiff(cxxopts::Options& options){

	//////////////////////////////////
	//Program options
	bool help = false;
	if(!options.count("input") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Calculating pvalue for selection using output of mode Frequency." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Calculating frequency changes for " + options["input"].as<std::string>() + ".\n";

	//TODO: make this for multiple chromosomes

	std::string line_freq, line_lin, read1, read2;

	std::vector<std::vector<double>> mean, sd;
	std::vector<std::vector<int>> freq_count;

	std::string line;
	std::vector<std::string> chromosomes;
	std::vector<std::string> filenames_in, filenames_out;
	if(options.count("chr") > 0){

		igzstream is_chr(options["chr"].as<std::string>());
		if(is_chr.fail()){
			std::cerr << "Error while opening file " << options["chr"].as<std::string>() << std::endl;
		}
		while(getline(is_chr, line)){
			chromosomes.push_back(line);
			filenames_in.push_back(options["input"].as<std::string>() + "_chr" + line);
			filenames_out.push_back(options["output"].as<std::string>() + "_chr" + line);
		}
		is_chr.close();

	}else{
		chromosomes.resize(1);
		chromosomes[0] = "NA"; 
		filenames_in.push_back(options["input"].as<std::string>());
		filenames_out.push_back(options["output"].as<std::string>());
	}

	int N, k;
	float fN, fk;
	std::vector<float> num_lin, num_freq;
	double diff;

	//read in .freq and .lin files
	//obtain N, fN, k, fk for every SNP.
	//Calculate pvalue

	for(int chr = 0; chr < chromosomes.size(); chr++){

		igzstream is_freq(filenames_in[chr] + ".freq");
		if(is_freq.fail()) is_freq.open(filenames_in[chr] + ".freq.gz");
		if(is_freq.fail()){
			std::cerr << "Error while opening file " << filenames_in[chr] + ".freq(.gz)" << std::endl;
			exit(1);
		}
		igzstream is_lin(filenames_in[chr] + ".lin");
		if(is_lin.fail()) is_lin.open(filenames_in[chr] + ".freq.gz");
		if(is_lin.fail()){
			std::cerr << "Error while opening file " << filenames_in[chr] + ".lin(.gz)" << std::endl;
			exit(1);
		}
		std::ofstream os(filenames_out[chr] + ".freqdiff");
		if(os.fail()){
			std::cerr << "Error while opening file " << filenames_out[chr] + ".freqdiff" << std::endl;
			exit(1);
		}

		//skip headers
		getline(is_freq, line_freq);
		getline(is_lin, line_lin);

		os << line_freq.substr(0, line_freq.size() - 9) << "\n";

		//read line by line
		while(getline(is_freq, line_freq)){
			getline(is_lin, line_lin);

			//get fN, N, k, fk from the line
			std::stringstream s_freq(line_freq);
			std::stringstream s_lin(line_lin);

			s_freq >> read1;
			os << read1 << " ";
			s_freq >> read2;
			os << read2 << " ";
			s_lin >> read1;
			s_lin >> read2;

			int add_entries = 2;

			//read in k from s_lin and fk from s_freq
			if(num_lin.size() == 0){

				float foo;
				while(s_lin >> foo) num_lin.push_back(foo);
				num_lin.resize(num_lin.size()-add_entries);
				std::reverse(num_lin.begin(), num_lin.end());
				num_freq.resize(num_lin.size());
				for(int i = num_freq.size()-1; i >= 0; i--){
					s_freq >> num_freq[i];
				}
				N = (int)num_lin[0];
				mean.resize(N);
				sd.resize(N);
				freq_count.resize(N);
				for(int i = 0; i < N; i++){
					mean[i].resize(num_freq.size()-1);
					std::fill(mean[i].begin(), mean[i].end(), 0);
					sd[i].resize(num_freq.size()-1);
					std::fill(sd[i].begin(), sd[i].end(), 0);
					freq_count[i].resize(num_freq.size()-1);
					std::fill(freq_count[i].begin(), freq_count[i].end(), 0);
				}

			}else{

				for(int i = num_lin.size()-1; i >= 0; i--){
					s_lin >> num_lin[i];
				}
				for(int i = num_freq.size() - 1; i >= 0; i--){
					s_freq >> num_freq[i];
				}

			}

			fN = num_freq[0]; //frequency when N lineages are remaining


			if(0){
				for(int i = 0; i < num_freq.size(); i++){
					std::cerr << num_freq[i] << "," << num_lin[i] << " ";
				}
				std::cerr << std::endl;
				exit(1);
			}

			//num_freq and num_lin
			//1. calculate frequency difference between time periods
			//2. While doing this, calculate average and sd for each present-day frequency and each time
			//3. Reread all files, standardise differences.

			for(int i = num_lin.size()-2; i >= 0; i--){
				if(num_freq[i+1] > 0 && num_lin[i+1] > 0.1*N){
					diff = ((double)num_freq[i+1])/num_lin[i+1] - ((double)num_freq[i])/num_lin[i];
				}else{
					diff = -10;
				}
				//exclude once num_freq[i+1] == 0
				if(num_freq[i+1] > 0 && num_lin[i+1] > 0.1*N){
					mean[fN][i] += diff;
					sd[fN][i]   += diff*diff;
					freq_count[fN][i]++;
				}
				os << diff << " ";
			}
			os << fN << "\n";

		}

		os.close();
		is_freq.close();
		is_lin.close();

	}

	for(int f = 0; f < N; f++){
		//std::cerr << f << std::endl;
		for(int i = 0; i < freq_count[f].size(); i++){
			if(freq_count[f][i] > 0){
				mean[f][i] /= (double)freq_count[f][i];
				sd[f][i]    = sqrt( (sd[f][i] - freq_count[f][i]*mean[f][i]*mean[f][i])/((double)freq_count[f][i] - 1.0) );
			}else{
				mean[f][i] = 0.0;
				sd[f][i]   = 0.0;
			}
		}
	}

	for(int chr = 0; chr < chromosomes.size(); chr++){

		igzstream is_freq(filenames_out[chr] + ".freqdiff");
		if(is_freq.fail()) is_freq.open(filenames_out[chr] + ".freqdiff.gz");
		if(is_freq.fail()){
			std::cerr << "Error while opening file " << filenames_out[chr] + ".freqdiff(.gz)" << std::endl;
			exit(1);
		}
		std::ofstream os(filenames_out[chr] + ".zfreqdiff");
		if(os.fail()){
			std::cerr << "Error while opening file " << filenames_out[chr] + ".zfreqdiff" << std::endl;
			exit(1);
		}

		//skip headers
		getline(is_freq, line_freq);

		os << line_freq << "\n";

		//read line by line
		while(getline(is_freq, line_freq)){

			//get fN, N, k, fk from the line
			std::stringstream s_freq(line_freq);

			s_freq >> read1;
			s_freq >> read2;

			int add_entries = 2;

			//read in k from s_lin and fk from s_freq
			std::fill(num_freq.begin(), num_freq.end(), 0);
			for(int i = num_freq.size()-2; i >= 0; i--){
				s_freq >> num_freq[i];
			}
			s_freq >> fN;

			if(fN > 1){

				os << read1 << " " << read2 << " ";
				for(int i = num_freq.size()-2; i >= 0; i--){
					//std::cerr << num_freq[i] << " " << fN << " " << mean[fN][i] << " " << sd[fN][i] << " " << (num_freq[i] - mean[fN][i])/sd[fN][i] << std::endl;
					if(num_freq[i] != -10){
						os << (num_freq[i] - mean[fN][i])/sd[fN][i] << " ";
					}else{
						os << "NA ";
					}
					//os << num_freq[i]/sd[fN][i] << " ";
				}
				os << fN << "\n";
			}

		}

		os.close();
		is_freq.close();

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



void 
Quality(cxxopts::Options& options){

	//////////////////////////////////
	//Program options
	bool help = false;
	if(!options.count("input") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: input, output. Optional: first_snp, last_snp." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "..." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Annotating quality of SNPs for " + options["input"].as<std::string>() + ".\n";

	std::string line, read;

	////////// PARSE DATA //////////

	int N;
	igzstream is_N;
	is_N.open(options["input"].as<std::string>() + ".anc");
	if(is_N.fail()) is_N.open(options["input"].as<std::string>() + ".anc.gz");
	if(is_N.fail()){
		std::cerr << "Error while opening .anc file." << std::endl;
		exit(1);
	} 
	is_N.ignore(256, ' ');
	is_N >> N;
	is_N.close();

	int L = 0;
	igzstream is_L;
	is_L.open(options["input"].as<std::string>() + ".mut");
	if(is_L.fail()) is_L.open(options["input"].as<std::string>() + ".mut.gz"); 
	if(is_L.fail()){
		std::cerr << "Error while opening .mut file." << std::endl;
		exit(1);
	} 
	std::string unused;
	std::getline(is_L, unused); 
	while ( std::getline(is_L, unused) ){
		++L;
	}


	Data data(N,L);
	int N_total = 2*data.N-1;

	int first_snp, last_snp;
	if(!options.count("first_snp")){
		first_snp = 0;
	}else{
		first_snp = options["first_snp"].as<int>();
	}
	if(!options.count("last_snp")){
		last_snp = data.L-1;
	}else{
		last_snp = options["last_snp"].as<int>();
	}


	//read mutations file
	Mutations mutations(data);
	mutations.Read(options["input"].as<std::string>() + ".mut");

	std::vector<int> SNPmapping(mutations.info.size()), tree_index(mutations.info.size());
	std::vector<int>::iterator it_SNPmapping = SNPmapping.begin(), it_tree_index = tree_index.begin();
	for(std::vector<SNPInfo>::iterator it_info = mutations.info.begin(); it_info != mutations.info.end(); it_info++){
		*it_SNPmapping = ((*it_info).branch.size() > 1);
		*it_tree_index = (*it_info).tree;
		it_tree_index++;
		it_SNPmapping++;
	}


	//open anc file
	MarginalTree mtr;
	int root = N_total-1;
	int snp_of_next_tree;

	igzstream is_anc(options["input"].as<std::string>() + ".anc");
	if(is_anc.fail()) is_anc.open(options["input"].as<std::string>() + ".anc.gz");
	if(is_anc.fail()){
		std::cerr << "Error while opening .anc file." << std::endl;
		exit(1);
	}
	std::ofstream os(options["output"].as<std::string>() + ".qual");
	if(is_anc.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	}

	getline(is_anc,line);
	getline(is_anc,line);
	getline(is_anc,line);

	//read tree
	mtr.Read(line, N);
	float num_snps_on_tree = 0.0;
	float frac_branches_with_snp = 0.0;
	int num_snps_not_mapping_to_tree = 0;
	int i = first_snp;
	while(tree_index[i] == tree_index[first_snp]){
		num_snps_not_mapping_to_tree += SNPmapping[i];
		i++;
	}
	int num_snps_tree_persisting = i - first_snp;

	for(std::vector<Node>::iterator it_node = std::next(mtr.tree.nodes.begin(), data.N); it_node != mtr.tree.nodes.end(); it_node++){
		if((*it_node).num_events >= 1.0) frac_branches_with_snp += 1.0;
		num_snps_on_tree += (*it_node).num_events;
	}
	frac_branches_with_snp /= data.N - 1.0;

	SNPInfo snp_info;
	int freq;
	int count_tree = 0;
	int num_snps_not_mapping = 0;

	os << "ID pos frac_branches_with_snp num_snps_on_tree fraction_snps_not_mapping\n";

	if(last_snp - first_snp < 1000){
		std::cerr << "Need at least 1000 SNPs." << std::endl;
		exit(1);
	}

	for(int snp = first_snp; snp < first_snp + 500; snp++){
		num_snps_not_mapping += SNPmapping[snp];
	}

	for(int snp = first_snp; snp <= last_snp; snp++){

		snp_info = mutations.info[snp];

		if(snp - first_snp < 500){
			num_snps_not_mapping += SNPmapping[snp+500];
		}else if(last_snp - snp < 500){
			num_snps_not_mapping -= SNPmapping[snp-500];
		}else{
			num_snps_not_mapping += SNPmapping[snp+500] - SNPmapping[snp - 500];
		}

		if(count_tree < snp_info.tree){
			while(count_tree < snp_info.tree){
				if(!getline(is_anc,line)){
					break; 
				};
				count_tree++;
			} 
			assert(count_tree == snp_info.tree);

			//read tree
			mtr.Read(line, N);

			i = snp; 
			num_snps_not_mapping_to_tree = 0;
			while(tree_index[i] == tree_index[snp]){
				num_snps_not_mapping_to_tree += SNPmapping[i];
				i++;
			}
			num_snps_tree_persisting = i - snp;

			num_snps_on_tree = 0.0;
			frac_branches_with_snp = 0.0;
			for(std::vector<Node>::iterator it_node = std::next(mtr.tree.nodes.begin(), data.N); it_node != mtr.tree.nodes.end(); it_node++){
				if((*it_node).num_events >= 1.0) frac_branches_with_snp += 1.0;
				num_snps_on_tree += (*it_node).num_events;
			}
			frac_branches_with_snp /= data.N - 1.0;

		}

		int b = *snp_info.branch.begin();
		int carriers = 0;
		for(std::vector<int>::iterator it_freq = snp_info.freq.begin(); it_freq != snp_info.freq.end(); it_freq++){
			carriers += *it_freq;
		}
		os << snp_info.rs_id << " " << snp_info.pos << " ";
		os << frac_branches_with_snp << " " << num_snps_on_tree << " ";
		//os << num_snps_not_mapping_to_tree/((float) num_snps_tree_persisting) << " " << num_snps_tree_persisting << " ";

		if(snp - first_snp < 500){
			os << num_snps_not_mapping/((float) (500 + snp - first_snp + 1));
		}else if(last_snp - snp < 500){
			os << num_snps_not_mapping/((float) (500 + last_snp - snp));
		}else{
			os << num_snps_not_mapping/1000.0;
		}

		os << "\n";


	}

	is_anc.close();
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


int main(int argc, char* argv[]){

	//////////////////////////////////
	//Program options
	cxxopts::Options options("RelateSelection");
	options.add_options()
		("help", "Print help")
		("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>())
		("first_snp", "Index of first SNP. Optional.", cxxopts::value<int>())
		("last_snp", "Index of last SNP. Optional.", cxxopts::value<int>())
		("threshold", "Optional: Threshold for number of mutations that trees need for inclusion. Default = 0.", cxxopts::value<int>())
		("years_per_gen", "Optional: Years per generation (float). Default: 28.", cxxopts::value<float>())
		("bins", "Specify epoch bins. Format: lower, upper, stepsize for function c(0,10^seq(lower, upper, stepsize)).", cxxopts::value<std::string>())
		("chr", "Optional: File specifying chromosomes to use. Overrides first_chr, last_chr.", cxxopts::value<std::string>()) 
		("i,input", "Filename of .anc and .mut file without file extension", cxxopts::value<std::string>())
		("o,output", "Output file", cxxopts::value<std::string>());

	options.parse(argc, argv);
	std::string mode = options["mode"].as<std::string>();

	if(!mode.compare("Selection")){

		Selection(options);

	}else if(!mode.compare("Frequency")){

		Frequency(options);

	}else if(!mode.compare("Quality")){

		Quality(options);

	}else if(!mode.compare("SDS")){

		SDS(options);

	}else if(!mode.compare("FreqDiff")){

		FreqDiff(options);

	}else{

		std::cout << "####### error #######" << std::endl;
		std::cout << "Invalid or missing mode." << std::endl;
		std::cout << "Options for --mode are:" << std::endl;
		std::cout << "Frequency, Selection, Quality, SDS, FreqDiff." << std::endl;

	}

	bool help = false;
	if(!options.count("mode")){
		std::cout << "Not enough arguments supplied." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		exit(0);
	}

}

