#include <iostream>
#include "plot.hpp"
#include "sample.hpp"
#include "anc.hpp"
#include "anc_builder.hpp"
#include "tree_comparer.hpp"
#include "cxxopts.hpp"

#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

int FinalizePopulationSize(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Estimate population size using coalescence rate." << std::endl;
		exit(0);
	} 

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Finalizing coalescence rate..." << std::endl;

	////////// read labels of sequences //////////

	std::string line, seq;
	std::string read;

	/////////////////// INITIALIZE //////////////////////

	int num_epochs;
	std::vector<float> epoch;

	FILE* fp = fopen((options["output"].as<std::string>() + ".bin").c_str(),"rb");
	assert(fp != NULL);

	fread(&num_epochs, sizeof(int), 1, fp);
	epoch.resize(num_epochs);
	fread(&epoch[0], sizeof(float), num_epochs, fp);
	std::vector<CollapsedMatrix<float>> coalescent_rate_data(num_epochs);
	for(int e = 0; e < num_epochs; e++){
		coalescent_rate_data[e].ReadFromFile(fp);
	}
	fclose(fp);
	int N = coalescent_rate_data[0].size();

	std::vector<CollapsedMatrix<float>> coalescent_rate_num(num_epochs);
	for(std::vector<CollapsedMatrix<float>>::iterator it_c = coalescent_rate_num.begin(); it_c != coalescent_rate_num.end(); it_c++){
		(*it_c).resize(1,1);
		std::fill((*it_c).vbegin(), (*it_c).vend(), 0.0);
	}

	std::vector<CollapsedMatrix<float>> coalescent_rate_denom(num_epochs);
	for(std::vector<CollapsedMatrix<float>>::iterator it_c = coalescent_rate_denom.begin(); it_c != coalescent_rate_denom.end(); it_c++){
		(*it_c).resize(1,1);
		std::fill((*it_c).vbegin(), (*it_c).vend(), 0.0);
	}

	for(int i = 0; i < N; i++){
		for(int j = i+1; j < N; j++){

			//int i = 0, j = 1;
			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_data  = coalescent_rate_data.begin();
			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_num   = coalescent_rate_num.begin();
			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_denom = coalescent_rate_denom.begin();

			for(; it_coalescent_rate_data != std::prev(coalescent_rate_data.end(),1);){
				//i < j always holds
				//if((*it_coalescent_rate_data)[i][j] > 0.0){
				//  (*it_coalescent_rate)[0][0] += (*it_coalescent_rate_data)[i][j]/(*it_coalescent_rate_data)[j][i];
				//}
				(*it_coalescent_rate_num)[0][0] += (*it_coalescent_rate_data)[i][j];
				(*it_coalescent_rate_denom)[0][0] += (*it_coalescent_rate_data)[j][i];

				it_coalescent_rate_num++;
				it_coalescent_rate_denom++;
				it_coalescent_rate_data++;
			}    

		}
	}

	for(int e = 0; e < num_epochs; e++){
		if(coalescent_rate_denom[e][0][0] == 0){
			coalescent_rate_num[e][0][0] = 0.0;
			coalescent_rate_denom[e][0][0] = 0.0;
		}
	}

	std::ofstream os(options["output"].as<std::string>() + ".coal");
	if(os.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	} 

	os << "group1\n";

	for(int e = 0; e < num_epochs; e++){
		os << epoch[e]  << " "; 
	}
	os << "\n";

	std::vector<double> coal(epoch.size());

	os << 0 << " " << 0 << " ";
	for(int e = 0; e < num_epochs; e++){
		//coal[e] =  2.0 * coalescent_rate[e][0][0]/((float) (N*(N-1.0)));
		coal[e] = coalescent_rate_num[e][0][0]/coalescent_rate_denom[e][0][0];
		os << coal[e] << " ";
	}
	os << "\n"; 

	for(int e = 0; e < num_epochs; e++){
		if(coal[e] != 0.0) coal[e] = 0.5/coal[e]; 
	}

	plot p(60,10);
	p.draw(epoch, coal);

	os.close();

	/////////////////////////////////////////////
	//Resource Usage

	rusage usage;
	getrusage(RUSAGE_SELF, &usage);

	std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6) << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
	std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

	return 0;

}

int FinalizePopulationSizeByGroup(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("poplabels") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: poplabels, output. (Optional: bins - recommended with aDNA)" << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Estimate population size using coalescence rate." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Finalizing coalescence rate..." << std::endl;

	////////// read labels of sequences //////////

	Sample sample;
	sample.Read(options["poplabels"].as<std::string>());
	int N = sample.group_of_haplotype.size();

	/////////////////// INITIALIZE EPOCHES //////////////////////

	int num_epochs;
	std::vector<float> epoch;

	FILE* fp = fopen((options["output"].as<std::string>() + ".bin").c_str(),"rb");
	assert(fp != NULL);

	fread(&num_epochs, sizeof(int), 1, fp);
	epoch.resize(num_epochs);
	fread(&epoch[0], sizeof(float), num_epochs, fp);
	
	std::vector<CollapsedMatrix<float>> coalescent_rate_data(num_epochs);

	for(int e = 0; e < num_epochs; e++){
		coalescent_rate_data[e].ReadFromFile(fp);
		if(coalescent_rate_data[e].size() != N || coalescent_rate_data[e].subVectorSize(0) != N){
			std::cerr << N << " " << coalescent_rate_data[e].size() << std::endl;
			std::cerr << "Error: number of haplotypes in anc/mut does not match number of samples in .poplabels file" << std::endl;
			std::cerr << "You can just rerun this step using:" << std::endl;
			std::cerr << "PATH_TO_RELATE/bin/RelateCoalescentRate --mode FinalizePopulationSize -o example --poplabels example.poplabels" << std::endl;
			exit(1);
		}
	}
	fclose(fp);


	//I would need to check which epochs are proper, what is the minimum age per group, and then aggregate over improper epochs.
	std::vector<int> proper_epochs(num_epochs, 1);
	std::vector<float> group_min_age(sample.groups.size(), 0);
	std::vector<float> sample_ages(N, 0);
	if(options.count("bins")){

      float years_per_gen = 28.0;
      if(options.count("years_per_gen")){
        years_per_gen = options["years_per_gen"].as<float>();
      }

      int num_epochs_tmp;
      std::vector<float> epochs_tmp;
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
        epochs_tmp.resize(1);
        epochs_tmp[ep] = 0.0;
        ep++; 
        double epoch_boundary = 0.0;
        if(log_age < epoch_lower && age != 0.0){
          epochs_tmp.push_back(age);
          ep++;
        }
        epoch_boundary = epoch_lower;
        while(epoch_boundary < epoch_upper){
          if(log_age < epoch_boundary){
            if(ep == 1 && age != 0.0) epochs_tmp.push_back(age);
            epochs_tmp.push_back( std::exp(log_10 * epoch_boundary)/years_per_gen );
            ep++;
          }
          epoch_boundary += epoch_step;
        }
        epochs_tmp.push_back( std::exp(log_10 * epoch_upper)/years_per_gen );
        epochs_tmp.push_back( std::max(1e8, 10.0*epochs_tmp[epochs_tmp.size()-1])/years_per_gen );
        num_epochs_tmp = epochs_tmp.size();

      }else{
        num_epochs_tmp = 31;
        epochs_tmp.resize(num_epochs_tmp);
        epochs_tmp[0] = 0.0;
        epochs_tmp[1] = 1e3/years_per_gen;
        for(int e = 2; e < num_epochs_tmp-1; e++){
          epochs_tmp[e] = std::exp( log_10 * ( 3.0 + 4.0 * (e-1.0)/(num_epochs_tmp-3.0) ))/years_per_gen;
        }
        epochs_tmp[num_epochs_tmp-1] = 1e8/years_per_gen;
      }

			int e = 0, f = 0;
			for(; f < num_epochs_tmp; f++){
        while(epoch[e] < epochs_tmp[f]){
          e++;
					if(e == epoch.size()) break;
				}
	      if(e < epoch.size()){
				  if(epoch[e] != epochs_tmp[f]) break;
				}
			}
			if(f < num_epochs_tmp){
				std::cerr << "Warning: not all bin edges are present in the original epochs (different --bin?)" << std::endl;
			}else{
				std::fill(proper_epochs.begin(), proper_epochs.end(), 0);
				int e_tmp = 0;
				for(int e = 0; e < num_epochs; e++){
					while(epochs_tmp[e_tmp] < epoch[e]){
						e_tmp++;
						if(e_tmp == epochs_tmp.size()) break;
					}
					if(e_tmp == epochs_tmp.size()) break;
					if(epoch[e] == epochs_tmp[e_tmp]){
						proper_epochs[e] = 1;
					}
				}

				for(int i = 0; i < N; i++){
					float age = epoch[0];
					bool found = false;
					for(int e = 0; e < num_epochs; e++){
						for(int j = 0; j < N; j++){
							if(coalescent_rate_data[e][i][j] > 0){
								age = epoch[e];
								found = true;
								sample_ages[i] = age;
								break;
							} 
						}
						if(found) break;
					}
				}

				std::fill(group_min_age.begin(), group_min_age.end(), std::numeric_limits<int>::max());
				for(int i = 0 ; i < N; i++){
					if(sample_ages[i] < group_min_age[sample.group_of_haplotype[i]]){
						group_min_age[sample.group_of_haplotype[i]] = sample_ages[i];
					}
				}

				for (int g = 0; g < static_cast<int>(sample.groups.size()); ++g) {
					float t = group_min_age[g];

					auto it = std::lower_bound(epoch.begin(), epoch.end(), t);
					int idx;

					if (it == epoch.begin()) {
							idx = 0;
					} else if (it == epoch.end()) {
							idx = static_cast<int>(epoch.size() - 1);
					} else {
							int hi = static_cast<int>(it - epoch.begin());
							int lo = hi - 1;
							// Tie-break toward the lower epoch with <=
							idx = (t - epoch[lo] <= epoch[hi] - t) ? lo : hi;
					}

					// Store the closest epoch index (consider using a separate int vector)
					group_min_age[g] = static_cast<float>(idx);
				}
			}
	}

	/////////////////// SUMMARIZE //////////////////////

	std::vector<CollapsedMatrix<float>> coalescent_rate_num(num_epochs), coalescent_rate_denom(num_epochs);
	for(std::vector<CollapsedMatrix<float>>::iterator it_c = coalescent_rate_num.begin(); it_c != coalescent_rate_num.end(); it_c++){
		(*it_c).resize(sample.groups.size(), sample.groups.size());
		std::fill((*it_c).vbegin(), (*it_c).vend(), 0.0);
	}
	for(std::vector<CollapsedMatrix<float>>::iterator it_c = coalescent_rate_denom.begin(); it_c != coalescent_rate_denom.end(); it_c++){
		(*it_c).resize(sample.groups.size(), sample.groups.size());
		std::fill((*it_c).vbegin(), (*it_c).vend(), 0.0);
	}

	for(int i = 0; i < N; i++){
		for(int j = i+1; j < N; j++){

			//choose entry for groups i, j
			int group_i = sample.group_of_haplotype[i];
			int group_j = sample.group_of_haplotype[j];
			//make index of group_i < group_j
			if(group_i > group_j){
				int foo = group_i;
				group_i = group_j;
				group_j = foo;
			}

			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_data   = coalescent_rate_data.begin();
			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_num    = coalescent_rate_num.begin();
			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_denom  = coalescent_rate_denom.begin();

			for(; it_coalescent_rate_data != std::prev(coalescent_rate_data.end(),1);){
				//i < j always holds 

				//if((*it_coalescent_rate_data)[i][j] == 0.0){
				//  (*it_coalescent_rate)[group_i][group_j] += 0.0;
				//}else{ 
				//  (*it_coalescent_rate)[group_i][group_j] += (*it_coalescent_rate_data)[i][j]/(*it_coalescent_rate_data)[j][i];
				//}
				(*it_coalescent_rate_num)[group_i][group_j]   += (*it_coalescent_rate_data)[i][j];
				(*it_coalescent_rate_denom)[group_i][group_j] += (*it_coalescent_rate_data)[j][i];

				it_coalescent_rate_num++;
				it_coalescent_rate_denom++;
				it_coalescent_rate_data++;
			}    

		}
	}

	if(sample_ages.size() > 0){

		int proper_e = 0;
		//add the non-proper epochs to the earlier proper epoch.
		//however, if min sample age is > proper_epoch, add it to the min epoch
		for(int e = 1; e < num_epochs; e++){
			if(proper_epochs[e] == 1){
				proper_e = e;
			}else{
				for(int i = 0; i < sample.groups.size(); i++){
					for(int j = 0; j < sample.groups.size(); j++){
						coalescent_rate_num[proper_e][i][j]   += coalescent_rate_num[e][i][j];
						coalescent_rate_denom[proper_e][i][j] += coalescent_rate_denom[e][i][j];
					}
				}
			}
		}

		proper_e = 0;
		for(int e = 1; e < num_epochs; e++){
			if(proper_epochs[e] == 1){
				proper_e = e;
			}else{
				for(int i = 0; i < sample.groups.size(); i++){
					for(int j = 0; j < sample.groups.size(); j++){
						coalescent_rate_num[e][i][j]   = coalescent_rate_num[proper_e][i][j];
						coalescent_rate_denom[e][i][j] = coalescent_rate_denom[proper_e][i][j];
					}
				}
			}
		}

		for(int e = 0; e < num_epochs; e++){
			for(int i = 0; i < sample.groups.size(); i++){
				for(int j = 0; j < sample.groups.size(); j++){
					if(group_min_age[i] > e || group_min_age[j] > e){
						coalescent_rate_num[e][i][j]   = 0;
						coalescent_rate_denom[e][i][j] = 0;
					}
				}
			}
		}

	}

	for(int e = 0; e < num_epochs; e++){
		for(int i = 0; i < sample.groups.size(); i++){
			for(int j = 0; j < sample.groups.size(); j++){
				if(coalescent_rate_denom[e][i][j] == 0){
					coalescent_rate_num[e][i][j] = 0.0;
					coalescent_rate_denom[e][i][j] = 0.0;
				}
			}
		}
	}

	/////////////////// OUTPUT //////////////////////

	std::ofstream os(options["output"].as<std::string>() + ".coal");
	if(os.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	} 

	for(std::vector<std::string>::iterator it_groups = sample.groups.begin(); it_groups != sample.groups.end(); it_groups++){
		os << *it_groups << " ";
	}
	os << "\n";

	for(int e = 0; e < num_epochs; e++){
		os << epoch[e]  << " "; 
	}
	os << "\n";

	for(int i = 0; i < sample.groups.size(); i++){
		for(int j = 0; j < sample.groups.size(); j++){
			os << i << " " << j << " ";
			for(int e = 0; e < num_epochs; e++){
				if(i == j){
					double rate = coalescent_rate_num[e][i][j]/coalescent_rate_denom[e][i][j];
					os << rate << " ";
					//os << 2.0 * rate/((float) (sample.group_sizes[i] * (sample.group_sizes[j]-1.0))) << " ";
				}else if(i < j){
					double rate = coalescent_rate_num[e][i][j]/coalescent_rate_denom[e][i][j];
					os << rate << " ";
					//os << rate/((float) (sample.group_sizes[i] * sample.group_sizes[j])) << " ";
				}else{
					double rate = coalescent_rate_num[e][j][i]/coalescent_rate_denom[e][j][i];
					os << rate << " ";
					//os << rate/((float) (sample.group_sizes[j] * sample.group_sizes[i])) << " ";
				}
			}
			os << "\n"; 
		}
	}

	os.close();

	/////////////////////////////////////////////
	//Resource Usage

	rusage usage;
	getrusage(RUSAGE_SELF, &usage);

	std::cerr << "CPU Time spent: " << usage.ru_utime.tv_sec << "." << std::setfill('0') << std::setw(6) << usage.ru_utime.tv_usec << "s; Max Memory usage: " << usage.ru_maxrss/1000.0 << "Mb." << std::endl;
	std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

	return 0;


}

int FinalizePopulationSizeByHaplotype(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Estimate population size using coalescence rate." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Finalizing coalescence rate..." << std::endl;



	////////// read labels of sequences //////////

	/////////////////// INITIALIZE //////////////////////

	int num_epochs;
	std::vector<float> epoch;

	FILE* fp = fopen((options["output"].as<std::string>() + ".bin").c_str(),"rb");
	assert(fp != NULL);

	fread(&num_epochs, sizeof(int), 1, fp);
	epoch.resize(num_epochs);
	fread(&epoch[0], sizeof(float), num_epochs, fp);
	std::vector<CollapsedMatrix<float>> coalescent_rate_data(num_epochs);
	for(int e = 0; e < num_epochs; e++){
		coalescent_rate_data[e].ReadFromFile(fp);
	}
	fclose(fp);

	int N = coalescent_rate_data[0].size();

	std::vector<CollapsedMatrix<float>> coalescent_rate(num_epochs);
	for(std::vector<CollapsedMatrix<float>>::iterator it_c = coalescent_rate.begin(); it_c != coalescent_rate.end(); it_c++){
		(*it_c).resize(N, N);
		std::fill((*it_c).vbegin(), (*it_c).vend(), 0.0);
	}

	for(int i = 0; i < N; i++){
		for(int j = i+1; j < N; j++){

			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_data = coalescent_rate_data.begin();
			std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate      = coalescent_rate.begin();

			for(; it_coalescent_rate_data != std::prev(coalescent_rate_data.end(),1);){
				//i < j always holds 
				if((*it_coalescent_rate_data)[i][j] == 0.0){
					(*it_coalescent_rate)[i][j] += 0.0;
				}else{ 
					(*it_coalescent_rate)[i][j] += (*it_coalescent_rate_data)[i][j]/(*it_coalescent_rate_data)[j][i];
				}

				it_coalescent_rate++;
				it_coalescent_rate_data++;
			}    

		}
	}


	std::ofstream os(options["output"].as<std::string>() + ".coal");
	if(os.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	} 

	for(int i = 0; i < N; i++){
		os << i << " ";
	}
	os << "\n";

	for(int e = 0; e < num_epochs; e++){
		os << epoch[e]  << " "; 
	}
	os << "\n";

	for(int i = 0; i < N; i++){
		for(int j = i+1; j < N; j++){
			os << i << " " << j << " ";
			for(int e = 0; e < num_epochs; e++){
				os << coalescent_rate[e][i][j] << " ";
			}
			os << "\n"; 
		}
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

	return 0;


}

int FinalizeCoalescenceCount(cxxopts::Options& options){

	//////////////////////////////////
	//Program options

	bool help = false;
	if(!options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: output." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Estimate population size using coalescence rate." << std::endl;
		exit(0);
	}  

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Finalizing coalescence count..." << std::endl;



	////////// read labels of sequences //////////

	/////////////////// INITIALIZE //////////////////////

	int num_epochs;
	std::vector<float> epoch;

	FILE* fp = fopen((options["output"].as<std::string>() + ".bin").c_str(),"rb");
	assert(fp != NULL);

	fread(&num_epochs, sizeof(int), 1, fp);
	epoch.resize(num_epochs);
	fread(&epoch[0], sizeof(float), num_epochs, fp);
	std::vector<CollapsedMatrix<float>> coalescent_rate_data(num_epochs);
	for(int e = 0; e < num_epochs; e++){
		coalescent_rate_data[e].ReadFromFile(fp);
	}
	fclose(fp);

	int N = coalescent_rate_data[0].size();

	std::vector<CollapsedMatrix<float>> coalescent_rate(num_epochs);
	for(std::vector<CollapsedMatrix<float>>::iterator it_c = coalescent_rate.begin(); it_c != coalescent_rate.end(); it_c++){
		(*it_c).resize(N, N);
		std::fill((*it_c).vbegin(), (*it_c).vend(), 0.0);
	}

	std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate_data = coalescent_rate_data.begin();
	std::vector<CollapsedMatrix<float>>::iterator it_coalescent_rate      = coalescent_rate.begin();

	int tree_index = 0;
	int snp = 0;
	int block_size = 1e6;
	int chr = 1;
	Mutations mut;
	mut.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");

	for(; it_coalescent_rate_data != std::prev(coalescent_rate_data.end(),1);){

		//TODO: fix 
		//get proportion of 1Mv that each tree is persisting for
		float prop = 0.0;
		while(mut.info[snp].tree == tree_index){
			prop += mut.info[snp].dist;
			snp++;
			if(snp == mut.info.size()) break;
		}
		prop /= block_size;
		if(prop > 1) std::cerr << prop << std::endl;

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){ 
				(*it_coalescent_rate)[i][j] += (*it_coalescent_rate_data)[i][j] * prop;
			}
		}
		tree_index++;
		it_coalescent_rate++;
		it_coalescent_rate_data++;

		if(chr <= 22 && mut.info.size()== snp){
			chr++;
			snp = 0;
			tree_index = 0;
			mut.Read(options["input"].as<std::string>() + "_chr" + std::to_string(chr) + ".mut");
		}
	}    


	std::ofstream os(options["output"].as<std::string>() + ".coal");
	if(os.fail()){
		std::cerr << "Error while opening file." << std::endl;
		exit(1);
	} 

	for(int i = 0; i < N; i++){
		os << i << " ";
	}
	os << "\n";

	for(int e = 0; e < num_epochs; e++){
		os << epoch[e]  << " "; 
	}
	os << "\n";

	for(int i = 0; i < N; i++){
		for(int j = i+1; j < N; j++){
			os << i << " " << j << " ";
			for(int e = 0; e < num_epochs; e++){
				os << coalescent_rate[e][i][j] << " ";
			}
			os << "\n"; 
		}
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

	return 0;


}

//new Finalise function

