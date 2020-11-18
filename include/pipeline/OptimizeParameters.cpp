/*! \file OptimizeParameters.cpp 
 *  \brief Optimize parameters for painting
 *
 *  Input: data
 *  Output: optimal parameters for painting (theta and recombination_factor)
 */

#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>

#include "cxxopts.hpp"
#include "data.hpp"
#include "anc.hpp"
#include "fast_painting.hpp"
#include "anc_builder.hpp"

#include "MakeChunks.cpp"
#include "Paint.cpp"
#include "Clean.cpp"

int OptimizeParameters(cxxopts::Options& options){

	bool help = false;
	if(!options.count("haps") || !options.count("sample") || !options.count("map") || !options.count("output")){
		std::cout << "Not enough arguments supplied." << std::endl;
		std::cout << "Needed: haps, sample, map, output. Optional: dist." << std::endl;
		help = true;
	}
	if(options.count("help") || help){
		std::cout << options.help({""}) << std::endl;
		std::cout << "Use to make smaller chunks from the data." << std::endl;
		exit(0);
	}

	std::cerr << "############" << std::endl;
	std::cerr << "Optimizing Parameters..." << std::endl;

	int N, L;
	double memory_size;
	int start_chunk, end_chunk;
	std::string file_out = options["output"].as<std::string>() + "/";

	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Using:" << std::endl;
	std::cerr << "  " << options["haps"].as<std::string>() << std::endl;    
	std::cerr << "  " << options["sample"].as<std::string>() << std::endl;
	std::cerr << "  " << options["map"].as<std::string>() << std::endl;
	if(!options.count("coal")){
		std::cerr << "with mu = " << options["mutation_rate"].as<float>() << " and 2Ne = " << options["effectiveN"].as<float>() << "." << std::endl;
	}else{
		std::cerr << "with mu = " << options["mutation_rate"].as<float>() << " and coal = " << options["coal"].as<std::string>() << "." << std::endl;
	}

	std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

	MakeChunks(options);

	FILE* fp = fopen((file_out + "parameters.bin").c_str(), "r");
	assert(fp != NULL);
	fread(&N, sizeof(int), 1, fp);
	fread(&L, sizeof(int), 1, fp);
	fread(&end_chunk, sizeof(int), 1, fp);
	fread(&memory_size, sizeof(double), 1, fp);
	fclose(fp);

	end_chunk--;
	start_chunk = 0;

	Data data(N,L);
	std::cerr << "---------------------------------------------------------" << std::endl;
	std::cerr << "Read " << data.N << " haplotypes with " << data.L << " SNPs per haplotype." << std::endl;
	std::cerr << "Expected minimum memory usage: " << memory_size << "Gb." << std::endl;
	std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

	std::vector<float> theta;//      = {1e-4, 1e-3, 1e-2, 1e-1};
	std::vector<float> rec_factor;// = {0.001, 0.1, 1, 10, 100};

	float val;
	std::string line;
  std::ifstream is(options["input"].as<std::string>());
  
	getline(is, line);
	std::istringstream itheta;
	itheta.str(line);
  while(itheta >> val){
		if(val >= 1.0 || val <= 0){
      std::cerr << "Error: theta value has to be in (0,1)" << std::endl;
			exit(1);
		}
    theta.push_back(val);
	}

	getline(is, line);
	std::istringstream irec;
	irec.str(line);
	while(irec >> val){
		if(val <= 0){
      std::cerr << "Error: rho value has to be positive" << std::endl;
			exit(1);
		}
		rec_factor.push_back(val);
	}

	is.close();

  std::vector<std::vector<int>> num_notmapping(theta.size());
	for(int i = 0; i < num_notmapping.size(); i++){
    num_notmapping[i].resize(rec_factor.size());
		std::fill(num_notmapping[i].begin(), num_notmapping[i].end(), 0);
	}

	for(int c = start_chunk; c <= end_chunk; c++){

		std::cerr << "---------------------------------------------------------" << std::endl;
		std::cerr << "Starting chunk " << c << " of " << end_chunk << "." << std::endl;
		std::cerr << "---------------------------------------------------------" << std::endl << std::endl;

		std::string file_out = options["output"].as<std::string>() + "/";

		int N, L, num_windows;
		std::vector<int> window_boundaries;
		FILE* fp = fopen((file_out + "parameters_c" + std::to_string(c) + ".bin").c_str(), "r");
		assert(fp != NULL);
		fread(&N, sizeof(int), 1, fp);
		fread(&L, sizeof(int), 1, fp);
		fread(&num_windows, sizeof(int), 1, fp);
		window_boundaries.resize(num_windows);
		fread(&window_boundaries[0], sizeof(int), num_windows, fp);
		fclose(fp);
		num_windows--;

		Data data((file_out + "chunk_" + std::to_string(c) + ".hap").c_str(), (file_out + "chunk_" + std::to_string(c) + ".bp").c_str(), (file_out + "chunk_" + std::to_string(c) + ".r").c_str(), (file_out + "chunk_" + std::to_string(c) + ".rpos").c_str()); //struct data is defined in data.hpp
		data.name = (file_out + "chunk_" + std::to_string(c) + "/paint/relate");
		const std::string dirname = file_out + "chunk_" + std::to_string(c) + "/";

		std::vector<double> rec_rate = data.r;
		for(int theta_index = 0; theta_index < (int) theta.size(); theta_index++){
			for(int rec_index = 0; rec_index < (int) rec_factor.size(); rec_index++){

				float mean_rec = 0.0;
				data.theta     = theta[theta_index];
				data.ntheta    = 1.0 - data.theta;

				for(int l = 0; l < (int)data.r.size(); l++){
					data.r[l] = rec_rate[l] * rec_factor[rec_index];
					mean_rec  += data.r[l];
				}

				Paint(options,c);
				for(int section = 0; section < num_windows; section++){

					AncesTree anc;
					AncesTreeBuilder ancbuilder(data);

					int section_startpos = window_boundaries[section];
					int section_endpos   = window_boundaries[section+1]-1;
					if(section_endpos >= data.L) section_endpos = data.L-1;
			
					int seed = c + section + std::time(0) + getpid();
					bool ancestral_state = true;
					//need to sort out data
					num_notmapping[theta_index][rec_index] += ancbuilder.OptimizeParameters(section, section_startpos, section_endpos, data, seed);

				}
			}

		}

	}

	Clean(options);

  std::ofstream os(options["output"].as<std::string>() + ".opt");
	for(int theta_index = 0; theta_index < (int) theta.size(); theta_index++){
		for(int rec_index = 0; rec_index < (int) rec_factor.size(); rec_index++){
      os << theta[theta_index] << " " << rec_factor[rec_index] << " " << num_notmapping[theta_index][rec_index] << std::endl;
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
