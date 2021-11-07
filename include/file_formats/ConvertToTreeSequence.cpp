#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>
#include <bitset>
#include <err.h>

#include "gzstream.hpp"
#include "collapsed_matrix.hpp"
#include "data.hpp"
#include "sample.hpp"
#include "mutations.hpp"
#include "cxxopts.hpp"
#include "tskit.h"

#define check_tsk_error(val) if (val < 0) {\
  errx(EXIT_FAILURE, "line %d: %s", __LINE__, tsk_strerror(val));\
}

void
ConvertToTreeSequenceTxt(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("output") ){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Converts anc/mut file format (Relate) to tree sequence file format (tskit) as txt files.." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Converting anc/mut file format (Relate) to tree sequence file format (tskit) as txt files.." << std::endl;

  ////////////////////////
  //read in anc file
  AncMutIterators ancmut(options["input"].as<std::string>() + ".anc", options["mut"].as<std::string>() + ".mut");
  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0;
  Data data(N,L);
  int root = 2*data.N - 2;

  Mutations mut;
  mut.Read(options["input"].as<std::string>() + ".mut");

  //Individuals table
  std::ofstream os_indiv_table(options["output"].as<std::string>() + ".indiv_table");
  if(os_indiv_table.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".indiv_table" << std::endl;
  }
  os_indiv_table << "flags\tlocation\n";
  for(int i = 0; i < data.N; i++){
    os_indiv_table << "0\t0.0,0.0\n";
  }
  os_indiv_table.close();

  //Site table
  std::ofstream os_site_table(options["output"].as<std::string>() + ".site_table");
  if(os_site_table.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".site_table" << std::endl;
  }
  os_site_table << "position\tancestral_state\n";
  for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
    os_site_table << (*it_mut).pos << "\t" << (*it_mut).mutation_type[0] << std::endl;
  }
  os_site_table.close();

  //Population table

  std::ofstream os_population_table(options["output"].as<std::string>() + ".population_table");
  if(os_population_table.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".population_table" << std::endl;
  }
  os_population_table << "id\tmetadata\n";
  for(int i = 0; i < data.N; i++){
    os_population_table << "0\t\n";
  }
  os_population_table.close();

  ////////////////////////////////////////////
  //The remaining tables actually concern the genealogy

  //Node table

  std::ofstream os_node_table(options["output"].as<std::string>() + ".node_table");
  if(os_node_table.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".node_table" << std::endl;
  }
  os_node_table << "is_sample\tindividual\ttime\tmetadata" << std::endl;

  if(ancmut.sample_ages.size() > 0){
    for(int i = 0; i < data.N; i++){
      os_node_table << "1\t" << i << "\t" << ancmut.sample_ages[i] << "\n";
    }
  }else{
    for(int i = 0; i < data.N; i++){
      os_node_table << "1\t" << i << "\t0.0\n";
    }
  }

  //Edge table

  std::ofstream os_edge_table(options["output"].as<std::string>() + ".edge_table");
  if(os_edge_table.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".edge_table" << std::endl;
  }
  os_edge_table << "left\tright\tparent\tchild" << std::endl;

  //Mutation table

  std::ofstream os_mut_table(options["output"].as<std::string>() + ".mut_table");
  if(os_mut_table.fail()){
    std::cerr << "Error while opening file " << options["output"].as<std::string>() + ".mut_table" << std::endl;
  }
  os_mut_table << "site\tnode\tderived_state" << std::endl;

  //If I treat trees separately, I can need to change labels of nodes, record edges, and a new mutation in the sites file

  std::vector<float> coordinates(2*data.N-1,0.0);
 
  int pos, snp, pos_end, snp_end, tree_count = 0, node, node_const, site_count = 0, count = 0;
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  while(num_bases_tree_persists >= 0.0){

    mtr.tree.GetCoordinates(coordinates);

    for(int i = 0; i < mtr.tree.nodes.size()-1; i++){
      if(!(coordinates[(*mtr.tree.nodes[i].parent).label] - coordinates[i] > 0.0)){
        int parent = (*mtr.tree.nodes[i].parent).label, child = i;
        while(coordinates[parent] - coordinates[child] < 1e-5){
          coordinates[parent] = coordinates[child] + 1e-5;
          if(parent == root) break;
          child  = parent;
          parent = (*mtr.tree.nodes[parent].parent).label;
        }
      }
    }


    pos = mut.info[mtr.pos].pos;
    if(mtr.pos == 0) pos = 0;
    snp = mtr.pos;

    tree_count = mut.info[snp].tree;
    node_const = count * (data.N - 1);

    //Mutation table
    int l = snp;
    while(mut.info[l].tree == tree_count){
      if(mut.info[l].branch.size() == 1){
        node = (*mut.info[l].branch.begin());
        if(node < data.N){
          os_mut_table << l << "\t" << node << "\t" << mut.info[l].mutation_type[2] << "\n"; //<< "\t"<< (*(mtr.tree.nodes[*mut.info[l].branch.begin()]).parent).label + node_const << "\n";
        }else{
          os_mut_table << l << "\t" << node + node_const << "\t" << mut.info[l].mutation_type[2] << "\n"; //<< "\t" << (*(mtr.tree.nodes[*mut.info[l].branch.begin()]).parent).label + node_const << "\n"; 
        }
        site_count++;
      }

      l++; 
      if(l == data.L) break;
    }
    snp_end = l;
    if(snp_end < data.L){
      pos_end = mut.info[snp_end].pos;
    }else{
      pos_end = mut.info[data.L-1].pos + 1;
    }

    //Node table
    std::vector<Node>::iterator it_node = std::next(mtr.tree.nodes.begin(), data.N);
    for(std::vector<float>::iterator it_coords = std::next(coordinates.begin(), data.N); it_coords != coordinates.end(); it_coords++){
      //(*(*it_node)child_left).SNP_begin, (*(*it_node)child_left).SNP_end, (*(*it_node)child_right).SNP_begin, (*(*it_node)child_right).SNP_end  
      os_node_table << "0\t-1\t" << *it_coords << "\n";
    }

    //Edge table
    for(it_node = mtr.tree.nodes.begin(); it_node != std::prev(mtr.tree.nodes.end(),1); it_node++){
      node = (*it_node).label;
      if(node >= data.N) node += node_const;
      os_edge_table << pos << "\t" << pos_end << "\t" << (*(*it_node).parent).label + node_const << "\t" << node << "\n";
    }

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
    count++;

  } 
  os_mut_table.close();
  os_node_table.close();
  os_edge_table.close();

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
ConvertToTreeSequence(cxxopts::Options& options){

  //////////////////////////////////
  //Program options

  bool help = false;
  if( !options.count("input") || !options.count("output") ){
    std::cout << "Not enough arguments supplied." << std::endl;
    std::cout << "Needed: input, output." << std::endl;
    help = true;
  }
  if(options.count("help") || help){
    std::cout << options.help({""}) << std::endl;
    std::cout << "Converts anc/mut file format (Relate) to tree sequence file format (tskit).." << std::endl;
    exit(0);
  }  

  std::cerr << "---------------------------------------------------------" << std::endl;
  std::cerr << "Converting anc/mut file format (Relate) to tree sequence file format (tskit).." << std::endl;

 
	////////////////////////
	//read in anc file

	MarginalTree mtr, prev_mtr; //stores marginal trees. mtr.pos is SNP position at which tree starts, mtr.tree stores the tree
	Muts::iterator it_mut, it_mut_first, it_mut_tmp; //iterator for mut file
	float num_bases_tree_persists = 0.0;

	////////// 1. Read one tree at a time /////////

	//We open anc file and read one tree at a time. File remains open until all trees have been read OR ancmut.CloseFiles() is called.
	//The mut file is read once, file is closed after constructor is called.
	AncMutIterators ancmut(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");

	num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
	it_mut_first = it_mut;
	int N = (mtr.tree.nodes.size() + 1)/2.0, root = 2*N - 2, L = ancmut.NumSnps();
	Data data(N,L);
	std::vector<float> coordinates(2*data.N-1,0.0);

	Mutations mut;
	mut.Read(options["input"].as<std::string>() + ".mut");
	
	//........................................................................
	//Populate ts tables

	int ret;
	tsk_table_collection_t tables;
	ret = tsk_table_collection_init(&tables, 0);
	check_tsk_error(ret);

	tables.sequence_length = (*std::prev(ancmut.mut_end(),1)).pos + 1;
	for(int i = 0; i < N; i++){
		tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0 , NULL, 0);
	}




	//population table

	//sites table
	char ancestral_allele[1];
	double pos, pos_begin, pos_end;
	std::vector<double> bps(L);
	int bps_index = 0;
	//tsk_site_table_add_row(&tables.sites, 1, ancestral_allele, sizeof(ancestral_allele), NULL, 0);
	for(; it_mut != ancmut.mut_end();){
		ancestral_allele[0] = (*it_mut).mutation_type[0];
		pos = (*it_mut).pos;
		int count = 0;

		it_mut_tmp = it_mut;
		while((*it_mut_tmp).pos == pos){
			it_mut_tmp++;
			count++;
			if(it_mut_tmp == ancmut.mut_end()) break;
		}
		assert(count > 0);

		if(count == 1){
			ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
			bps[bps_index] = (*it_mut).pos;
			bps_index++;
			it_mut++;
		}else{

			if(it_mut_tmp != ancmut.mut_end()){
				pos_end = ((*it_mut_tmp).pos + (*std::prev(it_mut_tmp)).pos)/2.0;
			}else{
				pos_end = (*std::prev(it_mut_tmp)).pos;
			}
			it_mut_tmp = it_mut;
			if(it_mut_tmp != it_mut_first){
				pos_begin = ((*it_mut_tmp).pos + (*std::prev(it_mut_tmp)).pos)/2.0;
			}else{
				pos_begin = pos;
			}
			int i = 0;
			while((*it_mut_tmp).pos == pos){
				ret = tsk_site_table_add_row(&tables.sites, ((i+1.0)/(count+1.0))*(pos_end - pos_begin) + pos_begin, ancestral_allele, 1, NULL, 0);
				bps[bps_index] = ((i+1.0)/(count+1.0))*(pos_end - pos_begin) + pos_begin;
				bps_index++;
				it_mut_tmp++;
				i++;
				if(it_mut_tmp == ancmut.mut_end()) break;
			}
			it_mut = it_mut_tmp;

		}

		//std::cerr << (*it_mut).pos << " " << count << " " << bps_index-1 << " " << bps[bps_index-1] << std::endl;
		check_tsk_error(ret);
	}
	assert(bps_index == L);

	if(ancmut.sample_ages.size() > 0){
		for(int i = 0; i < data.N; i++){
			ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, ancmut.sample_ages[i], TSK_NULL, i, NULL, 0);
			check_tsk_error(ret);
		}
	}else{
		for(int i = 0; i < data.N; i++){
			ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, i, NULL, 0);
			check_tsk_error(ret);
		}
	}

	///////////////////////////////////////////////////////////////////////////// 

	it_mut = it_mut_first; 
	int snp, snp_end, tree_count = 0, node, node_const, site_count = 0, count = 0;
	int node_count = data.N, edge_count = 0;

	char derived_allele[1];
	while(num_bases_tree_persists >= 0.0){

		mtr.tree.GetCoordinates(coordinates);
		for(int i = 0; i < mtr.tree.nodes.size()-1; i++){
			if(!(coordinates[(*mtr.tree.nodes[i].parent).label] - coordinates[i] > 0.0)){
				int parent = (*mtr.tree.nodes[i].parent).label, child = i;
				while(coordinates[parent] <= coordinates[child] + std::nextafter(coordinates[child], coordinates[child] + 1)){
					coordinates[parent] = coordinates[child] + std::nextafter(coordinates[child], coordinates[child] + 1);
					if(parent == root) break;
					child  = parent;
					parent = (*mtr.tree.nodes[parent].parent).label;
				}
			}
		}

		for(int i = 0; i < mtr.tree.nodes.size()-1; i++){	
			assert(coordinates[i] < coordinates[(*mtr.tree.nodes[i].parent).label]);
		}

		snp = mtr.pos;
		if(snp == 0){
			pos = 0;
		}else{
			pos = (bps[snp] + bps[snp-1])/2.0;
		}

		tree_count = (*it_mut).tree;
		node_const = tree_count * (data.N - 1);

		//Mutation table
		int l = snp;
		while((*it_mut).tree == tree_count){
			if((*it_mut).branch.size() == 1){
				node = *(*it_mut).branch.begin();
				if(node < N){
					derived_allele[0] = (*it_mut).mutation_type[2];
					ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}else{
					derived_allele[0] = (*it_mut).mutation_type[2];
					ret = tsk_mutation_table_add_row(&tables.mutations, l, node + node_const, TSK_NULL, derived_allele, 1, NULL, 0);
					check_tsk_error(ret);
				}
				site_count++;
			}

			l++;
			it_mut++; 
			if(l == L) break;
		}

		snp_end = l;
		if(snp_end < L){
			pos_end = (bps[snp_end-1] + bps[snp_end])/2.0;
		}else{
			pos_end = bps[L-1] + 1;
		}

		assert(pos != pos_end);
		assert(pos <= bps[snp]);
		assert(pos_end >= bps[snp]);

		//Node table
		std::vector<Node>::iterator it_node = std::next(mtr.tree.nodes.begin(), data.N);
		int n = N;
		for(std::vector<float>::iterator it_coords = std::next(coordinates.begin(), data.N); it_coords != coordinates.end(); it_coords++){   
			ret = tsk_node_table_add_row(&tables.nodes, 0, *it_coords, TSK_NULL, TSK_NULL, NULL, 0);   
			check_tsk_error(ret);
			n++;
			node_count++;
		}

		//Edge table
		for(it_node = mtr.tree.nodes.begin(); it_node != std::prev(mtr.tree.nodes.end(),1); it_node++){
			node = (*it_node).label;
			if(node >= data.N) node += node_const;
			ret = tsk_edge_table_add_row(&tables.edges, pos, pos_end, (*(*it_node).parent).label + node_const, node);    
			check_tsk_error(ret);
			edge_count++;
		}

		num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
		count++;

	} 

	tsk_table_collection_sort(&tables, NULL, 0);
	check_tsk_error(ret);

  //////////////////////////

  // Write out the tree sequence
  ret = tsk_table_collection_dump(&tables, (options["output"].as<std::string>() + ".trees").c_str(), 0);        
  check_tsk_error(ret);
  tsk_table_collection_free(&tables); 

  /*
  int iter;
  tsk_treeseq_t ts;
  tsk_tree_t tree;

  ret = tsk_treeseq_load(&ts, (options["output"].as<std::string>() + ".trees").c_str(), 0);
  check_tsk_error(ret);
  ret = tsk_tree_init(&tree, &ts, 0);
  check_tsk_error(ret);

  printf("Iterate forwards\n");
  FILE* fp;
  fp = fopen("state.txt", "w");
  for (iter = tsk_tree_first(&tree); iter == 1; iter = tsk_tree_next(&tree)) {
    tsk_tree_print_state(&tree, fp);
    printf("\ttree %d has %d roots\n", tsk_tree_get_index(&tree), tsk_tree_get_num_roots(&tree));
  }
  fclose(fp);
  check_tsk_error(iter);
  */

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

