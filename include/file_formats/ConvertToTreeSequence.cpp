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

  for(int i = 0; i < data.N; i++){
    os_node_table << "1\t" << i << "\t0.0\n";
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
 
  int pos, snp, pos_end, snp_end, tree_count = 0, node, node_const, site_count = 0;
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  while(num_bases_tree_persists >= 0.0){

    mtr.tree.GetCoordinates(coordinates);
    pos = mut.info[mtr.pos].pos;
    if(mtr.pos == 0) pos = 0;
    snp = mtr.pos;

    tree_count = mut.info[snp].tree;
    node_const = tree_count * (data.N - 1);

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
  AncMutIterators ancmut(options["input"].as<std::string>() + ".anc", options["input"].as<std::string>() + ".mut");
  int N = ancmut.NumTips();
  int L = ancmut.NumSnps();
  MarginalTree mtr;
  Muts::iterator it_mut;
  float num_bases_tree_persists = 0.0;
  Data data(N,L);
  int root = 2*data.N - 2;

  Mutations mut;
  mut.Read(options["input"].as<std::string>() + ".mut");

  //////////////////////////////////////////////////////////////////////////////////////////////

  int ret;
  tsk_table_collection_t tables;
  ret = tsk_table_collection_init(&tables, 0);
  check_tsk_error(ret);

  tables.sequence_length = mut.info[data.L-1].pos + 1;
  //tables.individuals.num_rows = data.N; //individuals table
  for(int i = 0; i < data.N; i++){
    tsk_individual_table_add_row(&tables.individuals, 0, NULL, 0 , NULL, 0);
  }

  //population table

  //sites table
  char ancestral_allele[1];
  //tsk_site_table_add_row(&tables.sites, 1, ancestral_allele, sizeof(ancestral_allele), NULL, 0);
  for(std::vector<SNPInfo>::iterator it_mut = mut.info.begin(); it_mut != mut.info.end(); it_mut++){
    ancestral_allele[0] = (*it_mut).mutation_type[0];
    ret = tsk_site_table_add_row(&tables.sites, (*it_mut).pos, ancestral_allele, 1, NULL, 0);
    check_tsk_error(ret);
  }

  for(int i = 0; i < data.N; i++){
    ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, 0, TSK_NULL, i, NULL, 0);
    check_tsk_error(ret);
  }

  ///////////////////////////////////////////////////////////////////////////// 

  std::vector<float> coordinates(2*data.N-1,0.0);
  int pos, snp, pos_end, snp_end, tree_count = 0, node, node_const, site_count = 0;

	int count = 0;
  char derived_allele[1];
  num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);
  while(num_bases_tree_persists >= 0.0){

    mtr.tree.GetCoordinates(coordinates);

    pos = mut.info[mtr.pos].pos;
    if(mtr.pos == 0) pos = 0;
    snp = mtr.pos;

    tree_count = mut.info[snp].tree;
    node_const = tree_count * (data.N - 1);

    //Mutation table
    int l = snp;
    while(mut.info[l].tree == tree_count){
      if(mut.info[l].branch.size() == 1){
        node = (*mut.info[l].branch.begin());
        if(node < data.N){
          derived_allele[0] = mut.info[l].mutation_type[2];
          ret = tsk_mutation_table_add_row(&tables.mutations, l, node, TSK_NULL, derived_allele, 1, NULL, 0);
          check_tsk_error(ret);
        }else{
          derived_allele[0] = mut.info[l].mutation_type[2];
          ret = tsk_mutation_table_add_row(&tables.mutations, l, node + node_const, TSK_NULL, derived_allele, 1, NULL, 0);
          check_tsk_error(ret);
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
    int n = N;
    for(std::vector<float>::iterator it_coords = std::next(coordinates.begin(), data.N); it_coords != coordinates.end(); it_coords++){   
      ret = tsk_node_table_add_row(&tables.nodes, 0, *it_coords, TSK_NULL, TSK_NULL, NULL, 0);   
      //ret = tsk_node_table_add_row(&tables.nodes, TSK_NODE_IS_SAMPLE, leaves[n].num_leaves, TSK_NULL, TSK_NULL, NULL, 0);    
      check_tsk_error(ret);
      n++;
    }

    //Edge table
    for(it_node = mtr.tree.nodes.begin(); it_node != std::prev(mtr.tree.nodes.end(),1); it_node++){
      node = (*it_node).label;
      if(node >= data.N) node += node_const;
      ret = tsk_edge_table_add_row(&tables.edges, pos, pos_end, (*(*it_node).parent).label + node_const, node);    
      check_tsk_error(ret);
    }

    num_bases_tree_persists = ancmut.NextTree(mtr, it_mut);

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

