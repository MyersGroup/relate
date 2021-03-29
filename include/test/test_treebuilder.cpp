#include "catch_amalgamated.hpp"
#include <algorithm>

#include "data.hpp"
#include "anc.hpp"
#include "tree_builder.hpp"


TEST_CASE( "Testing tree builder" ){

  SECTION("Testing tree builder with good distance matrix");

  int N = 5;
  int L = 1;
  Data data(N,L); //struct data is defined in data.hpp
  data.theta = 0.025;

  Tree tree;
  MinMatch tb(data);

  CollapsedMatrix<float> d;
  d.resize(N,N);

  //build tree with empty distance matrix

  for(int n = 0; n < N; n++){
    for(int m = 0; m < N; m++){
      d[n][m] = 0.0;
    }
  }
  std::vector<double> sample_ages(N,0);
  tb.QuickBuild(d,tree, sample_ages);

  //built tree in case where I know the correct tree

  d[0][0] = 0; 
  d[0][1] = 0;
  d[0][2] = 1;
  d[0][3] = 2;
  d[0][4] = 2;
  d[1][0] = 2; 
  d[1][1] = 0;
  d[1][2] = 3;
  d[1][3] = 4;
  d[1][4] = 4;
  d[2][0] = 0; 
  d[2][1] = 0;
  d[2][2] = 0;
  d[2][3] = 1;
  d[2][4] = 1;
  d[3][0] = 1; 
  d[3][1] = 1;
  d[3][2] = 1;
  d[3][3] = 0;
  d[3][4] = 0;
  d[4][0] = 1; 
  d[4][1] = 1;
  d[4][2] = 1;
  d[4][3] = 0;
  d[4][4] = 0;

  //build tree

  tb.QuickBuild(d,tree,sample_ages);

  //compare to true tree

  Tree true_tree;
  true_tree.nodes.resize(2*N-1);

  true_tree.nodes[0].parent = &true_tree.nodes[6]; 
  true_tree.nodes[1].parent = &true_tree.nodes[6]; 
  true_tree.nodes[2].parent = &true_tree.nodes[7]; 
  true_tree.nodes[3].parent = &true_tree.nodes[5]; 
  true_tree.nodes[4].parent = &true_tree.nodes[5]; 
  true_tree.nodes[5].parent = &true_tree.nodes[8]; 
  true_tree.nodes[6].parent = &true_tree.nodes[7]; 
  true_tree.nodes[7].parent = &true_tree.nodes[8]; 

  for(int n = 0; n < 2*N-1; n++){
    true_tree.nodes[n].label = n;
  }

  for(int n = 0; n < 2*N-2; n++){
    REQUIRE( (*tree.nodes[n].parent).label == (*true_tree.nodes[n].parent).label );
  }

  SECTION("Testing tree builder with bad distance matrix");

  N = 4;
  Data data2(N,L); //struct data is defined in data.hpp
  data2.theta = 0.025;
  
  MinMatch tb2(data2);

  CollapsedMatrix<float> d2;
  d2.resize(N,N);

  //built tree in case where I know the correct tree

  d2[0][0] = 0; 
  d2[0][1] = 1;
  d2[0][2] = 2;
  d2[0][3] = 2;
  d2[1][0] = 3; 
  d2[1][1] = 0;
  d2[1][2] = 1;
  d2[1][3] = 1;
  d2[2][0] = 0; 
  d2[2][1] = 1;
  d2[2][2] = 0;
  d2[2][3] = 1;
  d2[3][0] = 1; 
  d2[3][1] = 1;
  d2[3][2] = 0;
  d2[3][3] = 0;

  //build tree
  tb2.QuickBuild(d2,tree, sample_ages);

  //compare to true tree
  true_tree.nodes.resize(2*N-1);

  true_tree.nodes[0].parent = &true_tree.nodes[6]; 
  true_tree.nodes[1].parent = &true_tree.nodes[5]; 
  true_tree.nodes[2].parent = &true_tree.nodes[4]; 
  true_tree.nodes[3].parent = &true_tree.nodes[4]; 
  true_tree.nodes[4].parent = &true_tree.nodes[5]; 
  true_tree.nodes[5].parent = &true_tree.nodes[6]; 

  for(int n = 0; n < 2*N-1; n++){
    true_tree.nodes[n].label = n;
  }

  for(int n = 0; n < 2*N-2; n++){
    REQUIRE( (*tree.nodes[n].parent).label == (*true_tree.nodes[n].parent).label );
  }

}

void
GetCoordinates(Node& n, std::vector<float>& coordinates){

  if(n.child_left != NULL){
    GetCoordinates(*n.child_left, coordinates);
    GetCoordinates(*n.child_right, coordinates);

    coordinates[n.label] = coordinates[(*n.child_left).label] + (*n.child_left).branch_length;

  }else{
    coordinates[n.label] = 0.0;
  }

}

TEST_CASE( "Testing inference of branch lengths" ){

  SECTION( "Testing inference of branch lengths" );

  //first test that branch lengths are determined by prior if no information is given

  int N = 5;
  int L = 1;
  Data data(N,L); //struct data is defined in data.hpp
  data.theta = 0.025;
  data.pos.resize(L);
  data.rpos.resize(L);
  data.pos[0] = 1;
  data.rpos[0] = 0;

  Tree tree;
  MinMatch tb(data);

  CollapsedMatrix<float> d;
  d.resize(N,N);

  d[0][0] = 0; 
  d[0][1] = 1;
  d[0][2] = 2;
  d[0][3] = 3;
  d[0][4] = 4;
  d[1][0] = 0; 
  d[1][1] = 0;
  d[1][2] = 1;
  d[1][3] = 2;
  d[1][4] = 3;
  d[2][0] = 0; 
  d[2][1] = 0;
  d[2][2] = 0;
  d[2][3] = 1;
  d[2][4] = 2;
  d[3][0] = 0; 
  d[3][1] = 0;
  d[3][2] = 0;
  d[3][3] = 0;
  d[3][4] = 1;
  d[4][0] = 0; 
  d[4][1] = 0;
  d[4][2] = 0;
  d[4][3] = 0;
  d[4][4] = 0;

  // build tree  
  std::vector<double> sample_ages(N,0);
  tb.QuickBuild(d,tree, sample_ages);

  InferBranchLengths bl(data);
  bl.MCMC(data,tree); //tree is an empty tree, so branch lengths should be given by prior

  //Sort branches
  std::vector<float> coordinates(2*N-1); 
  GetCoordinates(tree.nodes[2*N-2], coordinates); 
  std::sort(coordinates.begin(), coordinates.end());

  for(int n = N; n < 2*N-1; n++){
    int num_lineages = 2*N - n;
    //std::cerr << std::fabs(coordinates[n] - coordinates[n-1]) << " " <<  2.0/(num_lineages * (num_lineages - 1.0)) * data.Ne  << std::endl;
    //REQUIRE( std::fabs(coordinates[n] - coordinates[n-1] - 2.0/(num_lineages * (num_lineages - 1.0)) * data.Ne) < 1e-3 );
  }
  //std::cerr << coordinates[2*N-2] << std::endl;

  //next test the order of branch lengths when there is information

  tree.nodes[3].num_events = 1;
  tree.nodes[4].num_events = 2;
  tree.nodes[6].num_events = 3;

  bl.EM(data, tree,true);

  std::vector<int> true_order = {0,1,2,3,4,6,5,7,8}, inferred_order(2*N-1);
  GetCoordinates(tree.nodes[2*N-2], coordinates); 

  std::size_t n(0);
  std::generate(std::begin(inferred_order) + N, std::end(inferred_order), [&]{ return n++; });
  std::sort(std::begin(inferred_order) + N, std::end(inferred_order), [&](int i1, int i2) { return coordinates[i1 + N] < coordinates[i2 + N]; } );

  for(int n = N; n < N; n++){
    REQUIRE( true_order[n] == inferred_order[n] );
  }

}

TEST_CASE( "Testing MCMC of branch lengths" ){

  SECTION( "Testing MCMC of branch lengths" );

  //first test that branch lengths are determined by prior if no information is given

  int N = 5;
  int L = 1;
  Data data(N,L); //struct data is defined in data.hpp
  data.theta = 0.025;
  data.pos.resize(L);
  data.rpos.resize(L);
  data.pos[0] = 0;
  data.rpos[0] = 0;


  Tree tree;
  MinMatch tb(data);
  InferBranchLengths bl(data);

  CollapsedMatrix<float> d;
  d.resize(N,N);

  //do MCMC on empty tree
  //tb.QuickBuild(d,tree);
  //bl.MCMC(tree);

  d[0][0] = 0; 
  d[0][1] = 0;
  d[0][2] = 1;
  d[0][3] = 2;
  d[0][4] = 2;
  d[1][0] = 2; 
  d[1][1] = 0;
  d[1][2] = 3;
  d[1][3] = 4;
  d[1][4] = 4;
  d[2][0] = 0; 
  d[2][1] = 0;
  d[2][2] = 0;
  d[2][3] = 1;
  d[2][4] = 1;
  d[3][0] = 1; 
  d[3][1] = 1;
  d[3][2] = 1;
  d[3][3] = 0;
  d[3][4] = 0;
  d[4][0] = 1; 
  d[4][1] = 1;
  d[4][2] = 1;
  d[4][3] = 0;
  d[4][4] = 0;

  // build tree
  std::vector<double> sample_ages(N,0);
  tb.QuickBuild(d,tree,sample_ages);
  bl.EM(data,tree,true); //tree is an empty tree, so branch lengths should be given by prior
  bl.MCMC(data,tree);



}

TEST_CASE( "Testing MCMC of branch lengths (2)" ){

  SECTION( "Testing MCMC of branch lengths" );

  //first test that branch lengths are determined by prior if no information is given

  int N = 4;
  int L = 2;
  Data data(N,L); //struct data is defined in data.hpp
  data.theta = 0.025;
  data.mu = 1e-8;
  data.pos.resize(L);
  data.rpos.resize(L);
  data.pos[0] = 0;
  data.rpos[0] = 0;
  data.pos[1] = 1;
  data.rpos[1] = 1;

  Tree tree;
  MinMatch tb(data);
  InferBranchLengths bl(data);

  CollapsedMatrix<float> d;
  d.resize(N,N);

  //do MCMC on empty tree
  //tb.QuickBuild(d,tree);
  //bl.MCMC(tree);

  d[0][0] = 0; 
  d[0][1] = 0;
  d[0][2] = 1;
  d[0][3] = 1;
  d[1][0] = 0; 
  d[1][1] = 0;
  d[1][2] = 1;
  d[1][3] = 1;
  d[2][0] = 1; 
  d[2][1] = 1;
  d[2][2] = 0;
  d[2][3] = 0;
  d[3][0] = 1; 
  d[3][1] = 1;
  d[3][2] = 0;
  d[3][3] = 0;

  // build tree

  std::vector<double> sample_ages(N,0);
  tb.QuickBuild(d,tree, sample_ages);
  bl.MCMC(data,tree);

  for(int i = 0; i < 2*N-1; i++){
    tree.nodes[i].num_events = 1.0*data.mu*tree.nodes[i].branch_length;
    tree.nodes[i].SNP_begin = 0;
    tree.nodes[i].SNP_end   = 1;
  }

  bl.MCMC(data,tree);

  //for(int i = 0; i < 2*N-1; i++){
  //  std::cerr << i << " " << tree.nodes[i].branch_length << " " << tree.nodes[i].num_events/(1.0*data.mu) << std::endl;
  //}
}

