#include "catch_amalgamated.hpp"

#include "data.hpp"
#include "anc.hpp"
#include "tree_builder.hpp"
#include "anc_builder.hpp"

TEST_CASE( "Testing Pearson Correlation "){

 Correlation cor(10);

 Leaves set1, set2;

 set1.num_leaves = 2;
 set1.member.resize(2);
 set1.member[0] = 1;
 set1.member[1] = 5;

 set2.num_leaves = 2;
 set2.member.resize(2);
 set2.member[0] = 1;
 set2.member[1] = 9;

 REQUIRE( std::fabs( cor.Pearson(set1, set2) - 0.375 ) < 1e-5 );

 set1.num_leaves = 6;
 set2.num_leaves = 6;
 set1.member.resize(6);
 set2.member.resize(6);
 for(int i = 0; i < 6; i++){
   set1.member[i] = i;
   set2.member[i] = i;
 }

 REQUIRE( std::fabs( cor.Pearson(set1, set2) - 1.0 ) < 1e-5 );

 Correlation cor_lance(6000);

 set1.num_leaves = 5000;
 set2.num_leaves = 5000;
 set1.member.resize(5000);
 set2.member.resize(5000);
 for(int i = 0; i < 5000; i++){
   set1.member[i] = i;
   set2.member[i] = i;
 }

 REQUIRE( std::fabs( cor_lance.Pearson(set1, set2) - 1.0 ) < 1e-5 );

}

TEST_CASE( "Testing code that finds equivalent branches" ){

  SECTION( "Testing code that finds equivalent branches" );

  int N = 5;
  int L = 1;
  Data data(N,L); //struct data is defined in data.hpp
  data.theta = 0.025;

  Tree tree;
  MinMatch tb(data);

  CollapsedMatrix<float> d;
  d.resize(N,N);

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
  tb.QuickBuild(d,tree, sample_ages);

  // I will find equivalent branches between the same tree.

  Tree same_tree = tree;

  AncesTreeBuilder ancbuilder(data, sample_ages);
  ancbuilder.PreCalcPotentialBranches(); 

  std::vector<int> equivalent_branches;
  ancbuilder.BranchAssociation(tree, same_tree, equivalent_branches); //O(N^2)

  for(int n = 0; n < 2*N-2; n++){
    REQUIRE(equivalent_branches[n] == n); //branch n is associated with branch n, since the trees are the same
  }

}

/*
TEST_CASE( "Testing optimize parameters" ){

  // populate data

  int N = 5;
  int L = 2*1000 + 10;
  Data data(N,L); //struct data is defined in data.hpp

  //populate data.sequence data.pos data.r

  data.sequence.resize(L,N);
  data.pos.resize(L);
  data.r.resize(L);
  data.rpos.resize(L);
  data.rpos[0] = 0.0;
  data.pos[0]  = 0;
  data.r[0]    = 0.0;
  for(int snp = 1; snp < L; snp++){
    data.pos[snp]  = snp;
    data.r[snp]    = 0.001;
    data.rpos[snp] = data.rpos[snp-1] + data.r[snp-1];
  }

  std::fill(data.sequence.vbegin(), data.sequence.vend(), 0);
  std::vector<char> row = {'0','1','1','0','0','0','0','0','0','0'};
  for(int snp = data.window_length; snp < data.window_length + 10; snp++){
    data.sequence[snp][0] = row[snp-data.window_length]; 
  }
  row = {'0','1','1','0','0','1','0','1','0','0'};
  for(int snp = data.window_length; snp < data.window_length + 10; snp++){
    data.sequence[snp][1] = row[snp-data.window_length]; 
  }
  row = {'0','1','1','0','0','0','0','0','0','0'};
  for(int snp = data.window_length; snp < data.window_length + 10; snp++){
    data.sequence[snp][2] = row[snp-data.window_length]; 
  }
  row = {'0','0','0','0','1','0','0','0','0','0'};
  for(int snp = data.window_length; snp < data.window_length + 10; snp++){
    data.sequence[snp][3] = row[snp-data.window_length]; 
  }
  row = {'0','0','0','0','1','0','0','0','0','0'};
  for(int snp = data.window_length; snp < data.window_length + 10; snp++){
    data.sequence[snp][4] = row[snp-data.window_length]; 
  }


  data.theta = 0.025;
  data.ntheta = 1.0 - data.theta;
  int section_startpos = 0;
  int section_endpos   = L-1;

  AncesTreeBuilder ab(data);
  REQUIRE(ab.OptimizeParameters(section_startpos, section_endpos, data).first == 0); 

  //now check a matrix explicitely
  
  //prepare by painting section
  const float log_theta_term = log((1.0 - data.theta)/data.theta);
  float min;
  DistanceMeasure d(data);

  std::vector<CollapsedMatrix<float>> topology;
  std::vector<std::vector<float>> logscales;
  CollapsedMatrix<float> alpha_begin, beta_end;

  topology.resize(data.N);
  logscales.resize(data.N);
  alpha_begin.resize(1,data.N);
  beta_end.resize(1,data.N);
  float logscale_alpha = 0.0, logscale_beta = 0.0;
  int boundarySNP_begin = section_startpos, boundarySNP_end = section_endpos;

  // populate variables that are independent of k 
  float prior_theta  = data.theta/(data.N-1.0) - data.ntheta/(data.N-1.0);
  float prior_ntheta = data.ntheta/(data.N-1.0);
  int derived;
  std::fill(beta_end.vbegin(), beta_end.vend(), 1.0);

  // paint
  FastPainting painter(data);
  float rescale = fast_log(data.theta/(1.0-data.theta));

  for(int k = 0; k < data.N; k++){
    for(int n = 0; n < data.N; n++){
      derived                   = (double) (data.sequence[0][k] > data.sequence[0][n]);
      alpha_begin[0][n]         = derived * prior_theta + prior_ntheta;
    }
    painter.RePaintSection(data, topology[k], logscales[k], alpha_begin, beta_end, boundarySNP_begin, boundarySNP_end, logscale_alpha, logscale_beta, k);

    float normalizing_constant = fast_log(N-1.0) - topology.size() * fast_log(data.ntheta); //topology.size() needs to be derived mutations of section
    for(int l = 0; l < (int) logscales[k].size(); l++){
      logscales[k][l] += normalizing_constant;
    }

  }
  int snp = section_startpos + data.window_length;
  d.AssignTopology(topology, logscales, section_startpos, section_endpos, snp);

  //look at snp 501

  d.GetMatrix(data.window_length * 1); //calculate d

  CollapsedMatrix<float> d_hiddenSNP;
  d_hiddenSNP.resize(N,N);

  for(int i = 0; i < data.N; i++){
    for(int j = 0; j < data.N; j++){
      d_hiddenSNP[i][j] = d.matrix[i][j];
    }
  }

  //modify distance matrix so that current SNP is cancelled
  for(int i = 0; i < data.N; i++){
    if(data.sequence[data.window_length + 1][i] == '1'){
      min = std::numeric_limits<float>::infinity();
      for(int j = 0; j < data.N; j++){
        if( data.sequence[data.window_length + 1][j] == '0' && i != j) d.matrix[i][j] -= log_theta_term;
        if(min > d.matrix[i][j]) min = d.matrix[i][j];
        assert(d.matrix[i][j] < std::numeric_limits<float>::infinity());
      }
    }
  }

  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 3; j++){
      REQUIRE( std::fabs( d.matrix[i][j] - d_hiddenSNP[i][j]) < 1e-4 );
    }
    for(int j = 3; j < 5; j++){
      REQUIRE( std::fabs( d.matrix[i][j] - d_hiddenSNP[i][j] - log(data.theta/(1.0-data.theta)) ) < 1e-4 );
    }
  }
  for(int i = 3; i < 5; i++){
    for(int j = 0; j < 5; j++){
      REQUIRE( std::fabs( d.matrix[i][j] - d_hiddenSNP[i][j]) < 1e-4 );
    }
  }

}
*/
