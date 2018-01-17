#include "catch.hpp"

#include "data.hpp"
#include "fast_painting.hpp"


TEST_CASE( "Testing painting" ){

  SECTION("Testing painting");

  //test two cases: recombination rate = 0, recombination_rate = infinity
  //Use RePaintSection

  // populate data

  int N = 5;
  int L = 10;
  Data data(N,L); //struct data is defined in data.hpp
  data.theta = 0.025;
  data.ntheta = 1.0 - data.theta;

  //populate data.sequence data.pos data.r

  data.sequence.resize(L,N);
  data.pos.resize(L);
  data.r.resize(L);
  for(int snp = 0; snp < L; snp++){
    data.pos[snp] = snp;
    data.r[snp]   = 0.0;
  }

  std::vector<char> row = {'0','1','1','0','0','0','0','0','0','0'};
  for(int snp = 0; snp < L; snp++){
    data.sequence[snp][0] = row[snp]; 
  }
  row = {'0','1','1','0','0','1','0','1','0','0'};
  for(int snp = 0; snp < L; snp++){
    data.sequence[snp][1] = row[snp]; 
  }
  row = {'0','1','0','0','0','0','0','0','0','0'};
  for(int snp = 0; snp < L; snp++){
    data.sequence[snp][2] = row[snp]; 
  }
  row = {'0','0','0','0','1','0','0','0','0','0'};
  for(int snp = 0; snp < L; snp++){
    data.sequence[snp][3] = row[snp]; 
  }
  row = {'0','0','0','0','1','0','0','0','0','0'};
  for(int snp = 0; snp < L; snp++){
    data.sequence[snp][4] = row[snp]; 
  }

  //the matrix I should be obtaining when data.r = 0
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


  /////////////////

  // declare painting variables

  CollapsedMatrix<float> topology, alpha_begin, beta_end;
  alpha_begin.resize(1,N);
  beta_end.resize(1,N);
  std::vector<float> logscales;
  float logscale_alpha = 0.0, logscale_beta = 0.0;
  int boundarySNP_begin = 0, boundarySNP_end = L-1;

  // populate variables that are independent of k 

  float prior_theta  = data.theta/(N-1.0) - data.ntheta/(N-1.0);
  float prior_ntheta = data.ntheta/(N-1.0);
  int derived;
  std::fill(beta_end.vbegin(), beta_end.vend(), 1.0);

  // paint

  FastPainting painter(data);
  float rescale              = fast_log(data.theta/(1.0-data.theta));

  for(int k = 0; k < N; k++){
    for(int n = 0; n < N; n++){
      derived                   = (double) (data.sequence[0][k] > data.sequence[0][n]);
      alpha_begin[0][n]         = derived * prior_theta + prior_ntheta;
    }
    painter.RePaintSection(data, topology, logscales, alpha_begin, beta_end, boundarySNP_begin, boundarySNP_end, logscale_alpha, logscale_beta, k);

    float normalizing_constant = fast_log(N-1.0) - topology.size() * fast_log(data.ntheta);
    
    //topology containts the posterior
    std::vector<float> ref_topology(N); //this is to check that values stay the same accross sequence (when data.r = 0)
    float ref_logscales = logscales[0];
    for(int n = 0; n < N; n++){
      ref_topology[n] = topology[0][n];
    }
    
    for(int l = 0; l < topology.size(); l++){
      REQUIRE(std::fabs( ref_logscales - logscales[l] ) < 1e-5);
      for(int n = 0; n < N; n++){
        REQUIRE(std::fabs( topology[l][n] - ref_topology[n] ) < 1e-5);
        if(k != n){
          //I need to round the output because I am using fast_log and some small numerical errors might be present
          REQUIRE( d[k][n] == std::round( (fast_log(topology[l][n]) + logscales[l] + normalizing_constant)/rescale ) );
        }
      }
    }
  }

}

