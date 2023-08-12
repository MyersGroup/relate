#include "fast_painting.hpp"

/****
 * Input
 * data: contains all the data
 * window_boundaries: SNP at which window begins (e.g. first window goes from window_boundaries[0] to window_boundaries[1]-1) 
 * pfiles: File pointers to files where painting gets dumped
 * k: haplotype to be painted
 *
 * Output
 * CollapsedMatrix<float> alpha, beta; contains the alpha and betas (forward & backward probabilities) at the window boundaries
 * std::vector<float> logscales_alpha, logscales_beta; : contains the logscales at the window boundaries
 * forward probability = exp(log(alpha[i]) + logscales_alpha[i])
 * boundarySNP_begin and boundarySNP_end are the SNP indices at which a window begins ends. May be different to window_boundaries because these might not be derived SNPs.
 */

void 
FastPainting::PaintSteppingStones(const Data& data, std::vector<int>& window_boundaries, std::vector<FILE*> pfiles, const int k){

  //additional vector that records boundary SNPs
  int num_windows = window_boundaries.size() - 1;
  std::vector<int> boundarySNP_begin(num_windows), boundarySNP_end(num_windows);
  std::vector<int>::iterator it_boundarySNP_begin = boundarySNP_begin.begin();
  std::vector<int>::iterator it_boundarySNP_end = boundarySNP_end.begin();
  std::vector<int>::reverse_iterator rit_boundarySNP_end = boundarySNP_end.rbegin();
  int window_index = 1;
  int window_end   = window_boundaries[1];
  //std::cerr << window_boundaries[num_windows] << " " << data.L << " " << window_boundaries[num_windows] - window_boundaries[0] << std::endl;;
  assert(window_boundaries[num_windows] == data.L);

  *it_boundarySNP_begin = 0;
  it_boundarySNP_begin++;

  //k is the sequence that is being painted
  CollapsedMatrix<float> alpha, beta; 
  std::vector<float> logscales_alpha, logscales_beta;

  /////////////////
  //precalculate quantities 

  std::vector<double> r_prob(data.L), nor_x_theta(data.L);
  std::vector<int> derived_k(data.L);    
  std::vector<double>::iterator it_r_prob = r_prob.begin(), it_r_prob_prev, it_nor_x_theta = nor_x_theta.begin();
  std::vector<int>::iterator it_derived_k = derived_k.begin();    

  char seq_k;
  double derived;
  int snp, snp_next; //previous SNP
  int last_snp = data.L - 1; //last SNP
  double tmp; //just a temporary variable for intermediate calculations

  *it_derived_k   = 0;
  *it_r_prob      = data.r[0];
  snp = 1;
  int num_derived_sites = 1;
  while(data.sequence[snp][k] != '1' && snp != last_snp){
    *it_r_prob   += data.r[snp];
    snp++;
  }
  if(snp >= window_end && (*it_derived_k) < window_end){ 
    while(window_end <= snp){
      *it_boundarySNP_end = snp;
      it_boundarySNP_end++;
      *it_boundarySNP_begin = *it_derived_k;
      it_boundarySNP_begin++;
      window_index++;
      window_end = window_boundaries[window_index];
    }
  }
  it_derived_k++;

  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));
  //if(*it_r_prob > 0.999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.99){
    *it_r_prob = 0.99;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  it_r_prob++;
  it_nor_x_theta++;

  *it_r_prob      = data.r[snp];
  *it_derived_k   = snp;

  snp++;
  num_derived_sites++;
  while(snp < data.L){
    //skip all non-derived sites
    while(data.sequence[snp][k] != '1' && snp != last_snp){
      *it_r_prob += data.r[snp];
      snp++;
    }

    if(snp >= window_end && (*it_derived_k) < window_end){ 
      while(window_end <= snp){
        *it_boundarySNP_end = snp;
        it_boundarySNP_end++;
        *it_boundarySNP_begin = *it_derived_k;
        it_boundarySNP_begin++;
        window_index++;
        window_end = window_boundaries[window_index];
      }
    }

    //this is a derived site
    it_derived_k++;

    *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
    *it_r_prob      = 1.0 - exp(-(*it_r_prob));
    //if(*it_r_prob > 0.999){
    //  *it_r_prob = 1.0;
    //  *it_nor_x_theta = 0.0;
    //}
    if(*it_r_prob > 0.99){
      *it_r_prob = 0.99;
      *it_nor_x_theta = log_small + log_ntheta;
    }

    it_r_prob++;
    it_nor_x_theta++;   

    *it_derived_k   = snp;
    *it_r_prob      = data.r[snp];

    snp++;
    num_derived_sites++;
  }
  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));
  //if(*it_r_prob > 0.999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.99){
    *it_r_prob = 0.99;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  *it_r_prob++;
  *it_r_prob      = 1.0; //just a technicality

  derived_k.resize(num_derived_sites);
  r_prob.resize(num_derived_sites+1);
  nor_x_theta.resize(num_derived_sites);

  *it_boundarySNP_end = last_snp;
  it_boundarySNP_end++;
  assert(it_boundarySNP_begin == boundarySNP_begin.end());
  assert(it_boundarySNP_end == boundarySNP_end.end());

  it_derived_k = derived_k.begin() + 1;
  it_r_prob = r_prob.begin();
  it_nor_x_theta = nor_x_theta.begin();

  /////////
  //Allocate memory
  alpha.resize(num_windows, data.N);
  beta.resize(num_windows, data.N);
  logscales_alpha.resize(num_windows);
  logscales_beta.resize(num_windows);

  auto it1_alpha     = alpha.ibegin();
  auto rit1_beta     = beta.irbegin();

  std::vector<float>::iterator it2_alpha;
  std::vector<float>::iterator it2_beta;
  std::vector<float>::iterator it_logscale_alpha = logscales_alpha.begin();
  std::vector<float>::reverse_iterator rit_logscale_beta = logscales_beta.rbegin();

  //alpha_aux and beta_aux contain the alpha and beta values along the sequence. I am alternating between two rows, to keep the previous and the current values
  CollapsedMatrix<double> alpha_aux, beta_aux;
  alpha_aux.resize(2,data.N);
  beta_aux.resize(2,data.N);
  std::vector<std::vector<double>::iterator> alpha_aux_rowbegin(2),alpha_aux_rowend(2);
  alpha_aux_rowbegin[0] = alpha_aux.vbegin();
  alpha_aux_rowbegin[1] = std::next(alpha_aux.vbegin(),data.N);
  alpha_aux_rowend[0]   = alpha_aux_rowbegin[1];
  alpha_aux_rowend[1]   = alpha_aux.vend();
  std::vector<std::vector<double>::iterator> beta_aux_rowbegin(2),beta_aux_rowend(2);
  beta_aux_rowbegin[0]  = beta_aux.vbegin();
  beta_aux_rowbegin[1]  = std::next(beta_aux.vbegin(),data.N);
  beta_aux_rowend[0]    = beta_aux_rowbegin[1];
  beta_aux_rowend[1]    = beta_aux.vend();
  int aux_index = 0, aux_index_prev = 1;

  double alpha_sum;
  double beta_sum;
  std::vector<double> logscale(2);

  std::vector<double>::iterator it2_alpha_aux;
  std::vector<double>::iterator it2_beta_aux;
  std::vector<double>::iterator it2_alpha_aux_prev;
  std::vector<double>::iterator it2_beta_aux_next;
  std::vector<char>::const_iterator it2_sequence;


  /////////
  //Forward Algorithm

  ////
  //SNP 0
  snp                         = *derived_k.begin();
  aux_index                   = 0;
  logscale[aux_index]         = 0.0;
  alpha_sum                   = 0.0; 

  it2_sequence                = data.sequence.rowbegin(snp);
  seq_k                       = *std::next(it2_sequence, k);

  it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
  //this loop vectorizes
  for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
    derived                   = (double) (seq_k > *it2_sequence);
    *it2_alpha_aux            = derived * prior_theta + prior_ntheta;

    it2_sequence++;
    it2_alpha_aux++;
  }

  it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
  *std::next(alpha_aux_rowbegin[aux_index],k) = 0.0;
  //this loop does not vectorize
  for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
    alpha_sum                += *it2_alpha_aux;
    it2_alpha_aux++;
  }


  it_boundarySNP_begin = boundarySNP_begin.begin();
  while(*it_boundarySNP_begin == snp){
    //store the start boundary of the current chunk  

    //copy to alpha, logscale
    it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
    it2_alpha                   = alpha.rowbegin(it1_alpha);
    for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
      *it2_alpha = *it2_alpha_aux;
      it2_alpha_aux++;
      it2_alpha++;
    }
    *it_logscale_alpha          = logscale[aux_index];
    it1_alpha++;
    it_logscale_alpha++;

    it_boundarySNP_begin++;
    if(it_boundarySNP_begin == boundarySNP_begin.end()) break;

  }


  ////
  //SNP > 0

  double r_x_alpha_sum;
  if(*it_r_prob < 1.0){
    r_x_alpha_sum = (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone) * alpha_sum;
  }else{
    r_x_alpha_sum = alpha_sum;
  }
  for(; it_derived_k != derived_k.end();){

    /////////////////
    //precalculated quantities 
    snp                     = *it_derived_k;
    it2_sequence            = data.sequence.rowbegin(snp);
    seq_k                   = *std::next(it2_sequence, k);
    aux_index               = (aux_index + 1) % 2;
    aux_index_prev          = 1 - aux_index;

    //////////////////
    //inner loop of forward algorithm

    if(*it_r_prob < 1){

      logscale[aux_index]         = logscale[aux_index_prev];
      logscale[aux_index]        += (*it_nor_x_theta);
      alpha_sum                   = 0.0; 

      it2_sequence                = data.sequence.rowbegin(snp);
      it2_alpha_aux_prev          = alpha_aux_rowbegin[aux_index_prev];

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      //this loop vectorizes
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        *it2_alpha_aux            = *it2_alpha_aux_prev + r_x_alpha_sum;
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha_aux           *= derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha_aux++;
        it2_alpha_aux_prev++;
        it2_sequence++;
      }

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      *std::next(alpha_aux_rowbegin[aux_index],k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        alpha_sum                += *it2_alpha_aux;
        it2_alpha_aux++;
      } 

    }else{

      logscale[aux_index_prev]    = logscale[aux_index];
      logscale[aux_index]        += log(data.ntheta/(Nminusone)*(alpha_sum));
      alpha_sum                   = 0.0;

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      it2_sequence                = data.sequence.rowbegin(snp);
      //this loop vectorizes
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha_aux            = derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha_aux++;
        it2_sequence++;
      }

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      *std::next(alpha_aux_rowbegin[aux_index],k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        alpha_sum                += *it2_alpha_aux;
        it2_alpha_aux++;
      } 

    }

    r_x_alpha_sum = alpha_sum;

    //check if alpha_sums get too small, if they do, rescale
    if(r_x_alpha_sum < lower_rescaling_threshold || r_x_alpha_sum > upper_rescaling_threshold){
      tmp           = r_x_alpha_sum; 
      it2_alpha_aux = alpha_aux_rowbegin[aux_index];
      //r_x_alpha_sum = 0.0;
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        *it2_alpha_aux /= tmp;
        //r_x_alpha_sum += *it2_alpha_aux;
        it2_alpha_aux++;
      }
      logscale[aux_index]   += log(tmp);
      r_x_alpha_sum          = 1.0;
      //std::cerr << r_x_alpha_sum << std::endl;
      assert(logscale[aux_index] < std::numeric_limits<double>::infinity());
    } 

    it_r_prob++; 
    if(*it_r_prob < 1.0){
      r_x_alpha_sum *= (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone);
    }

    if(it_boundarySNP_begin != boundarySNP_begin.end()){
      while(*it_boundarySNP_begin == snp){
        //store first the end boundary of the current chunk and then the start boundary of the next chunk  

        //copy to alpha, logscale
        it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
        it2_alpha                   = alpha.rowbegin(it1_alpha);
        for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
          *it2_alpha = *it2_alpha_aux;
          it2_alpha_aux++;
          it2_alpha++;
        }
        *it_logscale_alpha          = logscale[aux_index];
        it1_alpha++;
        it_logscale_alpha++;

        it_boundarySNP_begin++;
        if(it_boundarySNP_begin == boundarySNP_begin.end()) break;

      }
    }

    it_derived_k++;
    it_nor_x_theta++;
  }

  /* 
     std::cerr << *it_boundarySNP_begin << std::endl; 
     for(it_boundarySNP_begin = boundarySNP_begin.begin(); it_boundarySNP_begin != boundarySNP_begin.end(); it_boundarySNP_begin++){
     std::cerr << *it_boundarySNP_begin << " ";
     }
     std::cerr << std::endl;
     for(it_boundarySNP_end = boundarySNP_end.begin(); it_boundarySNP_end != boundarySNP_end.end(); it_boundarySNP_end++){
     std::cerr << *it_boundarySNP_end << " ";
     }
     std::cerr << std::endl;
     */

  assert(it_boundarySNP_begin == boundarySNP_begin.end());
  assert(it1_alpha == alpha.iend());
  assert(it_logscale_alpha == logscales_alpha.end());

  ///////////////
  //Backward algorithm

  normalizing_constant = (double) log(Nminusone) - num_derived_sites * log_ntheta;

  /////
  //SNP L-1
  it_derived_k--;
  assert(last_snp == *it_derived_k);

  logscale[aux_index]              = normalizing_constant;
  logscale[aux_index_prev]         = normalizing_constant;
  beta_sum                         = 0.0;

  it2_sequence                     = data.sequence.rowbegin(last_snp);
  seq_k                            = *std::next(it2_sequence, k);

  it2_beta_aux                     = beta_aux_rowbegin[aux_index];
  //this loop vectorizs
  for(; it2_beta_aux != beta_aux_rowend[aux_index];){
    *it2_beta_aux                  = 1.0;
    it2_beta_aux++;
  }

  it2_beta_aux                     = beta_aux_rowbegin[aux_index];
  //this loop does not vectorize
  for(; it2_beta_aux != beta_aux_rowend[aux_index];){
    if(seq_k  > *it2_sequence){
      beta_sum                    += data.theta;
    }else{
      beta_sum                    += data.ntheta;       
    }
    it2_sequence++;
    it2_beta_aux++;
  }
  beta_sum                        -= data.ntheta;

  rit_boundarySNP_end = boundarySNP_end.rbegin();
  while(*rit_boundarySNP_end == last_snp){
    //copy beta
    it2_beta                          = beta.rowbegin(rit1_beta);
    it2_beta_aux                      = beta_aux_rowbegin[aux_index];
    for(; it2_beta_aux != beta_aux_rowend[aux_index];){ 
      *it2_beta                       = *it2_beta_aux;
      it2_beta++;
      it2_beta_aux++;
    }
    *rit_logscale_beta                = logscale[aux_index];
    rit1_beta++;
    rit_logscale_beta++;
    rit_boundarySNP_end++;
    if(rit_boundarySNP_end == boundarySNP_end.rend()) break;
  }

  /////
  //SNP < L-1
  double beta_sum_oneminustheta, beta_sum_theta;
  double r_x_beta_sum;
  if(*it_r_prob < 1.0){
    r_x_beta_sum = (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone) * beta_sum;
  }else{
    r_x_beta_sum = beta_sum;
  } 
  snp_next = last_snp; 

  while(snp > 0){

    it_derived_k--;
    snp             = *it_derived_k;
    aux_index       = (aux_index + 1) %2;
    aux_index_prev  = 1 - aux_index;

    if(*it_r_prob < 1.0){

      //inner loop of backwards algorithm
      logscale[aux_index]                    = logscale[aux_index_prev];
      logscale[aux_index]                   += *it_nor_x_theta;
      beta_sum                               = 0.0;
      beta_sum_oneminustheta                 = r_x_beta_sum / data.ntheta;
      beta_sum_theta                         = r_x_beta_sum / data.theta - beta_sum_oneminustheta;
      it2_sequence                           = data.sequence.rowbegin(snp_next);

      it2_beta_aux_next                      = beta_aux_rowbegin[aux_index_prev];
      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      //this loop vectorizes
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){
        derived                              = (double) (seq_k > *it2_sequence);
        *it2_beta_aux                        = *it2_beta_aux_next + derived * beta_sum_theta + beta_sum_oneminustheta;
        *it2_beta_aux                       *= derived * theta_ratio + 1.0; 
        it2_sequence++;
        it2_beta_aux++;
        it2_beta_aux_next++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      *std::next(beta_aux_rowbegin[aux_index], k) = 0.0;
      //this loop does not vectorize
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){      
        if(seq_k > *it2_sequence){
          beta_sum                           += data.theta * (*it2_beta_aux); 
        }else{
          beta_sum                           += data.ntheta * (*it2_beta_aux); 
        }
        it2_sequence++;
        it2_beta_aux++;
      }

    }else{

      logscale[aux_index_prev]               = logscale[aux_index];
      logscale[aux_index]                   += log((data.ntheta/Nminusone) * alpha_sum);
      beta_sum                               = 0.0; 

      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      //this loop vectorizes
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){
        *it2_beta_aux                        = 1.0; 
        it2_beta_aux++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      *std::next(beta_aux_rowbegin[aux_index], k) = 0.0;
      //this loop does not vectorize
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){      
        if(seq_k > *it2_sequence){
          beta_sum                          += data.theta * (*it2_beta_aux); 
        }else{ 
          beta_sum                          += data.ntheta * (*it2_beta_aux); 
        }
        it2_sequence++;
        it2_beta_aux++;
      }

    }

    //std::cerr << beta_sum << std::endl;
    r_x_beta_sum = beta_sum;

    if(r_x_beta_sum < lower_rescaling_threshold || r_x_beta_sum > upper_rescaling_threshold){
      //if they get too small, rescale
      tmp          = r_x_beta_sum;
      it2_beta_aux = beta_aux_rowbegin[aux_index];
      //r_x_beta_sum = 0.0;
      for(; it2_beta_aux < beta_aux_rowend[aux_index];){
        *it2_beta_aux /= tmp;
        //r_x_beta_sum += *it2_beta_aux;
        it2_beta_aux++;
      }
      logscale[aux_index]      += fast_log(tmp);
      r_x_beta_sum              = 1.0;
      assert(logscale[aux_index] < std::numeric_limits<double>::infinity());
    }

    it_r_prob--;
    if(*it_r_prob < 1.0){
      r_x_beta_sum *= (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone);
    }

    //update beta and topology
    if(rit_boundarySNP_end != boundarySNP_end.rend()){ 
      while(*rit_boundarySNP_end == snp){
        //store first the start boundary of the current chunk and then the end boundary of the next chunk

        it2_beta                          = beta.rowbegin(rit1_beta);
        it2_beta_aux                      = beta_aux_rowbegin[aux_index];
        for(; it2_beta_aux != beta_aux_rowend[aux_index];){
          *it2_beta                       = *it2_beta_aux;
          //if(snp == 8493) std::cerr << *it2_beta << " " << fast_log(*it2_beta) + logscale[aux_index] << std::endl;
          it2_beta++;
          it2_beta_aux++;
        }
        *rit_logscale_beta                = logscale[aux_index];
        rit1_beta++;
        rit_logscale_beta++;

        rit_boundarySNP_end++;
        if(rit_boundarySNP_end == boundarySNP_end.rend()) break;
      }
    }

    snp_next = snp;
    it_nor_x_theta--; //I want this to be pointing at snp_next
  }

  assert(*it_derived_k == 0);
  assert(rit_boundarySNP_end == boundarySNP_end.rend());

  //Dump to file

  for(int i = 0; i < num_windows; i++){

    int startinterval = window_boundaries[i];
    int endinterval   = window_boundaries[i+1]-1;
    fwrite(&startinterval, sizeof(int), 1, pfiles[i]);
    fwrite(&endinterval, sizeof(int), 1, pfiles[i]);

    //dump alpha
    alpha.DumpToFile(pfiles[i], i, boundarySNP_begin, logscales_alpha); 
    //dump beta
    beta.DumpToFile(pfiles[i], i, boundarySNP_end, logscales_beta); 

  }

  //debug
  /*
     if(snp == 1422){
     for(int l = 0; l < num_derived_sites; l++){
     for(int n = 0; n < data.N; n++){
     if(n!=k){
     if( topology[l][n] == std::numeric_limits<double>::infinity() || std::isnan(topology[l][n]) ){
     std::cerr << k << " " << l << " " << n << " " << topology[l][n] << " " << alpha[l][n] << " " << beta[l][n] << std::endl;
     }
     }
     }
     }
     }
     */

}

void 
FastPainting::RePaintSection(const Data& data, CollapsedMatrix<float>& topology, std::vector<float>& logscales, CollapsedMatrix<float>& alpha_begin, CollapsedMatrix<float>& beta_end, int boundarySNP_begin, int boundarySNP_end, float logscale_alpha, float logscale_beta, const int k){

  /////////////////
  //precalculate quantities 

  int max_snps = boundarySNP_end - boundarySNP_begin + 2;
  CollapsedMatrix<double> alpha, beta; 
  std::vector<double> r_prob(max_snps), nor_x_theta(max_snps);
  std::vector<int> derived_k(max_snps);    
  std::vector<double>::iterator it_r_prob = r_prob.begin(), it_r_prob_prev, it_nor_x_theta = nor_x_theta.begin();
  std::vector<int>::iterator it_derived_k = derived_k.begin();    

  char seq_k;
  double derived;
  int snp, snp_next; //previous SNP
  int first_snp = boundarySNP_begin; 
  int last_snp  = boundarySNP_end; //last SNP
  double tmp; //just a temporary variable for intermediate calculations

  *it_derived_k   = first_snp;
  *it_r_prob      = data.r[first_snp];
  snp = first_snp+1;
  int num_derived_sites = 1;
  while(data.sequence[snp][k] != '1' && snp != last_snp){
    *it_r_prob   += data.r[snp];
    snp++;
  }

  it_derived_k++;

  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));
  //if(*it_r_prob > 0.999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.99){
    *it_r_prob = 0.99;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  it_r_prob++;
  it_nor_x_theta++;

  *it_r_prob      = data.r[snp];
  *it_derived_k   = snp;

  snp++;
  num_derived_sites++;
  while(snp <= last_snp){
    //skip all non-derived sites
    while(data.sequence[snp][k] != '1' && snp != last_snp){
      *it_r_prob += data.r[snp];
      snp++;
    }

    //this is a derived site
    it_derived_k++;

    *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
    *it_r_prob      = 1.0 - exp(-(*it_r_prob));
    //if(*it_r_prob > 0.999){
    //  *it_r_prob = 1.0;
    //  *it_nor_x_theta = 0.0;
    //}
    if(*it_r_prob > 0.99){
      *it_r_prob = 0.99;
      *it_nor_x_theta = log_small + log_ntheta;
    }

    it_r_prob++;
    it_nor_x_theta++;   
    *it_r_prob      = data.r[snp];

    *it_derived_k   = snp;

    snp++;
    num_derived_sites++;
  }
  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));
  //if(*it_r_prob > 0.999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.99){
    *it_r_prob = 0.99;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  *it_r_prob++;
  *it_r_prob      = 1.0; //just a technicality

  derived_k.resize(num_derived_sites);
  r_prob.resize(num_derived_sites+1);
  nor_x_theta.resize(num_derived_sites);

  it_derived_k = derived_k.begin() + 1;
  it_r_prob = r_prob.begin();
  it_nor_x_theta = nor_x_theta.begin();

  /////////
  //Allocate memory

  alpha.resize(num_derived_sites, data.N);
  beta.resize(num_derived_sites, data.N);
  topology.resize(num_derived_sites, data.N);
  logscales.resize(num_derived_sites);
  std::fill(logscales.begin(), logscales.end(), 0.0);

  double alpha_sum;
  double beta_sum;

  auto it1_alpha       = alpha.ibegin();
  auto rit1_beta       = beta.irbegin();
  auto rit1_alpha      = alpha.irbegin();
  auto rit1_topology   = topology.irbegin();

  std::vector<double>::iterator it2_alpha;
  std::vector<double>::iterator it2_beta;
  std::vector<double>::iterator it2_alpha_prev;
  std::vector<double>::iterator it2_beta_next;
  std::vector<float>::iterator it2_topology;

  std::vector<double>::iterator it2_alpha_rowbegin;
  std::vector<double>::iterator it2_beta_rowbegin;
  std::vector<char>::const_iterator it2_sequence;
  std::vector<float>::iterator it_logscale = logscales.begin(); 
  std::vector<float>::reverse_iterator rit_logscale = logscales.rbegin();


  /////////
  //Forward Algorithm

  ////
  //SNP 0
  *it_logscale                = logscale_alpha;
  alpha_sum                   = 0.0; 

  it2_sequence                = data.sequence.rowbegin(first_snp);
  seq_k                       = *std::next(it2_sequence, k);

  it2_alpha_rowbegin          = alpha.rowbegin(it1_alpha);
  it2_alpha                   = it2_alpha_rowbegin;
  std::vector<float>::iterator it2_alpha_begin = alpha_begin.vbegin(); 
  //this loop vectorizes
  for(; it2_alpha != alpha.rowend(it1_alpha);){
    *it2_alpha                = *it2_alpha_begin;
    it2_alpha_begin++;
    it2_alpha++;
  }

  it2_alpha                        = it2_alpha_rowbegin;
  *std::next(it2_alpha_rowbegin,k) = 0.0;
  //this loop does not vectorize
  for(; it2_alpha != alpha.rowend(it1_alpha);){
    alpha_sum                += *it2_alpha;
    it2_alpha++;
  }

  ////
  //SNP > 0

  double r_x_alpha_sum;
  if(*it_r_prob < 1.0){
    r_x_alpha_sum = (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone) * alpha_sum;
  }else{
    r_x_alpha_sum = alpha_sum;
  }
  double prev_logscale = *it_logscale;
  for(; it_derived_k != derived_k.end();){

    /////////////////
    //precalculated quantities 
    snp                     = *it_derived_k;
    it2_sequence            = data.sequence.rowbegin(snp);
    seq_k                   = *std::next(it2_sequence, k);


    //////////////////
    //inner loop of forward algorithm

    if(*it_r_prob < 1){

      it_logscale++;
      prev_logscale              += (*it_nor_x_theta);
      *it_logscale                = prev_logscale;
      alpha_sum                   = 0.0; 

      it2_sequence                = data.sequence.rowbegin(snp);
      it2_alpha_prev              = alpha.rowbegin(it1_alpha);

      it1_alpha++;
      it2_alpha_rowbegin          = alpha.rowbegin(it1_alpha);
      it2_alpha                   = it2_alpha_rowbegin;
      //this loop vectorizes
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        *it2_alpha                = *it2_alpha_prev + r_x_alpha_sum;
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha               *= derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha++;
        it2_alpha_prev++;
        it2_sequence++;
      }

      it2_alpha                        = it2_alpha_rowbegin;
      *std::next(it2_alpha_rowbegin,k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        alpha_sum                += (*it2_alpha);
        it2_alpha++;
      } 

    }else{

      it_logscale++;
      prev_logscale              += log((data.ntheta/Nminusone) * alpha_sum);
      *it_logscale                = prev_logscale;
      alpha_sum                   = 0.0;

      it1_alpha++;
      it2_alpha_rowbegin          = alpha.rowbegin(it1_alpha);
      it2_alpha                   = it2_alpha_rowbegin;
      it2_sequence                = data.sequence.rowbegin(snp);
      //this loop vectorizes
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha                = derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha++;
        it2_sequence++;
      }

      it2_alpha                        = it2_alpha_rowbegin;
      *std::next(it2_alpha_rowbegin,k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        alpha_sum                += (*it2_alpha);
        it2_alpha++;
      } 

    }
    r_x_alpha_sum = alpha_sum;

    //check if alpha_sums get too small, if they do, rescale
    if(r_x_alpha_sum < lower_rescaling_threshold || r_x_alpha_sum > upper_rescaling_threshold){
      tmp       = r_x_alpha_sum; 
      it2_alpha = it2_alpha_rowbegin;
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        *it2_alpha /= tmp;
        it2_alpha++;
      }
      prev_logscale         += log(tmp);
      *it_logscale          += log(tmp);
      r_x_alpha_sum          = 1.0;
      assert(*it_logscale < std::numeric_limits<double>::infinity());
    } 

    it_r_prob++; 
    if(*it_r_prob < 1.0){
      r_x_alpha_sum *= (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone);
    }

    it_derived_k++;
    it_nor_x_theta++;
  }

  ///////////////
  //Backward algorithm

  /////
  //SNP L-1
  it_derived_k--;
  assert(last_snp == *it_derived_k);

  *it_logscale                    += logscale_beta;
  beta_sum                         = 0.0;

  it2_sequence                     = data.sequence.rowbegin(last_snp);
  seq_k                            = *std::next(it2_sequence, k);

  it2_beta_rowbegin                = beta.rowbegin(rit1_beta);
  it2_beta                         = it2_beta_rowbegin;
  std::vector<float>::iterator it2_beta_end = beta_end.vbegin();
  //this loop vectorizes
  for(; it2_beta != beta.rowend(rit1_beta);){
    *it2_beta                      = *it2_beta_end;
    it2_beta_end++;
    it2_beta++;
  }

  it2_beta                         = it2_beta_rowbegin;
  *std::next(it2_beta_rowbegin, k)       = 0.0;
  //this loop does not vectorize
  for(; it2_beta != beta.rowend(rit1_beta);){      
    if(seq_k > *it2_sequence){
      beta_sum                    += data.theta * (*it2_beta); 
    }else{
      beta_sum                    += data.ntheta * (*it2_beta); 
    }
    it2_sequence++;
    it2_beta++;
  }

  //calculate topology matrix
  it2_topology                     = topology.rowbegin(rit1_topology);
  it2_alpha                        = alpha.rowbegin(rit1_alpha);
  it2_beta                         = it2_beta_rowbegin;
  for(; it2_topology != topology.rowend(rit1_topology);){ 
    //*it2_topology                   = std::log(*it2_alpha) + std::log(*it2_beta);
    *it2_topology                  = (*it2_alpha) * (*it2_beta);
    //std::cout << *it2_topology << " ";
    it2_topology++;
    it2_alpha++;
    it2_beta++;
  }
  //std::cout << std::endl;

  rit1_alpha++;
  rit1_topology++;

  /////
  //SNP < L-1
  double beta_sum_oneminustheta, beta_sum_theta;
  double r_x_beta_sum;
  if(*it_r_prob < 1.0){
    r_x_beta_sum = (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone) * beta_sum;
  }else{
    r_x_beta_sum = beta_sum;
  } 
  snp_next = last_snp;
  prev_logscale = logscale_beta;
  while(rit1_topology != topology.irend()){

    it_derived_k--;
    snp          = *it_derived_k;

    if(*it_r_prob < 1.0){

      //inner loop of backwards algorithm
      it_logscale--;

      prev_logscale                         += *it_nor_x_theta;
      *it_logscale                          += prev_logscale;  
      beta_sum                               = 0.0;
      beta_sum_oneminustheta                 = r_x_beta_sum / data.ntheta;
      beta_sum_theta                         = r_x_beta_sum / data.theta - beta_sum_oneminustheta;
      it2_sequence                           = data.sequence.rowbegin(snp_next);
      it2_beta_next                          = beta.rowbegin(rit1_beta);

      rit1_beta++;
      it2_beta_rowbegin                      = beta.rowbegin(rit1_beta);
      it2_beta                               = it2_beta_rowbegin;
      //this loop vectorizes
      for(; it2_beta != beta.rowend(rit1_beta);){
        derived                              = (double) (seq_k > *it2_sequence);
        *it2_beta                            = *it2_beta_next + derived * beta_sum_theta + beta_sum_oneminustheta;
        *it2_beta                           *= derived * theta_ratio + 1.0; 
        it2_sequence++;
        it2_beta++;
        it2_beta_next++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta                               = it2_beta_rowbegin;
      *std::next(it2_beta_rowbegin, k)       = 0.0;
      //this loop does not vectorize
      for(; it2_beta != beta.rowend(rit1_beta);){      
        if(seq_k > *it2_sequence){
          beta_sum                          += data.theta * (*it2_beta); 
        }else{
          beta_sum                          += data.ntheta * (*it2_beta); 
        }
        it2_sequence++;
        it2_beta++;
      }

    }else{

      it_logscale--;
      prev_logscale                         += log(beta_sum/Nminusone);
      *it_logscale                          += prev_logscale;  
      beta_sum                               = 0.0; 

      rit1_beta++;
      it2_beta_rowbegin                      = beta.rowbegin(rit1_beta);
      it2_beta                               = it2_beta_rowbegin;
      //this loop vectorizes
      for(; it2_beta != beta.rowend(rit1_beta);){
        *it2_beta                            = 1.0; 
        it2_beta++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta                               = it2_beta_rowbegin;
      *std::next(it2_beta_rowbegin, k)       = 0.0;
      //this loop does not vectorize
      for(; it2_beta != beta.rowend(rit1_beta);){      
        if(seq_k > *it2_sequence){
          beta_sum                          += data.theta; 
        }else{  
          beta_sum                          += data.ntheta; 
        }
        it2_sequence++;
        it2_beta++;
      }

    }

    r_x_beta_sum = beta_sum;

    //calculate topology matrix
    it2_topology                           = topology.rowbegin(rit1_topology);
    it2_alpha                              = alpha.rowbegin(rit1_alpha);
    it2_beta                               = it2_beta_rowbegin;
    for(; it2_topology != topology.rowend(rit1_topology);){  
      //*it2_topology                        = std::log(*it2_alpha) + std::log(*it2_beta);
      *it2_topology                        = (*it2_alpha) * (*it2_beta);
      //std::cout << *it2_topology << " ";
      it2_topology++;
      it2_alpha++;
      it2_beta++;
    }
    //std::cout << std::endl;

    if(r_x_beta_sum < lower_rescaling_threshold || r_x_beta_sum > upper_rescaling_threshold){
      //if they get too small, rescale
      tmp        = r_x_beta_sum;
      it2_beta   = it2_beta_rowbegin;
      //r_x_beta_sum = 0.0;
      for(; it2_beta < beta.rowend(rit1_beta);){
        *it2_beta /= tmp;
        //r_x_beta_sum += *it2_beta;
        it2_beta++;
      }
      prev_logscale            += log(tmp);
      *it_logscale             += log(tmp);
      r_x_beta_sum              = 1.0;
      assert(*it_logscale < std::numeric_limits<double>::infinity());
    }

    it_r_prob--;
    if(*it_r_prob < 1.0){
      r_x_beta_sum *= (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone);
    }

    snp_next = snp;
    //I want these to be pointing at snp_next
    it_nor_x_theta--;
    rit1_alpha++;
    rit1_topology++;
  }

  assert(*it_derived_k == first_snp);

  //debug
  /*
     if(snp == 1422){
     for(int l = 0; l < num_derived_sites; l++){
     for(int n = 0; n < data.N; n++){
     if(n!=k){
     if( topology[l][n] == std::numeric_limits<double>::infinity() || std::isnan(topology[l][n]) ){
     std::cerr << k << " " << l << " " << n << " " << topology[l][n] << " " << alpha[l][n] << " " << beta[l][n] << std::endl;
     }
     }
     }
     }
     }
     */

}



void 
FastPainting::PaintSteppingStonesNormal(const Data& data, std::vector<int>& window_boundaries, std::vector<FILE*> pfiles, const int k){

  //additional vector that records boundary SNPs
  int num_windows = window_boundaries.size() - 1;
  std::vector<int> boundarySNP_begin(num_windows), boundarySNP_end(num_windows);
  std::vector<int>::iterator it_boundarySNP_begin = boundarySNP_begin.begin();
  std::vector<int>::iterator it_boundarySNP_end = boundarySNP_end.begin();
  std::vector<int>::reverse_iterator rit_boundarySNP_end = boundarySNP_end.rbegin();
  int window_index = 1;
  int window_end   = window_boundaries[1];
  //std::cerr << window_boundaries[num_windows] << " " << data.L << " " << window_boundaries[num_windows] - window_boundaries[0] << std::endl;;
  assert(window_boundaries[num_windows] == data.L);

  *it_boundarySNP_begin = 0;
  it_boundarySNP_begin++;

  //k is the sequence that is being painted
  CollapsedMatrix<float> alpha, beta; 
  std::vector<float> logscales_alpha, logscales_beta;

  /////////////////
  //precalculate quantities 

  std::vector<double> r_prob(data.L+10), nor_x_theta(data.L+10);
  std::vector<int> derived_k(data.L+10);    
  std::vector<double>::iterator it_r_prob = r_prob.begin(), it_r_prob_prev, it_nor_x_theta = nor_x_theta.begin();
  std::vector<int>::iterator it_derived_k = derived_k.begin();    

  char seq_k;
  double derived;
  int snp, snp_next; //previous SNP
  int last_snp = data.L - 1; //last SNP
  double tmp; //just a temporary variable for intermediate calculations

  *it_derived_k   = 0;
  *it_r_prob      = data.r[0];
  snp = 1;
  int num_derived_sites = 1;

  if(snp >= window_end && (*it_derived_k) < window_end){ 
    while(window_end <= snp){
      *it_boundarySNP_end = snp;
      it_boundarySNP_end++;
      *it_boundarySNP_begin = *it_derived_k;
      it_boundarySNP_begin++;
      window_index++;
      window_end = window_boundaries[window_index];
    }
  }
  it_derived_k++;

  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));

  //if(*it_r_prob > 0.9999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.999){
    *it_r_prob = 0.999;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  it_r_prob++;
  it_nor_x_theta++;

  *it_r_prob      = data.r[snp];
  *it_derived_k   = snp;

  snp++;
  num_derived_sites++;
  while(snp < data.L){

    if(snp >= window_end && (*it_derived_k) < window_end){ 
      while(window_end <= snp){
        *it_boundarySNP_end = snp;
        it_boundarySNP_end++;
        *it_boundarySNP_begin = *it_derived_k;
        it_boundarySNP_begin++;
        window_index++;
        window_end = window_boundaries[window_index];
      }
    }

    //this is a derived site
    it_derived_k++;

    *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
    *it_r_prob      = 1.0 - exp(-(*it_r_prob));

    //if(*it_r_prob > 0.9999){
    //  *it_r_prob = 1.0;
    //  *it_nor_x_theta = 0.0;
    //}
    if(*it_r_prob > 0.999){
      *it_r_prob = 0.999;
      *it_nor_x_theta = log_small + log_ntheta;
    }

    it_r_prob++;
    it_nor_x_theta++;   

    *it_derived_k   = snp;
    *it_r_prob      = data.r[snp];

    snp++;
    num_derived_sites++;
  }

  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));
  //if(*it_r_prob > 0.9999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.999){
    *it_r_prob = 0.999;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  *it_r_prob++;
  *it_r_prob      = 1.0; //just a technicality

  derived_k.resize(num_derived_sites);
  r_prob.resize(num_derived_sites+1);
  nor_x_theta.resize(num_derived_sites);

  *it_boundarySNP_end = last_snp;
  it_boundarySNP_end++;
  assert(it_boundarySNP_begin == boundarySNP_begin.end());
  assert(it_boundarySNP_end == boundarySNP_end.end());

  it_derived_k = derived_k.begin() + 1;
  it_r_prob = r_prob.begin();
  it_nor_x_theta = nor_x_theta.begin();

  /////////
  //Allocate memory
  alpha.resize(num_windows, data.N);
  beta.resize(num_windows, data.N);
  logscales_alpha.resize(num_windows);
  logscales_beta.resize(num_windows);

  auto it1_alpha     = alpha.ibegin();
  auto rit1_beta     = beta.irbegin();

  std::vector<float>::iterator it2_alpha;
  std::vector<float>::iterator it2_beta;
  std::vector<float>::iterator it_logscale_alpha = logscales_alpha.begin();
  std::vector<float>::reverse_iterator rit_logscale_beta = logscales_beta.rbegin();

  //alpha_aux and beta_aux contain the alpha and beta values along the sequence. I am alternating between two rows, to keep the previous and the current values
  CollapsedMatrix<double> alpha_aux, beta_aux;
  alpha_aux.resize(2,data.N);
  beta_aux.resize(2,data.N);
  std::vector<std::vector<double>::iterator> alpha_aux_rowbegin(2),alpha_aux_rowend(2);
  alpha_aux_rowbegin[0] = alpha_aux.vbegin();
  alpha_aux_rowbegin[1] = std::next(alpha_aux.vbegin(),data.N);
  alpha_aux_rowend[0]   = alpha_aux_rowbegin[1];
  alpha_aux_rowend[1]   = alpha_aux.vend();
  std::vector<std::vector<double>::iterator> beta_aux_rowbegin(2),beta_aux_rowend(2);
  beta_aux_rowbegin[0]  = beta_aux.vbegin();
  beta_aux_rowbegin[1]  = std::next(beta_aux.vbegin(),data.N);
  beta_aux_rowend[0]    = beta_aux_rowbegin[1];
  beta_aux_rowend[1]    = beta_aux.vend();
  int aux_index = 0, aux_index_prev = 1;

  double alpha_sum;
  double beta_sum;
  std::vector<double> logscale(2);

  std::vector<double>::iterator it2_alpha_aux;
  std::vector<double>::iterator it2_beta_aux;
  std::vector<double>::iterator it2_alpha_aux_prev;
  std::vector<double>::iterator it2_beta_aux_next;
  std::vector<char>::const_iterator it2_sequence;

  /////////
  //Forward Algorithm

  ////
  //SNP 0
  snp                         = *derived_k.begin();
  aux_index                   = 0;
  logscale[aux_index]         = 0.0;
  alpha_sum                   = 0.0; 

  it2_sequence                = data.sequence.rowbegin(snp);
  seq_k                       = *std::next(it2_sequence, k);

  it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
  //this loop vectorizes
  for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
    derived                   = (double) (seq_k > *it2_sequence);
    *it2_alpha_aux            = derived * prior_theta + prior_ntheta;

    it2_sequence++;
    it2_alpha_aux++;
  }

  it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
  *std::next(alpha_aux_rowbegin[aux_index],k) = 0.0;
  //this loop does not vectorize
  for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
    alpha_sum                += *it2_alpha_aux;
    it2_alpha_aux++;
  }

  it_boundarySNP_begin = boundarySNP_begin.begin();
  while(*it_boundarySNP_begin == snp){
    //store the start boundary of the current chunk  

    //copy to alpha, logscale
    it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
    it2_alpha                   = alpha.rowbegin(it1_alpha);
    for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
      *it2_alpha = *it2_alpha_aux;
      it2_alpha_aux++;
      it2_alpha++;
    }
    *it_logscale_alpha          = logscale[aux_index];
    it1_alpha++;
    it_logscale_alpha++;

    it_boundarySNP_begin++;
    if(it_boundarySNP_begin == boundarySNP_begin.end()) break;

  }

  ////
  //SNP > 0

  double r_x_alpha_sum;
  if(*it_r_prob < 1.0){
    r_x_alpha_sum = (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone) * alpha_sum;
  }else{
    r_x_alpha_sum = alpha_sum;
  }
  for(; it_derived_k != derived_k.end();){

    /////////////////
    //precalculated quantities 
    snp                     = *it_derived_k;
    it2_sequence            = data.sequence.rowbegin(snp);
    seq_k                   = *std::next(it2_sequence, k);
    aux_index               = (aux_index + 1) % 2;
    aux_index_prev          = 1 - aux_index;

    //////////////////
    //inner loop of forward algorithm

    if(*it_r_prob < 1){

      logscale[aux_index]         = logscale[aux_index_prev];
      logscale[aux_index]        += (*it_nor_x_theta);
      alpha_sum                   = 0.0; 

      it2_sequence                = data.sequence.rowbegin(snp);
      it2_alpha_aux_prev          = alpha_aux_rowbegin[aux_index_prev];

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      //this loop vectorizes
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        *it2_alpha_aux            = *it2_alpha_aux_prev + r_x_alpha_sum;
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha_aux           *= derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha_aux++;
        it2_alpha_aux_prev++;
        it2_sequence++;
      }

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      *std::next(alpha_aux_rowbegin[aux_index],k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        alpha_sum                += *it2_alpha_aux;
        it2_alpha_aux++;
      } 

    }else{

      logscale[aux_index_prev]    = logscale[aux_index];
      logscale[aux_index]        += log(data.ntheta/(Nminusone)*(alpha_sum));
      alpha_sum                   = 0.0;

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      it2_sequence                = data.sequence.rowbegin(snp);
      //this loop vectorizes
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha_aux            = derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha_aux++;
        it2_sequence++;
      }

      it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
      *std::next(alpha_aux_rowbegin[aux_index],k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        alpha_sum                += *it2_alpha_aux;
        it2_alpha_aux++;
      } 

    }

    r_x_alpha_sum = alpha_sum;

    //check if alpha_sums get too small, if they do, rescale
    if(r_x_alpha_sum < lower_rescaling_threshold || r_x_alpha_sum > upper_rescaling_threshold){
      tmp           = r_x_alpha_sum; 
      it2_alpha_aux = alpha_aux_rowbegin[aux_index];
      //r_x_alpha_sum = 0.0;
      for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
        *it2_alpha_aux /= tmp;
        //r_x_alpha_sum += *it2_alpha_aux;
        it2_alpha_aux++;
      }
      logscale[aux_index]   += log(tmp);
      r_x_alpha_sum          = 1.0;
      //std::cerr << r_x_alpha_sum << std::endl;
      assert(logscale[aux_index] < std::numeric_limits<double>::infinity());
    } 

    it_r_prob++; 
    if(*it_r_prob < 1.0){
      r_x_alpha_sum *= (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone);
    }

    if(it_boundarySNP_begin != boundarySNP_begin.end()){
      while(*it_boundarySNP_begin == snp){
        //store first the end boundary of the current chunk and then the start boundary of the next chunk  

        //copy to alpha, logscale
        it2_alpha_aux               = alpha_aux_rowbegin[aux_index];
        it2_alpha                   = alpha.rowbegin(it1_alpha);
        for(; it2_alpha_aux != alpha_aux_rowend[aux_index];){
          *it2_alpha = *it2_alpha_aux;
          it2_alpha_aux++;
          it2_alpha++;
        }
        *it_logscale_alpha          = logscale[aux_index];
        it1_alpha++;
        it_logscale_alpha++;

        it_boundarySNP_begin++;
        if(it_boundarySNP_begin == boundarySNP_begin.end()) break;

      }
    }

    it_derived_k++;
    it_nor_x_theta++;
  }

  /* 
     std::cerr << *it_boundarySNP_begin << std::endl; 
     for(it_boundarySNP_begin = boundarySNP_begin.begin(); it_boundarySNP_begin != boundarySNP_begin.end(); it_boundarySNP_begin++){
     std::cerr << *it_boundarySNP_begin << " ";
     }
     std::cerr << std::endl;
     for(it_boundarySNP_end = boundarySNP_end.begin(); it_boundarySNP_end != boundarySNP_end.end(); it_boundarySNP_end++){
     std::cerr << *it_boundarySNP_end << " ";
     }
     std::cerr << std::endl;
     */

  assert(it_boundarySNP_begin == boundarySNP_begin.end());
  assert(it1_alpha == alpha.iend());
  assert(it_logscale_alpha == logscales_alpha.end());

  ///////////////
  //Backward algorithm

  normalizing_constant = (double) log(Nminusone) - num_derived_sites * log_ntheta;

  /////
  //SNP L-1
  it_derived_k--;
  assert(last_snp == *it_derived_k);

  logscale[aux_index]              = normalizing_constant;
  logscale[aux_index_prev]         = normalizing_constant;
  beta_sum                         = 0.0;

  it2_sequence                     = data.sequence.rowbegin(last_snp);
  seq_k                            = *std::next(it2_sequence, k);

  it2_beta_aux                     = beta_aux_rowbegin[aux_index];
  //this loop vectorizs
  for(; it2_beta_aux != beta_aux_rowend[aux_index];){
    *it2_beta_aux                  = 1.0;
    it2_beta_aux++;
  }

  it2_beta_aux                     = beta_aux_rowbegin[aux_index];
  //this loop does not vectorize
  for(; it2_beta_aux != beta_aux_rowend[aux_index];){
    if(seq_k  > *it2_sequence){
      beta_sum                    += data.theta;
    }else{
      beta_sum                    += data.ntheta;       
    }
    it2_sequence++;
    it2_beta_aux++;
  }
  beta_sum                        -= data.ntheta;

  rit_boundarySNP_end = boundarySNP_end.rbegin();
  while(*rit_boundarySNP_end == last_snp){
    //copy beta
    it2_beta                          = beta.rowbegin(rit1_beta);
    it2_beta_aux                      = beta_aux_rowbegin[aux_index];
    for(; it2_beta_aux != beta_aux_rowend[aux_index];){ 
      *it2_beta                       = *it2_beta_aux;
      it2_beta++;
      it2_beta_aux++;
    }
    *rit_logscale_beta                = logscale[aux_index];
    rit1_beta++;
    rit_logscale_beta++;
    rit_boundarySNP_end++;
    if(rit_boundarySNP_end == boundarySNP_end.rend()) break;
  }

  /////
  //SNP < L-1
  double beta_sum_oneminustheta, beta_sum_theta;
  double r_x_beta_sum;
  if(*it_r_prob < 1.0){
    r_x_beta_sum = (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone) * beta_sum;
  }else{
    r_x_beta_sum = beta_sum;
  } 
  snp_next = last_snp; 

  while(snp > 0){

    it_derived_k--;
    snp             = *it_derived_k;
    aux_index       = (aux_index + 1) %2;
    aux_index_prev  = 1 - aux_index;

    if(*it_r_prob < 1.0){

      //inner loop of backwards algorithm
      logscale[aux_index]                    = logscale[aux_index_prev];
      logscale[aux_index]                   += *it_nor_x_theta;
      beta_sum                               = 0.0;
      beta_sum_oneminustheta                 = r_x_beta_sum / data.ntheta;
      beta_sum_theta                         = r_x_beta_sum / data.theta - beta_sum_oneminustheta;
      it2_sequence                           = data.sequence.rowbegin(snp_next);

      it2_beta_aux_next                      = beta_aux_rowbegin[aux_index_prev];
      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      //this loop vectorizes
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){
        derived                              = (double) (seq_k > *it2_sequence);
        *it2_beta_aux                        = *it2_beta_aux_next + derived * beta_sum_theta + beta_sum_oneminustheta;
        *it2_beta_aux                       *= derived * theta_ratio + 1.0; 
        it2_sequence++;
        it2_beta_aux++;
        it2_beta_aux_next++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      *std::next(beta_aux_rowbegin[aux_index], k) = 0.0;
      //this loop does not vectorize
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){      
        if(seq_k > *it2_sequence){
          beta_sum                           += data.theta * (*it2_beta_aux); 
        }else{
          beta_sum                           += data.ntheta * (*it2_beta_aux); 
        }
        it2_sequence++;
        it2_beta_aux++;
      }

    }else{

      logscale[aux_index_prev]               = logscale[aux_index];
      logscale[aux_index]                   += log((data.ntheta/Nminusone) * alpha_sum);
      beta_sum                               = 0.0; 

      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      //this loop vectorizes
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){
        *it2_beta_aux                        = 1.0; 
        it2_beta_aux++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta_aux                           = beta_aux_rowbegin[aux_index];
      *std::next(beta_aux_rowbegin[aux_index], k) = 0.0;
      //this loop does not vectorize
      for(; it2_beta_aux != beta_aux_rowend[aux_index];){      
        if(seq_k > *it2_sequence){
          beta_sum                          += data.theta * (*it2_beta_aux); 
        }else{ 
          beta_sum                          += data.ntheta * (*it2_beta_aux); 
        }
        it2_sequence++;
        it2_beta_aux++;
      }

    }

    //std::cerr << beta_sum << std::endl;
    r_x_beta_sum = beta_sum;

    if(r_x_beta_sum < lower_rescaling_threshold || r_x_beta_sum > upper_rescaling_threshold){
      //if they get too small, rescale
      tmp          = r_x_beta_sum;
      it2_beta_aux = beta_aux_rowbegin[aux_index];
      //r_x_beta_sum = 0.0;
      for(; it2_beta_aux < beta_aux_rowend[aux_index];){
        *it2_beta_aux /= tmp;
        //r_x_beta_sum += *it2_beta_aux;
        it2_beta_aux++;
      }
      logscale[aux_index]      += fast_log(tmp);
      r_x_beta_sum              = 1.0;
      assert(logscale[aux_index] < std::numeric_limits<double>::infinity());
    }

    it_r_prob--;
    if(*it_r_prob < 1.0){
      r_x_beta_sum *= (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone);
    }

    //update beta and topology
    if(rit_boundarySNP_end != boundarySNP_end.rend()){ 
      while(*rit_boundarySNP_end == snp){
        //store first the start boundary of the current chunk and then the end boundary of the next chunk

        it2_beta                          = beta.rowbegin(rit1_beta);
        it2_beta_aux                      = beta_aux_rowbegin[aux_index];
        for(; it2_beta_aux != beta_aux_rowend[aux_index];){
          *it2_beta                       = *it2_beta_aux;
          //if(snp == 8493) std::cerr << *it2_beta << " " << fast_log(*it2_beta) + logscale[aux_index] << std::endl;
          it2_beta++;
          it2_beta_aux++;
        }
        *rit_logscale_beta                = logscale[aux_index];
        rit1_beta++;
        rit_logscale_beta++;

        rit_boundarySNP_end++;
        if(rit_boundarySNP_end == boundarySNP_end.rend()) break;
      }
    }

    snp_next = snp;
    it_nor_x_theta--; //I want this to be pointing at snp_next
  }

  assert(*it_derived_k == 0);
  assert(rit_boundarySNP_end == boundarySNP_end.rend());

  //Dump to file

  for(int i = 0; i < num_windows; i++){

    int startinterval = window_boundaries[i];
    int endinterval   = window_boundaries[i+1]-1;
    fwrite(&startinterval, sizeof(int), 1, pfiles[i]);
    fwrite(&endinterval, sizeof(int), 1, pfiles[i]);

    //dump alpha
    alpha.DumpToFile(pfiles[i], i, boundarySNP_begin, logscales_alpha); 
    //dump beta
    beta.DumpToFile(pfiles[i], i, boundarySNP_end, logscales_beta); 

  }


  //debug
  /*
     if(snp == 1422){
     for(int l = 0; l < num_derived_sites; l++){
     for(int n = 0; n < data.N; n++){
     if(n!=k){
     if( topology[l][n] == std::numeric_limits<double>::infinity() || std::isnan(topology[l][n]) ){
     std::cerr << k << " " << l << " " << n << " " << topology[l][n] << " " << alpha[l][n] << " " << beta[l][n] << std::endl;
     }
     }
     }
     }
     }
     */

}

void 
FastPainting::RePaintSectionNormal(const Data& data, CollapsedMatrix<float>& topology, std::vector<float>& logscales, CollapsedMatrix<float>& alpha_begin, CollapsedMatrix<float>& beta_end, int boundarySNP_begin, int boundarySNP_end, float logscale_alpha, float logscale_beta, const int k){

  /////////////////
  //precalculate quantities 

  int max_snps = boundarySNP_end - boundarySNP_begin + 2;
  CollapsedMatrix<double> alpha, beta; 
  std::vector<double> r_prob(max_snps), nor_x_theta(max_snps);
  std::vector<int> derived_k(max_snps);    
  std::vector<double>::iterator it_r_prob = r_prob.begin(), it_r_prob_prev, it_nor_x_theta = nor_x_theta.begin();
  std::vector<int>::iterator it_derived_k = derived_k.begin();    

  char seq_k;
  double derived;
  int snp, snp_next; //previous SNP
  int first_snp = boundarySNP_begin; 
  int last_snp  = boundarySNP_end; //last SNP
  double tmp; //just a temporary variable for intermediate calculations

  *it_derived_k   = first_snp;
  *it_r_prob      = data.r[first_snp];
  snp = first_snp+1;
  int num_derived_sites = 1;

  it_derived_k++;

  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));

  //if(*it_r_prob > 0.9999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.999){
    *it_r_prob = 0.999;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  it_r_prob++;
  it_nor_x_theta++;

  *it_r_prob      = data.r[snp];
  *it_derived_k   = snp;

  snp++;
  num_derived_sites++;
  while(snp <= last_snp){

    //this is a derived site
    it_derived_k++;

    *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
    *it_r_prob      = 1.0 - exp(-(*it_r_prob));

    //if(*it_r_prob > 0.9999){
    //  *it_r_prob = 1.0;
    //  *it_nor_x_theta = 0.0;
    //}
    if(*it_r_prob > 0.999){
      *it_r_prob = 0.999;
      *it_nor_x_theta = log_small + log_ntheta;
    }

    it_r_prob++;
    it_nor_x_theta++;   
    *it_r_prob      = data.r[snp];

    *it_derived_k   = snp;

    snp++;
    num_derived_sites++;
  }
  *it_nor_x_theta = -(*it_r_prob) + log_ntheta;
  *it_r_prob      = 1.0 - exp(-(*it_r_prob));
  //if(*it_r_prob > 0.9999){
  //  *it_r_prob = 1.0;
  //  *it_nor_x_theta = 0.0;
  //}
  if(*it_r_prob > 0.999){
    *it_r_prob = 0.999;
    *it_nor_x_theta = log_small + log_ntheta;
  }

  *it_r_prob++;
  *it_r_prob      = 1.0; //just a technicality

  derived_k.resize(num_derived_sites);
  r_prob.resize(num_derived_sites+1);
  nor_x_theta.resize(num_derived_sites);

  it_derived_k = derived_k.begin() + 1;
  it_r_prob = r_prob.begin();
  it_nor_x_theta = nor_x_theta.begin();

  /////////
  //Allocate memory

  alpha.resize(num_derived_sites, data.N);
  beta.resize(num_derived_sites, data.N);
  topology.resize(num_derived_sites, data.N);
  logscales.resize(num_derived_sites);
  std::fill(logscales.begin(), logscales.end(), 0.0);

  double alpha_sum;
  double beta_sum;

  auto it1_alpha       = alpha.ibegin();
  auto rit1_beta       = beta.irbegin();
  auto rit1_alpha      = alpha.irbegin();
  auto rit1_topology   = topology.irbegin();

  std::vector<double>::iterator it2_alpha;
  std::vector<double>::iterator it2_beta;
  std::vector<double>::iterator it2_alpha_prev;
  std::vector<double>::iterator it2_beta_next;
  std::vector<float>::iterator it2_topology;

  std::vector<double>::iterator it2_alpha_rowbegin;
  std::vector<double>::iterator it2_beta_rowbegin;
  std::vector<char>::const_iterator it2_sequence;
  std::vector<float>::iterator it_logscale = logscales.begin(); 
  std::vector<float>::reverse_iterator rit_logscale = logscales.rbegin();


  /////////
  //Forward Algorithm

  ////
  //SNP 0
  *it_logscale                = logscale_alpha;
  alpha_sum                   = 0.0; 

  it2_sequence                = data.sequence.rowbegin(first_snp);
  seq_k                       = *std::next(it2_sequence, k);

  it2_alpha_rowbegin          = alpha.rowbegin(it1_alpha);
  it2_alpha                   = it2_alpha_rowbegin;
  std::vector<float>::iterator it2_alpha_begin = alpha_begin.vbegin(); 
  //this loop vectorizes
  for(; it2_alpha != alpha.rowend(it1_alpha);){
    *it2_alpha                = *it2_alpha_begin;
    it2_alpha_begin++;
    it2_alpha++;
  }

  it2_alpha                        = it2_alpha_rowbegin;
  *std::next(it2_alpha_rowbegin,k) = 0.0;
  //this loop does not vectorize
  for(; it2_alpha != alpha.rowend(it1_alpha);){
    alpha_sum                += *it2_alpha;
    it2_alpha++;
  }

  ////
  //SNP > 0

  double r_x_alpha_sum;
  if(*it_r_prob < 1.0){
    r_x_alpha_sum = (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone) * alpha_sum;
  }else{
    r_x_alpha_sum = alpha_sum;
  }
  double prev_logscale = *it_logscale;
  for(; it_derived_k != derived_k.end();){

    /////////////////
    //precalculated quantities 
    snp                     = *it_derived_k;
    it2_sequence            = data.sequence.rowbegin(snp);
    seq_k                   = *std::next(it2_sequence, k);


    //////////////////
    //inner loop of forward algorithm

    if(*it_r_prob < 1){

      it_logscale++;
      prev_logscale              += (*it_nor_x_theta);
      *it_logscale                = prev_logscale;
      alpha_sum                   = 0.0; 

      it2_sequence                = data.sequence.rowbegin(snp);
      it2_alpha_prev              = alpha.rowbegin(it1_alpha);

      it1_alpha++;
      it2_alpha_rowbegin          = alpha.rowbegin(it1_alpha);
      it2_alpha                   = it2_alpha_rowbegin;
      //this loop vectorizes
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        *it2_alpha                = *it2_alpha_prev + r_x_alpha_sum;
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha               *= derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha++;
        it2_alpha_prev++;
        it2_sequence++;
      }

      it2_alpha                        = it2_alpha_rowbegin;
      *std::next(it2_alpha_rowbegin,k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        alpha_sum                += (*it2_alpha);
        it2_alpha++;
      } 

    }else{

      it_logscale++;
      prev_logscale              += log((data.ntheta/Nminusone) * alpha_sum);
      *it_logscale                = prev_logscale;
      alpha_sum                   = 0.0;

      it1_alpha++;
      it2_alpha_rowbegin          = alpha.rowbegin(it1_alpha);
      it2_alpha                   = it2_alpha_rowbegin;
      it2_sequence                = data.sequence.rowbegin(snp);
      //this loop vectorizes
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        derived                   = (double) (seq_k > *it2_sequence);
        *it2_alpha                = derived * theta_ratio + 1.0; //trick to get around if statement (slightly inefficient but don't know any alternative)
        it2_alpha++;
        it2_sequence++;
      }

      it2_alpha                        = it2_alpha_rowbegin;
      *std::next(it2_alpha_rowbegin,k) = 0.0;
      //this for loop for some reason does not vectorize
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        alpha_sum                += (*it2_alpha);
        it2_alpha++;
      } 

    }
    r_x_alpha_sum = alpha_sum;

    //check if alpha_sums get too small, if they do, rescale
    if(r_x_alpha_sum < lower_rescaling_threshold || r_x_alpha_sum > upper_rescaling_threshold){
      tmp       = r_x_alpha_sum; 
      it2_alpha = it2_alpha_rowbegin;
      for(; it2_alpha != alpha.rowend(it1_alpha);){
        *it2_alpha /= tmp;
        it2_alpha++;
      }
      prev_logscale         += log(tmp);
      *it_logscale          += log(tmp);
      r_x_alpha_sum          = 1.0;
      assert(*it_logscale < std::numeric_limits<double>::infinity());
    } 

    it_r_prob++; 
    if(*it_r_prob < 1.0){
      r_x_alpha_sum *= (*it_r_prob)/((1.0 - (*it_r_prob)) * Nminusone);
    }

    it_derived_k++;
    it_nor_x_theta++;
  }

  ///////////////
  //Backward algorithm

  /////
  //SNP L-1
  it_derived_k--;
  assert(last_snp == *it_derived_k);

  *it_logscale                    += logscale_beta;
  beta_sum                         = 0.0;

  it2_sequence                     = data.sequence.rowbegin(last_snp);
  seq_k                            = *std::next(it2_sequence, k);

  it2_beta_rowbegin                = beta.rowbegin(rit1_beta);
  it2_beta                         = it2_beta_rowbegin;
  std::vector<float>::iterator it2_beta_end = beta_end.vbegin();
  //this loop vectorizes
  for(; it2_beta != beta.rowend(rit1_beta);){
    *it2_beta                      = *it2_beta_end;
    it2_beta_end++;
    it2_beta++;
  }

  it2_beta                         = it2_beta_rowbegin;
  *std::next(it2_beta_rowbegin, k)       = 0.0;
  //this loop does not vectorize
  for(; it2_beta != beta.rowend(rit1_beta);){      
    if(seq_k > *it2_sequence){
      beta_sum                    += data.theta * (*it2_beta); 
    }else{
      beta_sum                    += data.ntheta * (*it2_beta); 
    }
    it2_sequence++;
    it2_beta++;
  }

  //calculate topology matrix
  it2_topology                     = topology.rowbegin(rit1_topology);
  it2_alpha                        = alpha.rowbegin(rit1_alpha);
  it2_beta                         = it2_beta_rowbegin;
  for(; it2_topology != topology.rowend(rit1_topology);){ 
    //*it2_topology                   = std::log(*it2_alpha) + std::log(*it2_beta);
    *it2_topology                  = (*it2_alpha) * (*it2_beta);
    //std::cout << *it2_topology << " ";
    it2_topology++;
    it2_alpha++;
    it2_beta++;
  }
  //std::cout << std::endl;

  rit1_alpha++;
  rit1_topology++;

  /////
  //SNP < L-1
  double beta_sum_oneminustheta, beta_sum_theta;
  double r_x_beta_sum;
  if(*it_r_prob < 1.0){
    r_x_beta_sum = (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone) * beta_sum;
  }else{
    r_x_beta_sum = beta_sum;
  } 
  snp_next = last_snp;
  prev_logscale = logscale_beta;
  while(rit1_topology != topology.irend()){

    it_derived_k--;
    snp          = *it_derived_k;

    if(*it_r_prob < 1.0){

      //inner loop of backwards algorithm
      it_logscale--;

      prev_logscale                         += *it_nor_x_theta;
      *it_logscale                          += prev_logscale;  
      beta_sum                               = 0.0;
      beta_sum_oneminustheta                 = r_x_beta_sum / data.ntheta;
      beta_sum_theta                         = r_x_beta_sum / data.theta - beta_sum_oneminustheta;
      it2_sequence                           = data.sequence.rowbegin(snp_next);
      it2_beta_next                          = beta.rowbegin(rit1_beta);

      rit1_beta++;
      it2_beta_rowbegin                      = beta.rowbegin(rit1_beta);
      it2_beta                               = it2_beta_rowbegin;
      //this loop vectorizes
      for(; it2_beta != beta.rowend(rit1_beta);){
        derived                              = (double) (seq_k > *it2_sequence);
        *it2_beta                            = *it2_beta_next + derived * beta_sum_theta + beta_sum_oneminustheta;
        *it2_beta                           *= derived * theta_ratio + 1.0; 
        it2_sequence++;
        it2_beta++;
        it2_beta_next++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta                               = it2_beta_rowbegin;
      *std::next(it2_beta_rowbegin, k)       = 0.0;
      //this loop does not vectorize
      for(; it2_beta != beta.rowend(rit1_beta);){      
        if(seq_k > *it2_sequence){
          beta_sum                          += data.theta * (*it2_beta); 
        }else{
          beta_sum                          += data.ntheta * (*it2_beta); 
        }
        it2_sequence++;
        it2_beta++;
      }

    }else{

      it_logscale--;
      prev_logscale                         += log(beta_sum/Nminusone);
      *it_logscale                          += prev_logscale;  
      beta_sum                               = 0.0; 

      rit1_beta++;
      it2_beta_rowbegin                      = beta.rowbegin(rit1_beta);
      it2_beta                               = it2_beta_rowbegin;
      //this loop vectorizes
      for(; it2_beta != beta.rowend(rit1_beta);){
        *it2_beta                            = 1.0; 
        it2_beta++;
      }    

      it2_sequence                           = data.sequence.rowbegin(snp);
      seq_k                                  = *std::next(it2_sequence, k);
      it2_beta                               = it2_beta_rowbegin;
      *std::next(it2_beta_rowbegin, k)       = 0.0;
      //this loop does not vectorize
      for(; it2_beta != beta.rowend(rit1_beta);){      
        if(seq_k > *it2_sequence){
          beta_sum                          += data.theta; 
        }else{  
          beta_sum                          += data.ntheta; 
        }
        it2_sequence++;
        it2_beta++;
      }

    }

    r_x_beta_sum = beta_sum;

    //calculate topology matrix
    it2_topology                           = topology.rowbegin(rit1_topology);
    it2_alpha                              = alpha.rowbegin(rit1_alpha);
    it2_beta                               = it2_beta_rowbegin;
    for(; it2_topology != topology.rowend(rit1_topology);){  
      //*it2_topology                        = std::log(*it2_alpha) + std::log(*it2_beta);
      *it2_topology                        = (*it2_alpha) * (*it2_beta);
      //std::cout << *it2_topology << " ";
      it2_topology++;
      it2_alpha++;
      it2_beta++;
    }
    //std::cout << std::endl;

    if(r_x_beta_sum < lower_rescaling_threshold || r_x_beta_sum > upper_rescaling_threshold){
      //if they get too small, rescale
      tmp        = r_x_beta_sum;
      it2_beta   = it2_beta_rowbegin;
      //r_x_beta_sum = 0.0;
      for(; it2_beta < beta.rowend(rit1_beta);){
        *it2_beta /= tmp;
        //r_x_beta_sum += *it2_beta;
        it2_beta++;
      }
      prev_logscale            += log(tmp);
      *it_logscale             += log(tmp);
      r_x_beta_sum              = 1.0;
      assert(*it_logscale < std::numeric_limits<double>::infinity());
    }

    it_r_prob--;
    if(*it_r_prob < 1.0){
      r_x_beta_sum *= (*it_r_prob)/((1.0 - (*it_r_prob))*Nminusone);
    }

    snp_next = snp;
    //I want these to be pointing at snp_next
    it_nor_x_theta--;
    rit1_alpha++;
    rit1_topology++;
  }

  assert(*it_derived_k == first_snp);


  /*
  int size = std::max(1e5, num_derived_sites * data.N * 0.1);

  int total = unique.size();
  unique.resize(total + size);
  times.resize(total + size);

  auto it1_topology = topology.ibegin();
  auto it_unique = std::next(unique.begin(), total);
  auto it_times = std::next(times.begin(), total);

  float current_unique;
  for(; it1_topology != topology.iend(); it1_topology++){ 
  
    it2_topology = topology.rowbegin(it1_topology);
    current_unique = *it2_topology;
    *it_unique = current_unique;
    *it_times = 1;
    it2_topology++;
    for(; it2_topology != topology.rowend(it1_topology); it2_topology++){
      if( abs( current_unique - *it2_topology ) < 0.01 * std::min(current_unique,*it2_topology) ){
      //if( current_unique == *it2_topology ){
        *it_times += 1;
      }else{
        total++;
        it_unique++;
        it_times++;
        if(it_unique == unique.end()){
          unique.resize(unique.size() + size);
          times.resize(unique.size() + size);
          it_unique = std::next(unique.begin(), total);
          it_times  = std::next(times.begin(), total);
        }
        *it_unique = *it2_topology;
        current_unique = *it_unique;
        *it_times = 1;
      }
    }
    total++;
    it_unique++;
    it_times++;
    if(it_unique == unique.end()){
      unique.resize(unique.size() + size);
      times.resize(unique.size() + size);
      it_unique = std::next(unique.begin(), total);
      it_times  = std::next(times.begin(), total);
    }

  }

  unique.resize(total);
  times.resize(total);
  */

  //debug
  /*
     if(snp == 1422){
     for(int l = 0; l < num_derived_sites; l++){
     for(int n = 0; n < data.N; n++){
     if(n!=k){
     if( topology[l][n] == std::numeric_limits<double>::infinity() || std::isnan(topology[l][n]) ){
     std::cerr << k << " " << l << " " << n << " " << topology[l][n] << " " << alpha[l][n] << " " << beta[l][n] << std::endl;
     }
     }
     }
     }
     }
     */

}




