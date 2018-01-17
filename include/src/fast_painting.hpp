#ifndef FAST_PAINTING_HPP
#define FAST_PAINTING_HPP

#include "fast_log.hpp"
#include "data.hpp"
#include "collapsed_matrix.hpp"

#include <iostream>
#include <list>
#include <vector>
#include <limits>
#include <tgmath.h>

#include <cassert>

class FastPainting{

  private:

    double lower_rescaling_threshold, upper_rescaling_threshold;
    double Nminusone, prior_theta, prior_ntheta, theta_ratio, log_ntheta, log_small;
    double normalizing_constant, distance_rescaling_constant;  

  public:

    FastPainting(Data& data){

      lower_rescaling_threshold = 1e-10;
      upper_rescaling_threshold = 1.0/lower_rescaling_threshold;

      assert(data.theta < 1.0);
      Nminusone = data.N - 1.0;
      prior_theta  = data.theta/Nminusone - data.ntheta/Nminusone;
      prior_ntheta = data.ntheta/Nminusone;
      theta_ratio = data.theta/(1.0 - data.theta) - 1.0;
      log_ntheta = log(data.ntheta);
      log_small  = log(0.01);

      //normalizing_constant = (double) log(Nminusone) - data.L * log_ntheta;
      //distance_rescaling_constant = log(data.theta/data.ntheta);  

    }

    /* 
    void Paint(const Data& data, CollapsedMatrix<float>& topology, std::vector<float>& logscales, const int k);
    void PaintToFile(const Data& data, CollapsedMatrix<float>& topology, std::vector<float>& logscales, CollapsedMatrix<int>& boundarySNP, const int k);
    void PaintAll(const Data& data, CollapsedMatrix<float>& topology, std::vector<float>& logscales, const int k);
    void PaintSteppingStones2(const Data& data, std::vector<FILE*> pfiles, const int k);
    */
    void PaintSteppingStones(const Data& data, std::vector<int>& window_boundaries, std::vector<FILE*> pfiles, const int k);
    void RePaintSection(const Data& data, CollapsedMatrix<float>& topology, std::vector<float>& logscales, CollapsedMatrix<float>& alpha_begin, CollapsedMatrix<float>& beta_end, int boundarySNP_begin, int boundarySNP_end, float logscale_alpha, float logscale_beta, const int k);

};

#endif //FAST_PAINTING_HPP
