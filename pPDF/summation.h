/*
  SUPPORT FOR CONVERTING INPUT FILEDB DATA INTO PERSONAL CORRELATORS TYPE
*/
#ifndef __summation_h__
#define __summation_h__

/* #include<vector> */
#include "personal_fns.h"
#include "fit_util.h"
#include<stdlib.h>
#include<map>
#include<vector>

#include<gsl/gsl_matrix.h>


namespace Summation
{
  // Big class to hold correlation function information
  class corrFunc
  {
  public:
    // Containers of the actual correlation functions
    correlators    data;
    correlators    dataJk;
    correlators    dataJkEnsAvg,dataSummedJkEnsAvg;
    
    
    // Constructors
    corrFunc() {}; // default 
    corrFunc(int g,int t) // parametrized
      {
	gauge_configs = g; num_tslices = t;
      };
    /* // Destructor */
    /* virtual ~corrFunc(); */
    
    void initialize();
    void initializeJks();
    void initializeJkEnsAvgs();
    void initializeSummedJkEnsAvgs();
    
    int getT() { return num_tslices; }
    int getCfgs() { return gauge_configs; }

  private:
    int gauge_configs, num_tslices;
  };


  class Ratios
  {
  public:

    struct NtCorr_t
    {
      double avgR, avgI, errR, errI;
      correlators tData;
    NtCorr_t(int cfgs) : tData(cfgs) {};
    };

    std::map<int, NtCorr_t> ratio;

    std::map<int, FitRes_t> ratioFit;

    struct covariances
    {
      gsl_matrix * cov, * inv;
      int svsR, svsI; // singular values removed in inverse computation
    } covR, covI;

    int tmin, tstep, tmax;
    std::vector<int> tseries;
    

    // constructors
    Ratios() {};
  Ratios(int cfgs, int _tmin, int _tstep, int _tmax) : tmin(_tmin), tstep(_tstep), tmax(_tmax)
    {
      /* tmin = _tmin; tstep = _tstep; tmax = _tmax; */
      for ( int t = tmin; t <= tmax; t+=tstep )
	{
	  tseries.push_back( t );
	  NtCorr_t dum(cfgs);
	  std::pair<int, NtCorr_t> tsepR(t, dum);
	  ratio.insert(tsepR);
	}
    };

    void makeSR(std::vector<corrFunc> &_3ptFuncs, corrFunc &_2ptFunc)
    {
      for ( auto r = ratio.begin(); r != ratio.end(); ++r )
	{
	  int t_idx = (r->first - tmin)/tstep;
	  summedRatio( _2ptFunc.dataJkEnsAvg, _3ptFuncs[t_idx].dataSummedJkEnsAvg, r->second.tData, r->first );
	}
    }


    // Get a central value estimate of summed ratios
    void mean();
    // Determine the data covariance matrices
    void makeCovs();
    // Determine the inverses of data covariances
    void makeInvCovs();

    // Perform a fit to jackknife samples
    void fit(std::string &s);
    
  }; // Ratios


} // Summation
#endif
