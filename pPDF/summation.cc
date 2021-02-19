/*
  Implentation of corrFunc class defined methods
*/
#include "summation.h"
#include "fit_util.h"
// #include<iostream>

/*
  Helpful utilities in adat/lib/hadron/irrep_util.cc

  Array<int> canonicalOrder(const Array<int>& mom)
  std::string shortMom(const Array<int>& mom)
  Array<int> Mom3d(int p1, int p2, int p3) - helper function; make an Array<int> of momenta

  Mult replace each Array<T>& element: *=val
*/

using namespace FIT;

namespace Summation
{
  void corrFunc::initialize()
  {
    data.ncor.resize(1);
    data.ncor[0].real.resize(gauge_configs); data.ncor[0].imag.resize(gauge_configs);
    for ( int j = 0; j < data.ncor[0].real.size(); ++j )
      {
	data.ncor[0].real[j].resize(num_tslices); data.ncor[0].imag[j].resize(num_tslices);
      }
  }

  void corrFunc::initializeJks()
  {
    dataJk.ncor.resize(gauge_configs);
    initializeJkArrays(dataJk,num_tslices); // external call
  }

  void corrFunc::initializeJkEnsAvgs()
  {
    dataJkEnsAvg.ncor.resize(gauge_configs);
  }

  void corrFunc::initializeSummedJkEnsAvgs()
  {
    dataSummedJkEnsAvg.ncor.resize(gauge_configs);
  }


  
  // Get a central value estimate of summed ratios
  void Ratios::mean()
  {
    for ( auto r = ratio.begin(); r != ratio.end(); ++r )
      {
	double _r(0.0), _i(0.0);
	for ( auto g = r->second.tData.ncor.begin(); g < r->second.tData.ncor.end(); g++ )
	  {
	    _r += g->real[0][0];
	    _i += g->imag[0][0];
	  }
	r->second.avgR = _r / r->second.tData.ncor.size();
	r->second.avgI = _i / r->second.tData.ncor.size();

	// _r = 0.0; _i = 0.0;
	// for ( auto g = r->second.tData.begin(); g < 
      }	    
  }

  // Determine the data covariance matrices
  void Ratios::makeCovs()
  {
    covR.cov = gsl_matrix_alloc(ratio.size(), ratio.size());
    covI.cov = gsl_matrix_alloc(ratio.size(), ratio.size());
    for ( auto i = ratio.begin(); i != ratio.end(); ++i )
      {
	for ( auto j = ratio.begin(); j != ratio.end(); ++j )
	  {
	    
	    double _r(0.0), _i(0.0);
	    int cfgs = i->second.tData.ncor.size();
	    for ( int g = 0; g < cfgs; g++ )
	      {
		_r += ( i->second.tData.ncor[g].real[0][0] - i->second.avgR )*
		  ( j->second.tData.ncor[g].real[0][0] - j->second.avgR );
		_i += ( i->second.tData.ncor[g].imag[0][0] - i->second.avgI )*
		  ( j->second.tData.ncor[g].imag[0][0] - j->second.avgI );
	      } // end g

	    // Set covariance entries and proceed
	    gsl_matrix_set(covR.cov, std::distance(ratio.begin(), i),
			   std::distance(ratio.begin(), j), (( cfgs - 1 )/(1.0*cfgs))*_r );
	    gsl_matrix_set(covI.cov, std::distance(ratio.begin(), i),
			   std::distance(ratio.begin(), j), (( cfgs - 1 )/(1.0*cfgs))*_i );

	  } // end j
      } // end i
  }

  // Determine the inverses of data covariances
  void Ratios::makeInvCovs()
  {
    covR.inv = gsl_matrix_alloc(ratio.size(), ratio.size());
    covI.inv = gsl_matrix_alloc(ratio.size(), ratio.size());

    covR.svsR = matrixInv(covR.cov, covR.inv);
    covI.svsI = matrixInv(covI.cov, covI.inv);
  }



} // Summation
