/*
  Implentation of corrFunc class defined methods
*/
#include "summation.h"
#include "fit_util.h"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h> // multidimensional minimization
#include <gsl/gsl_multifit_nlin.h>


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



  // Correct for bias in jackknifing
  void Ratios::antiJkRatio()
  {
    int num_jk_samples = ratio.begin()->second.tData.ncor.size();
    
    for ( auto r = ratio.begin(); r != ratio.end(); ++r )
      {
	for ( auto g = r->second.tData.ncor.begin(); g != r->second.tData.ncor.end(); ++g )
	  {
	    // std::cout << "BEFORE = " << g->real[0][0] << std::endl;
	    // std::cout << "NUM JK = " << num_jk_samples << std::endl;
	    // std::cout << "AVG EST = " << r->second.avgR << std::endl;

	    double tmp = g->real[0][0]; 

	    // std::cout << "TMP = " << tmp << std::endl;
	    // std::cout << "(NJK-1)*tmp = " << ( num_jk_samples - 1 )*tmp << std::endl;
	    // std::cout << "BIAS CORRECTED = " << num_jk_samples*r->second.avgR - ( num_jk_samples - 1 )*tmp
	    // 	      << std::endl;
	    g->real[0][0] = ( num_jk_samples*r->second.avgR - ( num_jk_samples - 1 )*tmp );
	    // std::cout << "AFTER = " << g->real[0][0] << std::endl;
	    // exit(80);
	    tmp = g->imag[0][0];
	    g->imag[0][0] = ( num_jk_samples*r->second.avgI - ( num_jk_samples - 1 )*tmp );
	  }
      }
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

    covR.svs = matrixInv(covR.cov, covR.inv);
    covI.svs = matrixInv(covI.cov, covI.inv);
  }



  /*
    Write out the SR data
  */
  void Ratios::writeSR(std::string out)
  {
    // Stream to a file
    std::ofstream os;
    os.open(out+".dat");
    os << ratio.begin()->second.tData.ncor.size() << " " << ratio.size() << " 1 0 1\n";

    for ( int J = 0; J < ratio.begin()->second.tData.ncor.size(); J++ )
      {
	for ( auto ir = ratio.begin(); ir != ratio.end(); ++ir )
	  {
	    int idx = std::distance(ratio.begin(), ir); // need 0-based index for time series written to file
	    os << std::setprecision(10) << idx << " "
	       << ir->second.tData.ncor[J].real[0][0] << " "
	       << ir->second.tData.ncor[J].imag[0][0] << "\n";
	  }
      }
    os.close();
  }


  
  /*
    HERE FOR FIT STUFF
  */
  // Perform a fit to jackknife samples
  void Ratios::fit(std::string s, std::string outName)
  {
    size_t order;
    gsl_vector *_ini, *_iniSteps;
    if ( s == "LINEAR" )
      {
	order = 2;
	_ini      = gsl_vector_alloc(order);
	_iniSteps = gsl_vector_alloc(order);
	for ( size_t p = 0; p < order; p++ ) { gsl_vector_set(_ini,p,0.4); gsl_vector_set(_iniSteps,p,0.2); }
      }

    std::ofstream out; // the out file stream

    // Fit both the real/imag data
    for ( int COMP = 1; COMP < 3; COMP++ )
      {
	std::string component;
	if ( COMP == 1 )
	  component = "RE";
	if ( COMP == 2 )
	  component = "IM";
	// Open the file to contain fit results
	out.open(outName+"."+component+"_fit.dat");

	// Initialize the solver here
	const gsl_multimin_fminimizer_type *minimizer = gsl_multimin_fminimizer_nmsimplex2rand;
	gsl_multimin_fminimizer * fmin = gsl_multimin_fminimizer_alloc(minimizer,_ini->size);
	
	
        // Fit per jackknife sample
	for ( int J = 1; J < ratio.begin()->second.tData.ncor.size(); J++ )
	  {
	    
	    // Initialize struct to hold this jack's data
	    linFit_t * jfit;
	    std::vector<double> dum;
	    for ( auto it = ratio.begin(); it != ratio.end(); ++it )
	      {
		if ( COMP == 1 )
		  dum.push_back( it->second.tData.ncor[J].real[0][0] );
		if ( COMP == 2 )
		  dum.push_back( it->second.tData.ncor[J].imag[0][0] );
	      }

	    if ( COMP == 1 )
	      jfit = new linFit_t(dum, tseries, covR.inv);
	    if ( COMP == 2 )
	      jfit = new linFit_t(dum, tseries, covI.inv);

	    // for ( auto a = jfit->m.begin(); a != jfit->m.end(); ++a )
	    //   std::cout << *a << " ";
	    // std::cout << "\n";
	    
	    
	    // Define the gsl_multimin_function
	    gsl_multimin_function Chi2;
	    // Dimension of the system
	    Chi2.n = _ini->size;
	    // Function to minimize
	    Chi2.f = &chi2Linear;
	    Chi2.params = jfit;
	    
	    
	    std::cout << "Establishing initial state for minimizer..." << std::endl;
	    int status = gsl_multimin_fminimizer_set(fmin,&Chi2,_ini,_iniSteps);
	    std::cout << "Minimizer established..." << std::endl;
	    
	    // Iteration count
	    int k = 1;
	    double tolerance = 0.00000000000001; // 0.0001
	    int maxIters = 10000;                // 1000
	    
	    while ( gsl_multimin_test_size( gsl_multimin_fminimizer_size(fmin), tolerance) == GSL_CONTINUE )
	      {
		// End after maxIters
		if ( k > maxIters ) { break; }
		// Iterate
		gsl_multimin_fminimizer_iterate(fmin);

		// std::cout << "Current params = " << gsl_vector_get(gsl_multimin_fminimizer_x(fmin),0)
		// 	  << " " << gsl_vector_get(gsl_multimin_fminimizer_x(fmin),1) << std::endl;
		
		k++;
	      }
	    
	    // Return the best fit parameters
	    gsl_vector *bestFitParams = gsl_vector_alloc(_ini->size);
	    bestFitParams = gsl_multimin_fminimizer_x(fmin);
	    
	    // Return the best correlated Chi2
	    double chiSq = gsl_multimin_fminimizer_minimum(fmin);
	    // Determine the reduced chi2
	    double reducedChiSq;
	    if ( COMP == 1 )
	      reducedChiSq = chiSq / (ratio.size() - _ini->size - covR.svs);
	    if ( COMP == 2 )
	      reducedChiSq = chiSq / (ratio.size() - _ini->size - covI.svs);

	    
	    out << std::setprecision(12) << "A " << gsl_vector_get(bestFitParams,0) << "\n"
		<< "B " << gsl_vector_get(bestFitParams,1) << "\n"
		<< "rChi2 " << reducedChiSq << std::endl;

	    // exit(90);

	    delete jfit;

	  } // end J

	out.close();
      } // end COMP
  } // end fit


} // Summation
