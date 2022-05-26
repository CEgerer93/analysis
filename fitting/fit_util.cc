/*
  Define classes/structs/methods needed to handle matrix element extraction
  From adat based correlators; fitting with gsl
*/
#include "fit_util.h"
#include "cov_utils.h"
#include "varpro.h"
// #include "summation.h"

#include<iostream>
#include<gsl/gsl_permutation.h> // permutation header for matrix inversions
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h> // linear algebra
#include<gsl/gsl_multimin.h>
#include<gsl/gsl_multifit_nlinear.h>


namespace NFIT
{

  struct bundle_t
  {
    NCOR::dat_t jkDat;
    NCOR::fitFunc_t fit;

    std::string comp;
    gsl_matrix covInv;

    bool VARPRO;

    // Parameterized
    bundle_t(NCOR::dat_t _j, NCOR::fitFunc_t _f, std::string &_comp, gsl_matrix _cv, bool _vpBit)
    {
      jkDat = _j; fit = _f; comp = _comp; covInv = _cv; VARPRO = _vpBit;
    };

    // // Destructor
    // ~bundle_t()
    // {
    //   for ( auto a = fit.fitCov.dat.begin(); a != fit.fitCov.dat.end(); ++a )
    // 	gsl_matrix_free(a->second);
    //   for ( auto a = fit.fitCov.inv.begin(); a != fit.fitCov.inv.end(); ++a )
    // 	gsl_matrix_free(a->second);

    //   gsl_matrix_free(&covInv);
    // }
  };

  /*
    Evaluate correlated chi2
  */
  double chi2(const gsl_vector * x, void *data)
  {
    
    // Get a pointer to the void data structure
    bundle_t * b = (bundle_t *)data;


    // Covenience
    int _min, _max, _step;
    _min = b->fit.theFit.range.min;
    _max = b->fit.theFit.range.max;
    _step = b->fit.theFit.range.step;
    std::vector<int> subDomain = b->fit.theFit.range.makeDomain();



    // Initialize any vectors gsl_vectors need in chi2 evaluation (including in VarPro format)
    std::vector<double> predict;
    gsl_vector *iDiffVec = gsl_vector_calloc(subDomain.size());
    gsl_vector *jDiffVec = gsl_vector_calloc(subDomain.size());
    gsl_vector *dataVec  = gsl_vector_calloc(subDomain.size()); // vector of data to be fit
    // Initialize any constants needed in evaluation
    double chi2(0.0);
    int buff = ( _min - b->jkDat.T[0]) / ( b->jkDat.T[1] - b->jkDat.T[0] );
    
    /*
      VarPro doesn't need explicit diff btwn data & fit!
    */
    if ( ! b->VARPRO )
      {
	for ( auto t = subDomain.begin(); t != subDomain.end(); ++t )
	  // Times are set correctly, and predict vector set correctly
	  predict.push_back( b->fit.func(*t, x) );
      }
    //--------------------------------------------------------------

    // More DEBUGS
#if 0
    std::cout << "Predict = ";
    for ( auto p = predict.begin(); p != predict.end(); ++p )
      std::cout << *p;
    std::cout << "\n";
#endif

    // Populate dataVec & iDiffVec vectors
    for ( int l = 0; l < subDomain.size(); ++l )
      {
	if ( b->comp == "real" )
	  {
	    gsl_vector_set(dataVec, l, b->jkDat.avg[buff+l].real());
	    if ( ! b->VARPRO )
	      gsl_vector_set(iDiffVec, l, gsl_vector_get(dataVec,l) - predict[l]);
	  }
	if ( b->comp == "imag" )
	  {
	    gsl_vector_set(dataVec, l, b->jkDat.avg[buff+l].imag());
	    if ( ! b->VARPRO )
	      gsl_vector_set(iDiffVec, l, gsl_vector_get(dataVec,l) - predict[l]);
	  }
      }

    // More DEBUGS
#if 0
    std::cout << "DATAVEC = " << dataVec << std::endl;
#endif


    /*
      If VARPRO bit is set, compute chi2 using variable projection
    */
    if ( b->VARPRO )
      {
	// Create a VarPro instance
	VarPro::varPro vp(x->size, subDomain.size(), b->fit.theFit.bayesianFit);
	vp.makeBasis(x, subDomain);
	vp.makeY(dataVec, &(b->covInv), b->fit.theFit.priors.prior, b->fit.theFit.priors.width);
	vp.makePhi(&(b->covInv), b->fit.theFit.priors.width);
	vp.getInvPhi();

	// Collect results of data vecs sandwiched btwn inv of data cov
	// & result of varPro mat/vec ops
	double dataSum(0.0), varProSum(0.0);

	// Identity
	gsl_matrix *id = gsl_matrix_alloc(vp.rank,vp.rank);
	gsl_matrix_set_identity(id);

	/*
	  Determine VarPro solution for chi^2
	*/
	gsl_matrix *product = gsl_matrix_alloc(vp.invPhi->size1,vp.invPhi->size2);
	gsl_matrix_memcpy(product, vp.invPhi);

	gsl_blas_dgemm(CblasTrans,CblasNoTrans,-1.0,vp.invPhi,id,2.0,product);
	gsl_vector *rightMult = gsl_vector_alloc(vp.Y->size);
	gsl_blas_dgemv(CblasNoTrans,1.0,product,vp.Y,0.0,rightMult);
	gsl_blas_ddot(vp.Y,rightMult,&varProSum);
	varProSum *= -1;

	// If doing a Bayesian fit, don't forget to include constant contributions
	// ...coming from priors/widths of LINEAR PARAMS!
	// ...as well as standard priors/widths of NON-LINEAR PARAMS!
	if ( b->fit.theFit.bayesianFit )
	  {
	    for ( auto it = b->fit.theFit.priors.prior.begin();
		  it != b->fit.theFit.priors.prior.end(); ++it )
	      {
		int idx = std::distance(b->fit.theFit.priors.prior.begin(),it);
		
		// Standard priors of non-linear params...
		if ( idx < x->size )
		  varProSum += pow(gsl_vector_get(x,idx)-(*it),2)/pow(b->fit.theFit.priors.width[idx],2);
		// constant contributions of linear params...
		else
		  varProSum += pow( (*it)/b->fit.theFit.priors.width[idx],2 );
	      }
	  }

    
	/*
	  Now compute (data)^T x Cov^-1 x (data)
	*/
	gsl_vector *dataRMult = gsl_vector_alloc(b->covInv.size1);
	gsl_blas_dgemv(CblasNoTrans,1.0,&(b->covInv),dataVec,0.0,dataRMult);
	gsl_blas_ddot(dataVec,dataRMult,&dataSum);

	// Final Chi2 from VarPro
	chi2 = dataSum + varProSum;

	// Free memory
	gsl_vector_free(dataRMult);
	gsl_vector_free(rightMult);
	gsl_matrix_free(product);
	gsl_matrix_free(id);
      }
    else
      {
#if 0
	// Some more debugs
	std::cout << "Init chi2 = " << chi2 << std::endl;
	std::cout << "IDIFFVEC = " << iDiffVec << std::endl;
	std::cout << "CovInv = ";
	LinAlg::printMat(&(b->covInv));
#endif
	// Copy the difference vector
	gsl_vector_memcpy(jDiffVec, iDiffVec);

	// Initialize cov^-1 right multiplying jDiffVec
	gsl_vector *invCovRightMult = gsl_vector_alloc(subDomain.size());
	gsl_blas_dgemv(CblasNoTrans,1.0,&(b->covInv),jDiffVec,0.0,invCovRightMult);

	// Form the scalar dot product of iDiffVec & result of invCov x jDiffVec
	gsl_blas_ddot(iDiffVec,invCovRightMult,&chi2);

	// Since VarPro wasn't used, add Gaussian priors if they exist
	if ( b->fit.theFit.bayesianFit )
	  {
	    for ( auto it = b->fit.theFit.priors.prior.begin();
		  it != b->fit.theFit.priors.prior.end(); ++it )
	      {
		int idx = std::distance(b->fit.theFit.priors.prior.begin(), it);
		
		chi2 += pow( (gsl_vector_get(x,idx)-(*it))/b->fit.theFit.priors.width[idx], 2);
	      }
	  }
	      

	// Free some memory
	gsl_vector_free(iDiffVec);
	gsl_vector_free(jDiffVec);
	gsl_vector_free(dataVec);
	gsl_vector_free(invCovRightMult);
      } // VARPRO
    //-------------------------------------------------------------------------------------------

    return chi2;
  } // chi2

  /*
    Big Driver for Fits
  */
  void driver(NCOR::correlator *c, std::string &comp, bool varPro = false)
  {
    int numP;
    // If we are performing fit w/ varpro, drop c->fit.num to # non-lin. params
    if ( varPro )
      numP = c->fit.theFit.getNumNonLin();
    else
      numP = c->fit.num;
    //-------------------------------------------------------------------------

    // Set starting values & step sizes of each fitted parameter
    gsl_vector *ini = gsl_vector_calloc(numP);
    gsl_vector *step = gsl_vector_alloc(numP);
    for ( int i = 0; i < step->size; ++i )
      {
	gsl_vector_set(ini, i, c->fit.theFit.initFitParams.start[i]); 
	gsl_vector_set(step,i, c->fit.theFit.initFitParams.step[i]);
      }

    // Initialize solver
    const gsl_multimin_fminimizer_type *mini = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *fmin = gsl_multimin_fminimizer_alloc(mini,numP);
    
    // Fit each jackknife sample
    for ( std::vector<NCOR::dat_t>::iterator it = c->jack.begin(); it != c->jack.end(); ++it )
      {

	bundle_t * jkFit = new bundle_t(*it, c->fit, comp, *(c->fit.fitCov.inv[comp]), varPro);

	
	
	// Define the gsl_multimin_function
	gsl_multimin_function F;
	// Set properties of gsl_multimin_function
	F.n = numP;       // Dimension - num fitted params
	F.f = &chi2;      // Function to minimize
	F.params = jkFit; // The data - this is a dat_t reference
	
	// Establish initial state for minimizer
	int stat = gsl_multimin_fminimizer_set(fmin,&F,ini,step);
	
	
	// Tolerance, etc
	int k = 1;
	double tolerance = c->fit.theFit.initFitParams.tolerance;
	int maxIters     = c->fit.theFit.initFitParams.maxIters;
	
	// Iterate
	while ( gsl_multimin_test_size( gsl_multimin_fminimizer_size(fmin), tolerance) == GSL_CONTINUE )
	  {
	    // End after maxIters
	    if ( k > maxIters ) { break; }

	    gsl_multimin_fminimizer_iterate(fmin);
	    k++;
	  }
	
	// Best fit params for this iteration
	gsl_vector *fminBest = gsl_vector_alloc(numP);
	// Vector to contain non-linear params & constants determined by varpro
	gsl_vector *finParams = gsl_vector_calloc(c->fit.num);

	gsl_vector_memcpy(fminBest,gsl_multimin_fminimizer_x(fmin));


	/*
	  If varpro has been used in optimization, form a fresh varPro instance
	  to determine vector of constants
	*/
	if ( varPro )
	  {
	    // Must pack fminBest into finParams to print non-linear/constant params of fit
	    double dum;
	    for ( size_t b = 0; b < fminBest->size; ++b )
	      {
		dum = gsl_vector_get(fminBest,b);
		gsl_vector_set(finParams,b, dum);
	      }

	    // Remake data vec - this is ugly! Can I make this better??
	    gsl_vector *dataVec  = gsl_vector_alloc(c->fit.theFit.range.numT());
	    int buff = ( c->fit.theFit.range.min - it->T[0]) / ( it->T[1] - it->T[0] );
	    for ( int l = 0; l < c->fit.theFit.range.numT(); ++l )
	      {
		if ( comp == "real" )
		  gsl_vector_set(dataVec, l, it->avg[buff+l].real());
		if ( comp == "imag" )
		  gsl_vector_set(dataVec, l, it->avg[buff+l].imag());
	      }

	    // Here is the varPro instance - to get solution vec of constants
	    VarPro::varPro soln(numP,c->fit.theFit.range.numT(),c->fit.theFit.bayesianFit);
	    soln.makeBasis(fminBest, c->fit.theFit.range.makeDomain());
	    soln.makeY(dataVec, c->fit.fitCov.inv[comp], c->fit.theFit.priors.prior, c->fit.theFit.priors.width);
	    soln.makePhi(c->fit.fitCov.inv[comp], c->fit.theFit.priors.width);
	    soln.getInvPhi();
	    soln.getSoln();


	    // Append constants determined via VarPro onto finParams
	    for ( size_t b = fminBest->size; b < finParams->size; ++b )
	      {
		dum = gsl_vector_get(soln.soln, b - fminBest->size);
		gsl_vector_set(finParams,b,dum);
	      }

	    gsl_vector_free(dataVec);
	  }
	else
	  {
	    // If we haven't used VarPro, then fitParams is just fminBest
	    gsl_vector_memcpy(finParams, fminBest);
	  }



	// Access the cost
	double cost = gsl_multimin_fminimizer_minimum(fmin);
	int ndof    = c->fit.theFit.range.numT() - c->fit.num - c->fit.fitCov.svs[comp];

	double ellSq = cost;
	double chiSq = cost;

	
	if ( c->fit.theFit.bayesianFit )
	  {
	    // Compute the L^2 value
	    ellSq -= 2*log(1.0/sqrt(2*M_PI));
	    for ( auto w = c->fit.theFit.priors.width.begin(); w != c->fit.theFit.priors.width.end(); ++w )
	      ellSq -= 2*log(1.0/(sqrt(2*M_PI)*(*w)));
	    
	    // Compute the true Chi^2 value (absent any potential priors)
	    for ( auto w = c->fit.theFit.priors.width.begin(); w != c->fit.theFit.priors.width.end(); ++w )
	      {
		int idx = std::distance(c->fit.theFit.priors.width.begin(),w);
		chiSq -= pow( (gsl_vector_get(finParams, idx) - c->fit.theFit.priors.prior[idx])/(*w), 2);
	      }
	  }


	// Return the associated reduced L^2, Chi^2
	ellSq /= ndof; chiSq /= ndof;

	std::cout << "(L2,Chi2,ndof) w/ finparams = "
		  << "(" << ellSq << "," << chiSq << "," << ndof << ")" << std::endl;



	
	// Print the results of this fit - map each value to correct char param
	std::map<std::string, double> aFitMap = jkFit->fit.printFit(finParams);
	
	
	for ( auto it = aFitMap.begin(); it != aFitMap.end(); ++it )
	  c->res.params[it->first].push_back(it->second);

	// Whether Chi2 or L^2 are reported
	if ( c->fit.theFit.bayesianFit )
	  c->res.chi2.push_back(ellSq);
	else
	  c->res.chi2.push_back(chiSq);


	gsl_vector_free(fminBest);
	gsl_vector_free(finParams);
	
	// gsl_vector_free(best);
	delete jkFit;
      } // iterator it (jackknife)

    gsl_vector_free(ini);
    gsl_vector_free(step);
    gsl_multimin_fminimizer_free(fmin);

  } // driver


}

