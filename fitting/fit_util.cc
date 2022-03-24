/*
  Define classes/structs/methods needed to handle matrix element extraction
 65;6003;1c From adat based correlators; fitting with gsl
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
	// // Grab all the fitted parameters
	// std::vector<double> nonLinParam = b->fit.getNonLinParam(x);

#warning "FIX ME! - Priors and widths set to zero in VarPro"
	// Priors/Widths
	std::vector<double> prior(x->size,0.0), width(x->size,0.4);
	prior[0] = 0.5; prior[1] = 1.5;
	// std::vector<double> prior(nonLinParam.size(),0.0), width(nonLinParam.size(),1.0);

	// Create a VarPro instance
	VarPro::varPro vp(x->size, subDomain.size());
	vp.makeBasis(x, subDomain);
	vp.makeY(dataVec, &(b->covInv), prior, width);
	vp.makePhi(&(b->covInv), prior);
	vp.getInvPhi();
#if 0
	std::cout << "Tmp solution from varpro" << std::endl;
	vp.getSoln();
	std::cout << vp.soln << std::endl;
#endif
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

	for ( int l = 0; l < vp.rank; ++l )
	  varProSum += pow(prior[l]/width[l],2);
    
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


    // Collect priors & widths
#warning "FIX ME! - Priors & Widths hard coded"
    std::vector<double> prior(numP,0.0), width(numP,0.4);
    prior[0] = 0.5; prior[1] = 1.5;

    gsl_vector *ini = gsl_vector_calloc(numP);
    gsl_vector *step = gsl_vector_alloc(numP);
    for ( int i = 0; i < step->size; ++i )
      {
	// gsl_vector_set(ini, i, 0.3*(i+1));
	gsl_vector_set(ini, i, prior[i]);
	gsl_vector_set(step,i,0.15);
      }
    // gsl_vector_set(ini,0,0.5);
    // gsl_vector_set(ini,3,2);


    // Initialize solver
    const gsl_multimin_fminimizer_type *mini = gsl_multimin_fminimizer_nmsimplex2rand;
    gsl_multimin_fminimizer *fmin = gsl_multimin_fminimizer_alloc(mini,numP);


    
    // Fit each jackknife sample
    for ( std::vector<NCOR::dat_t>::iterator it = c->jack.begin(); it != c->jack.end(); ++it )
      {

	bundle_t * jkFit = new bundle_t(*it, c->fit, comp, *(c->fit.fitCov.inv[comp]), varPro);

	
	
	// Define the gsl_multimin_function
	gsl_multimin_function F;
	// Dimension
	F.n = numP; // num fitted params
	// Function to minimize
	F.f = &chi2;
	F.params = jkFit; // this is a dat_t reference
	
	// Establish initial state for minimizer
	int stat = gsl_multimin_fminimizer_set(fmin,&F,ini,step);
	
	
	// Tolerance, etc
	int k = 1;
	double tolerance = 0.0000001;
	int maxIters = 10000;
	
	// Iterate
	while ( gsl_multimin_test_size( gsl_multimin_fminimizer_size(fmin), tolerance) == GSL_CONTINUE )
	  {
	    // End after maxIters
	    if ( k > maxIters ) { break; }
	    
	    // Iterate
	    gsl_multimin_fminimizer_iterate(fmin);
	    k++;
	  }
	
	// Best fit params for this iteration
	gsl_vector *fminBest = gsl_vector_alloc(numP);
	// Vector to contain non-linear params & constants determined by varpro
	gsl_vector *finParams = gsl_vector_calloc(c->fit.num);
	// gsl_vector_set(finParams,0,*&(gsl_vector_get(fminBest,0)));
	// gsl_vector_set(finParams,finParams->size-1,gsl_vector_get(fminBest,fminBest->size-1));

	gsl_vector_memcpy(fminBest,gsl_multimin_fminimizer_x(fmin));
	// Return the associated reduced chi2	
	double chiSq = gsl_multimin_fminimizer_minimum(fmin);
	chiSq /= ( c->fit.theFit.range.numT() - c->fit.num - c->fit.fitCov.svs[comp] );
	std::cout << "Chi2 w/ finparams = " << chiSq << std::endl;

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
	    VarPro::varPro soln(numP,c->fit.theFit.range.numT());
	    soln.makeBasis(fminBest, c->fit.theFit.range.makeDomain());
	    soln.makeY(dataVec, c->fit.fitCov.inv[comp], prior, width);
	    soln.makePhi(c->fit.fitCov.inv[comp], prior);
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

	
	// Print the results of this fit - map each value to correct char param
	std::map<std::string, double> aFitMap = jkFit->fit.printFit(finParams);
	
	
	for ( auto it = aFitMap.begin(); it != aFitMap.end(); ++it )
	  c->res.params[it->first].push_back(it->second);
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

