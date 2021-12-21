/*
  Define classes/structs/methods needed for variable projection
*/
#include "varpro.h"
#include "cov_utils.h"
#include <gsl/gsl_blas.h>

namespace VarPro
{
  /*
    Populate the non-linear basis of functions
  */
  void varPro::makeBasis(const gsl_vector * nlParams, const std::vector<int> &T)
    // void varPro::makeBasis(std::vector<double> &nlParams, const std::vector<int> &T)
  {
    for ( int l = 0; l < rank; ++l )
      {
	// Evaluate the l^th basis function at each time slice
	for ( auto t = T.begin(); t != T.end(); ++t )
	  gsl_matrix_set(basis, l, std::distance(T.begin(), t),
			 exp(-gsl_vector_get(nlParams,l)* (*t) ));
      } // l
  }

  /*
    Populate Y Solution vector
  */
  void varPro::makeY(gsl_vector *data, gsl_matrix *invCov,
		     std::vector<double> &prior, std::vector<double> &width)
  {
    for ( int l = 0; l < rank; ++l )
      {
	double sum(0.0), priorsSum(0.0);

	gsl_vector *rMult = gsl_vector_alloc(invCov->size1);

	// Pull out l^th row of basis function matrix
	gsl_vector_view slice = gsl_matrix_row(basis, l);


	// Perform (covinv) x basis(l)
	gsl_blas_dgemv(CblasNoTrans,1.0,invCov,&slice.vector,0.0,rMult);
	// Perform (data)^T \cdot (rMult)
	gsl_blas_ddot(data,rMult,&sum);

	// Collect the contributions from any priors
	priorsSum += prior[l]/pow(width[l],2);

	// Now set the l^th entry of Y
	gsl_vector_set(Y, l, sum+priorsSum);
      }
  }

  /*
    Populate Phi matrix
  */
  void varPro::makePhi(gsl_matrix *invCov, std::vector<double> prior)
  {
    for ( int k = 0; k < rank; ++k )
      {
	// Pull out k^th row of basis function matrix
	gsl_vector_view k_slice = gsl_matrix_row(basis, k);

	for ( int l = 0; l < rank; ++l )
	  {
	    double sum(0.0);
	    // Pull out l^th row of basis function matrix
	    gsl_vector_view l_slice = gsl_matrix_row(basis, l);

	    gsl_vector * rMult = gsl_vector_alloc(invCov->size1);
	    // Perform (covinv) x basis(l)
	    gsl_blas_dgemv(CblasNoTrans,1.0,invCov,&l_slice.vector,0.0,rMult);
	    // Perform basis(k)^T \cdot (rMult)
	    gsl_blas_ddot(&k_slice.vector,rMult,&sum);


	    // Include contributions from priors!!!!!! - only appear on diagonal of Phi matrix
#warning "Priors off in Phi construction!"
#if 0
	    if ( k == l )
	      sum += 1.0/pow(prior[k],2);
#endif

	    // Insert this value into the Phi matrix
	    gsl_matrix_set(Phi, k, l, sum);

	  } // l
      } // k
  }

  /*
    Get the inverse of Phi matrix
  */
  void varPro::getInvPhi()
  {
    int catchSVs = LinAlg::matrixInv(Phi, invPhi);
  }

  /*
    Solution of VarPro - Soln = \Phi^{-1} x Y
  */
  void varPro::getSoln()
  {
    gsl_blas_dgemv(CblasNoTrans,1.0,invPhi,Y,0.0,soln);
  }
}
