/*
  Define classes/structs/methods needed to handle any linear algebra
*/
//#include "fit_util.h"
//#include "summation.h"

#include<iostream>
#include<gsl/gsl_permutation.h> // permutation header for matrix inversions
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h> // linear algebra

#include "cov_utils.h"

namespace LinAlg
{
  // // Define operator to easily print contents of a gsl_vector
  // std::ostream& operator<<(std::ostream& os, const gsl_vector *v)
  // {
  //   if ( v->size > 0 )
  //     {
  // 	os << gsl_vector_get(v,0);
  // 	for ( int i = 1; i < v->size; ++i )
  // 	  os << "," << gsl_vector_get(v,i);
  //     }
  //   return os;
  // }

  // double operator[](const gsl_vector *v, int i)
  // {
  //   return gsl_vector_get(v,i);
  // }

  // Generic gsl_matrix viewer
  void printMat(gsl_matrix *g)
  {
    std::cout << "{";
    for ( size_t i = 0; i < g->size1; i++ ) {
      std::cout << "{";
      for ( size_t j = 0; j < g->size2; j++ )
	std::cout << gsl_matrix_get(g,i,j) << ",";
      std::cout << "},\n";
    }    
  }
  
  // Generic gsl_matrix_complex viewer
  void printMat(gsl_matrix_complex *g)
  {
    std::cout << "{";
    for ( size_t i = 0; i < g->size1; i++ ) { 
      std::cout << "{";
      for ( size_t j = 0; j < g->size2; j++ )
	{
	  std::complex<double> dum(0.0,0.0);
	  gsl_complex gc = gsl_matrix_complex_get(g,i,j);
	  dum.real(gc.dat[0]);
	  dum.imag(gc.dat[1]);
	  std::cout << dum << ", ";
	}
      std::cout << "},\n";
    }
  }


  // Form tensor product of two matrices
  gsl_matrix_complex* tensorProd(gsl_matrix_complex *g, gsl_matrix_complex *h)
  {
    gsl_matrix_complex * t = gsl_matrix_complex_calloc(g->size1 * h->size1, g->size2 * h->size2);

    for ( int gi = 0; gi < g->size1; ++gi )
      {
	for ( int gj = 0; gj < g->size2; ++gj )
	  {
	    for ( int hi = 0; hi < h->size1; ++hi )
	      {
		for ( int hj = 0; hj < h->size2; ++hj )
		  {
		    gsl_complex prod = gsl_complex_mul(gsl_matrix_complex_get(g,gi,gj),gsl_matrix_complex_get(h,hi,hj));
		    gsl_matrix_complex_set(t,gi*h->size1+hi,gj*h->size2+hj,prod);
		  }
	      }

	  }
      }
    return t;
  }


  /*
    PERFORM INVERSION OF PASSED MATRIX - RETURN # OF SVs REMOVED
  */
  int matrixInv(gsl_matrix * M, gsl_matrix * MInv)
  {
#if VERBOSITY>2
    std::cout << "Inside LinAlg::matrixInv" << std::endl;
#endif
    size_t dataDim = M->size1;
    gsl_matrix * toInvert = gsl_matrix_alloc(dataDim,dataDim);
    gsl_matrix_memcpy(toInvert,M); // make a copy of M, as toInvert is modified below
    gsl_matrix * V = gsl_matrix_alloc(dataDim,dataDim);
    gsl_vector *S = gsl_vector_alloc(dataDim);
    gsl_vector *work = gsl_vector_alloc(dataDim);

    /*
      PERFORM THE SINGULAR VALUE DECOMPOSITION OF DATA COVARIANCE MATRIX (A)
      
      A = USV^T
          ---> A is an MxN matrix
    	  ---> S is the singular value matrix (diagonal w/ singular values along diag - descending)
    	  ---> V is returned in an untransposed form
    */
    gsl_linalg_SV_decomp(toInvert,V,S,work); // SVD decomp;  'toInvert' replaced w/ U on return

    // Define an svd cut
    double svdCut = 1e-16;
    // Initialize the inverse of the S diag
    gsl_vector *pseudoSInv = gsl_vector_alloc(dataDim);
    gsl_vector_set_all(pseudoSInv,0.0);

    // Vector of singular values that are larger than specified cut
    std::vector<double> aboveCutVals;
    
#if VERBOSITY>1
    std::cout << "The singular values above SVD Cut = " << svdCut << " are..." << std::endl;
#endif
    for ( int s = 0; s < dataDim; s++ )
      {
    	double dum = gsl_vector_get(S,s);
    	if ( dum >= svdCut )
	  {
	    aboveCutVals.push_back(dum);
#if VERBOSITY>1
	    std::cout << dum << " ";
#endif
	  }
#if VERBOSITY>1
	std::cout << "\n";
#endif
      }
    
    // Assign the inverse of aboveCutVals to the pseudoSInv vector
    for ( std::vector<double>::iterator it = aboveCutVals.begin(); it != aboveCutVals.end(); ++it )
      {
    	gsl_vector_set(pseudoSInv,it-aboveCutVals.begin(),1.0/(*it));
      }

    // Promote this pseudoSInv vector to a matrix, where entries are placed along diagonal
    gsl_matrix * pseudoSInvMat = gsl_matrix_alloc(dataDim,dataDim);
    gsl_matrix_set_zero(pseudoSInvMat);
    for ( int m = 0; m < dataDim; m++ )
      {
    	gsl_matrix_set(pseudoSInvMat,m,m,gsl_vector_get(pseudoSInv,m));
    	// std::cout << gsl_vector_get(pseudoSInv,m) << std::endl;
      }

    /*
      With singular values that are zero, best we can do is construct a pseudo-inverse

      In general, the inverse we are after is
                               VS^(-1)U^T
		     with S the diagonal matrix of singular values
    */
    // In place construct the transpose of U ('toInvert' was modified in place to U above)
    gsl_matrix_transpose(toInvert);

    gsl_matrix * SinvUT = gsl_matrix_alloc(dataDim,dataDim); gsl_matrix_set_zero(SinvUT);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pseudoSInvMat,toInvert,0.0,SinvUT);

    // Now make the inverse of 'toInvert'
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,SinvUT,0.0,MInv);

#if VERBOSITY>2
    std::cout << "Leaving LinAlg::matrixInv" << std::endl;
#endif
    return pseudoSInv->size - aboveCutVals.size();
  }
	       

#if 0
  /*
    PERFORM INVERSION OF PASSED COMPLEX MATRIX - RETURN # OF SVs REMOVED
  */
  int matrixInv(gsl_matrix_complex * M, gsl_matrix_complex * MInv)
  {
    size_t dataDim = M->size1;
    gsl_matrix_complex * toInvert = gsl_matrix_complex_alloc(dataDim,dataDim);
    gsl_matrix_complex_memcpy(toInvert,M); // make a copy of M, as toInvert is modified below
    gsl_matrix_complex * V = gsl_matrix_complex_alloc(dataDim,dataDim);
    gsl_vector_complex *S = gsl_vector_complex_alloc(dataDim);
    gsl_vector_complex *work = gsl_vector_complex_alloc(dataDim);

    /*
      PERFORM THE SINGULAR VALUE DECOMPOSITION OF DATA COVARIANCE MATRIX (A)
      
      A = USV^T
          ---> A is an MxN matrix
          ---> S is the singular value matrix (diagonal w/ singular values along diag - descending)
          ---> V is returned in an untransposed form
    */
    gsl_linalg_SV_decomp(toInvert,V,S,work); // SVD decomp;  'toInvert' replaced w/ U on return
#if 0
    // Define an svd cut
    double svdCut = 1e-16;
    // Initialize the inverse of the S diag
    gsl_vector *pseudoSInv = gsl_vector_alloc(dataDim);
    gsl_vector_set_all(pseudoSInv,0.0);

    // Vector of singular values that are larger than specified cut
    std::vector<double> aboveCutVals;
    
    std::cout << "The singular values above SVD Cut = " << svdCut << " are..." << std::endl;
    for ( int s = 0; s < dataDim; s++ )
      {
        double dum = gsl_vector_get(S,s);
        if ( dum >= svdCut ) { aboveCutVals.push_back(dum); std::cout << dum << " ";}
	std::cout << "\n";
      }
    
    // Assign the inverse of aboveCutVals to the pseudoSInv vector
    for ( std::vector<double>::iterator it = aboveCutVals.begin(); it != aboveCutVals.end(); ++it )
      {
        gsl_vector_set(pseudoSInv,it-aboveCutVals.begin(),1.0/(*it));
      }

    // Promote this pseudoSInv vector to a matrix, where entries are placed along diagonal
    gsl_matrix * pseudoSInvMat = gsl_matrix_alloc(dataDim,dataDim);
    gsl_matrix_set_zero(pseudoSInvMat);
    for ( int m = 0; m < dataDim; m++ )
      {
        gsl_matrix_set(pseudoSInvMat,m,m,gsl_vector_get(pseudoSInv,m));
        // std::cout << gsl_vector_get(pseudoSInv,m) << std::endl;
      }

    /*
      With singular values that are zero, best we can do is construct a pseudo-inverse

      In general, the inverse we are after is
                               VS^(-1)U^T
                     with S the diagonal matrix of singular values
    */
    // In place construct the transpose of U ('toInvert' was modified in place to U above)
    gsl_matrix_transpose(toInvert);

    gsl_matrix * SinvUT = gsl_matrix_alloc(dataDim,dataDim); gsl_matrix_set_zero(SinvUT);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pseudoSInvMat,toInvert,0.0,SinvUT);

    // Now make the inverse of 'toInvert'
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,SinvUT,0.0,MInv);


    std::cout << "Leaving LinAlg::matrixInv" << std::endl;
    return pseudoSInv->size - aboveCutVals.size();
#endif
    return 0;
  }
#endif
}

