/*
  Define classes/structs/methods needed to handle any linear algebra
*/
#include<iostream>
#include<gsl/gsl_permutation.h> // permutation header for matrix inversions
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h> // linear algebra

#include "cov_utils.h"

namespace LinAlg
{
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

  // std::ostream& operator<<(std::ostream& os, const gsl_matrix_complex *g)
  // {
  //   os << "{";
  //   for ( size_t i = 0; i < g->size1; i++ ) {
  //     os << "{";
  //     for ( size_t j = 0; j < g->size2; j++ )
  // 	{
  // 	  std::complex<double> dum(0.0,0.0);
  // 	  gsl_complex gc = gsl_matrix_complex_get(g,i,j);
  // 	  dum.real(gc.dat[0]); dum.imag(gc.dat[1]);
  // 	  os << dum << ", ";
  // 	}
  //     os << "},\n";
  //   }
  //   return os;
  // }


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

    int pseudoSInvSize = pseudoSInv->size;
    // Free memory
    gsl_matrix_free(toInvert);
    gsl_matrix_free(V);
    gsl_vector_free(S);
    gsl_vector_free(work);
    gsl_vector_free(pseudoSInv);
    gsl_matrix_free(pseudoSInvMat);
    gsl_matrix_free(SinvUT);

    return pseudoSInvSize - aboveCutVals.size();
  }
	       


  // Get singular values of a complex-valued matrix (Eigen)
  void getSVs(Eigen::MatrixXcd *h)
  {
    Eigen::JacobiSVD<Eigen::MatrixXcd> s(*h,Eigen::ComputeFullU | Eigen::ComputeFullV);
    std::cout << s.singularValues() << std::endl;
  }

  // Get left/right singular vectors of a complex-valued matrix (Eigen)
  void getSingVecs(Eigen::MatrixXcd *h)
  {
    Eigen::JacobiSVD<Eigen::MatrixXcd> s(*h,Eigen::ComputeFullU | Eigen::ComputeFullV);
    std::cout << "Left singular vectors: " << std::endl << s.matrixU() << std::endl;
    std::cout << "Right singular vectors:" << std::endl << s.matrixV() << std::endl;
  }


#if 0
  void extAmplitudes(Eigen::Matrix<std::complex<double>, 4, 2> * K)
  {
    // Eigen::MatrixXcd M = Eigen::MatrixXcd::Random(4,4);
    // Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(4,4);
    // M(0,1) = std::complex<double>(2,1);
    // std::cout << "RANDOM EIGEN MATRIX = " << std::endl << M << std::endl;

    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(*K, Eigen::ComputeFullU | Eigen::ComputeFullV);
    std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
    std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
    std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
    Eigen::VectorXcd rhs = Eigen::VectorXcd(4);
    rhs(0) = 1; rhs(1) = 1; rhs(2) = 1; rhs(3) = -1;
    std::cout << "Now consider this rhs vector:" << std::endl << rhs << std::endl;
    std::cout << "A least-squares solution of m*x = rhs is:" << std::endl << svd.solve(rhs) << std::endl;

  }
#endif



#if 0
  /*
    EXTRACT INVARIANT AMPLTUDES USING (IN GENERAL) AN SVD DECOMPOSITION
       --> K = kinematic matrix
  */
  void extAmplitudes(gmc * K, gvc * corrVec, gvc * amplitudes)
  {
    size_t nRow = K->size1; size_t nCol = K->size2;
    gsl_matrix_complex * U    = gsl_matrix_complex_alloc(nRow,nCol);
    gsl_matrix_complex * V    = gsl_matrix_complex_alloc(nCol,nCol);
    gsl_vector_complex * S    = gsl_vector_complex_alloc(nCol);
    gsl_vector_complex * work = gsl_vector_complex_alloc(nCol);

    gsl_matrix_complex_memcpy(U,K); // make a copy of K, as U is modified below

    /*
      PERFORM THE SINGULAR VALUE DECOMPOSITION OF KINEMATIC MATRIX 'K'
      
      K = USV^T
          ---> K an MxN matrix
          ---> S an NxN singular value matrix (diagonal w/ singular values along diag - descending)
          ---> V an NxN matrix (returned from 'gsl_linalg_SV_decomp' in an untransposed form)
    */
    gsl_linalg_complex_SV_decomp(U,V,S,work);
    gsl_complex_linalg_SV_decomp(U,V,S,work);
    gsl_linalg_SV_complex_decomp(U,V,S,work);
    gsl_linalg_complex_SV_decomp_complex(U,V,S,work);

    /*
      Solve the linear system using the SVD decomposition above
    */
    gsl_linalg_complex_SV_solve(U,V,S,corrVec,amplitudes);
  }
#endif
}

