/*
  Define classes/structs/methods needed to handle any linear algebra
*/

#include<iostream>
#include<vector>
#include<complex>
#include<gsl/gsl_permutation.h> // permutation header for matrix inversions
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h> // linear algebra
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

/* #include "shortcuts_gsl.h" */

#include <Eigen/Dense>

namespace LinAlg
{
  // Generic gsl_matrix viewer
  void printMat(gsl_matrix *g);

  // Generic gsl_matrix_complex viewer
  void printMat(gsl_matrix_complex *g);

  // Form tensor product of two matrices
  gsl_matrix_complex* tensorProd(gsl_matrix_complex *g, gsl_matrix_complex *h);


  /*
    PERFORM INVERSION OF PASSED DOUBLE MATRIX - RETURN # OF SVs REMOVED
  */
  int matrixInv(gsl_matrix * M, gsl_matrix * MInv);


  // TRY EIGEN SINCE IT HAS SUPPORT FOR SVD OF COMPLEX MATRICES
  /*
    EXTRACT INVARIANT AMPLTUDES USING (IN GENERAL) AN SVD DECOMPOSITION
  */
#if 1
  void extAmplitudes(Eigen::Matrix<std::complex<double>, 4, 2> * K);
#endif

#if 0
  void extAmplitudes(gmc * K, gvc * corrVec, gvc * amplitudes);
#endif
}

