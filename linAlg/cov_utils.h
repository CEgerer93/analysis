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

namespace LinAlg
{
  // Define operator to easily print contents of a gsl_vector
  /* std::ostream& operator<<(std::ostream& os, const gsl_vector *v); */
  /* // Accessor gsl_vector */
  /* double gsl_vector::operator[](const gsl_vector *v, int i ); */

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
#if 0
  /*
    PERFORM INVERSION OF PASSED COMPLEX MATRIX - RETURN # OF SVs REMOVED
  */
  int matrixInv(gsl_matrix_complex * M, gsl_matrix_complex * MInv);
#endif

  
}

