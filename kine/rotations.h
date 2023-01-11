/*
  Establish Euler rotations
*/
#ifndef _rotations_h_
#define _rotations_h_

#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

#include "shortcuts_gsl.h"

namespace Rotations
{
  // Rotation matrix for spin = 1/2 vectors using Euler angles in zyz-convention
  /*
    U[R(a,b,g)] = exp(-i a/2 \sigma_z) * exp(-i b/2 \sigma_y) * exp(-i g/2 \sigma_z)
  */
  gmc * eulerRotMat2(double alpha, double beta, double gamma);  

  // Rotation Matrix for Lorentz 4-vectors using euler angles z-y-z convention
  /*
    R(a,b,g) = Rz(a) Ry(b) Rz(g) acting on Cartesian four-vectors
  */
  gmc * eulerRotMat(double alpha, double beta, double gamma);
}

#endif
