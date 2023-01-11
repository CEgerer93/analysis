/*
  Establish Euler rotations
*/
#include<iostream>
#include "rotations.h"
#include<math.h>

namespace Rotations
{
  // Rotation matrix for spin = 1/2 vectors using Euler angles in zyz-convention
  /*
    U[R(a,b,g)] = exp(-i a/2 \sigma_z) * exp(-i b/2 \sigma_y) * exp(-i g/2 \sigma_z)
  */
  gmc * eulerRotMat2(double alpha, double beta, double gamma)
  {    
    gmc * Rza = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set(Rza,0,0,gc_rect(cos(alpha/2),-sin(alpha/2)));
    gsl_matrix_complex_set(Rza,1,1,gc_rect(cos(alpha/2),sin(alpha/2)));

    gmc * Ryb = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set(Ryb,0,0,gc_rect(cos(beta/2),0));
    gsl_matrix_complex_set(Ryb,0,1,gc_rect(-sin(beta/2),0));
    gsl_matrix_complex_set(Ryb,1,0,gc_rect(sin(beta/2),0));
    gsl_matrix_complex_set(Ryb,1,1,gc_rect(cos(beta/2),0));

    gmc * Rzg = gsl_matrix_complex_calloc(2,2);
    gsl_matrix_complex_set(Rzg,0,0,gc_rect(cos(gamma/2),-sin(gamma/2)));
    gsl_matrix_complex_set(Rzg,1,1,gc_rect(cos(gamma/2),sin(gamma/2)));

    gmc * mat = gsl_matrix_complex_calloc(2,2);
    gmc * res = gsl_matrix_complex_calloc(2,2);
    
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,Rza,Ryb,zero,mat);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,mat,Rzg,zero,res);

    gsl_matrix_complex_free(Rza);
    gsl_matrix_complex_free(Ryb);
    gsl_matrix_complex_free(Rzg);
    gsl_matrix_complex_free(mat);

    return res;
  }

  // Rotation Matrix for Lorentz 4-vectors using euler angles z-y-z convention
  /*
    R(a,b,g) = Rz(a) Ry(b) Rz(g) acting on Cartesian four-vectors
  */
  gmc * eulerRotMat(double alpha, double beta, double gamma)
  {
    gmc * Rza = gsl_matrix_complex_calloc(4,4);
    gsl_matrix_complex_set(Rza,0,0,gc_rect(1,0));
    gsl_matrix_complex_set(Rza,1,1,gc_rect(cos(alpha),0));
    gsl_matrix_complex_set(Rza,1,2,gc_rect(-sin(alpha),0));
    gsl_matrix_complex_set(Rza,2,1,gc_rect(sin(alpha),0));
    gsl_matrix_complex_set(Rza,2,2,gc_rect(cos(alpha),0));
    gsl_matrix_complex_set(Rza,3,3,gc_rect(1,0));

    gmc * Ryb = gsl_matrix_complex_calloc(4,4);
    gsl_matrix_complex_set(Ryb,0,0,gc_rect(1,0));
    gsl_matrix_complex_set(Ryb,1,1,gc_rect(cos(beta),0));
    gsl_matrix_complex_set(Ryb,1,3,gc_rect(sin(beta),0));
    gsl_matrix_complex_set(Ryb,2,2,gc_rect(1,0));
    gsl_matrix_complex_set(Ryb,3,1,gc_rect(-sin(beta),0));
    gsl_matrix_complex_set(Ryb,3,3,gc_rect(cos(beta),0));

    gmc * Rzg = gsl_matrix_complex_calloc(4,4);
    gsl_matrix_complex_set(Rzg,0,0,gc_rect(1,0));
    gsl_matrix_complex_set(Rzg,1,1,gc_rect(cos(gamma),0));
    gsl_matrix_complex_set(Rzg,1,2,gc_rect(-sin(gamma),0));
    gsl_matrix_complex_set(Rzg,2,1,gc_rect(sin(gamma),0));
    gsl_matrix_complex_set(Rzg,2,2,gc_rect(cos(gamma),0));
    gsl_matrix_complex_set(Rzg,3,3,gc_rect(1,0));

    gmc * mat = gsl_matrix_complex_calloc(4,4);
    gmc * res = gsl_matrix_complex_calloc(4,4);
    
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,Rza,Ryb,zero,mat);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,mat,Rzg,zero,res);

    gsl_matrix_complex_free(Rza);
    gsl_matrix_complex_free(Ryb);
    gsl_matrix_complex_free(Rzg);
    gsl_matrix_complex_free(mat);

    return res;
  }
}
