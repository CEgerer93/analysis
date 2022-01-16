/*
  Header to declare/define a structure to manage kinematic
  factors originating from threept function traces
*/

#ifndef _threept_tr_h_
#define _threept_tr_h_

#include<string>
#include<vector>
#include<complex>
#include<math.h>
#include "io/adat_xmlio.h"
#include "operators.h"
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h> // linear algebra
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

using namespace Pseudo;

namespace projections
{
  enum projSelect { UNPOL = 1, POL = 2 };
}


struct projector_t
{
  projections::projSelect projType;
  int         gamma;

  // Evaluate the trace
  std::complex<double> eval(XMLArray::Array<int> pf, XMLArray::Array<int> pi,
			    double Ef, double Ei, double m, int L)
  {
    std::complex<double> trace(0.0,0.0);
    switch(projType)
      {
      case projections::UNPOL:
	{
	  switch(gamma)
	    {
	    case 8:
	      trace.real( 2*((Ef+m)*(Ei+m)+pf*pi*pow(2*M_PI/L,2)) ); break;
	    case 4:
	      trace.imag( -2*(pf[2]*(2*M_PI/L)*(Ei+m)+pi[2]*(2*M_PI/L)*(Ef+m)) ); break;
	    case 2:
	      trace.imag( -2*(pf[1]*(2*M_PI/L)*(Ei+m)+pi[1]*(2*M_PI/L)*(Ef+m)) ); break;
	    case 1:
	      trace.imag( -2*(pf[0]*(2*M_PI/L)*(Ei+m)+pi[0]*(2*M_PI/L)*(Ef+m)) ); break;
	    default:
	      std::cerr << "Kinematic factor unknown for Gamma = " << gamma << std::endl;
	    }
	  break;
	}
      case projections::POL:
	{
	  std::cout << "You selected polarized" << std::endl;
	  break;
	}
      default:
	{
	  std::cerr << "What is your projector?" << std::endl;
	  exit(2);
	}
      } // switch

    return trace;
  }
  
  // Parameterized constructor
  projector_t(projections::projSelect p, int g)
  {
    projType = p;
    gamma    = g;
  }
};



struct pauli_t
{
  gsl_matrix_complex * m = gsl_matrix_complex_calloc(2,2);

  pauli_t(int n)
  {
    gsl_complex c;

    switch(n)
      {
      case 1:
	c = gsl_complex_rect(1.0,0);
	gsl_matrix_complex_set(m,0,1,c); gsl_matrix_complex_set(m,1,0,c); break;
      case 2:
	c = gsl_complex_rect(0,1.0);
	gsl_matrix_complex_set(m,1,0,c);
	c = gsl_complex_rect(0,-1.0);
	gsl_matrix_complex_set(m,0,1,c); break;
      case 3:
	c = gsl_complex_rect(1.0,0);
	gsl_matrix_complex_set(m,0,0,c);
	c = gsl_complex_rect(-1.0,0);
	gsl_matrix_complex_set(m,1,1,c); break;
      // Return 2x2 null matrix otherwise
      }	
  }
};


struct diracMat_t
{
  gsl_matrix_complex * gamma = gsl_matrix_complex_calloc(4,4);
  gsl_complex one = gsl_complex_rect(1.0,0.0);
  gsl_complex mone = gsl_complex_rect(-1.0,0.0);

  // Default
  diracMat_t() {}

  // Parameterized constructor
  diracMat_t(int n)
  {
    pauli_t s(n);
    switch(n)
      {
      case 1:
	for ( size_t i = 0; i < s.m->size1; ++i )
	  {
	    for ( size_t j = 0; j < s.m->size2; ++j )
	      {
		gsl_matrix_complex_set(gamma,i,2+j,gsl_matrix_complex_get(s.m,i,j));
		gsl_matrix_complex_set(gamma,i+2,j,
				       gsl_complex_mul(mone,gsl_matrix_complex_get(s.m,i,j)));
	      }
	  }
	break;
      case 2:
	for ( size_t i = 0; i < s.m->size1; ++i )
	  {
	    for ( size_t j = 0; j < s.m->size2; ++j )
	      {
		gsl_matrix_complex_set(gamma,i,2+j,gsl_matrix_complex_get(s.m,i,j));
		gsl_matrix_complex_set(gamma,i+2,j,
				       gsl_complex_mul(mone,gsl_matrix_complex_get(s.m,i,j)));
	      }
	  }
	break;
      case 3:
	for ( size_t i = 0; i < s.m->size1; ++i )
	  {
	    for ( size_t j = 0; j < s.m->size2; ++j )
	      {
		gsl_matrix_complex_set(gamma,i,2+j,gsl_matrix_complex_get(s.m,i,j));
		gsl_matrix_complex_set(gamma,i+2,j,
				       gsl_complex_mul(mone,gsl_matrix_complex_get(s.m,i,j)));
	      }
	  }
	break;
      case 4:
	for ( size_t i = 0; i < 2; ++i )
	  gsl_matrix_complex_set(gamma,i,i,one);
	for ( size_t i = 2; i < 4; ++i )
	  gsl_matrix_complex_set(gamma,i,i,mone);
	break;
      case 5:
	for ( int i = 3; i > -1; --i ) // int here because size_t decrements aren't intuitive at i = 0
	  gsl_matrix_complex_set(gamma,-i+3,i,one);
	break;
      default:
	std::cout << "Not a valid Dirac matrix!" << std::endl;
	std::cout << "You entered n = " << n << std::endl;
	exit(3);
      }
  }
};


/*
  Polarization Vector - S^\mu = (1/2m)\bar{u}(p,s) \gamma^\mu \gamma^5 u(p,s)
*/
struct polVec_t
{
  // Common stuff
  gsl_complex zero   = gsl_complex_rect(0.0,0.0);
  gsl_complex one    = gsl_complex_rect(1.0,0.0);
  
  gsl_vector_complex * spinUp, * spinDown, * matVec;
  gsl_matrix_complex * ID2d, *ID4d;


  // Pauli matrices
  pauli_t sx = pauli_t(1); pauli_t sy = pauli_t(2); pauli_t sz = pauli_t(3);
  // \gamma^4 & \gamma^5 are always apart of S^\mu evaluation
  diracMat_t g4 = diracMat_t(4); diracMat_t g5 = diracMat_t(5);

  // Dirac matrix \gamma^\mu that defines S^\mu
  diracMat_t d;



  // Evaluate the polarization vector
  std::complex<double> eval(XMLArray::Array<int>& p, double E, double m, int twoJz, int L)
  {
    std::complex<double> val(0.0,0.0);
    
    gsl_vector_complex * u = gsl_vector_complex_calloc(4);

    gsl_matrix_complex * matProd1 = gsl_matrix_complex_calloc(4,4);
    gsl_matrix_complex * matProd2 = gsl_matrix_complex_calloc(4,4);


    // Split 'p' into vec of gsl_complexes
    std::vector<gsl_complex> pc(3);
    for ( auto a = pc.begin(); a != pc.end(); ++a )
      *a = gsl_complex_rect(((2*M_PI/L)*p[std::distance(pc.begin(),a)])/(E+m),0);



    /*
      Form sigma \dot p
    */
    gsl_matrix_complex * sDotP = gsl_matrix_complex_calloc(2,2);

    // px \dot sigma_x OVER (E+m)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc[0],sx.m,ID2d,zero,sDotP);
    // PLUS py \dot sigma_y OVER (E+m)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc[1],sy.m,ID2d,one,sDotP);
    // PLUS pz \dot sigma_z OVER (E+m)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc[2],sz.m,ID2d,one,sDotP);
    //------------------------------------------------------------------------

#if 0
    std::cout << "SDOTP!" << std::endl;
    std::cout << gsl_matrix_complex_get(sDotP,0,0).dat[0] << " " << gsl_matrix_complex_get(sDotP,0,0).dat[1] << std::endl;
    std::cout << gsl_matrix_complex_get(sDotP,0,1).dat[0] << " " << gsl_matrix_complex_get(sDotP,0,1).dat[1] << std::endl;
    std::cout << gsl_matrix_complex_get(sDotP,1,0).dat[0] << " " << gsl_matrix_complex_get(sDotP,1,0).dat[1] << std::endl; 
    std::cout << gsl_matrix_complex_get(sDotP,1,1).dat[0] << " " << gsl_matrix_complex_get(sDotP,1,1).dat[1] << std::endl;
    exit(8);
#endif

    // Access either spinUp or spinDown basis vector
    gsl_vector_complex * lower = gsl_vector_complex_calloc(2);
    gsl_vector_complex * tmp = gsl_vector_complex_alloc(2);
    switch(twoJz)
      {
      case 1:
	gsl_vector_complex_memcpy(tmp,spinUp); break;
      case 2:
	gsl_vector_complex_memcpy(tmp,spinDown); break;
      }


    /*
      Build up u-spinor
    */
    for ( int s = 0; s < 2; ++s )
      gsl_vector_complex_set(u,s,gsl_vector_complex_get(tmp,s));

    // Mult 2-comp spinor by sDotP
    gsl_blas_zgemv(CblasNoTrans,one,sDotP,tmp,zero,lower);

    // Modified 2-comp spinor now becomes lower 2 comps of u-spinor
    for ( int s = 0; s < 2; ++s )
      gsl_vector_complex_set(u,2+s,gsl_vector_complex_get(lower,s));

    // Rescale u-spinor by \sqrt( (E+m)/(2m) )
    gsl_blas_zdscal(sqrt((E+m)/(2*m)),u);
    //---------------- Now we have u-spinor ---------------------------------

    //    gsl_vector_complex_free(tmp);


    // U-spinor seems correct!
#if 0
    for ( int s = 0; s < 4; ++s )
      std::cout << "u[" << s << "] = " << gsl_vector_complex_get(u,s).dat[0] << " " << gsl_vector_complex_get(u,s).dat[1] << std::endl;
    exit(8);
#endif


#if 1
    // More debugs
    gsl_blas_zgemv(CblasNoTrans,one,g4.gamma,u,zero,matVec);
    for ( int s = 0; s < 4; ++s )
      std::cout << "r[" << s << "] = " << gsl_vector_complex_get(matVec,s).dat[0] << " " << gsl_vector_complex_get(matVec,s).dat[1] << std::endl;
    exit(8);
#endif


    // Multiply \gamma^4\gamma^\mu
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,g4.gamma,d.gamma,zero,matProd1);
    // Right multiply by \gamma^5
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,matProd1,g5.gamma,zero,matProd2);
    // Matrix vector product (gamma's on u-spinor)
    gsl_blas_zgemv(CblasNoTrans,one,matProd2,u,zero,matVec);
    // Inner product of u^\dagger & (\gamma's \times u)  --  result placed in 'zero'
    gsl_blas_zdotc(u,matVec,&zero);


    // Push real/imag components of result "zero" into final result 'val'
    val.real(GSL_REAL(zero)); val.imag(GSL_IMAG(zero));
    // Normalize val
    val *= (1.0/(2*m));

    zero = gsl_complex_rect(0.0,0.0); //reset

    return val;
  }

  // Default
  polVec_t(int mu)
  {
    spinUp = gsl_vector_complex_calloc(2); gsl_vector_complex_set(spinUp,0,one);
    spinDown = gsl_vector_complex_calloc(2); gsl_vector_complex_set(spinDown,1,one);
    matVec = gsl_vector_complex_calloc(4);
    ID2d = gsl_matrix_complex_alloc(2,2); gsl_matrix_complex_set_identity(ID2d);
    ID4d = gsl_matrix_complex_alloc(4,4); gsl_matrix_complex_set_identity(ID4d);
    d = diracMat_t(mu);
  };
};


#endif
