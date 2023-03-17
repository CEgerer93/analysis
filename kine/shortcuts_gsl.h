/*
  Shortcuts for accessing gsl
*/

#ifndef _shortcuts_gsl_h_
#define _shortcuts_gsl_h_

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<complex>
/*
  Shortcuts
*/
typedef gsl_matrix_complex gmc;
typedef gsl_vector_complex gvc;
typedef gsl_complex        gc;

/*
  Common
*/
const gc zero = gsl_complex_rect(0.0,0.0);
const gc one  = gsl_complex_rect(1.0,0.0);
const gc half = gsl_complex_rect(0.5,0.0);
const gc mone = gsl_complex_rect(-1.0,0.0);
const gc I    = gsl_complex_rect(0.0,1.0);
const gc mI   = gsl_complex_rect(0.0,-1.0);

const double PI = 3.1415926535897931;
const std::complex<double> _I_(0.0,1.0);
const std::complex<double> _mI_(0.0,-1.0);
const std::complex<double> _ZERO_(0.0,0.0);

const double aLat  = 0.094; // fm
const double hbarc = 0.1973269804; // GeV x fm

/*
  Macros
*/
#define gc_rect(r,i) gsl_complex_rect(r,i)
#define gc_div(n,d)

#endif
