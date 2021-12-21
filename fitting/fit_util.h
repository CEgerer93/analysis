/*
  Define classes/structs/methods needed to handle the reduced pseudo-ioffe-time distributions
*/

#ifndef __fit_util_h__
#define __fit_util_h__

#include<vector>
#include<string>
#include<math.h>
#include<map>

#include<gsl/gsl_vector.h>
#include "corr_utils.h"

/* enum TYPES {LIN, LIN_tExp, TwoState}; */
/* enum TYPES {1 = "LIN", 2= "LIN_tExp", 3 = "2-state"}; */

namespace NFIT
{
  /*
    Big Driver for Fits
  */
  void driver(NCOR::correlator *c, std::string &comp, bool varPro);
  /* void driver(NCOR::correlator *c, const char * comp); */

  /*
    Evaluate correlated chi2
  */
  double chi2(const gsl_vector * x, void *data);
}
#endif
