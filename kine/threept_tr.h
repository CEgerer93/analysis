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
	      trace.real( 2*((Ef+m)*(Ei+m)+pf*pi*(2*M_PI/L)) ); break;
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

#endif
