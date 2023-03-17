/*
  For cases where momenta used irrep_util version that did not split Rref from Rlatt

  The below is pulled from hadron/irrep_util.h at https://github.com/JeffersonLab/adat on commit 182c0295b86e9c7b40399e47c5c6edc08fcafa95
*/
#include "old_irrep_util_angles.h"

namespace SingleRotations
{
  //----------------------------------------------------------------------------------
  //! Return rotation angles for a vector to a canonical direction
  angles_t momSingleRotationAngles(const XMLArray::Array<int>& mom)
  {
    const double pi(3.14159265359);

    // Lists of momentum directions
    const int momListD2[][3] = {{1, 1, 0}, {0, 1, 1}, {1, 0, 1}, {1, -1, 0}, {0, 1, -1}, {-1, 0, 1}, {-1, 1, 0}, {0, -1, 1}, {1, 0, -1}, {-1, -1, 0}, {0, -1, -1}, {-1, 0, -1}};
    const int dimMomListD2 = 12;
    const double rotAnglesListD2[][3] = {{1.0/4.0, 1.0/2.0, 0}, {1.0/2.0, 1.0/4.0, -(1.0/2.0)}, {0, 1.0/4.0, -(1.0/2.0)}, {-(1.0/4.0), 1.0/2.0, 0}, {-(1.0/2.0), -(3.0/4.0), -(1.0/2.0)}, {-1.0, 1.0/4.0, -(1.0/2.0)}, {-(1.0/4.0), -(1.0/2.0), 0}, {-(1.0/2.0), 1.0/4.0, -(1.0/2.0)}, {0, 3.0/4.0, -(1.0/2.0)}, {1.0/4.0, -(1.0/2.0), 0}, {1.0/2.0, -(3.0/4.0), -(1.0/2.0)}, {-1.0, 3.0/4.0, -(1.0/2.0)}};   // {phi, theta, psi} in units of pi

    const int momListD3[][3] = {{1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, -1}};
    const int dimMomListD3 = 8;
    const double rotAnglesListD3[][3] = {{pi/4., acos(1/sqrt(3.0)), 0.0}, {(3*pi)/4., acos(1/sqrt(3.0)), 0.0}, {-pi/4., acos(1/sqrt(3.0)), 0.0}, {pi/4., acos(-(1/sqrt(3.0))), -pi/3.}, {(-3*pi)/4., acos(1/sqrt(3.0)), 0.0}, {(3*pi)/4., -acos(-(1/sqrt(3.0))), 0.0}, {-pi/4., -acos(-(1/sqrt(3.0))), 0.0}, {pi/4., -acos(-(1/sqrt(3.0))), 0.0}};

    const int momListC4nm0[][3] = {{0, 1, 2}, {1, 2, 0}, {2, 0, 1}, {0, 2, 1}, {2, 1, 0}, {1, 0, 2}, {0, 1, -2}, {1, -2, 0}, {-2, 0, 1}, {0, -2, 1}, {-2, 1, 0}, {1, 0, -2}, {0, -1, 2}, {-1, 2, 0}, {2, 0, -1}, {0, 2, -1}, {2, -1, 0}, {-1, 0, 2}, {0, -1, -2}, {-1, -2, 0}, {-2, 0, -1}, {0, -2, -1}, {-2, -1, 0}, {-1, 0, -2}};
    const int dimMomListC4nm0 = 24;
    const double rotAnglesListC4nm0[][3] = {{pi/2.0, acos(2.0/sqrt(5.0)), 0.0}, {-acos(-(1.0/sqrt(5.0))), -pi/2.0, pi/2.0}, {-pi, -acos(1.0/sqrt(5.0)), 0.0}, {-pi/2.0, -acos(1.0/sqrt(5.0)), 0.0}, {acos(2.0/sqrt(5.0)), pi/2.0, pi/2.0}, {0.0, acos(2.0/sqrt(5.0)), 0.0}, {-pi/2.0, -acos(-2.0/sqrt(5.0)), 0.0}, {-acos(1.0/sqrt(5.0)), pi/2.0, pi/2.0}, {0.0, -acos(1.0/sqrt(5.0)), 0.0}, {pi/2.0, -acos(1.0/sqrt(5.0)), 0.0}, {-acos(2.0/sqrt(5.0)), -pi/2.0, pi/2.0}, {-pi, -acos(-2.0/sqrt(5.0)), 0.0}, {-pi/2.0, acos(2.0/sqrt(5.0)), 0.0}, {acos(-(1.0/sqrt(5.0))), pi/2.0, pi/2.0}, {0.0, acos(-(1.0/sqrt(5.0))), 0.0}, {pi/2.0, acos(-(1.0/sqrt(5.0))), 0.0}, {acos(-2.0/sqrt(5.0)), -pi/2.0, pi/2.0}, {-pi, acos(2.0/sqrt(5.0)), 0.0}, {pi/2.0, -acos(-2.0/sqrt(5.0)), 0.0}, {acos(1.0/sqrt(5.0)), -pi/2.0, pi/2.0}, {-pi, acos(-(1.0/sqrt(5.0))), 0.0}, {-pi/2.0, acos(-(1.0/sqrt(5.0))), 0.0}, {-acos(-2.0/sqrt(5.0)), pi/2.0, pi/2.0}, {0.0, -acos(-2.0/sqrt(5.0)), 0.0}};

    const int momListC4nnm[][3] = {{1, 1, 2}, {1, 2, 1}, {2, 1, 1}, {-1, 1, 2}, {-1, 2, 1}, {1, -1, 2}, {1, 2, -1}, {2, -1, 1}, {2, 1, -1}, {1, 1, -2}, {1, -2, 1}, {-2, 1, 1}, {-1, -1, 2}, {-1, 2, -1}, {2, -1, -1}, {-1, 1, -2}, {-1, -2, 1}, {1, -1, -2}, {1, -2, -1}, {-2, -1, 1}, {-2, 1, -1}, {-1, -1, -2}, {-1, -2, -1}, {-2, -1, -1}};
    const int dimMomListC4nnm = 24;
    const double rotAnglesListC4nnm[][3] = {{(-3.0*pi)/4.0, -acos(sqrt(2.0/3.0)), 0.0}, {-acos(-(1.0/sqrt(5.0))), -acos(1.0/sqrt(6.0)), pi/2.0 + acos(-sqrt(3.0/5.0))}, {acos(2.0/sqrt(5.0)), acos(1.0/sqrt(6.0)), pi/2.0 - acos(-sqrt(3.0/5.0))}, {-pi/4.0, -acos(sqrt(2.0/3.0)), 0.0}, {acos(-(1.0/sqrt(5.0))), acos(1.0/sqrt(6.0)), pi/2.0 - acos(-sqrt(3.0/5.0))}, {(3.0*pi)/4.0, -acos(sqrt(2.0/3.0)), 0.0}, {-acos(-(1.0/sqrt(5.0))), -acos(-(1.0/sqrt(6.0))), pi/2.0 - acos(-sqrt(3.0/5.0))}, {acos(-2.0/sqrt(5.0)), -acos(1.0/sqrt(6.0)), pi/2.0 + acos(-sqrt(3.0/5.0))}, {acos(2.0/sqrt(5.0)), acos(-(1.0/sqrt(6.0))), pi/2.0 + acos(-sqrt(3.0/5.0))}, {pi/4.0, acos(-sqrt(2.0/3.0)), 0.0}, {-acos(1.0/sqrt(5.0)), acos(1.0/sqrt(6.0)), pi/2.0 - acos(-sqrt(3.0/5.0))}, {-acos(2.0/sqrt(5.0)), -acos(1.0/sqrt(6.0)), pi/2.0 + acos(-sqrt(3.0/5.0))}, {pi/4.0, -acos(sqrt(2.0/3.0)), 0.0}, {acos(-(1.0/sqrt(5.0))), acos(-(1.0/sqrt(6.0))), pi/2.0 + acos(-sqrt(3.0/5.0))}, {acos(-2.0/sqrt(5.0)), -acos(-(1.0/sqrt(6.0))), pi/2.0 - acos(-sqrt(3.0/5.0))}, {(3.0*pi)/4.0, acos(-sqrt(2.0/3.0)), 0.0}, {acos(1.0/sqrt(5.0)), -acos(1.0/sqrt(6.0)), pi/2.0 + acos(-sqrt(3.0/5.0))}, {-pi/4.0, acos(-sqrt(2.0/3.0)), 0.0}, {-acos(1.0/sqrt(5.0)), acos(-(1.0/sqrt(6.0))), pi/2.0 + acos(-sqrt(3.0/5.0))}, {-acos(-2.0/sqrt(5.0)), acos(1.0/sqrt(6.0)), pi/2.0 - acos(-sqrt(3.0/5.0))}, {-acos(2.0/sqrt(5.0)), -acos(-(1.0/sqrt(6.0))), pi/2.0 - acos(-sqrt(3.0/5.0))}, {(-3.0*pi)/4.0, acos(-sqrt(2.0/3.0)), 0.0}, {acos(1.0/sqrt(5.0)), -acos(-(1.0/sqrt(6.0))), pi/2.0 - acos(-sqrt(3.0/5.0))}, {-acos(-2.0/sqrt(5.0)), acos(-(1.0/sqrt(6.0))), pi/2.0 + acos(-sqrt(3.0/5.0))}};

    const int momListC4nnm113[][3] = {{1, 1, 3}, {1, 3, 1}, {3, 1, 1}, {-1, 1, 3}, {-1, 3, 1}, {1, -1, 3}, {1, 3, -1}, {3, -1, 1}, {3, 1, -1}, {1, 1, -3}, {1, -3, 1}, {-3, 1, 1}, {-1, -1, 3}, {-1, 3, -1}, {3, -1, -1}, {-1, 1, -3}, {-1, -3, 1}, {1, -1, -3}, {1, -3, -1}, {-3, -1, 1}, {-3, 1, -1}, {-1, -1, -3}, {-1, -3, -1}, {-3, -1, -1}};
    const double rotAnglesListC4nnm113[][3] = {{(-3*pi)/4., -atan(sqrt(2)/3.), 0}, {atan(3), atan(sqrt(10)), atan(sqrt(11)/3.)}, {atan(1./3.),   atan(sqrt(10)), -atan(sqrt(11)/3.)}, {-pi/4., -atan(sqrt(2)/3.), 0}, {-atan(3), -atan(sqrt(10)), pi - atan(sqrt(11)/3.)}, {-pi/4., atan(sqrt(2)/3.), pi}, {atan(3), pi - atan(sqrt(10)), pi - atan(sqrt(11)/3.)}, {-atan(1./3.), atan(sqrt(10)), atan(sqrt(11)/3.)}, {atan(1./3.), pi - atan(sqrt(10)), -pi + atan(sqrt(11)/3.)}, {(-3*pi)/4., -pi + atan(sqrt(2)/3.), pi}, {-atan(3), atan(sqrt(10)), -atan(sqrt(11)/3.)}, {-atan(1./3.), -atan(sqrt(10)), -pi + atan(sqrt(11)/3.)}, {(-3*pi)/4., atan(sqrt(2)/3.), pi}, {-atan(3), -pi + atan(sqrt(10)), atan(sqrt(11)/3.)}, {-atan(1./3.), pi - atan(sqrt(10)), pi - atan(sqrt(11)/3.)}, {-pi/4., -pi + atan(sqrt(2)/3.), pi}, {atan(3), -atan(sqrt(10)), -pi + atan(sqrt(11)/3.)}, {-pi/4., pi - atan(sqrt(2)/3.), 0}, {-atan(3), pi - atan(sqrt(10)), -pi + atan(sqrt(11)/3.)}, {atan(1./3.), -atan(sqrt(10)), pi - atan(sqrt(11)/3.)}, {-atan(1./3.), -pi + atan(sqrt(10)), -atan(sqrt(11)/3.)}, {(-3*pi)/4., pi - atan(sqrt(2)/3.), 0}, {atan(3), -pi + atan(sqrt(10)), -atan(sqrt(11)/3.)}, {atan(1./3.), -pi + atan(sqrt(10)), atan(sqrt(11)/3.)}};

    const int momListC4nnm122[][3] = {{2,2,1},{2,1,2},{1,2,2},{-2,2,1},{-2,1,2},{2,-2,1},{2,1,-2},{1,-2,2},{1,2,-2},{2,2,-1},{2,-1,2},{-1,2,2},{-2,-2,1},{-2,1,-2},{1,-2,-2},{-2,2,-1},{-2,-1,2},{2,-2,-1},{2,-1,-2},{-1,-2,2},{-1,2,-2},{-2,-2,-1},{-2,-1,-2},{-1,-2,-2}};
    const double rotAnglesListC4nnm122[][3] = {{(-3*pi)/4.,-atan(2*sqrt(2)),0},{atan(0.5),atan(sqrt(5)/2.),atan(3)},{atan(2),atan(sqrt(5)/2.),-atan(3)},{-pi/4.,-atan(2*sqrt(2)),0},{-atan(0.5),-atan(sqrt(5)/2.),pi - atan(3)},{-pi/4.,atan(2*sqrt(2)),pi},{atan(0.5),pi - atan(sqrt(5)/2.),pi - atan(3)},{-atan(2),atan(sqrt(5)/2.),atan(3)},{atan(2),pi - atan(sqrt(5)/2.),-pi + atan(3)},{(-3*pi)/4.,-pi + atan(2*sqrt(2)),pi},{-atan(0.5),atan(sqrt(5)/2.),-atan(3)},{-atan(2),-atan(sqrt(5)/2.),-pi + atan(3)},{(-3*pi)/4.,atan(2*sqrt(2)),pi},{-atan(0.5),-pi + atan(sqrt(5)/2.),atan(3)},{-atan(2),pi - atan(sqrt(5)/2.),pi - atan(3)},{-pi/4.,-pi + atan(2*sqrt(2)),pi},{atan(0.5),-atan(sqrt(5)/2.),-pi + atan(3)},{-pi/4.,pi - atan(2*sqrt(2)),0},{-atan(0.5),pi - atan(sqrt(5)/2.),-pi + atan(3)},{atan(2),-atan(sqrt(5)/2.),pi - atan(3)},{-atan(2),-pi + atan(sqrt(5)/2.),-atan(3)},{(-3*pi)/4.,pi - atan(2*sqrt(2)),0},{atan(0.5),-pi + atan(sqrt(5)/2.),-atan(3)},{atan(2),-pi + atan(sqrt(5)/2.),atan(3)}};

    double momMag = sqrt( pow(double(mom[0]),2) + pow(double(mom[1]),2) + pow(double(mom[2]),2) );
    double alpha = 0.0;
    double beta = 0.0;
    double gamma = 0.0;

    std::string littleGroup = Hadron::generateLittleGroup(mom);
    if (littleGroup == "D4")   // Momenta proportional to 001 etc
      {
	alpha = atan2(double(mom[1]),double(mom[0]));
	beta = acos(double(mom[2])/momMag);
	gamma = 0.0;
      }
    else if (littleGroup == "D2")
      {
	int momIndex = -1;
	for (int i = 0; i < dimMomListD2; i++)
	  {
	    if ( (momListD2[i][0] == int(double(mom[0])*sqrt(2.0)/momMag))
		 && (momListD2[i][1] == int(double(mom[1])*sqrt(2.0)/momMag))
		 && (momListD2[i][2] == int(double(mom[2])*sqrt(2.0)/momMag)) )   // Momenta proportional to 011 etc
	      {
		alpha = rotAnglesListD2[i][0] * pi;
		beta = rotAnglesListD2[i][1] * pi;
		gamma = rotAnglesListD2[i][2] * pi;
		momIndex = i;
		break;
	      }
	  }
	if (momIndex < 0)
	  {
	    std::cerr << __func__ << ": ERROR: can't find match for LG= " << littleGroup << "  and momentum " << mom[0] << " " << mom[1] << " " << mom[2] << std::endl;

	    exit(1);
	  }
      }
    else if (littleGroup == "D3")
      {
	int momIndex = -1;
	for (int i = 0; i < dimMomListD3; i++)
	  {
	    if ( (momListD3[i][0] == int(double(mom[0])*sqrt(3.0)/momMag))
		 && (momListD3[i][1] == int(double(mom[1])*sqrt(3.0)/momMag))
		 && (momListD3[i][2] == int(double(mom[2])*sqrt(3.0)/momMag)) )   // Momenta proportional to 111 etc
	      {
		//std::cout << __func__ << ": momentum direction " << momListD3[i][0] << momListD3[i][1] << momListD3[i][2] << std::endl;
		alpha = rotAnglesListD3[i][0];
		beta = rotAnglesListD3[i][1];
		gamma = rotAnglesListD3[i][2];
		momIndex = i;
		break;
	      }
	  }
	if (momIndex < 0)
	  {
	    std::cerr << __func__ << ": ERROR: can't find match for LG= " << littleGroup << "  and momentum " << mom[0] << " " << mom[1] << " " << mom[2] << std::endl;

	    exit(1);
	  }
      }
    else if (littleGroup == "C4nm0")
      {
	int momIndex = -1;
	for (int i = 0; i < dimMomListC4nm0; i++)
	  {
	    if ( (momListC4nm0[i][0] == int(double(mom[0])*sqrt(5.0)/momMag))
		 && (momListC4nm0[i][1] == int(double(mom[1])*sqrt(5.0)/momMag))
		 && (momListC4nm0[i][2] == int(double(mom[2])*sqrt(5.0)/momMag)) )   // Momenta proportional to 012 etc
	      {
		//std::cout << __func__ << ": momentum direction " << momListC4nm0[i][0] << momListC4nm0[i][1] << momListC4nm0[i][2] << std::endl;
		alpha = rotAnglesListC4nm0[i][0];
		beta = rotAnglesListC4nm0[i][1];
		gamma = rotAnglesListC4nm0[i][2];
		momIndex = i;
		break;
	      }
	  }
	if (momIndex < 0)
	  {
	    std::cerr << __func__ << ": ERROR: can't find match for LG= " << littleGroup << "  and momentum " << mom[0] << " " << mom[1] << " " << mom[2] << std::endl;

	    exit(1);
	  }
      }
    else if (littleGroup == "C4nnm")
      {
	int momIndex = -1;
	for (int i = 0; i < dimMomListC4nnm; i++)
	  {
	    if ( (momListC4nnm[i][0] == int(double(mom[0])*sqrt(6.0)/momMag))
		 && (momListC4nnm[i][1] == int(double(mom[1])*sqrt(6.0)/momMag))
		 && (momListC4nnm[i][2] == int(double(mom[2])*sqrt(6.0)/momMag)) )   // Momenta proportional to 112 etc
	      {
		//std::cout << __func__ << ": momentum direction " << momListC4nnm[i][0] << momListC4nnm[i][1] << momListC4nnm[i][2] << std::endl;
		alpha = rotAnglesListC4nnm[i][0];
		beta = rotAnglesListC4nnm[i][1];
		gamma = rotAnglesListC4nnm[i][2];
		momIndex = i;
		break;
	      }
	    else if ( (momListC4nnm113[i][0] == int(double(mom[0])*sqrt(11.0)/momMag))
		      && (momListC4nnm113[i][1] == int(double(mom[1])*sqrt(11.0)/momMag))
		      && (momListC4nnm113[i][2] == int(double(mom[2])*sqrt(11.0)/momMag)) )    // Momenta proportional to 113 etc
	      {
		//std::cout << __func__ << ": momentum direction " << momListC4nnm113[i][0] << momListC4nnm113[i][1] << momListC4nnm113[i][2] << std::endl;
		alpha = rotAnglesListC4nnm113[i][0];
		beta = rotAnglesListC4nnm113[i][1];
		gamma = rotAnglesListC4nnm113[i][2];
		momIndex = i;
		break;
	      }
	    else if ( (momListC4nnm122[i][0] == int(double(mom[0])*sqrt(9.0)/momMag))
		      && (momListC4nnm122[i][1] == int(double(mom[1])*sqrt(9.0)/momMag))
		      && (momListC4nnm122[i][2] == int(double(mom[2])*sqrt(9.0)/momMag)) )    // Momenta proportional to 122 etc
	      {
		//std::cout << __func__ << ": momentum direction " << momListC4nnm122[i][0] << momListC4nnm122[i][1] << momListC4nnm122[i][2] << std::endl;
		alpha = rotAnglesListC4nnm122[i][0];
		beta = rotAnglesListC4nnm122[i][1];
		gamma = rotAnglesListC4nnm122[i][2];
		momIndex = i;
		break;
	      }
	  }
	if (momIndex < 0)
	  {
	    std::cerr << __func__ << ": ERROR: can't find match for LG= " << littleGroup << "  and momentum " << mom[0] << " " << mom[1] << " " << mom[2] << std::endl;

	    exit(1);
	  }
      }
    else
      {
	std::cerr << __func__ << ": ERROR: unsupported momentum " << mom[0] << " " << mom[1] << " " << mom[2] << std::endl;
	exit(1);
      }
    //    std::cout << "alpha= " << alpha << ", beta= " << beta << ", gamma= " << gamma << std::endl;

    angles_t rot;

    rot.alpha = alpha;
    rot.beta  = beta;
    rot.gamma = gamma;

    return rot;
  }
}
