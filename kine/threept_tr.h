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

#include "hadron/clebsch.h"
#include "hadron/irreps_cubic_factory.h"
#include "hadron/irreps_cubic_oct_factory.h"
#include "hadron/irreps_cubic_helicity_factory.h"
#include "hadron/subduce_tables_oct_factory.h"
#include "hadron/subduce_tables_lg_factory.h"
#include "hadron/subduce_tables.h"
#include "hadron/subduce_tables_factory.h"
#include "hadron/single_hadron_coeffs.h"

#include "adat/handle.h"

#include "cov_utils.h"
#include "pseudo_utils.h"

using namespace Pseudo;
using namespace std::complex_literals;

/*
  Shortcuts
*/
typedef gsl_matrix_complex gmc;
typedef gsl_vector_complex gvc;
typedef gsl_complex        gc;


std::complex<double> zeroComplex(const std::complex<double>& w);

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
	      std::cerr << "Unknown kinematic factor for unpolarized projector and Gamma = "
			<< gamma << std::endl;
	    }
	  break;
	}
      case projections::POL:
	{
#warning "Kine factors for pol. projection are correct ONLY in forward limit!"
	  // Only correct for forward case!
	  switch(gamma)
	    {
	    case 14: // \gamma_x\gamma_5   ---   n.b. need to check - factor may be imag instead
	      trace.real( 4*pf[2]*pi[1]*pow(2*M_PI/L,2) ); break;
	    case 13: // \gamma_y\gamma_5   ---   n.b. need to check - factor may be real instead
	      trace.imag( 4*pf[2]*pi[2]*pow(2*M_PI/L,2) ); break;
	    case 11: // \gamma_z\gamma_5   ---   n.b. redstar gzg5 has "i", so factor is real
	      trace.real( 4*(m*Ei+pow(m,2)+pi[2]*pi[2]*pow(2*M_PI/L,2)) ); break;
	      // Should this be 4\times X  or 2\times X ???  -- prolly 4 because forward v. off-forward
	    case 7: // \gamma_4\gamma_5
	      trace.real( -4*pi[2]*(2*M_PI/L)*(Ei+m) ); break;
	    default:
	      std::cerr << "Unknown kinematic factor for polarized projector and Gamma = "
			<< gamma << std::endl;
	    }
	  break;
	}
      default:
	{
	  std::cerr << "What is your projector?  -- Unpol. -or- Pol." << std::endl;
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
  gmc * m = gsl_matrix_complex_calloc(2,2);

  // Constructor
  pauli_t(int n);

  /* // Destructor */
  /* ~pauli_t() */
  /* { */
  /*   gsl_matrix_complex_free(m); */
  /* } */
};


struct diracMat_t
{
  gmc * gamma = gsl_matrix_complex_calloc(4,4);

  // Default
  diracMat_t() {}

  // Parameterized constructor
  diracMat_t(int n, bool MINK); // MINK = True for Minkowski construction

  /* // Destructor */
  /* ~diracMat_t() */
  /* { */
  /*   gsl_matrix_complex_free(gamma); */
  /* } */
};



struct spinor_t
{
  std::map<int, gsl_vector_complex> twoJz;

  // Default
  spinor_t() {}

  template<typename T> void build(XMLArray::Array<T>& p, double E, double m, int L);
};


/*
  Class to manage subduction info
*/
class subduceInfo
{
 public:
  // Public members
  std::string opHelIrrepLG, opIrrepLG, opIrrep, opIrrepNoP, opLG; // irrep info
  int irrep_dim;                                                  // irrep dim  
  std::string name;                                               // op name
  std::string contLGLabel;                                        // operator name from which
                                                                  // all other string members derived
  
  // Handle for how operator subduces
  ADAT::Handle< Hadron::SubduceTable > H;

  // Default
  subduceInfo() { name = "None"; }
  
  // Parametrized constructor
  subduceInfo(std::string s, XMLArray::Array<int> m)
    {
      name = s;
    
      opHelIrrepLG = Hadron::getSingleHadronOpIrrep(name);
      opIrrepLG    = Hadron::removeHelicity(opHelIrrepLG);
      opIrrep      = Hadron::removeIrrepLG(opIrrepLG);
      opIrrepNoP   = Hadron::getCubicRepNoParity(opHelIrrepLG);
      opLG         = Hadron::getIrrepLG(opIrrepLG);
    
      irrep_dim    = Hadron::getIrrepDim(opIrrepLG);

      // Set handle based on momentum
      if ( shortMom(m,"") == "000" )
	{
	  contLGLabel = "J1o2->" + opIrrepNoP + ",1";
	  H = Hadron::TheSubduceTableFactory::Instance().createObject(contLGLabel);
	}
      else {
	contLGLabel = "H1o2->H1o2" + opIrrepNoP + ",1";
	H = Hadron::TheSubduceTableFactory::Instance().createObject(contLGLabel);
      }
    }

  void printInfo()
  {
    std::cout << "  opHelIrrepLG = " << opHelIrrepLG << std::endl;
    std::cout << "  opIrrepLG    = " << opIrrepLG << std::endl;
    std::cout << "  opIrrep      = " << opIrrep << std::endl;
    std::cout << "  opIrrepNoP   = " << opIrrepNoP << std::endl;
    std::cout << "  opLG         = " << opLG << std::endl;
    std::cout << "  irrep_dim    = " << irrep_dim << std::endl;
    std::cout << "  contLGLabel  = " << contLGLabel << std::endl;
    std::cout << "  subduction table = " << std::endl;
    for ( int _i = 1; _i <= irrep_dim; ++_i )
      {
	std::cout << "                     ";
        for ( int _j = 1; _j <= irrep_dim; ++_j )
	  std::cout << (*H).operator()(_i,_j) << " ";
	std::cout << "\n";
      }
  }
};


/*
  Class to manage a spin-1/2 Spinor
*/
class Spinor
{
  const int            twoJ = 1;
  int                  L;       // Spatial volume
  double               E, m; // Energy and mass
  XMLArray::Array<int> mom;

  // Euler angles to rotate from |{p}|\hat{z} to \vec{p}
  Hadron::CubicCanonicalRotation_t rot;

  // Subduction information
  subduceInfo subductInfo;

  // Wigner-D
  gmc * wig;
  
  // Subduction coefficients
  gmc * coeffS;

 public:
  spinor_t canon, subduced, absolute;


  /*
    Accessors
  */
  int                  getTwoJ()     const { return twoJ; }
  int                  getL()        const { return L; }
  double               getE()        const { return E; }
  double               getM()        const { return m; }
  XMLArray::Array<int> getMom()      const { return mom; }
  double               absMom()      const { return sqrt(mom*mom); }
  int                  getIrrepDim() const { return subductInfo.irrep_dim; }

  /*
    Public Methods
  */
  void initSubduce(const std::string& s);
  void buildSpinors();

  /*
    Destructor
  */
  ~Spinor()
    {
      gsl_matrix_complex_free(wig);
      gsl_matrix_complex_free(coeffS);
    }

  /*
    Constructors
  */
  // Default
  Spinor() {}
  // Parameterized
 Spinor(const std::string& name, const XMLArray::Array<int>& _p,
	double _E, double _m, int _L): L(_L), E(_E), m(_m), mom(_p)
  {
    // Determine the rotation angles outright
    if ( shortMom(mom,"") != "000" )
      rot = Hadron::cubicCanonicalRotation(mom);
    else
      {
	rot.alpha=0; rot.beta=0; rot.gamma=0;
      }

    // Set the rank of Wigner-D matrix
    wig = gsl_matrix_complex_calloc(twoJ+1,twoJ+1);
    // Initialize subductions
    initSubduce(name);
  }
};

/*
  Polarization Vector S^\mu - Lorentz vector formed by (1/2m)\bar{u}(p,s) \gamma^\mu \gamma^5 u(p,s)
*/
struct polVec_t
{
  // Dirac matrices \gamma^4 & \gamma^5 that are always apart of polVect evaluation
  // + Dirac matrix d characterizing S^\mu
  diracMat_t g4, g5, d;

  // Evaluate the polarization vector
  std::complex<double> eval(gsl_vector_complex * left, gsl_vector_complex * right);
  /* std::complex<double> eval(XMLArray::Array<int>& p, double E, double m, int twoJz, int L); */

  polVec_t(int mu, bool MINK)
  {
    d  = diracMat_t(mu,MINK);
    g4 = diracMat_t(4,MINK);
    g5 = diracMat_t(5,MINK);
  }
};

/*
  Lorentz Scalar formed by \bar{u}\left(p_f,s_f\right) 1 u\left(p_i,s_i\right)
*/
struct u1u_t
{
  // Dirac matrix \gamma^4 needed for \bar{u}
  diracMat_t g4;

  // Evaluate scalar product of final/initial spinors
  std::complex<double> eval(gsl_vector_complex * left, gsl_vector_complex * right);

  // Default
  u1u_t(bool MINK)
  {
    g4 = diracMat_t(4,MINK);
  }

  // Destructor
  ~u1u_t()
  {
    gsl_matrix_complex_free(g4.gamma);
  }
};

/*
  Lorentz Vector formed by \bar{u}\left(p_f,s_f\right) \gamma^\mu u\left(p_i,s_i\right)
*/
struct ugu_t
{
  // Dirac matrix \gamma^4 needed for \bar{u}; Dirac matrix \gamma^\mu characterizing Lorentz vector
  diracMat_t g4, d;

  // Evaluate between given two spinors - adjoint of left spinor handled internally
  std::complex<double> eval(gsl_vector_complex * left, gsl_vector_complex * right);

  // Default
  ugu_t(int mu, bool MINK)
  {
    d = diracMat_t(mu,MINK);
    g4 = diracMat_t(4,MINK);
  }

  // Destructor
  ~ugu_t()
  {
    gsl_matrix_complex_free(g4.gamma);
    gsl_matrix_complex_free(d.gamma);
  }
};

/*
  Lorentz Tensor formed by \bar{u}\left(p_f,s_f\right) \sigma^{\mu\nu} u\left(p_i,s_i\right)
*/
struct utu_t
{
  // Dirac matrix \gamma^4 needed for \bar{u}
  // + Dirac matrices \gamma^\mu, \gamma^\nu characterizing Lorentz tensor
  diracMat_t g4, dl, dr; // dl == left Dirac matrix  &  dr == right Dirac matrix

  int mu, nu; // Convenience

  // Evaluate between given two spinors - adjoint of left spinor handled internally
  std::complex<double> eval(gsl_vector_complex * left, gsl_vector_complex * right);
  
  // Default
utu_t(int mu, int nu, bool MINK) : mu(mu), nu(nu)
  {
    dl = diracMat_t(mu,MINK); dr = diracMat_t(nu,MINK);
    g4 = diracMat_t(4,MINK);
  }

  // Destructor 
  ~utu_t()
  {
    gsl_matrix_complex_free(g4.gamma);
    gsl_matrix_complex_free(dl.gamma);
    gsl_matrix_complex_free(dr.gamma);
  }
};


/*
  Complex kinematic matrix
*/
struct kinMat_t
{
  /* gmc * mat; */
  /* Eigen::MatrixXcd mat; */
  Eigen::Matrix<std::complex<double>, 4, 2> mat;

  void assemble(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini);

  // Default
  kinMat_t() {}

  // Constructor
  kinMat_t(Spinor *fin, Spinor *ini)
  {
    /* /\* mat = gsl_matrix_complex_calloc(fin->getIrrepDim()*ini->getIrrepDim(),2); *\/ */
    /* mat = Eigen::MatrixXcd(fin->getIrrepDim()*ini->getIrrepDim(),2); */
  }
};




/*
  EXTRACT INVARIANT AMPLITUDES USING (IN GENERAL) AN SVD DECOMPOSITION
*/
void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * MAT,
		   std::vector<kinMat_t> * KIN,
		   std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * AMP);

/*
  Utilities to help manage three pt function traces
*/
void writePolVec(int npt, int mu, int cfgs, const XMLArray::Array<int> &mom,
		 std::vector<std::complex<double> > &polVec);

void writePrefactor(int npt, int mu, int cfgs, const XMLArray::Array<int> &mom,
		    std::string strRoot, std::vector<std::complex<double> > &P);

void writeSpinorContract(int npt, int mu, int cfgs, const XMLArray::Array<int>& pf,
			 const XMLArray::Array<int>& pi, int twoJz_f, int twoJz_i,
			 std::vector<std::complex<double> >& C);
#endif
