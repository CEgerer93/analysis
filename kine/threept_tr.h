/*
  Header to declare/define a structure to manage kinematic
  factors originating from threept function traces
*/

#ifndef _threept_tr_h_
#define _threept_tr_h_

#include<string>
#include<vector>
#include<math.h>
#include<complex>
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
#include "rotations.h"
#include "shortcuts_gsl.h"

using namespace Pseudo;
/* using namespace std::complex_literals; */


namespace projections
{
  enum projSelect { UNPOL = 1, POL = 2 };
}
enum current { VECTOR, AXIAL, TENSOR };


std::ostream& operator<<(std::ostream& os, const gsl_vector * v);
std::ostream& operator<<(std::ostream& os, const gsl_vector_complex * v);
std::ostream& operator<<(std::ostream& os, const gsl_matrix * m);
std::ostream& operator<<(std::ostream& os, const gsl_matrix_complex * m);
std::ostream& operator<<(std::ostream& os, const current c);

std::complex<double> zeroFuzz(const std::complex<double>& w);

/*
  Spinor sanity checks
*/
void iniFinSpinorsEqualL(const int l1, const int l2);


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



const Eigen::Matrix4i metric
{ {1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1} };

/* struct metric_t */
/* { */
/*   Eigen::Matrix4i m; */
/*   m(0,0) = 1; m(1,1) = -1; m(2,2) = -1; m(3,3) = -1; */
/* } metric; */

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
      opLG         = Hadron::getIrrepLG(opIrrepLG);

      // Set handle based on momentum
      if ( shortMom(m,"") == "000" )
	{
	  // Don't know why getIrrepDim/getCubicRepNoParity fails when given irrep G1g1
	  // - so do a hard code here
	  opIrrepNoP   = Hadron::getCubicRepNoParity("G1g");
	  irrep_dim    = Hadron::getIrrepDim("G1g");

	  contLGLabel = "J1o2->" + opIrrepNoP + ",1";
	  H = Hadron::TheSubduceTableFactory::Instance().createObject(contLGLabel);
	}
      else {
	opIrrepNoP   = Hadron::getCubicRepNoParity(opHelIrrepLG);
	irrep_dim    = Hadron::getIrrepDim(opIrrepLG);

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

  // Euler Rotation for spin=1/2
  gmc * eulerRot;
  
  // Subduction coefficients
  gmc * coeffS;

 public:
  spinor_t absolute, canon, helicity, subduced;


  /*
    Accessors
  */
  int                  getTwoJ()     const { return twoJ; }
  int                  getL()        const { return L; }
  double               getE()        const { return E; }
  double               getM()        const { return m; }
  XMLArray::Array<int> getMom()      const { return mom; }
  double               absMom()      const { return sqrt(1.0*(mom*mom)); }
  int                  getIrrepDim() const { return subductInfo.irrep_dim; }

  /*
    Public Methods
  */
  void initSubduce(const std::string& s);
  void buildSpinors();
  // Evaluate \sum_s u(p,s)\bar{u}(p,s)
  void projector(); 

  /*
    Destructor
  */
  ~Spinor()
    {
      gsl_matrix_complex_free(eulerRot);
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
    // Determine the rotation angles from |\vec{p}|\hat{z} to \vec{p} outright
    if ( shortMom(mom,"") != "000" )
      rot = Hadron::cubicCanonicalRotation(mom);
    else
      {
	rot.alpha=0; rot.beta=0; rot.gamma=0;
      }

    // Get the Euler rotation matrix
    eulerRot = Rotations::eulerRotMat2(rot.alpha, rot.beta, rot.gamma);

    // Build subductions
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
  /* std::complex<double> eval(gsl_vector_complex * left, gsl_vector_complex * right); */
  std::complex<double> eval(gsl_vector_complex * left, gsl_vector_complex * right, double mass);
  /* std::complex<double> eval(XMLArray::Array<int>& p, double E, double m, int twoJz, int L); */

  polVec_t(int mu, bool MINK)
  {
    d  = diracMat_t(mu,MINK);
    g4 = diracMat_t(4,MINK);
    g5 = diracMat_t(5,MINK);
  }

  // Destructor
  ~polVec_t()
  {
    gsl_matrix_complex_free(g4.gamma);
    gsl_matrix_complex_free(g5.gamma);
    gsl_matrix_complex_free(d.gamma);
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
  Generic class to evaluate spinor contractions
*/
class contractHandler
{
 public:
  void vectorContract(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini,
		      std::map<int, std::pair<int,int> > rowMap, Eigen::MatrixXcd &mat);
  void axialContract(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini,
		     const std::vector<int> &disp,
		     std::map<int, std::pair<int,int> > rowMap, Eigen::MatrixXcd &mat);
  void tensorContract(int mu, int nu, bool MINK, double mass, Spinor *fin, Spinor *ini,
		      const std::vector<int> &disp,
		      std::map<int, std::pair<int,int> > rowMap, Eigen::MatrixXcd &mat);

  /* // Parameterized contructor */
  /* contractHandler(current c) */
  /*   { */
  /*     switch(c) */
  /* 	{ */
  /* 	case VECTOR: contract = vectorContract; break; */
  /* 	/\* case AXIAL:  contract = axialContract;  break; *\/ */
  /* 	/\* case TENSOR: contract = tensorContract; break; *\/ */
  /* 	} */
  /*   } */

  contractHandler() {}
};


/*
  Complex kinematic matrix for PDFs
*/
struct kinMatPDF_t : contractHandler
{
  // Lorentz indices
  int mu = -1; int nu = -1;
  // Current type
  current TYPE;
  // Potential displacement
  std::vector<int> disp;
  // Dynamic matrix to hold entries of a PDF kinematic matrix
  Eigen::MatrixXcd mat;

  /* // Return a pointer to function that will handle contractions */
  /* /\* void * contract = &(contractHandler(current c).contract); *\/ */
  /* void * contract; */

  void assemble(bool MINK, double mass, Spinor *fin, Spinor *ini);
  void echoAction();

  // Ensure ini/fin spinor momenta are equal
  void spinorMomsEqual(bool truth);

  /*
    Constructors
  */
  // Default
  kinMatPDF_t() {}
  // Parameterized
 kinMatPDF_t(const int row, const int col, current c, std::vector<int> _d) : TYPE(c), disp(_d)
  {
    mat.resize(row, col);
    std::cout << "--> Init'd a Kinematic Matrix for PDFs of size " << row << " x "
	      << col << std::endl;
  }
 kinMatPDF_t(const int row, const int col, current c, int _m, std::vector<int> _d) : TYPE(c), mu(_m), disp(_d)
  {
    mat.resize(row, col);
    std::cout << "--> Init'd a Kinematic Matrix for PDFs of size " << row << " x "
	      << col << std::endl;
  }
 kinMatPDF_t(const int row, const int col, current c, int _m, int _n, std::vector<int> _d) : TYPE(c), mu(_m), nu(_n), disp(_d)
  {
    mat.resize(row, col);
    std::cout << "--> Init'd a Kinematic Matrix for PDFs of size " << row << " x "
	      << col << std::endl;
  }
};

/*
  Complex kinematic matrix
*/
struct kinMatGPD_t
{
  // Dynamic matrix to hold entries of a GPD kinematic matrix
  Eigen::MatrixXcd mat;

  void assemble(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini,
		const std::vector<int> &disp);

  // Ensure ini/fin spinor momenta are different
  void spinorMomsEqual(bool truth);

  /*
    Constructors
  */
  // Default
  kinMatGPD_t() {}
  // Parameterized
  kinMatGPD_t(const int row, const int col);
};

/*
  Complex kinematic matrix for helicity PDFs
*/
struct kinMat3_t
{
  Eigen::Matrix<std::complex<double>, 4, 3> mat;
  /* Eigen::Matrix<std::complex<double>, 4, 2> mat; */
  /* Eigen::Matrix<std::complex<double>, 2, 2> mat; */

  void assemble(int mu, bool MINK, double mass, Spinor *s, const std::vector<int> &disp);
  void assembleBig(int mu, bool MINK, double mass, Spinor *s, const std::vector<int> &disp);

  // Default
  kinMat3_t() {}
};

/*
  Complex-valued kinematic matric for form factor decompositions
*/
struct ffMat_t
{
  int                                       gamma; // Chroma indexing of gamma matrices
  Eigen::Matrix<std::complex<double>, 4, 2> mat;

  // Member assembler methods
  void assemble(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini);
  /* void assemble(int mu, int nu, bool MINK, double mass, Spinor *fin, Spinor *ini); */

  // Default
  ffMat_t() {}
  
  // Parameterized
ffMat_t(int _g) : gamma(_g) {}
};

/*
  EXTRACT INVARIANT AMPLITUDES USING (IN GENERAL) AN SVD DECOMPOSITION
*/
/* //--------- Unpol. PDF */
/* void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * MAT, */
/* 		   std::vector<kinMatPDF_t> * KIN, */
/* 		   std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * AMP); */
//--------- Unpol. GPD
void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * MAT,
		   std::vector<kinMatGPD_t> * KIN,
		   std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * AMP);
//---------
void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * MAT,
		   std::vector<kinMat3_t> * KIN,
		   std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * AMP);
/* void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * MAT, */
/* 		   std::vector<kinMat3_t> * KIN, */
/* 		   std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * AMP); */
void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * MAT,
		   std::vector<ffMat_t> * KIN,
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
