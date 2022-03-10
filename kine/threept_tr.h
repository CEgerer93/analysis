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


/*
  Shortcuts
*/
typedef gsl_matrix_complex gmc;
typedef gsl_vector_complex gvc;
typedef gsl_complex        gc;

/* /\* */
/*   Common */
/* *\/ */
/* static gsl_complex zero = gsl_complex_rect(0.0,0.0); */
/* static gsl_complex one  = gsl_complex_rect(1.0,0.0); */
/* static gsl_complex mone = gsl_complex_rect(-1.0,0.0); */
/* static gsl_complex I    = gsl_complex_rect(0.0,1.0); */
/* static gsl_complex mI   = gsl_complex_rect(0.0,-1.0); */


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
};


struct diracMat_t
{
  gmc * gamma = gsl_matrix_complex_calloc(4,4);

  // Default
  diracMat_t() {}

  // Parameterized constructor
  diracMat_t(int n, bool MINK); // MINK = True for Minkowski construction
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
  subduceInfo subduce;

  // Wigner-D
  gmc * wig;
  
  // Subduction coefficients
  gmc * coeffS;

 public:
  spinor_t canon, subduced, absolute;


  /*
    Accessors
  */
  int                        getL() { return L; }
  double                     getE() { return E; }
  double                     getM() { return m; }
  XMLArray::Array<int>       getMom() { return mom; }
  double                     absMom() { return sqrt(mom*mom); }

  /*
    Public Methods
  */
  void initSubduce(const std::string& s);
  void buildSpinors();

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



#if 0
/*
  Polarization Vector - S^\mu = (1/2m)\bar{u}(p,s) \gamma^\mu \gamma^5 u(p,s)
*/
struct polVec_t
{
  // Common stuff
  /* gsl_complex zero   = gsl_complex_rect(0.0,0.0); */
  /* gsl_complex one    = gsl_complex_rect(1.0,0.0); */
  
  gvc * spinUp, * spinDown, * matVec;
  gmc * ID2d, *ID4d;


  // Pauli matrices
  pauli_t sx = pauli_t(1); pauli_t sy = pauli_t(2); pauli_t sz = pauli_t(3);
  // \gamma^4 & \gamma^5 are always apart of S^\mu evaluation
  diracMat_t g4 = diracMat_t(4,true); diracMat_t g5 = diracMat_t(5,true);

  // Dirac matrix \gamma^\mu that defines S^\mu
  diracMat_t d;



  // Evaluate the polarization vector
  std::complex<double> eval(XMLArray::Array<int>& p, double E, double m, int twoJz, int L)
  {
    std::complex<double> val(0.0,0.0);
    
    gvc * u = gsl_vector_complex_calloc(4);

    gmc * matProd1 = gsl_matrix_complex_calloc(4,4);
    gmc * matProd2 = gsl_matrix_complex_calloc(4,4);


    // Split 'p' into vec of gsl_complexes
    std::vector<gsl_complex> pc(3);
    for ( auto a = pc.begin(); a != pc.end(); ++a )
      *a = gsl_complex_rect(((2*M_PI/L)*p[std::distance(pc.begin(),a)])/(E+m),0);



    /*
      Form sigma \dot p
    */
    gmc * sDotP = gsl_matrix_complex_calloc(2,2);

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
    gvc * lower = gsl_vector_complex_calloc(2);
    gvc * tmp = gsl_vector_complex_alloc(2);
    switch(twoJz)
      {
      case 1:
	gsl_vector_complex_memcpy(tmp,spinUp); break;
      case -1:
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


#if 0
    diracMat_t G = diracMat_t(5,false);
    // More debugs
    gsl_blas_zgemv(CblasNoTrans,one,G.gamma,u,zero,matVec);
    for ( int s = 0; s < 4; ++s )
      std::cout << "r[" << s << "] = " << gsl_vector_complex_get(matVec,s).dat[0] << " " << gsl_vector_complex_get(matVec,s).dat[1] << std::endl;
    exit(8);
#endif


    // Multiply \gamma^4\gamma^\mu
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,g4.gamma,d.gamma,zero,matProd1);
    // Right multiply by \gamma^5
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,matProd1,g5.gamma,zero,matProd2);


    /* std::cout << "GAMMA MATRICES PRODUCT 1:" << std::endl; */
    /* LinAlg::printMat(matProd1); */
    /* std::cout << "GAMMA5 :" << std::endl; */
    /* LinAlg::printMat(g5.gamma); */
    /* std::cout << "GAMMA MATRICES PRODUCT 2:" << std::endl; */
    /* LinAlg::printMat(matProd2); */
    /* exit(10); */

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
  polVec_t(int mu, bool MINK)
  {
    spinUp = gsl_vector_complex_calloc(2); gsl_vector_complex_set(spinUp,0,one);
    spinDown = gsl_vector_complex_calloc(2); gsl_vector_complex_set(spinDown,1,one);
    matVec = gsl_vector_complex_calloc(4);
    ID2d = gsl_matrix_complex_alloc(2,2); gsl_matrix_complex_set_identity(ID2d);
    ID4d = gsl_matrix_complex_alloc(4,4); gsl_matrix_complex_set_identity(ID4d);
    d = diracMat_t(mu,MINK);
  };
};
#endif


/*
  Utilities to help manage three pt function traces
*/
void writePolVec(int npt, int mu, int cfgs, const XMLArray::Array<int> &mom,
		 std::vector<std::complex<double> > &polVec);


/*
  Lorentz Vector formed by \bar{u}\left(p_f,s_f\right) \gamma^\mu u\left(p_i,s_i\right)
*/
struct ugu_t
{
  // Dirac matrix \gamma^4 needed for \bar{u}; Dirac matrix \gamma^\mu characterizing Lorentz vector
  diracMat_t g4, d;

  // Default
  ugu_t(int mu, bool MINK)
  {
    d = diracMat_t(mu,MINK);
    g4 = diracMat_t(4,MINK);
  }

  // Evaluate given two spinors - adjoint of left spinor handled internally
  std::complex<double> eval(gsl_vector_complex * left, gsl_vector_complex * right);
};


#if 0
struct ugu_t
{
  // Common stuff
  /* gsl_complex zero = gsl_complex_rect(0.0,0.0); */
  /* gsl_complex one  = gsl_complex_rect(1.0,0.0); */

  gvc * spinUp, * spinDown, *matVec;
  gmc * ID2d, * ID4d;

  // Local Pauli matrices
  pauli_t sx = pauli_t(1); pauli_t sy = pauli_t(2); pauli_t sz = pauli_t(3);

  // \gamma^4 always needed for \bar{u}
  diracMat_t g4;

  // Dirac matrix \gamma^\mu characterizing Lorentz vector
  diracMat_t d;

  
  /*
    Evaluate
  */
  std::complex<double> eval(XMLArray::Array<int>& pf, XMLArray::Array<int>& pi, double Ef, double Ei,
			    double m, int twoJz_f, int twoJz_i, int L)
  {
    std::complex<double> val(0.0,0.0);

    // Initial/final state spinors
    gvc * ui = gsl_vector_complex_calloc(4);
    gvc * uf = gsl_vector_complex_calloc(4);

    // To hold \gamma^4 \times \gamma^\mu
    gmc * matProd = gsl_matrix_complex_calloc(4,4);

    // Split 'pf' & 'pi' into vec of gsl_complexes
    std::vector<gsl_complex> pc_i(3), pc_f(3);
    for ( auto a = pc_i.begin(); a != pc_i.end(); ++a )
      *a = gsl_complex_rect(((2*M_PI/L)*pi[std::distance(pc_i.begin(),a)])/(Ei+m),0);
    for ( auto a = pc_f.begin(); a != pc_f.end(); ++a )
      *a = gsl_complex_rect(((2*M_PI/L)*pf[std::distance(pc_f.begin(),a)])/(Ef+m),0);

#if 0
    for ( int i = 0; i < 3; ++i )
      {
	std::cout << "pc_f[" << i << "] = " << pc_f[i].dat[0] << " " << pc_f[i].dat[1] << std::endl;
	std::cout << "pc_i[" << i << "] = " << pc_i[i].dat[0] << " " << pc_i[i].dat[1] << std::endl;
      }
#endif

    std::cout << "v set" << std::endl;
    /*
      Form sigma \dot pi  &  sigma \dot pf
    */
    gmc * sDotPi = gsl_matrix_complex_calloc(2,2);
    gmc * sDotPf = gsl_matrix_complex_calloc(2,2);
    // px \dot sigma_x OVER (E+m)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc_i[0],sx.m,ID2d,zero,sDotPi);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc_f[0],sx.m,ID2d,zero,sDotPf);
    // PLUS py \dot sigma_y OVER (E+m)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc_i[1],sy.m,ID2d,one,sDotPi);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc_f[1],sy.m,ID2d,one,sDotPf);
    // PLUS pz \dot sigma_z OVER (E+m)
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc_i[2],sz.m,ID2d,one,sDotPi);
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc_f[2],sz.m,ID2d,one,sDotPf);
    //------------------------------------------------------------------------
#if 0
    std::cout << "SDOTPF" << std::endl;
    std::cout << gsl_matrix_complex_get(sDotPf,0,0).dat[0] << " " << gsl_matrix_complex_get(sDotPf,0,0).dat[1] << std::endl;
    std::cout << gsl_matrix_complex_get(sDotPf,0,1).dat[0] << " " << gsl_matrix_complex_get(sDotPf,0,1).dat[1] << std::endl;
    std::cout << gsl_matrix_complex_get(sDotPf,1,0).dat[0] << " " << gsl_matrix_complex_get(sDotPf,1,0).dat[1] << std::endl;
    std::cout << gsl_matrix_complex_get(sDotPf,1,1).dat[0] << " " << gsl_matrix_complex_get(sDotPf,1,1).dat[1] << std::endl;
#endif


    // The basis vectors comprising 4-component spinor
    gvc * basis_i = gsl_vector_complex_alloc(2);
    gvc * lower_i = gsl_vector_complex_calloc(2);
    gvc * basis_f = gsl_vector_complex_alloc(2);
    gvc * lower_f = gsl_vector_complex_calloc(2);
    switch(twoJz_i)
      {
      case 1:
	gsl_vector_complex_memcpy(basis_i,spinUp); break;
      case -1:
	gsl_vector_complex_memcpy(basis_i,spinDown); break;
      }
    switch(twoJz_f)
      {
      case 1:
	gsl_vector_complex_memcpy(basis_f,spinUp); break;
      case -1:
	gsl_vector_complex_memcpy(basis_f,spinDown); break;
      }

    /*
      Build up ui/uf spinors
    */
    for ( int s = 0; s < 2; ++s )
      {
	gsl_vector_complex_set(ui,s,gsl_vector_complex_get(basis_i,s));
	gsl_vector_complex_set(uf,s,gsl_vector_complex_get(basis_f,s));
      }
    // Mult 2-comp spinors by respective sDotP's
    gsl_blas_zgemv(CblasNoTrans,one,sDotPi,basis_i,zero,lower_i);
    gsl_blas_zgemv(CblasNoTrans,one,sDotPf,basis_f,zero,lower_f);
    // Modified 2-comp spinors not become lower 2 comps of ui/uf spinors
    for ( int s = 0; s < 2; ++s )
      {
	gsl_vector_complex_set(ui,2+s,gsl_vector_complex_get(lower_i,s));
	gsl_vector_complex_set(uf,2+s,gsl_vector_complex_get(lower_f,s));
      }
    // Rescale each of ui/uf by \sqrt( (E+m)/(2m) )
    gsl_blas_zdscal(sqrt((Ei+m)/(2*m)),ui);
    gsl_blas_zdscal(sqrt((Ef+m)/(2*m)),uf);
    //---------------- Now we have ui & uf spinors ------------------------------
#if 0
    for ( int s = 0; s < 4; ++s )
      std::cout << "uf[" << s << "] = " << gsl_vector_complex_get(uf,s).dat[0] << " " << gsl_vector_complex_get(uf,s).dat[1] << std::endl;
    for ( int s = 0; s < 4; ++s )
      std::cout << "ui[" << s << "] = " << gsl_vector_complex_get(ui,s).dat[0] << " " << gsl_vector_complex_get(ui,s).dat[1] << std::endl;
#endif

    
    // Multiply \gamma^4\gamma^\mu
    gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,g4.gamma,d.gamma,zero,matProd);

    // Matrix vector product (gamma^\mu on ui-spinor)
    gsl_blas_zgemv(CblasNoTrans,one,matProd,ui,zero,matVec);
    // Inner product of u_f^\dagger & (gamma^\mu \times ui)  --  result placed in 'zero'
    gsl_blas_zdotc(uf,matVec,&zero);

    
    // Push real/imag components of result "zero" into final result 'val'
    val.real(GSL_REAL(zero)); val.imag(GSL_IMAG(zero));
    
    zero = gsl_complex_rect(0.0,0.0); // reset

    std::cout << val << std::endl;
    
    gsl_vector_complex_free(spinUp);
    gsl_vector_complex_free(spinDown);
    gsl_vector_complex_free(matVec);
    gsl_vector_complex_free(ui);
    gsl_vector_complex_free(uf);
    gsl_matrix_complex_free(g4.gamma);
    gsl_matrix_complex_free(matProd);
    gsl_matrix_complex_free(sDotPi);
    gsl_matrix_complex_free(sDotPf);
    gsl_vector_complex_free(basis_i); gsl_vector_complex_free(basis_f);
    gsl_vector_complex_free(lower_i); gsl_vector_complex_free(lower_f);
    gsl_matrix_complex_free(ID2d);
    gsl_matrix_complex_free(ID4d);

    return val;
  }

  // Default
  ugu_t(int mu, bool MINK)
  {
    spinUp = gsl_vector_complex_calloc(2); gsl_vector_complex_set(spinUp,0,one);
    spinDown = gsl_vector_complex_calloc(2); gsl_vector_complex_set(spinDown,1,one);
    matVec = gsl_vector_complex_calloc(4);
    ID2d = gsl_matrix_complex_alloc(2,2); gsl_matrix_complex_set_identity(ID2d);
    ID4d = gsl_matrix_complex_alloc(4,4); gsl_matrix_complex_set_identity(ID4d);
    d = diracMat_t(mu,MINK);
    g4 = diracMat_t(4,MINK);
  };
}; // ugu_t
#endif

/*
  Utilities to help manage three pt function traces
*/
void writePrefactor(int npt, int mu, int cfgs, const XMLArray::Array<int> &mom,
		    std::string strRoot, std::vector<std::complex<double> > &P);

void writeSpinorContract(int npt, int mu, int cfgs, const XMLArray::Array<int>& pf,
			 const XMLArray::Array<int>& pi, int twoJz_f, int twoJz_i,
			 std::vector<std::complex<double> >& C);
#endif
