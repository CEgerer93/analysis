#ifndef __corr_utils_h__
#define __corr_utils_h__

#include<vector>
#include<map>
#include<complex>
#include<string>
#include<iostream>
#include<math.h>

#include<gsl/gsl_vector.h>

#include "hdf5.h"
#include "H5Cpp.h"

#include "cov_utils.h"
#include "pseudo_structs.h"

#include "adat/map_obj.h"
#include "hadron/hadron_sun_npart_npt_corr.h"


using namespace LinAlg;
using namespace H5;


/*
    OPERATORS
*/
// Define operator to easily print contents of a gsl_vector
std::ostream& operator<<(std::ostream& os, const gsl_vector *v);

namespace NCOR {
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T>& v);

/* template<typename T> */
std::ostream& operator<<(std::ostream& os, std::vector<std::complex<double> >& v);
}


namespace NCOR
{
  /* template<typename T> */
  std::vector<std::complex<double> > operator*=(std::vector<std::complex<double> >& v, double d);

  /* template<typename T> */
  std::vector<std::complex<double> > operator*=(std::vector<std::complex<double> >& v,
						std::complex<double> c);

  std::vector<std::complex<double> > operator+=(std::vector<std::complex<double> >& v1,
						std::vector<std::complex<double> >& v2);

  template<typename T>
  std::vector<std::complex<T> > operator-=(std::vector<std::complex<T> >& v1,
					   std::vector<std::complex<T> >& v2);


  typedef std::vector<std::vector<std::complex<double> > > VVC;

  // Data hold cfgs x Nt
  struct dat_t
  {
    std::vector<int> T; // temporal info
    VVC ens;
    std::vector<std::complex<double> > avg;
  };

  /* /\* */
  /*   OPERATORS */
  /* *\/ */
  /* // Define operator to easily print contents of a gsl_vector */
  /* std::ostream& operator<<(std::ostream& os, const gsl_vector *v); */

  /* template<typename T> */
  /*   std::ostream& operator<<(std::ostream& os, const std::vector<std::complex<T> >& v); */
  
  /* template<typename T, typename F> */
  /*   std::vector<T>& operator*=(std::vector<T>& v, F factor); */

  /* template<typename T> */
  /* std::vector<std::complex<double> > operator*=(std::vector<std::complex<double> >& v, */
  /* 						std::complex<double> c); */

  dat_t operator*(dat_t& d, std::complex<double> c);
  
  dat_t operator*=(dat_t& d, std::complex<double> c);

  // Direct add of ensembles
  dat_t operator+=(dat_t& d1, dat_t& d2);

  // Direct subtract of ensembles
  dat_t operator-=(dat_t& d1, dat_t& d2);

  // Fit info struct
  struct fitInfo_t
  {
    char type;
    Pseudo::domain_t range;

    std::string verbose()
    {
      std::string s;
      switch(type)
	{
	case 'l':
	  s = "a+bT"; break;
	case 'L':
	  s = "a+bT+cT*exp{-dE*T}"; break;
	case 'o':
	  s = "a*exp{-E0*T}"; break;
	case 't':
	  s = "exp{-E0*T}*{a+b*exp{-{E1-E0}*T}}"; break;
	}
      return s;
    }

    // Return # of non-linear fitted parameters
    int getNumNonLin()
    {
      switch(type)
	{
	case 'l':
	  return 0; break;
	case 'L':
	  return 1; break;
	case 'o':
	  return 1; break;
	case 't':
	  return 2; break;
	}
    }
  };


  // Struct to hold covariance and its inverse
  struct cov_t
  {
    std::map<std::string, gsl_matrix *> dat;
    std::map<std::string, gsl_matrix *> inv;
    std::map<std::string, int>          svs; // singular values removed in inverse
  };
    

  
  struct fitFunc_t
  {
    int num           = 0; // total # fitted params
    fitInfo_t theFit;
    cov_t fitCov;          // covariance & inverse on domain of fit

    /* // Try holding the priors/widths */
    /* struct { */
    /*   gsl_vector * mean; */
    /*   gsl_vector * width; */
    /* } priors; */

    // Some useful methods
    char getType() { return theFit.type; }

    /*
      Build covariance & its inverse on domain of fit
    */
    // Need correlator's domain_t to get correct submatrix of covariance
    void scov(std::map<std::string, gsl_matrix *> &c, Pseudo::domain_t &T)
    {
      // Sub - Covariance
      for ( auto it = c.begin(); it != c.end(); ++it )
	{
	  int subSize            = theFit.range.numT();
	  gsl_matrix * subCov    = gsl_matrix_alloc(subSize,subSize);
	  gsl_matrix * subCovInv = gsl_matrix_alloc(subSize,subSize);

	  for ( int i = 0; i < subSize; ++i )
	    {
	      for ( int j = 0; j < subSize; ++j )
		{
		  gsl_matrix_set(subCov,i,j,
				 gsl_matrix_get(it->second, (theFit.range.min+theFit.range.step*i
							     -T.min)/T.step,
						(theFit.range.min+theFit.range.step*j
						 -T.min)/T.step));
		} // j
	    } // i

	  /* std::cout << "******************Full covariance:"; */
	  /* LinAlg::printMat(it->second); */
	  /* std::cout << "**********************************" << std::endl;; */


	  /* std::cout << "Subcov:    "; LinAlg::printMat(subCov); */
	  /* std::cout << "Subcovinv (B4) inversion   "; LinAlg::printMat(subCovInv); */
	  fitCov.svs[it->first] = LinAlg::matrixInv(subCov,subCovInv);
	  /* std::cout << "Subcovinv: "; LinAlg::printMat(subCovInv); */

#if 0
	  gsl_matrix * id = gsl_matrix_alloc(subSize,subSize);
	  gsl_matrix_set_zero(id);
	  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,subCov,subCovInv,0.0,id);
	  std::cout << "ID?: "; LinAlg::printMat(id);
#endif

	  std::pair<std::string, gsl_matrix *> p = std::make_pair(it->first,subCov);
	  fitCov.dat.insert(p);

	  /* // Its inverse */
	  /* fitCov.svs[it->first] = LinAlg::matrixInv(subCov,subCovInv); */
	  std::pair<std::string, gsl_matrix *> q = std::make_pair(it->first,subCovInv);
	  fitCov.inv.insert(q);

	  /* std::cout << "SubCov check" <<  std::endl; */
	  /* LinAlg::printMat(subCov); */
	  /* LinAlg::printMat(subCovInv); */
	  /* std::cout << "End SubCov check" << std::endl; */


	} // it
    }


    double func(int T, const gsl_vector *p)
    {
      switch(theFit.type)
        {
        case 'l':
	  return gsl_vector_get(p,0) + gsl_vector_get(p,1)*T; break;
          /* return a + b*T; */
        case 'L':
	  return gsl_vector_get(p,0) + gsl_vector_get(p,1)*T
	    + gsl_vector_get(p,2)*T*exp(-gsl_vector_get(p,3)*T); break;
          /* return a + b*T + c*T*exp(-d*T); */
	case 'o':
	  return exp(-gsl_vector_get(p,0)*T)*gsl_vector_get(p,1); break;
	  /* return exp(-E0*T)*a; */
        case 't':
	  return gsl_vector_get(p,2)*exp(-gsl_vector_get(p,0)*T)+
	    gsl_vector_get(p,3)*exp(-gsl_vector_get(p,1)*T);
	  break;
          /* return a*exp(-E0*T)+b*exp(-E1*T); */
        }
    }

    
    std::map<std::string, double> printFit(const gsl_vector *p)
    {
      std::map<std::string, double> paramValMap;
      std::pair<std::string, double> dum;
      switch(theFit.type)
	{
	case 'l':
	  std::cout << "(a,b) = (" << p << ")" << std::endl;
	  dum = std::make_pair("a", gsl_vector_get(p,0)); paramValMap.insert(dum);
	  dum = std::make_pair("b", gsl_vector_get(p,1)); paramValMap.insert(dum);
	  break;
	case 'L':
	  std::cout << "(a,b,c,dE) = (" << p << ")" << std::endl;
	  dum = std::make_pair("a", gsl_vector_get(p,0));  paramValMap.insert(dum);
	  dum = std::make_pair("b", gsl_vector_get(p,1));  paramValMap.insert(dum);
	  dum = std::make_pair("c", gsl_vector_get(p,2));  paramValMap.insert(dum);
	  dum = std::make_pair("dE", gsl_vector_get(p,3)); paramValMap.insert(dum);
	  break;
	case 'o':
	  std::cout << "(E0,a) = (" << p << ")" << std::endl;
	  dum = std::make_pair("E0", gsl_vector_get(p,0)); paramValMap.insert(dum);
	  dum = std::make_pair("a", gsl_vector_get(p,1));  paramValMap.insert(dum);
	  // Would be nice to make an operator[] here...but many attempts failed
	  break;  
	case 't':
	  std::cout << "(E0,E1,a,b) = (" << p << ")" << std::endl;
	  dum = std::make_pair("E0", gsl_vector_get(p,0)); paramValMap.insert(dum);
	  dum = std::make_pair("E1", gsl_vector_get(p,1)); paramValMap.insert(dum);
	  dum = std::make_pair("a", gsl_vector_get(p,2));  paramValMap.insert(dum);
	  dum = std::make_pair("b", gsl_vector_get(p,3));  paramValMap.insert(dum);
	  break;
	}

      return paramValMap;
    } // printFit

    // Very hacky way for now to get the non-linear parameters of a fit
    std::vector<double> getNonLinParam(const gsl_vector *p)
    {
      std::vector<double> nl;
      switch(theFit.type)
	{
	case 'o':
	  nl.push_back(gsl_vector_get(p,0)); break;
	case 't':
	  nl.push_back(gsl_vector_get(p,0));
	  nl.push_back(gsl_vector_get(p,3)); break;
	default:
	  std::cout << "No non-linear parameters for type = " << theFit.type << std::endl;
	}
      return nl;
    }


    // Default
    fitFunc_t() {};
    // Constructor w/ initializer list
    fitFunc_t(fitInfo_t& _I, std::map<std::string, gsl_matrix *> &c, Pseudo::domain_t &T)
    {
      theFit = _I;
      scov(c, T); // submatrix of full data covariance - domain_t ensures correct submatrix is obtained
      switch(theFit.type)
        {
        case 'l':
          num = 2; break;
        case 'L':
          num = 4; break;
	case 'o':
	  num = 2; break;
        case 't':
          num = 4; break;
        }
    };
    /* fitFunc_t(std::string _t, double _a = 0.0, double _b = 0.0, double _c = 0.0, double _d = 0.0, */
    /*       double _E0 = 0.0, double _dE = 0.0) : */
    /* a(_a), b(_b), c(_c), d(_d) { type = _t; } */
  };

    
  struct fitRes_t
  {
    std::map<std::string, std::vector<double> > params;
    std::vector<double> chi2;
  };


  class correlator
  {
    Pseudo::prop_t aspects;
  public:
    dat_t ensemble;
    std::vector<dat_t> jack;
    cov_t cov; // full data covariance & inverse

    
    fitFunc_t fit;
    fitRes_t res;

    // Default
    correlator() {}

    // Parameterized
    correlator(dat_t _d, int _g, int _t)
      {
	aspects.Nt=_t; aspects.cfgs=_g;
	ensemble=_d;
	ensemble.T.resize(aspects.Nt);
      }
    /* correlator(int _g, Pseudo::domain_t _d) */
    correlator(Pseudo::prop_t _p)
      {
	aspects = _p;
	ensemble.T.resize(aspects.Nt);
	ensemble.ens.resize(aspects.cfgs);
	for ( auto it = ensemble.ens.begin(); it != ensemble.ens.end(); ++it )
	    it->resize(aspects.Nt);

	for ( int _t = aspects.domain.min; _t <= aspects.domain.max; _t+=aspects.domain.step )
	  ensemble.T[ ( _t - aspects.domain.min ) / aspects.domain.step ]=_t;
      }
    // Destructor
    virtual ~correlator() {};


    /*
      Accessors
    */
    int                                getCfgs()  { return aspects.cfgs; }
    int                                getNt()    { return aspects.Nt; }
    int                                npt()      { return aspects.npt; }
    int                                getGamma() { return aspects.gamma; }
    Hadron::KeyHadronSUNNPartNPtCorr_t key()      { return aspects.key; }


    XMLArray::Array<int>    getPi()   { return aspects.key.npoint[aspects.npt].irrep.irrep_mom.mom; }
    XMLArray::Array<int>    getPf()   { return aspects.key.npoint[1].irrep.irrep_mom.mom; }
    std::vector<int>        getDisp() { return aspects.key.npoint[2].irrep.op.ops[1].disp_list; }

    std::pair<std::string, int> getSrc()
    {
      std::pair<std::string, int> p(aspects.key.npoint[aspects.npt].irrep.op.ops[1].name,
				    aspects.key.npoint[aspects.npt].irrep.irrep_mom.row);
      return p;
    }

    std::pair<std::string, int> getSnk()
    {
      std::pair<std::string, int> p(aspects.key.npoint[1].irrep.op.ops[1].name,
				    aspects.key.npoint[1].irrep.irrep_mom.row);
      return p;
    }


    /*
      Public Methods
    */

    // Jackknife the data
    void jackknife();

    // Compute data covariance & its inverse
    void Cov();

    // Ensemble average
    void ensAvg();

    // Sum up temporal (e.g. summation method)
    void summation();

    // Correct for bias in forming ratios of Npt/2pt funcs (per jk ensemble avg)
    void removeBias();

    std::ostream& printParams();
  };


  /*
    Define a big class to hold correlation functions belonging to the same equivalence class
            i.e.  of same operator, abs(z), abs(p)
  */
  class corrEquivalence
  {
  private:
    int cfgs, Nt;
    Hadron::KeyHadronSUNNPartNPtCorr_t key;

  public:
    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlator> keyCorrMap;
    /* XMLArray::Array<XMLArray:: */
    /* int num_tslices; */

    /* corrEquivalence() {}; */
  corrEquivalence(int _g = 1, int _t = 1) : cfgs(_g), Nt(_t) {};
  };





  /*
    Operators on classes above
  */
  std::ostream& operator<<(std::ostream& os, correlator &c);

  
  /*
    Other methods
  */
  void read( std::ifstream& input, std::string& file, dat_t& corr);

  correlator mergeCorrs(std::vector<VVC>& v);

  void conj(VVC* c);

  /*
    H5 Support
  */
  void H5Read(char *inH5, correlator *c, int gauge_configs, int zmin, int zmax, int pmin,
              int pmax, std::string dTypeName);

  void writeCorr(correlator *c);

  void fitResW(correlator *c, std::string& comp); /* const char *comp); */
}
#endif
