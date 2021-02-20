/*
  Define classes/structs/methods needed to handle the reduced pseudo-ioffe-time distributions
*/

#ifndef __pitd_util_h__
#define __pitd_util_h__

#include<vector>
#include<string>
/* #include<map> */
/* #include<iostream> */
/* #include<iomanip> */
/* #include<complex> */

#include<gsl/gsl_matrix.h>

/* #include "hdf5.h" */


namespace FIT
{
  struct linFit_t
  {
    std::vector<double> m; // the matrix element at T
    std::vector<double> T; // the T's
    gsl_matrix * covinv;   // inverse of data covariance
  linFit_t(std::vector<double> &_m, std::vector<int> &_T, gsl_matrix *cov_inv): m(_m), T(_T), covinv(cov_inv) {}
  };

  struct FitRes_t
  {
    double a, b, c, d;
    std::string type;

    double func(int T)
    {
      if ( type == "LINEAR" )
	return a + b*T;
      if ( type == "LIN_tExp" )
	return a + b*T + c*T*exp(-d*T);
    }
    
  FitRes_t(std::string &_type, double _a = 0.0, double _b = 0.0, double _c = 0.0, double _d = 0.0) :
    a(_a), b(_b), c(_c), d(_d) { type = _type; }
  };
  /* struct momVals */
  /* { */
  /*   double                             IT; */
  /*   std::vector<std::complex<double> > mat; */
  /*   std::complex<double>               matAvg; */
  /* momVals(int g_ = 1) : mat(g_) {} */
  /* }; */


  /* struct zvals */
  /* { */
  /*   std::map<std::string, momVals> moms; */
  /*   // Polynomial fit results to pITD */
  /*   std::vector<polyFitParams_t>   polyR; */
  /*   std::vector<polyFitParams_t>   polyI; */
  /* }; */


  /* // PITD info for each of Real/Imag components */
  /* struct pitd */
  /* { */
  /*   std::map<int, zvals> disps;                     // data for each zsep */
  /*   std::map<int, gsl_matrix *> covsR, covsI;       // covariances for each zsep */
  /*   std::map<int, gsl_matrix *> covsRInv, covsIInv; // inverse of covariances for each zsep */
  /*   std::map<int, int> svsR, svsI; */

  /*   gsl_matrix *covR, *covI;        // Full data covariances */
  /*   gsl_matrix *invCovR, *invCovI; // Inverses of full data covariances */
  /*   int svsFullR, svsFullI;        // # singular values removed from full data covariance */
  /* }; */


  /* /\* */
  /*   Master class for the reducedPITD */
  /* *\/ */
  /* class reducedPITD */
  /* { */
  /* public: */
  /*   // Hold real/imaginary components */
  /*   pitd data; */


  /*   // Default */
  /*   reducedPITD() {} */
  /*   // Parametrized */
  /*   reducedPITD(int g) { gauge_configs = g; } */
  /*   reducedPITD(int g, int zmin, int zmax, int pmin, int pmax) */
  /*     { */
  /* 	gauge_configs = g; */
  /* 	zminCut = zmin; zmaxCut = zmax; */
  /* 	pminCut = pmin; pmaxCut = pmax; */
  /*     } */
  /*   // Destructor */
  /*   virtual ~reducedPITD() {}; */


  /*   // Quickly return the number of configs/jackknife samples */
  /*   int getCfgs() { return gauge_configs; } */

  /*   // Print all ensemble avg data with the same z value */
  /*   void ensemPrintZ(int zf, int comp); */
  /*   // Print the polynomial fit coefficients for a specified z */
  /*   void polyFitPrint(int zf, int comp); */


  /*   // Determine the data covariance for each zsep */
  /*   void calcCovPerZ(); */
  /*   // Determine the inverse of data covariance for each zsep */
  /*   void calcInvCovPerZ(); */

  /*   // Determine the full data covariance */
  /*   void calcCov(); */
  /*   void calcInvCov(); */

  /*   // Cut on Z's and P's to exclude from fit */
  /*   void cutOnPZ(int minz, int maxz, int minp, int maxp); */

  /*   // View a covariance matrix */
  /*   void viewZCovMat(int zsep); */
  /*   // View inverse of a covariance matrix */
  /*   void viewZCovInvMat(int zsep); */


  /*   /\* */
  /*     Isolate a single jackknife sample */
  /*   *\/ */

    
  /* private: */
  /*   int gauge_configs; */
  /*   int zminCut, zmaxCut, pminCut, pmaxCut; */
  /* }; */


  /* /\* // Print real/imag (comp) ensemble average data for provided zsep *\/ */
  /* /\* void reducedPITD::ensemPrintZ(int zf, int c) {} *\/ */

  /* /\* // Print the real/imag (comp) polynomial fit parameters to ensemble average data for provided zsep *\/ */
  /* /\* void reducedPITD::polyFitPrint(int zf, int comp) {} *\/ */

  
  /* /\* */
  /*   READER FOR PASSED H5 FILES */
  /* *\/ */
  /* void H5Read(char *inH5, reducedPITD *dat, int gauge_configs, int zmin, int zmax, int pmin, */
  /* 	      int pmax, std::string dTypeName); */

  
  /* /\* */
  /*   WRITER FOR MAKING NEW H5 FILES - E.G. EVOLVED/MATCHED DATASETS */
  /* *\/ */
  /* void H5Write(char *outH5, reducedPITD *dat, int gauge_configs, int zmin, int zmax, int pmin, */
  /* 	       int pmax, std::string dTypeName); */

  void printMat(gsl_matrix *g);
  int matrixInv(gsl_matrix * M, gsl_matrix * MInv);
  double chi2Linear(const gsl_vector * x, void *data);
}
#endif
