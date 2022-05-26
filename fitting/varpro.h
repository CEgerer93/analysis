/*
  Define classes/structs/methods needed for variable projection
*/

#ifndef __varpro_h__
#define __varpro_h__

#include<vector>
#include<map>
#include<iostream>
#include<iomanip>
#include<complex>

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

#include "hdf5.h"

namespace VarPro
{
  /*
    Construct a model function from vector of coefficients and vector of basis functions

                     f(x_k)=\sum_j c_j \Phi_j(non-linear params, x_k)
  */
  class varPro
  {
  public:
    gsl_matrix * basis; // Non-linear Basis functions
    gsl_vector * Y;     // Y = \sum_i c_i*Phi_i(non-linear params)
    gsl_matrix * Phi;   // Outer product of basis functions;
    gsl_matrix * invPhi;

    gsl_vector * soln; // The variable projection solution for linear constants

    size_t rank;       // dim(Y) - # of fitted constants c_i removed by VarPro
    bool   bayesianFit;

    // Default
    varPro() {}
    // Parametrized
    varPro(size_t _rank, size_t numData, bool _bayesianFit)
      {
	bayesianFit = _bayesianFit;
	rank        = _rank;
	basis       = gsl_matrix_calloc(rank,numData);
	Y           = gsl_vector_alloc(rank);
	soln        = gsl_vector_alloc(rank);
	Phi         = gsl_matrix_calloc(rank,rank);
	invPhi      = gsl_matrix_calloc(rank,rank);
      }
 
    // Destructor
    ~varPro()
      {
	gsl_matrix_free(basis);
	gsl_vector_free(Y);
	gsl_matrix_free(Phi);
	gsl_matrix_free(invPhi);
	gsl_vector_free(soln);
      }

    // Populate the non-linear basis of functions
    void makeBasis(const gsl_vector * nlParams, const std::vector<int> &T);
    /* void makeBasis(std::vector<double> &nlParams, const std::vector<int> &T); */
    // Populate Y Solution vector
    void makeY(gsl_vector *data, gsl_matrix *invCov,
	       std::vector<double> &prior, std::vector<double> &width);
    // Populate Phi matrix
    void makePhi(gsl_matrix *invCov, std::vector<double> &width);
    // Get the inverse of Phi matrix
    void getInvPhi();
    // Solution of VarPro
    void getSoln();
    
  };



}
#endif
