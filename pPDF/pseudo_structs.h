/*
  Structures to aid in pseudo-distribution analyses
*/
#ifndef _pseudo_structs_h_
#define _pseudo_structs_h_

#include "pseudo_utils.h"
#include "/u/home/cegerer/src/personal_fns.h"

namespace Pseudo
{
  typedef XMLArray::Array<XMLArray::Array<int> > XMLInt2D;
  // Define a struct for global props
  struct global
  {
    std::string ensem, observable, state;
    int cfgs, t2pt, nvec;
    XMLArray::Array<int> pi, pf;
    
    // XMLArray::Array<XMLArray::Array<std::string> > opMomXML;
    // std::map<std::string,std::string> opMomMap;

    ADAT::MapObject<std::string,std::string> opMomXML;
  };

  // Define some structs to help determine what edbs to read
  struct t3pt
  {
    int tm3pt, ts3pt, tx3pt;
  };

  struct info3pt
  {
    XMLArray::Array<std::string> base;
    std::string momTag, tsnkTag, t0Tag, zTag;
    // How rows and real/imag components are to be combined
//    XMLArray::Array<XMLArray::Array<int> > rows, signs;

    ADAT::MapObject<XMLArray::Array<int>, double> rowWeights;
    
    /* ADAT::MapObject<ADAT::Array2dO<int>, ADAT::MapObject<XMLArray::Array<int>, double> rowWeights; */
  };
  
  struct info2pt
  {
    std::string base;
    std::string momTag, t0Tag;
//    XMLArray::Array<XMLArray::Array<int> > rows, signs;
    ADAT::MapObject<XMLArray::Array<int>, double> rowWeights;
  };



  /*
    Define a big class to hold correlation functions of same operator, abs(z), abs(p)
  */
  class typeCorr
  {
  private:
    int gauge_configs, num_tslices;
    Hadron::KeyHadronSUNNPartNPtCorr_t key;

  public:
    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> corrs;
    /* XMLArray::Array<XMLArray:: */
    /* int num_tslices; */

    typeCorr() {};
  /* typeCorr(int n = 1, int g = 1, int t = 1) : corrs(n), gauge_configs(g), num_tslices(t) {}; */
    


  /*   // Public members */
  /*   correlators corrs; */
    
  /*   // Constructors */
  /*   typeCorr(int,int,int) {}; */
  /*   /\* typeCorr(Hadron::KeyHadronSUNNPartNPtCorr_t k,int n,int g,int t) :  *\/ */
  /*   /\*   Hadron::KeyHadronSUNNPartNPtCorr_t key = k, corrs(n), gauge_configs(g), num_tslices(t) {}; *\/ */
    
  /* typeCorr(int n,int g,int t) : corrs(n), gauge_configs(g), num_tslices(t) {} */
    
 /* // some key handle */
 /*    correlators corrs(n); // # of correlators of type k */
 /*    gauge_configs = g; // configs of each  */
 /*    num_tslices = t; // time slices of each */


    /* void rowAvg(XMLArray::Array<XMLArray::Array<int> >& r, XMLArray::Array<XMLArray::Array<int> >& s); */

  /* private: */
  /*   Hadron::KeyHadronSUNNPartNPtCorr_t key; */
  /*   int gauge_configs, num_tslices; */
  }; // typeCorr

}
#endif
