/*
  Structures to aid in pseudo-distribution analyses
*/
#ifndef _pseudo_structs_h_
#define _pseudo_structs_h_

#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"
#include "adat/map_obj.h"
#include "hadron/hadron_sun_npart_npt_corr.h"
/* #include "/u/home/cegerer/src/personal_fns.h" */

namespace Pseudo
{
  typedef XMLArray::Array<XMLArray::Array<int> > XMLInt2D;
  // Define a struct for global props
  struct global_t
  {
    std::string ensem, observable, state;
    int cfgs, nvec, Lt, Lx;
    int t2ptRows; // 1, 2, -or- 0 (for both 1 & 2)
    XMLArray::Array<int> pi, pf, rest; // lazy me

    int chromaGamma; // gamma matrix appearing in 3pt trace

    bool momNegate;

    ADAT::MapObject<std::string,std::string> opMomXML;
  };


  /*
    Time domain info
  */
  struct domain_t
  {
    int min, step, max;
    /*
      Member functions
    */
    // Get number of T's
    int numT()    { return (max-min)/step + 1; }
    // Return vector, with each elem at Tsep
    std::vector<int> makeDomain()
    {
      std::vector<int> d;
      for ( int t = min; t <= max; t+= step )
	d.push_back(t);
      
      return d;
    }

    // Destructor
    ~domain_t() {};
    // Constructors
    domain_t() {};
    domain_t(int min_, int step_, int max_ )
    {
      min = min_; step = step_; max = max_;
    };
  };

  // Define some structs to help determine what edbs to read
  struct prop_t
  {
    int cfgs;
    int npt;
    int Nt, gamma;

    /* struct momenta_t */
    /* { */
    /*   XMLArray::Array<int> ini, fin; */
    /* } mom; */

    /* struct opt_t */
    /* { */
    /*   std::string name; */
    /*   int row; */
    /* } snk, src; */

    /* // Potential displacement in correlator */
    /* std::vector<int> disp;     */

    Hadron::KeyHadronSUNNPartNPtCorr_t key;

    domain_t domain;

    // Default
    prop_t() {};
    // Constructors
    prop_t(int c_, domain_t &d_)
    {
      cfgs = c_; domain = d_;
      Nt = domain.numT();
    };
    prop_t(int c_, domain_t &d_, Hadron::KeyHadronSUNNPartNPtCorr_t& k)
    {
      cfgs = c_; domain = d_; key = k;
      Nt = domain.numT();
    };
  };


  struct info3pt
  {
    XMLArray::Array<std::string> base;
    std::string momTag, tsnkTag, t0Tag, zTag;

    XMLArray::Array<XMLArray::Array<int> > rows;
  };
  
  struct info2pt
  {
    std::string base;
    std::string momTag, t0Tag;

    XMLArray::Array<XMLArray::Array<int> > rows;
  };
}
#endif
