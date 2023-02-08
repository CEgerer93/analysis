/*
  Structures to aid in pseudo-distribution analyses
*/
#ifndef _pseudo_structs_h_
#define _pseudo_structs_h_

#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"
#include "adat/map_obj.h"
#include "hadron/hadron_sun_npart_npt_corr.h"
#include "hadron/hadron_sun_npart_irrep.h"
/* #include "/u/home/cegerer/src/personal_fns.h" */

namespace Pseudo
{
  // Enum for types of current insertions
  enum currTypes { b_b0xDA__J0_A1pP, gamma_x, gamma_y };

  typedef XMLArray::Array<XMLArray::Array<int> > XMLInt2D;
  // Define a struct for global props
  struct global_t
  {
    std::string ensem, observable, state;
    int cfgs, nvec, Lt, Lx;
    int t2ptRows; // 1, 2, -or- 0 (for both 1 & 2)
    XMLArray::Array<int> pi, pf, rest; // lazy me
    std::vector<int> disp_list;


    /* projections::projSelect projector; // Don't want to mess w/ making enum reader */
    int projector; // whether unpol (1) or pol (2) projector is used in trace

    /* // Insertions properties */
    /* struct insProps_t */
    /* { */
    /*   KeyHadronSUNNPartIrrep_t key; */
    /* } insertions; */
    /* ADAT::MapObject<int, std::string> insertions; */

    // Basic operator structure for each i of Npt[i]
    Hadron::KeyHadronSUNNPartIrrep_t basic_op;

    // Npt properties
    struct nptProps_t
    {
      bool create, smear;
      Hadron::KeyCGCSU3_t cgc;
      ADAT::MapObject<std::string, int> cont_names;
      std::vector<int> disp_list;
    } snk, src, ins;

    bool dispNegate, momNegate;

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
    int Nt;

    currTypes adatCurrent;
    int gamma = -1;

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
    prop_t(int c_, domain_t &d_, const Hadron::KeyHadronSUNNPartNPtCorr_t& k)
    {
      cfgs = c_; domain = d_; key = k;
      Nt = domain.numT();

      if ( k.npoint[2].irrep.op.ops[1].name == "b_b0xDA__J0_A1pP" )
	adatCurrent = currTypes::b_b0xDA__J0_A1pP;
      if ( k.npoint[2].irrep.op.ops[1].name == "gamma_x" )
	adatCurrent = currTypes::gamma_x;
      if ( k.npoint[2].irrep.op.ops[1].name == "gamma_y" )
	adatCurrent = currTypes::gamma_y;
      
      // If we are dealing w/ a 3pt function, set the gamma int based on npoint[2] name
      if ( k.npoint.size() > 2 )
	{
	  switch(adatCurrent)
	    {
	    case b_b0xDA__J0_A1pP:
	      gamma = 8; break;
	    case gamma_x:
	      gamma = 1; break;
	    case gamma_y:
	      gamma = 2; break;
	    default:
	      std::cerr << "Expected b_b0xDA__J0_A1pP, gamma_x, or gamma_y insertions!" << std::endl;
	      std::cerr << "You provided op = " << k.npoint[2].irrep.op.ops[1].name << std::endl;
	      exit(2);
	    }
	} // if k.npoint.size() > 2

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
