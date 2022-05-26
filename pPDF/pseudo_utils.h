/*
  Utilities to help with any pseudo-distribution analyses
*/
#ifndef _pseudo_utils_h_
#define _pseudo_utils_h_

/* #include "io/adat_xmlio.h" */
/* #include "io/adat_xml_group_reader.h" */
/* #include "adat/map_obj.h" */
/* #include "hadron/hadron_sun_npart_npt_corr.h" */
#include "hadron/irrep_util.h"
#include "pseudo_structs.h"

namespace Pseudo
{
  // Ensure 3pt keys conserve 3-momentum
  void conserveMom3PtKey(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& k);
  // Set correct operator names and src/snk momenta based on global properties
  void setOpsMoms(Hadron::KeyHadronSUNNPartNPtCorr_t *k3, Hadron::KeyHadronSUNNPartNPtCorr_t& k2f,
		  Hadron::KeyHadronSUNNPartNPtCorr_t& k2i, Pseudo::global_t& g);
  void setOpsMoms(Hadron::KeyHadronSUNNPartNPtCorr_t *k3, Hadron::KeyHadronSUNNPartNPtCorr_t& k2f,
		  Hadron::KeyHadronSUNNPartNPtCorr_t& k2i, Hadron::KeyHadronSUNNPartNPtCorr_t& k2Rest,
		  Pseudo::global_t& g);

  // Negate all npoint[N] moms, accepting original full set of keys (k)
  void checkKeyMomNegate(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& k, global_t& g);

  // Build a short version of momentum - inspired by adat/lib/hadron/irrep_util.cc
  std::string shortMom(const XMLArray::Array<int>& mom, const std::string substr);

  // Build an intuitive displacement
  std::vector<int> shortZ(const std::vector<int>& d);

  // Ioffe time evaluation
  double ioffeTime(const XMLArray::Array<int>& mom, const std::vector<int>& disp, int L);

  /*
    Some debugging utilities
  */
  // print databases to be searched through
  template<typename T>
    void dumpDBs(std::vector<T>& vdb)
    {
      // n.b. typename needed in 'it' declaration, as 'it' depends on typename 'T'
      for ( typename std::vector<T>::const_iterator it = vdb.begin(); it != vdb.end(); ++it )
	std::cout << *it << std::endl;
    }
      

  // Print keys to be search for in databases
  void dumpKeys(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& vk, int n);
  
  // Print snk/src row combinations and signs
  void dumpRowInfo(XMLArray::Array<XMLArray::Array<int> >& r, XMLArray::Array<XMLArray::Array<int> >& s, int n);
}
#endif
