/*
  Definitions of utilties to support any pseudo-distribution analyses
*/
#include "pseudo_utils.h"
#include "pseudo_structs.h"

namespace Pseudo
{
  // Ensure 3pt keys conserve 3-momentum
  void conserveMom3PtKey(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& k)
  {
    // Iterate through each key, and enforce 3-momentum conservation
    for ( std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>::iterator i = k.begin(); i != k.end(); ++i )
      {
	XMLArray::Array<int> _pf, _pi;
	_pf = i->npoint[1].irrep.irrep_mom.mom;
	_pi = i->npoint[3].irrep.irrep_mom.mom;
	i->npoint[2].irrep.irrep_mom.mom = _pf - _pi; // q
	i->npoint[2].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(i->npoint[2].irrep.irrep_mom.mom);
      } // end i
  }

  // Set correct operator names and src/snk momenta based on global properties
  void setOpsMoms(Hadron::KeyHadronSUNNPartNPtCorr_t *k3, Hadron::KeyHadronSUNNPartNPtCorr_t& k2f,
		  Hadron::KeyHadronSUNNPartNPtCorr_t& k2i, global_t& g)
  {
    // Setting the 3pt ops/moms
    k3->npoint[1].irrep.irrep_mom.mom = g.pf; // snk mom
    k3->npoint[1].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pf,"")]; // snk op
    k3->npoint[1].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k3->npoint[1].irrep.irrep_mom.mom);
    k3->npoint[3].irrep.irrep_mom.mom = g.pi; // src mom
    k3->npoint[3].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pi,"")]; // src op
    k3->npoint[3].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k3->npoint[3].irrep.irrep_mom.mom);
    k3->npoint[2].irrep.irrep_mom.mom = g.pf - g.pi; // q
    k3->npoint[2].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k3->npoint[2].irrep.irrep_mom.mom);
    // Setting the 2pt ops/moms
    for ( int n = 1; n != k2f.npoint.size()+1; n++ )
      {
	// Setting the snk 2pt ops/moms
	k2f.npoint[n].irrep.irrep_mom.mom = g.pf;
	k2f.npoint[n].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pf,"")];
	k2f.npoint[n].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k2f.npoint[n].irrep.irrep_mom.mom);
	// Setting the src 2pt ops/moms
	k2i.npoint[n].irrep.irrep_mom.mom = g.pi;
	k2i.npoint[n].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pi,"")];
	k2i.npoint[n].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k2i.npoint[n].irrep.irrep_mom.mom);
      }
    std::cout << "  Set src/snk operators/momenta in template keys " << std::endl;
  }

  // Set correct operator names and src/snk momenta based on global properties
  void setOpsMoms(Hadron::KeyHadronSUNNPartNPtCorr_t *k3, Hadron::KeyHadronSUNNPartNPtCorr_t& k2f,
                  Hadron::KeyHadronSUNNPartNPtCorr_t& k2i,
		  Hadron::KeyHadronSUNNPartNPtCorr_t& k2Rest, global_t& g)
  {
    // Setting the 3pt ops/moms
    k3->npoint[1].irrep.irrep_mom.mom = g.pf; // snk mom
    k3->npoint[1].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pf,"")]; // snk op
    k3->npoint[1].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k3->npoint[1].irrep.irrep_mom.mom);
    k3->npoint[3].irrep.irrep_mom.mom = g.pi; // src mom
    k3->npoint[3].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pi,"")]; // src op
    k3->npoint[3].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k3->npoint[3].irrep.irrep_mom.mom);
    k3->npoint[2].irrep.irrep_mom.mom = g.pf - g.pi; // q
    k3->npoint[2].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k3->npoint[2].irrep.irrep_mom.mom);
    // Setting the 2pt ops/moms
    for ( int n = 1; n != k2f.npoint.size()+1; n++ )
      {
        // Setting the snk 2pt ops/moms
        k2f.npoint[n].irrep.irrep_mom.mom = g.pf;
        k2f.npoint[n].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pf,"")];
        k2f.npoint[n].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k2f.npoint[n].irrep.irrep_mom.mom);
        // Setting the src 2pt ops/moms
        k2i.npoint[n].irrep.irrep_mom.mom = g.pi;
        k2i.npoint[n].irrep.op.ops[1].name = g.opMomXML[shortMom(g.pi,"")];
        k2i.npoint[n].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k2i.npoint[n].irrep.irrep_mom.mom);

	// Setting the rest 2pt ops/moms
	k2Rest.npoint[n].irrep.irrep_mom.mom = g.rest;
	k2Rest.npoint[n].irrep.op.ops[1].name = g.opMomXML[shortMom(g.rest,"")];
	k2Rest.npoint[n].irrep.op.ops[1].mom_type = Hadron::canonicalOrder(k2Rest.npoint[n].irrep.irrep_mom.mom);
      }
    std::cout << "  Set src/snk operators/momenta in template keys " << std::endl;
  }

  // Negate all npoint[N] moms, accepting original full set of keys (k)
  void checkKeyMomNegate(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& k, global_t& g)
  {
    // Make a copy of the passed vector of keys, so they may be modified via
    // pointer, and thus not affect the contents of k
    std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> kCopy(k);
    for ( auto a = kCopy.begin(); a != kCopy.end(); ++a )
      {
	a->npoint[1].irrep.irrep_mom.mom *= -1;
	if ( a->npoint.size() == 3 )
	  {
	    a->npoint[3].irrep.irrep_mom.mom *= -1;
	    a->npoint[2].irrep.irrep_mom.mom = a->npoint[1].irrep.irrep_mom.mom - a->npoint[3].irrep.irrep_mom.mom;
	    // update the src op in 3pt
	    a->npoint[3].irrep.op.ops[1].name = g.opMomXML[shortMom(a->npoint[3].irrep.irrep_mom.mom,"")];
	  }
	else // caught a 2pt key instead
	  {
	    a->npoint[2].irrep.irrep_mom.mom *= -1;
	    // update the src op in 2pt
	    a->npoint[2].irrep.op.ops[1].name = g.opMomXML[shortMom(a->npoint[2].irrep.irrep_mom.mom,"")];
	  }	
	// Update the snk operator based on this changed momenta
	a->npoint[1].irrep.op.ops[1].name = g.opMomXML[shortMom(a->npoint[1].irrep.irrep_mom.mom,"")]; // snk op
	// Now add this negative momentum key to the std::vector
	// But only if src or snk mom != 0
	if ( shortMom(a->npoint[1].irrep.irrep_mom.mom,"") != "000" ||
	     shortMom(a->npoint[a->npoint.size()].irrep.irrep_mom.mom,"") != "000" )
	  {
	    k.push_back(*a);
	  }
      }
  }

  // Build a short version of momentum - inspired by adat/lib/hadron/irrep_util.cc
  std::string shortMom(const XMLArray::Array<int>& mom, const std::string substr)
  {
    std::ostringstream os;
    os << mom[0] << substr << mom[1] << substr << mom[2];
    // os << substr << mom[0] << substr << mom[1] << substr << mom[2];
    return os.str();
  }

  // Build an intuitive displacement
  std::vector<int> shortZ(const std::vector<int>& d)
  {
    std::vector<int> dsmart(3);
    int dx(0), dy(0), dz(0);

    for ( auto a = d.begin(); a != d.end(); ++a )
      {
	switch(*a) {
	case -1:
	  dx-=1; break;
	case 1:
	  dx+=1; break;
	case -2:
	  dy-=1; break;
	case 2:
	  dy+=1; break;
	case -3:
	  dz-=1; break;
	case 3:
	  dz+=1; break;
	}
      }
    dsmart[0] = dx; dsmart[1] = dy; dsmart[2] = dz;
    return dsmart;
  }

  // Ioffe time evaluation
  double ioffeTime(const XMLArray::Array<int>& mom, const std::vector<int>& disp, int L)
  {
    if ( mom.size() != disp.size() )
      {
	std::cerr << "Ioffe time Error: Displacement must be an intuitive 3-vector" << std::endl;
	exit(2);
      }

    double ioffe(0.0);
    for ( int i = 0; i < disp.size(); ++i )
      ioffe += (2*M_PI/L)*mom[i]*disp[i];
    return ioffe;
  }

  /*
    Some debugging utilities
  */
  // Print keys to be search for in databases
  void dumpKeys(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& vk, int n)
  {
    std::cout << "Number of " << std::to_string(n) << " pt keys = " << vk.size() << std::endl;
    for ( std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>::const_iterator it = vk.begin(); it != vk.end(); ++it )
      {
	std::cout << *it << std::endl;
      }
  }

  // Print snk/src row combinations and signs
  void dumpRowInfo(XMLArray::Array<XMLArray::Array<int> >& r, XMLArray::Array<XMLArray::Array<int> >& s, int n)
  {
    std::cout << std::to_string(n) << "PT ROW STUFF" << std::endl;
    for ( int ri = 0; ri != r.size(); ri++ )
      {
        std::cout << " R [ " << ri << " ] = " << s[ri][0] << " * " << r[ri][0] << " " << r[ri][1] << std::endl;
      }
  }


}
