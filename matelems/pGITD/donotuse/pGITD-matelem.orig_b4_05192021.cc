/*
  ROUTINE TO READ IN 2PT/3PT CORRELATORS COMPUTED IN DISTILLATION,
  WHERE 3PT FUNCTION ARE MATELEMS OF PSEUDO-GPD FORMALISM

  2PT & 3PT FUNCTIONS ARE COMBINED ACCORDING TO KINEMATICS/DISPLACEMENTS

  DATA WRITTEN TO FILE ACCORDING TO:
     + 2-STATE FITS: 2PT/3PT WRITTEN SEPARATELY FOR FITTING
     + SUMMATION: SUMMED 3PT/2PT RATIO WRITTEN TO FILE FOR FITTING
*/
#include<vector>
#include<fstream>
#include<math.h>
#include<string>
#include<tuple>
#include<sstream>
#include<stdlib.h>
#include<iomanip>
#include<iostream>
#include<algorithm>
#include<ctype.h>
#include<omp.h>

// Headers to manage xml input
#include<libxml/xmlreader.h>
#include<libxml/xpath.h>
#include<libxml/parser.h>

#include "summation.h"
#include "operators.h"
#include "pseudo_utils.h"
#include "pseudo_structs.h"

// Bring in neccesary adat headers
#include "AllConfStoreDB.h"
#include "DBFunc.h"
#include "adat/map_obj.h"
#include "hadron/hadron_sun_npart_npt_corr.h"
#include "hadron/ensem_filenames.h"
#include "hadron/clebsch.h"
#include "io/key_val_db.h"
#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"
#include "adat/handle.h"
#include "hadron/subduce_tables_oct_factory.h"
#include "hadron/subduce_tables_lg_factory.h"
#include "hadron/subduce_tables.h"
#include "hadron/subduce_tables_factory.h"
// #include "hadron/irreps_cubic_helicity_factory.h"


using namespace ADATXML;
using namespace Summation;
using namespace Pseudo;

typedef XMLArray::Array<XMLArray::Array<int> > XML2D;

// Define a struct for global props
global props;

// Define some structs to help determine what edbs to read
t3pt temporal3pt;
info3pt db3ptInfo;
info2pt db2ptInfo;


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
  };

  void printInfo()
  {
    std::cout << "  opHelIrrepLG = " << opHelIrrepLG << std::endl;
    std::cout << "  opIrrepLG    = " << opIrrepLG << std::endl;
    std::cout << "  opIrrep      = " << opIrrep << std::endl;
    std::cout << "  opIrrepNoP   = " << opIrrepNoP << std::endl;
    std::cout << "  opLG         = " << opLG << std::endl;
    std::cout << "  irrep_dim    = " << irrep_dim << std::endl;
    std::cout << "  contLGLabel  = " << contLGLabel << std::endl;
  }
  

};


/*
  Convert between lattice LG rows and helicity amplitudes
*/
std::complex<double> lgToHelicityAmps(Hadron::KeyHadronSUNNPartNPtCorr_t *k)
{
  std::complex<double> weight; // weight to be applied to < \mu_f| J | \mu_i > correlator

  
  // Track the sink subductions
  subduceInfo snkOp(k->npoint[1].irrep.op.ops[1].name,k->npoint[1].irrep.irrep_mom.mom);
  subduceInfo srcOp(k->npoint[3].irrep.op.ops[1].name,k->npoint[3].irrep.irrep_mom.mom);
  
  // snkOp.printInfo();
  // srcOp.printInfo();


  for ( int ki = 1; ki <= snkOp.irrep_dim; ki++ )
    {
      for ( int si = 1; si <= srcOp.irrep_dim; si++ )
	{
	  weight += (*snkOp.H)(k->npoint[1].irrep.irrep_mom.row, ki)*
	    (*srcOp.H)(k->npoint[3].irrep.irrep_mom.row, si);
	}
    }

  return weight;
}


/*
  Insertion is computed in Redstar from source interpolator (i.e. contact term) forward to
  tslice = Tsep - 1. So to neglect constant term in fits to summed ratio data, ignore the T=0 entry
 */
void summation(std::vector<corrFunc> &f)
{
  // #pragma omp parallel for
//   {
// #pragma omp for
  for ( std::vector<corrFunc>::iterator it = f.begin(); it != f.end(); ++it )
    {
      (*it).initializeSummedJkEnsAvgs();
      for ( int g = 0; g < (*it).dataJkEnsAvg.ncor.size(); ++g )
	{
	  it->dataSummedJkEnsAvg.ncor[g].real.resize(1);
	  it->dataSummedJkEnsAvg.ncor[g].imag.resize(1);

	  double dumR=0.0;
	  double dumI=0.0;
	  /*
	    Do the actual summation in t here
	    Start summation from T = 1 - neglecting contact term
	  */
	  for ( int t = 1; t < (*it).getT(); t++ )
	    {
	      dumR += (*it).dataJkEnsAvg.ncor[g].real[0][t];
	      dumI += (*it).dataJkEnsAvg.ncor[g].imag[0][t];
	    }
	  // Assignment
	  (*it).dataSummedJkEnsAvg.ncor[g].real[0].push_back(dumR);
	  (*it).dataSummedJkEnsAvg.ncor[g].imag[0].push_back(dumI);
	}
    }
  // } // #pragma omp parallel
}


/*
  SOME CALLS TO HELP ADAT READ
*/
// Read a correlator edb
// std::map<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators>
// ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators>
std::vector<typeCorr>
getCorrs(const std::vector<std::string>& dbases, std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& fetch)
{
  typedef Hadron::KeyHadronSUNNPartNPtCorr_t       K;
  typedef ENSEM::VectorComplex                     V;

  // Initialize map returned to main associating keys with correlators datatypes
  // n.b. cannot use a std::map here as I haven't created a comparator
  // ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> retMap;

  // vector of typeCorrs to return
  std::vector<typeCorr> ret( (temporal3pt.tx3pt-temporal3pt.tm3pt)/temporal3pt.ts3pt + 1 );


  // Make a vector large enough to hold just desired keys in fetch
  std::vector<ADATIO::SerialDBKey<K> > keyFetch(fetch.size());
  std::cout << "Will search for " << keyFetch.size() << " keys from databases" << std::endl;
  std::vector< std::vector<ADATIO::SerialDBData<V> > > dataFetch(keyFetch.size());

  // Now run over the keys
  for ( std::vector<ADATIO::SerialDBKey<K> >::iterator k = keyFetch.begin(); k != keyFetch.end(); ++k )
    {
      int cnt = k-keyFetch.begin();
      (*k).key() = fetch[cnt]; // set the key

      int ret;
      // Run over the databases
      for ( std::vector<std::string>::const_iterator dbase = dbases.begin(); dbase != dbases.end(); ++dbase)
	{

	  // Initialize arbitrary DB structure
	  FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<K>,  ADATIO::SerialDBData<V> > database;
	  
	  // Check the current db even exists
	  if (database.open(*dbase, O_RDONLY, 0400) != 0)
	    {
	      std::cerr << __func__ << ": error opening dbase= " << *dbase << std::endl;
	      exit(1);
	    }

	  ret = database.get(*k,dataFetch[cnt]);
	  // Attempt to access the key in the current dbase
	  if ( ret == 0 )
	    {
	      std::cout << "Successfully fetched key = " << k->key() << std::endl;
	      database.close();
	      break;
	    }
	  else {
	    continue;
	  }
	  database.close();
	} // end dbase


      // Abort if all dbs have been searched and current key has not been found
      if ( ret != 0 )
	{      
	  std::cerr << __func__ << ": key not found\n" << k->key();
	  exit(1);
	}
    } // end k


  // Hold the number of cnfgs in ensemble
  int enscfgs = dataFetch[0].size();
  // String stream for interpreting ENSEM::peekObs data
  std::stringstream sIn;

  // With all the keys found and SerialDBData obtained, parse into correlators form
  for ( int k = 0; k < keyFetch.size(); ++k )
    {
      // Make some array to associate with this key
      correlators ensAdat;

      ensAdat.ncor[0].real.resize(enscfgs);
      ensAdat.ncor[0].imag.resize(enscfgs);
      // Parse data for all configs for this key
      for ( int g = 0; g < dataFetch[k].size(); ++g )
	{
	  // Loop over the data elems for the config of this key
	  for ( int dd = 0; dd < dataFetch[k][g].data().numElem(); ++dd )
	    {
	      double r_, i_;
	      sIn << ENSEM::peekObs(ENSEM::real(dataFetch[k][g].data()),dd);
	      sIn >> r_; sIn.clear(); sIn.str(std::string());
	      sIn << ENSEM::peekObs(ENSEM::imag(dataFetch[k][g].data()),dd);
	      sIn >> i_; sIn.clear(); sIn.str(std::string());

	      ensAdat.ncor[0].real[g].push_back(r_);
	      ensAdat.ncor[0].imag[g].push_back(i_);
	    }
	}
      // Insert this key-correlator combo
      // retMap.insert(keyFetch[k].key(),ensAdat);

      ret[ (keyFetch[k].key().npoint[1].t_slice - temporal3pt.tm3pt)/temporal3pt.ts3pt ].
	corrs.insert(keyFetch[k].key(),ensAdat);
		

    }

  // // Read all the keys and values
  // std::vector< ADATIO::SerialDBKey<K> >                  keys;
  // std::vector< std::vector< ADATIO::SerialDBData<V> > >  vals;
  // database.keysAndData(keys, vals);
  // // Put into something more useful
  // ADAT::MapObject<K, std::vector<V> > dest;
  // for(int i = 0; i < keys.size(); ++i) {
  //   std::vector<V> v(vals[i].size());
  //   for(int n=0; n < v.size(); ++n) {
  //     v[n] = vals[i][n].data();
  //   }
  //   dest.insert(keys[i].key(), v);
  //   std::cout << keys[i].key() << std::endl;
  // }

  // return retMap;
  return ret;
}


// Reader for global properties
void read(XMLReader& xml, const std::string path, global& g)
{
  try {
    read(xml, path+"/ensem", g.ensem);
    read(xml, path+"/cfgs", g.cfgs);
    read(xml, path+"/t2pt", g.t2pt);
    read(xml, path+"/observable", g.observable);
    read(xml, path+"/state", g.state);
    read(xml, path+"/nvec", g.nvec);
    read(xml, path+"/pi", g.pi);
    read(xml, path+"/pf", g.pf);
    read(xml, path+"/OpMap", g.opMomXML);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse global properties " << e << std::endl;
  }
  // Now set the std::map of Moms + Ops
  
}

// Reader for temporal3pt struct
void read(XMLReader& xml, const std::string path, t3pt& t)
{
  try {
    read(xml, path+"/tmin", t.tm3pt);
    read(xml, path+"/tstep", t.ts3pt);
    read(xml, path+"/tmax", t.tx3pt);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse temporal3pt struct from ini xml " << e << std::endl;
  }
}

// Reader for db3ptInfo struct
void read(XMLReader& xml, const std::string path, info3pt& d)
{
  try {
    read(xml, path+"/base", d.base);
    read(xml, path+"/momTag", d.momTag);
    read(xml, path+"/tsnkTag", d.tsnkTag);
    read(xml, path+"/t0Tag", d.t0Tag);
    read(xml, path+"/zTag", d.zTag);
    read(xml, path+"/rowinfo", d.rowWeights);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse db3ptInfo struct from ini xml " << e << std::endl;
    exit(1);
  }
}

// Reader for db2ptInfo struct
void read(XMLReader& xml, const std::string path, info2pt& d)
{
  try {
    read(xml, path+"/base", d.base);
    read(xml, path+"/momTag", d.momTag);
    read(xml, path+"/t0Tag", d.t0Tag);
    read(xml, path+"/rowinfo", d.rowWeights);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse db2ptInfo struct from ini xml " << e << std::endl;
    exit(1);
  }
}


/*
  Modify the src/snk row structure of passed key to allow a second RowCombo (RC2) to be accessed
*/
void makeKeyRC2(Hadron::KeyHadronSUNNPartNPtCorr_t &kRC2, XMLArray::Array<int> &newRC)
{
  // newRC2 = db2ptInfo.rows[1]; i.e. an XML::Array<int> object
  kRC2.npoint[1].irrep.irrep_mom.row = newRC[0]; // new snk row
  kRC2.npoint[2].irrep.irrep_mom.row = newRC[1]; // new src row
}


/*
  Accept a std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> of keys, and
  row/sign modifiers and create new keys to append
*/
void makeKeyRCs(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& v,
		ADAT::MapObject<XMLArray::Array<int>, double>& r)
{
  // Make r.size local copies of the vector of keys
  std::vector<std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> > modV(r.size(),v);

  // Iterate over the *NEW* (starting at idx1) different row combinations
  auto ri = r.begin();
  while ( ri != r.end() )
    {
      // Perform a quick check on first key in &v to see if this row combo (ri) exists,
      // if so skip inputting this row combo
      if ( modV[0][0].npoint[1].irrep.irrep_mom.row == ri->first[0]
	   && modV[0][0].npoint[modV[0][0].npoint.size()].irrep.irrep_mom.row == ri->first[1] )
	{
	  continue;
	}
      else {
	int idx = std::distance(r.begin(),ri);
	for ( auto i = modV[idx].begin(); i != modV[idx].end(); ++i )
	  {
	    // Modify the 1st npoint struct
	    i->npoint[1].irrep.irrep_mom.row = ri->first[0];
	    // Modify the 2nd npoint struct, unless we are iterating through a 3pt key
	    if ( i->npoint.size() == 2 ) { i->npoint[2].irrep.irrep_mom.row = ri->first[1]; }
	    else { i->npoint[3].irrep.irrep_mom.row = ri->first[1]; }
	  }
	// Now tack these modified keys onto the end of v
	v.insert(v.end(), modV[idx].begin(), modV[idx].end());

      } // end check for existing row combo in first entry of &v

      ri++; // iterate
    }
}
  

/*
  Utility to make 2pt db list
*/
std::vector<std::string> makeDBList(global& g, info2pt& I, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<std::string> s;
  
  // Get the momentum of the template key
  Array<int> thisMom = kTemplate->npoint[1].irrep.irrep_mom.mom;

  // A temporary db to set
  std::string tmp_db_plus = I.base+"/"+g.ensem+"/"+I.t0Tag+"/momXYZ."+shortMom(thisMom,".") \
    +"/EDB/"+g.ensem+"."+g.state+"_p"+shortMom(thisMom,"")+".n"+g.nvec+"."+I.t0Tag+".edb";
      
  // Push to db list
  s.push_back(tmp_db_plus);
  // Make a db name for a potential negative momentum
  if ( shortMom(thisMom,"") != "000" )
    {
      Array<int> thisMomMinus = thisMom*-1;
      // A temporary db to set
      std::string tmp_db_minus = I.base+"/"+g.ensem+"/"+I.t0Tag+"/momXYZ."+shortMom(thisMomMinus,".") \
	+"/EDB/"+g.ensem+"."+g.state+"_p"+shortMom(thisMomMinus,"")+".n"+g.nvec+"."+I.t0Tag+".edb";
      // Push to db list
      s.push_back(tmp_db_minus);
    }
  return s;
}

/*
  Utility to make 3pt db list
*/
std::vector<std::string> makeDBList(global& g, info3pt& I, t3pt& t, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<std::string> s;
  for ( int ti = t.tm3pt; ti <= t.tx3pt; ti+=t.ts3pt )
    {
      // A temporary db to set
      std::string tmp_db_plus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/snk"+shortMom(g.pf,".") \
	+"_src"+shortMom(g.pi,".")+"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+shortMom(g.pf,"") \
	+"_pi"+shortMom(g.pi,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
	I.tsnkTag+ti+"."+I.zTag+".edb";
	  
      // Push to db list
      s.push_back(tmp_db_plus);
      // Check for any non-trivial 3pt momenta and make new db file names
      if ( shortMom(g.pf,"") != "000" || shortMom(g.pi,"") != "000" )
	{
	  // A temporary db to set
	  std::string tmp_db_minus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/snk"+shortMom(g.pf*-1,".") \
	    +"_src"+shortMom(g.pi*-1,".")+"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+shortMom(g.pf*-1,"") \
	    +"_pi"+shortMom(g.pi*-1,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
	    I.tsnkTag+ti+"."+I.zTag+".edb";
	  // Push to db list
	  s.push_back(tmp_db_minus);
	}
    }
  std::cout << "  Made 3pt db list" << std::endl;
  return s;
}


/*
  Utility to make 2pt key list
*/
std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> makeKeyList(Hadron::KeyHadronSUNNPartNPtCorr_t &kTemplate)
{
  std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> s;

  s.push_back(kTemplate);
  // Attempt to add the row2-row2 combo, but skip if at rest
  if ( shortMom(kTemplate.npoint[1].irrep.irrep_mom.mom,"") != "000" )
    {
      // Assuming only a single 2pt key populates s, and that it is a row1-row1 piece
      s.push_back(kTemplate);
      s[1].npoint[1].irrep.irrep_mom.row = 2;
      s[1].npoint[2].irrep.irrep_mom.row = 2;
    }
  checkKeyMomNegate(s,props);
  return s;
}


/*
  Utility to make 3pt key list
*/
std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> makeKeyList(t3pt& t, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> s;

  // With template 3pt key passed as pointer, iterate over tseps and make remaining 3pt keys
  for ( int ti = t.tm3pt; ti <= t.tx3pt; ti+=t.ts3pt )
    {
      // Make a new key from the template
      Hadron::KeyHadronSUNNPartNPtCorr_t tmp_k = *kTemplate;
      
      // Change sink t_slice
      tmp_k.npoint[1].t_slice=ti;
      // Push this tmp key
      s.push_back(tmp_k);

      // Make a key for opposing displacement, if disp_list != ''
      if ( tmp_k.npoint[2].irrep.op.ops[1].disp_list.size() > 0 )
	{
	  for ( int d = 0; d < tmp_k.npoint[2].irrep.op.ops[1].disp_list.size(); ++d )
	    {
	      tmp_k.npoint[2].irrep.op.ops[1].disp_list[d]*=-1;
	    }
	  // insert this key with opposite displacement
	  s.push_back(tmp_k);
	}
    }

  // Negate all momenta in keys
  checkKeyMomNegate(s,props);
  std::cout << "  Set all 3pt keys from template key" << std::endl;
  return s;
}

/*
  Utility to retrive & order correlators stored in keycorrMap MapObjects
*/
void pop3ptFuncs(std::vector<typeCorr> &F, ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators>& kc,
		 t3pt& t)
{
  // auto k = kc.begin();
  // std::vector<bool> truth(kc.size(),false); // track what keys have been touched
  // while ( k != kc.end() )
  //   {
  //     std::cout << "New 3pt k iterator location" << std::endl;
  //     int iterLoc = std::distance(kc.begin(),k); // location of iterator k
  //     if ( !truth[iterLoc] ) // if false, haven't yet touched key
  // 	{
  // 	  // Mark this key as touched
  // 	  truth[iterLoc] = true;
  // 	  int cntr = 0; // counter of # of correlators of same tsep
	  
  // 	  // Info on this correlator
  // 	  int thisT = k->first.npoint[1].t_slice;
  // 	  int fIdx = (thisT/2)-2; // index of F
  // 	  std::cout << " fIdx = " << fIdx << std::endl;
  // 	  std::cout << " thisT = " << thisT << std::endl;
	  
  // 	  // Set a corrs member of this typeCorr to kc correlator
  // 	  F[fIdx].corrs.ncor[cntr] = k->second.ncor[0];
  // 	  cntr++;
  // 	  std::cout << "Run over remaining keys" << std::endl;
  // 	  std::cout << "Size of this fIdx corrs = " << F[fIdx].corrs.ncor.size() << std::endl;
  // 	  // Run over remaining keys and find those with a matching tsep
  // 	  auto kk = kc.begin();
  // 	  std::advance(kk, std::distance(kc.begin(),k)+1 );
  // 	  while ( kk != kc.end() )
  // 	    {
  // 	      std::cout << "Got into kk loop : tslice = " << kk->first.npoint[1].t_slice << std::endl;
  // 	      if ( kk->first.npoint[1].t_slice == thisT )
  // 		{
  // 		  std::cout << "Attempting to set another" << std::endl;
  // 		  F[fIdx].corrs.ncor[cntr] = kk->second.ncor[0];
  // 		  std::cout << "Set the other" << std::endl;
  // 		  truth[std::distance(kc.begin(),kk)] = true;
  // 		  cntr++;
  // 		}
  // 	      ++kk;
  // 	    }
  // 	}
	  
  //     // Iterate k
  //     k++;
  //   }
  // std::cout << "Completed iteration through kc" << std::endl;


  // // Loop through elements of F (i.e. key/corrs of same tsep)
  // for ( std::vector<typeCorr>::iterator it = F.begin(); it != F.end(); ++it )
  //   {

  //     for ( int j = 0; j < it->corrs.ncor.size(); j++ )
  // 	{
  // 	  arr_print(it->corrs.ncor[j]);
  // 	}
  //     exit(9);
  //   }





  // for ( std::vector<corrFunc>::iterator f = F.begin(); f != F.end(); ++f )
  //   {
  //     // Local copy of this tsep & snk mom
  //     int t = f->getT();
  //     std::string mom = shortMom(h.npoint[1].irrep.irrep_mom.mom,"");
      
  //     // Set a key to the template key
  //     Hadron::KeyHadronSUNNPartNPtCorr_t sortKey = h;
  //     sortKey.npoint[1].t_slice=t;

  //     // Momentum combinations to lookup
  //     std::vector<Array<int> > momLookup(1,sortKey.npoint[1].irrep.irrep_mom.mom);
  //     if ( mom != "000" )
  // 	{
  // 	  Array<int> tmp = momLookup[0]; tmp *= -1;	  
  // 	  momLookup.push_back(tmp);
  // 	}

  //     // Displacement combinations to lookup
  //     std::vector<std::vector<int> > dispLookup(1,sortKey.npoint[2].irrep.op.ops[1].disp_list);
  //     if ( dispLookup[0].size() > 0 )
  // 	{
  // 	  std::vector<int> tmp;
  // 	  for ( auto dd = dispLookup[0].begin(); dd != dispLookup[0].end(); ++dd )
  // 	    {
  // 	      tmp.push_back( (*dd)*-1 );
  // 	    }
  // 	  dispLookup.push_back(tmp);
  // 	}
      

  //     // Now we know how many correlators need to be merged
  //     int num2Merge = momLookup.size()*dispLookup.size();
  //     correlators toMerge(num2Merge);



  //     int cnt = 0; // talley of correlator included in toMerge total
  //     for ( int p = 0; p < momLookup.size(); ++p )
  //       {
  //         for ( int z = 0; z < dispLookup.size(); ++z )
  //           {
  //             // Set remainder of sortKey
  //             sortKey.npoint[1].irrep.irrep_mom.mom = momLookup[p];
  //             sortKey.npoint[3].irrep.irrep_mom.mom = momLookup[p];
  //             sortKey.npoint[2].irrep.op.ops[1].disp_list = dispLookup[z];

  //             // Retrieve the correlators datatype associated with this key
  //             ensemble retV = kc[sortKey].ncor[0];


  // 	      // Sloppy capture of the sign of the ioffe time here!
  // 	      int ioffe;
  // 	      if ( dispLookup[z].empty() ) { ioffe = 0; }
  // 	      else { ioffe = momLookup[p]*dispLookup[z]; }

  // 	      // Conjugate if ioffe < 0
  //             if ( ioffe < 0 )
  //               {
  //                 cmplxConj(retV);
  //               }

  //             // Tack the retrieved correlators datatype onto toMerge object
  //             toMerge.ncor[cnt] = retV;
  //             cnt++;
  //           } // end z
  //       } // end p

  //     // Now we have all correlators at same tsep, and with same abs(pz) abs(z)
  //     // So now we merge them
  //     f->data = mergeCorrelators(toMerge);
  //   } // end F iterator
}


void rowAvg(std::vector<typeCorr> &tC, ADAT::MapObject<XMLArray::Array<int>, double> &r)
{
  std::cout << "New typeCorr to average" << std::endl;
  // Iterate over std::vec components - i.e. time slices
  for ( std::vector<typeCorr>::iterator it = tC.begin(); it != tC.end(); ++it )
    {
      // For this tsep, organize keys within each pz/z channel
      Hadron::KeyHadronSUNNPartNPtCorr_t kAvg = it->corrs.begin()->first;
      ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> cAvg;

      std::cout << "--- Running through this typeCorr" << std::endl;
      // Hold z's & p's as temporaries
      XMLArray::Array<int> _pi, _pf, _q;
      std::vector<int> _d;
      _pf = kAvg.npoint[1].irrep.irrep_mom.mom;
      _q  = kAvg.npoint[2].irrep.irrep_mom.mom;
      _pi = kAvg.npoint[3].irrep.irrep_mom.mom;
      _d  = kAvg.npoint[2].irrep.op.ops[1].disp_list;

      // Try to negate momentum & z's, and fetch new keys
      for ( int i = 1; i > -2; i-=2 )
	{
	  for ( int j = 1; j > -2; j-=2 )
	    {
	      std::cout << "------ Swap pair = ( " << i << " , " << j << " )" << std::endl;
	      kAvg.npoint[1].irrep.irrep_mom.mom = _pf * i; // allow swap of snk mom
	      kAvg.npoint[2].irrep.irrep_mom.mom = _q * i; // allow swap of ins mom
	      kAvg.npoint[3].irrep.irrep_mom.mom = _pi * i; // allow swap of src mom
	      // Since we are accessing a key from the tc[*].corrs map (unordered)
	      // We have no guarantee the operator name is correct after setting opposing momenta
	      // So we enforce the correct name here...
	      kAvg.npoint[1].irrep.op.ops[1].name =
		props.opMomXML[shortMom(kAvg.npoint[1].irrep.irrep_mom.mom,"")];
	      kAvg.npoint[3].irrep.op.ops[1].name =
		props.opMomXML[shortMom(kAvg.npoint[3].irrep.irrep_mom.mom,"")];
	      

	      // Allow swap of z
	      for ( int d = 0; d != _d.size(); ++d )
		{ 
		  kAvg.npoint[2].irrep.op.ops[1].disp_list[d] = _d[d] * j;
		}

	      
	      // Collect the correlators to average
	      correlators toAvg(r.size());

	      for ( auto rows = r.begin(); rows != r.end(); ++rows )
		{
		  kAvg.npoint[1].irrep.irrep_mom.row = rows->first[0];
		  kAvg.npoint[3].irrep.irrep_mom.row = rows->first[1];

		  std::cout << "*****Key = " << kAvg << std::endl;

		  // Map matelem btwn lattice LG rows to helicity amplitudes
		  std::complex<double> weight = lgToHelicityAmps(&kAvg);

		  std::cout << "*****Weight = " << weight << std::endl;
		  std::cout << "\n\n";

		  // // Fetch this correlator, applying weight in line
		  // toAvg.ncor[std::distance(r.begin(),rows)] = rows->second * it->corrs[kAvg].ncor[0];
		}
	      
	      // Set row values to "-60" to indicate averaging
	      kAvg.npoint[1].irrep.irrep_mom.row = -60;
	      kAvg.npoint[3].irrep.irrep_mom.row = -60;
	      // Insert this average key and result of merged correlators
	      cAvg.insert(kAvg, mergeCorrelators(toAvg));
	    } // for zswap (j)
	} // for 3-mom swap (i)

      // With all swapping done, clear std::vector<typeCorr>[*it].corrs and refill it with "cAvg"
      it->corrs.clear();
      it->corrs = cAvg;

    } // for std::vector<typeCorr> iterator
}


/*
  MUST INCLUDE THIS DUMB READER FOR INSTANCES WHERE 2PT AND 3PT CORRELATORS WERE MADE
  WITH DEVEL AND DEVEL-PDF BRANCHES OF ADAT, RESPECTIVELY!
*/
#ifndef HAVE_DISPLIST_2PT

bool checkExist(const std::string n)
{
  std::ifstream in(n);
  return in.good();
}


void Reads2pt(std::ifstream &inFile, corrFunc &f2pt, info2pt& I, Hadron::KeyHadronSUNNPartNPtCorr_t *h)
{
  // // Look and check for existence of potentially several files
  // for ( auto d = I.base.begin(); d != I.base.end(); ++d )
  //   {  --> *d
      std::string dumb = I.base+"/"+I.momTag+shortMom(h->npoint[1].irrep.irrep_mom.mom,".")
	+"/ENS/"+Hadron::ensemFileName(*h);
      // if ( checkExist(dumb) )
      // 	{
	  // READ HERE
	  read(inFile,dumb,f2pt.data.ncor[0]);
    // 	}
    //   else { continue; }
    // }
}
#endif


/*
  Utility to ensure data has been parsed correctly from dbs
*/
// void parseCheck(ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> A)
void parseCheck(std::vector<typeCorr> T)
{
  for ( auto tt = T.begin(); tt != T.end(); ++tt )
    {
      for ( auto kc = tt->corrs.begin(); kc != tt->corrs.end(); ++kc )
	{
	  std::cout << "Performing simple ensemble average check on key = \n"
		    << "     " << kc->first << std::endl;
	  
	  corrFunc funcDum;
	  funcDum.data = kc->second;
	  std::cout << "We have dimensions: " << funcDum.data.ncor.size() << " X "
		    << funcDum.data.ncor[0].real.size() << " X "
		    << funcDum.data.ncor[0].real[0].size() << std::endl;
	  // arr_print(funcDum.data.ncor[0]);
	  ensemble funcDumAvg = corrAvg(funcDum.data);
	  std::cout << " Now the average" << std::endl;
	  arr_print(funcDumAvg); 
	}
    }
}

int main(int argc, char *argv[])
{
  if ( argc != 2 )
    {
      std::cout << "\nUsage: $0 <ini xml>\n" << std::endl;
      exit(1);
    }

  // Instantiate a string streamer
  std::stringstream ss;

  XMLReader xmlSR(argv[1]);
  std::cout << "Reading from inixml = " << std::string(argv[1]) << std::endl;


  // Register all the necessary factories for relating lattice matrix elements and helicity amplitudes
  // Hadron::IrrepsCubicHelicityEnv::registerAll();
  Hadron::SubduceTablesOctEnv::registerAll();
  Hadron::SubduceTablesLgEnv::registerAll();


  /*
    Determine global properties
  */
  read(xmlSR, "/PGITD/global", props);
  /*
    Determine the db info - sets where/what edbs to search for keys
  */
  read(xmlSR, "/PGITD/dbInfo/threePt/tseries", temporal3pt);
  read(xmlSR, "/PGITD/dbInfo/threePt", db3ptInfo);
  read(xmlSR, "/PGITD/dbInfo/twoPt", db2ptInfo);
  // dumpRowInfo(db3ptInfo.rows,db3ptInfo.signs,3);
  // dumpRowInfo(db2ptInfo.rows,db2ptInfo.signs,2);

  // Set the number of tseps once and for all
  const int nTSeps = 1+(temporal3pt.tx3pt-temporal3pt.tm3pt)/temporal3pt.ts3pt;

  // Structure to hold all desired 3pt/2pt motion/normalizing keys
  struct keydb_t
  {
    std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> keys;
    std::vector<std::string>                        dbs;
    // Map of keys and correlators structs returned
    // std::map<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> keycorrMap; --> need a comparator
    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> keycorrMap;
  };

  struct nptKeysValues_t
  {
    keydb_t keydbs;
  } threePt, twoPi, twoPf;


  struct keyTemplate_t
  {
    Hadron::KeyHadronSUNNPartNPtCorr_t key;
  } tempKey3pt, tempKey2Pi, tempKey2Pf;


  /*
    Read the template keys from the ini xml
  */
  try {
    read(xmlSR, "/PGITD/Ops/threePt/keyTemplate", tempKey3pt.key);
    read(xmlSR, "/PGITD/Ops/twoPt/keyTemplate", tempKey2Pi.key);
    tempKey2Pf = tempKey2Pi;
  } catch ( std::string &e ) {
    std::cerr << "Unable to access key templates from ini xml " << e << std::endl;
  }
  std::cout << "READ THE TEMPLATE KEYS TO FETCH" << std::endl;

  /*
    Substitute correct momenta and src/snk operator names (based on passed global& props)
    into template 3pt & 2pt keys
  */
  setOpsMoms(&tempKey3pt.key, tempKey2Pf.key, tempKey2Pi.key, props);

  


  /*   Make the 3pt database structures to search   */
  threePt.keydbs.dbs    = makeDBList(props, db3ptInfo, temporal3pt, &tempKey3pt.key);
  // Make all 3pt keys from templates
  threePt.keydbs.keys   = makeKeyList(temporal3pt, &tempKey3pt.key);

  dumpKeys(threePt.keydbs.keys, 3);
  exit(9);

  makeKeyRCs(threePt.keydbs.keys, db3ptInfo.rowWeights); // create the other 3pt row combinations
  /* Ensure 3pt keys conserve 3-momentum */
  conserveMom3PtKey(threePt.keydbs.keys);


  // dumpDBs(threePt.keydbs.dbs);
  // dumpKeys(threePt.keydbs.keys, 3);
  // exit(8);
  // threePt.keydbs.keycorrMap    = getCorrs(threePt.keydbs.dbs,threePt.keydbs.keys);
  std::vector<typeCorr> funcs3pt = getCorrs(threePt.keydbs.dbs,threePt.keydbs.keys);

  // rowAvg(funcs3pt, db3ptInfo.rowWeights);

  for ( int i = 0; i != funcs3pt.size(); i++ )
    {
      std::cout << " This corrs map size = " << funcs3pt[i].corrs.size() << std::endl;
    }
  // HERE!!!

  /*   Make the 2pt database structures to search   */
  twoPi.keydbs.dbs      = makeDBList(props, db2ptInfo, &tempKey2Pi.key);
  twoPf.keydbs.dbs      = makeDBList(props, db2ptInfo, &tempKey2Pf.key);
  // dumpDBs(twoPi.keydbs.dbs);
  // dumpDBs(twoPf.keydbs.dbs);
  // Make all 2pt keys from templates
  std::cout << "Initial template 2pt key = " << tempKey2Pi.key << std::endl;
  std::cout << "Final template 2pt key = " << tempKey2Pf.key << std::endl;
  twoPi.keydbs.keys     = makeKeyList(tempKey2Pi.key);
  twoPf.keydbs.keys     = makeKeyList(tempKey2Pf.key);


  std::cout << "Initial 2pt keys" << std::endl;
  dumpKeys(twoPi.keydbs.keys,2);
  std::cout << "Final 2pt keys" << std::endl;
  dumpKeys(twoPf.keydbs.keys,2);

  std::cout << "Done setting 2pt shit" << std::endl;












  // Determine the snk/src rotation angles
  // Hadron::CubicCanonicalRotation_t snk_rot = Hadron::cubicCanonicalRotation(props.pf);
  // Hadron::CubicCanonicalRotation_t src_rot = Hadron::cubicCanonicalRotation(props.pi);

  // // From these angles, get the Wigner-Ds
  // std::complex<double> snkWignerD = Hadron::Wigner_D(int 2*J, int 2*m, int 2*n,
  // 						     snk_rot.alpha, snk_rot.beta, snk_rot.gamma);
  // std::complex<double> srcWignerD = Hadron::Wigner_D(int 2*J, int 2*m, int 2*n,
  // 						     src_rot.alpha, src_rot.beta, src_rot.gamma);

  // // Print the rotation angles
  // std::cout << "Snk mom = " << props.pf << " rotation angles ==> "
  // 	    << " Alpha = " << Hadron::cubicCanonicalRotation(props.pf).alpha
  // 	    << " Beta = " << Hadron::cubicCanonicalRotation(props.pf).beta
  // 	    << " Gamma = " << Hadron::cubicCanonicalRotation(props.pf).gamma << std::endl;
  // for ( int i = 1; i > -2; i-=2 )
  //   {
  //     for ( int j = 1; j > -2; j-=2 )
  // 	{
  // 	  std::cout << "Snk Wigner-D = "
  // 		    << Hadron::Wigner_D(1,i,j,
  // 					Hadron::cubicCanonicalRotation(props.pf).alpha,
  // 					Hadron::cubicCanonicalRotation(props.pf).beta,
  // 					Hadron::cubicCanonicalRotation(props.pf).gamma) << std::endl;
  // 	}
  //   }


  // std::complex<double> Hadron::zeroComplex(std::complex<double> w); // zero fuzz bits of a complex # (in single_hadron_coeffs.h)
  // int twoJ = Hadron::convertSpintoTwoJ(std::string J); // J1 -> twoJ=2 ( in sub_spin_hel_map.h )
  // std::string Hadron::convertTwoJtoHelicity(twoJ); // twoJ -> H

  // exit(9);

  // std::cout << "Src mom = " << props.pi << " rotation angles ==> "
  // 	    << " Alpha = " << Hadron::cubicCanonicalRotation(props.pi).alpha
  // 	    << " Beta = " << Hadron::cubicCanonicalRotation(props.pi).beta
  // 	    << " Gamma = " << Hadron::cubicCanonicalRotation(props.pi).gamma << std::endl;

  /*
    Access the 2pt correlators
  */
// #ifdef HAVE_DISPLIST_2PT
// #warning "Using getCorrs to read 2pt correlators"
//   // twoPi.keydbs.keycorrMap      = getCorrs(twoPi.keydbs.dbs,twoPi.keydbs.keys);
//   // twoPf.keydbs.keycorrMap      = getCorrs(twoPf.keydbs.dbs,twoPf.keydbs.keys);
//   std::vector<typeCorr> twoPtIni = getCorrs(twoPi.keydbs.dbs,twoPi.keydbs.keys);
//   std::vector<typeCorr> twoPtFin = getCorrs(twoPf.keydbs.dbs,twoPf.keydbs.keys);
// #else
//   #warning "   2pt correlators constructed w/o disp_list  -->  reverting to dumb 2pt reader"
//   std::ifstream inFile;
  
//   Reads2pt(inFile, func2ptRC1, db2ptInfo, &twoPt.keydbs.keys[0]);
//   if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
//     Reads2pt(inFile, func2ptRC2, db2ptInfo, &twoPtRC2.keydbs.keys[0]);
// #endif


#if 0
  // parseCheck(threePt.keydbs.keycorrMap); std::cout << "\n";
  // parseCheck(twoPi.keydbs.keycorrMap); std::cout << "\n";
  // parseCheck(twoPf.keydbs.keycorrMap); std::cout << "\n";

  parseCheck(funcs3pt); std::cout << "\n";
  // parseCheck(twoPtIni); std::cout << "\n";
  // parseCheck(twoPtFin); std::cout << "\n";
  exit(9);
#endif
  std::cout << "Successfully performed the reads!" << std::endl;




  /*
    Need to unpack the keys and push into corrFunc instances to reuse all code below

    funcs3pt:
          --> threePt.motion.keys  + momenta
	  --> threePt.motion.keys  - momenta
  */

  // Trying to use this to combine user provided row combinations
  // typeCorr b_b0(tempKey3pt.key,threePt.keydbs.keys.size()/nTSeps,props.cfgs,4);
  // b_b0.rowAvg(db3ptInfo.rows, db3ptInfo.signs);

  // A messy(ish) initialization of these std::vectors
  // std::cout << "The set size of funcs3pt.ncor = " << threePt.keydbs.keys.size()/nTSeps << std::endl;


  // HERE!!!
  // std::vector<typeCorr> funcs3pt(nTSeps,typeCorr(threePt.keydbs.keys.size()/nTSeps,props.cfgs,4));
  // pop3ptFuncs(funcs3pt, threePt.keydbs.keycorrMap, temporal3pt);
  

  // for ( int i = 0; i < funcs3pt.size(); ++i )
  //   {
  //     ensemble funcDumAvg1 = corrAvg(funcs3pt[i].data);
  //     arr_print(funcDumAvg1);
  //     std::cout << "------------------------" << std::endl;
  //     // arr_print(funcDumAvg2);
  //     // std::cout << "------------------------" << std::endl;
  //     // arr_print(funcDumAvg3);
  //     // std::cout << "------------------------" << std::endl;
  //     // arr_print(funcDumAvg4);
  //     // std::cout << "------------------------" << std::endl;
  //     // std::cout << "\n\n" << std::endl;
  //   }
  exit(8);
 

//   /*
//     Initialize the 2pt functions
//   */
//   std::cout << "Initializing the 2pt functions" << std::endl;

  
//   // Start by initializing single master 2pt correlator
//   corrFunc func2pt;
//   // Initialize a correlators database to hold all 2pt correlators to merge
//   correlators toMerge2pts; toMerge2pts.ncor.resize(2); // minimally 1-set to merge (e.g. rest w/ 1 RCs)
//                                                        // But will duplicate RC1 to not change functionality

//   // Create corrFunc objects for rest (or +p mom) RC 1 & 2
//   corrFunc func2ptRC1(props.cfgs,props.t2pt), func2ptRC2(props.cfgs,props.t2pt);
//   func2ptRC1.initialize(); func2ptRC2.initialize();
// #ifndef HAVE_DISPLIST_2PT
// #warning "   2pt correlators constructed w/o disp_list  -->  reverting to dumb 2pt reader"
//   std::ifstream inFile;
//   // Reads2pt(inFile, func2ptRC1, db2ptInfo, &twoPt.keydbs.keys[0]);
//   // if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
//   //   {
//   //     Reads2pt(inFile, func2ptRC2, db2ptInfo, &twoPtRC2.keydbs.keys[0]);
//   //   }
// #endif

//   // Push into toMerge2pts
//   toMerge2pts.ncor[0] = func2ptRC1.data.ncor[0];

//   /*   Now manage the RC2 piece for 2pt rest case   */
//   if ( shortMom(tempKey2Pi.key.npoint[1].irrep.irrep_mom.mom,"") == "000" )
//     {
//       // There is no RC2 for 2pt rest case, so copy 2pt RC1
//       toMerge2pts.ncor[1] = func2ptRC1.data.ncor[0];
//     }
//   else { toMerge2pts.ncor[1] = func2ptRC2.data.ncor[0]; }

  
//   // If shortMom != 000, then we must initialize two more corrFunc objects for -p mom RC 1 & 2
//   if ( shortMom(tempKey2Pi.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
//     {
//       corrFunc func2ptMinusRC1(props.cfgs,props.t2pt),func2ptMinusRC2(props.cfgs,props.t2pt);
//       func2ptMinusRC1.initialize(); func2ptMinusRC2.initialize();
// #ifndef HAVE_DISPLIST_2PT
// #warning "   2pt correlators constructed w/o disp_list  -->  reverting to dumb 2pt reader"
//       std::ifstream inFile;
//       Reads2pt(inFile, func2ptMinusRC1, db2ptInfo, &twoPi.keydbs.keys[1]);
//       // Reads2pt(inFile, func2ptMinusRC2, db2ptInfo, &twoPtRC2.keydbs.keys[1]);
// #endif
//       // Append these to the toMerge2pts object
//       toMerge2pts.ncor.push_back(func2ptMinusRC1.data.ncor[0]);
//       toMerge2pts.ncor.push_back(func2ptMinusRC2.data.ncor[0]);
      
// #if 0
//       std::cout << "2pt Minus Parse Check!" << std::endl;
//       ensemble dum = corrAvg(func2ptMinusRC1.data);
//       arr_print(dum);
//       dum = corrAvg(func2ptMinusRC2.data);
//       arr_print(dum);
// #endif
//     } // if shortMom(2ptkey"--") != 000

// #if 0
//   ensemble dum = corrAvg(func2ptRC1.data);
//   std::cout << "2pt Parse Check!" << std::endl;
//   arr_print(dum);
//   if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
//     { dum = corrAvg(func2ptRC2.data); arr_print(dum); }
//   // exit(9);
// #endif

//   // To get all info, reset func2pt equal to func2ptRC1
//   func2pt = func2ptRC1;
//   std::cout << "---> Averaging all rows/moms 2pt correlators..." << std::endl;
//   std::cout << toMerge2pts.ncor.size() << std::endl;
//   std::cout << toMerge2pts.ncor[0].real.size() << std::endl;
//   std::cout << toMerge2pts.ncor[0].real[0].size() << std::endl;
//   func2pt.data.ncor[0] = mergeCorrelators(toMerge2pts).ncor[0];
//   std::cout << "<--- Finished averaging all rows/moms 2pt correlators..." << std::endl;

//   /*
//     ###########################################################################################
//     BELOW HERE SHOULD BE ALL AVERAGED/CONJUGATED CORRELATORS...
//     ###########################################################################################
//   */

//   // Print averages of 2pt & 3pt correlator
// #if 0
//   ensemble C2pt_avg = corrAvg(func2pt.data);
//   std::vector<ensemble> C3pt_avgs(funcs3pt.size());
//   std::cout << "C2pt_avg = " << std::endl;
//   arr_print(C2pt_avg);
//   std::cout << "******************" << std::endl;
//   for ( int i = 0; i < funcs3pt.size(); ++i )
//     {
//       C3pt_avgs[i] = corrAvg(funcs3pt[i].data);
//       arr_print(C3pt_avgs[i]);
//       std::cout << "*************" << std::endl;
//     }
//   exit(8);
// #endif



//   /*
//     Make 2pt jackknife samples for each correlator
//   */
//   std::cout << "---> Making 2pt jackknife samples..." << std::endl;
//   func2pt.initializeJks();
//   fill_jk_arrays(func2pt.dataJk, func2pt.data);
//   std::cout << "---> Determining 2pt jackknife ensemble averages..." << std::endl;
//   func2pt.initializeJkEnsAvgs();
//   jkEnsAvg(func2pt.dataJk,func2pt.dataJkEnsAvg);

//   /*
//     Make 3pt jackknife samples for each correlator
//   */
//   std::cout << "---> Making 3pt jackknife samples for each correlator..." << std::endl;
//   std::cout << "+++> & Determining ensemble averages per jackknife sample..." << std::endl;
//   // Iterate over all 3pt data
//   //#pragma omp parallel for
// //   {
// // #pragma omp for
//   for ( std::vector<corrFunc>::iterator it = funcs3pt.begin(); it != funcs3pt.end(); ++it )
//     {
//       auto c = it - funcs3pt.begin();
//       //*****************************
//       // Initialize the 3pts
//       (*it).initializeJks();
//       fill_jk_arrays((*it).dataJk,(*it).data);
//       (*it).initializeJkEnsAvgs();
//       jkEnsAvg((*it).dataJk,(*it).dataJkEnsAvg);
//     }
//   // } // #pragma omp parallel

//   std::cout << "<--- Completed jackknifing/jk averaging of 3pt data..." << std::endl;
  

//   std::cout << "---> Performing the summation of current insertions per jackknife ensemble average"
// 	    << std::endl;
//   // Now sum up the insertions for each jackknife ensemble average
//   summation(funcs3pt);
//   std::cout << "<--- Finished the summation of current insertions per jackknife ensemble average"
// 	    << std::endl;



//   // For each summed 3pt function, divide by the 2pt function at the same tsep
//   // per jackknife ensemble average
//   std::vector<correlators> ratioSum(funcs3pt.size());
//   // Form the summed ratio for each jackknife ensemble average
//   //#pragma omp parallel for
// //   {
// // #pragma omp for
//     for ( int i = 0; i < funcs3pt.size(); ++i )
//       {
// 	ratioSum[i].ncor.resize(func2pt.getCfgs());
// 	// Ensure the corresponding tsep is passed
// 	summedRatio(func2pt.dataJkEnsAvg,funcs3pt[i].dataSummedJkEnsAvg,ratioSum[i],funcs3pt[i].getT());
//       }
//   // } // #pragma omp parallel
//   std::cout << "<-- Finished forming jackknife ensemble average summed ratios" << std::endl;



//   /*
//     Get a central value estimate of summed ratios
//   */
//   std::vector<ensemble> ratioSumAvgEst(funcs3pt.size());
//   //#pragma omp parallel for
//   for ( auto it = ratioSumAvgEst.begin(); it != ratioSumAvgEst.end(); ++it )
//     {
//       it->real.resize(1); it->imag.resize(1);
//       double sum_r(0.0), sum_i(0.0);
//       for ( int g = 0; g < func2pt.getCfgs(); ++g )
// 	{
// 	  auto c = it - ratioSumAvgEst.begin();
// 	  sum_r += ratioSum[c].ncor[g].real[0][0];
// 	  sum_i += ratioSum[c].ncor[g].imag[0][0];
// 	}
//       sum_r *= 1.0/func2pt.getCfgs();
//       sum_i *= 1.0/func2pt.getCfgs();
//       it->real[0].push_back( sum_r );
//       it->imag[0].push_back( sum_i );


//       double err_r(0.0), err_i(0.0);
//       for ( int g = 0; g < func2pt.getCfgs(); ++g )
// 	{
// 	  auto c = it - ratioSumAvgEst.begin();
// 	  err_r += pow( ratioSum[c].ncor[g].real[0][0] - sum_r, 2);
// 	  err_i += pow( ratioSum[c].ncor[g].imag[0][0] - sum_i, 2);
// 	}

//       err_r = sqrt( ((1.0*(func2pt.getCfgs()-1))/func2pt.getCfgs())*err_r);
//       err_i = sqrt( ((1.0*(func2pt.getCfgs()-1))/func2pt.getCfgs())*err_i);
      
//       std::cout << "   RATIO SUM AVG EST  " << std::endl;
//       std::cout << " R = " << it->real[0] << " +/- " << err_r
// 		<<" I = " << it->imag[0] << " +/- " << err_i << std::endl;
//     }



//   /*
//     Correct for bias in each jackknife ensemble average estimate of summed ratio
//     & concatenate corrected jackknife ensemble average estimates per tsep
//   */
//   std::vector<correlators> antiJkSummedRatio(funcs3pt.size());
//   std::vector<correlators> catAntiJkSummedRatio(funcs3pt.size());
//   //#pragma omp parallel for
// //   {
// // #pragma omp for
//   for ( std::vector<correlators>::iterator it = antiJkSummedRatio.begin(); it != antiJkSummedRatio.end(); ++it )
//     {
//       // Again catch the iterator value
//       auto c = it-antiJkSummedRatio.begin();
//       // Correct for the bias
//       *it = antiJkRatio(ratioSumAvgEst[c],ratioSum[c]);
//       // Concatenate bias corrected jackknife estimates
//       catAntiJkSummedRatio[c] = catAntiJkEnsRatio(*it);
//     }
//   // } // #pragma omp parallel


//   /*
//     Write out ensemble of bias corrected jackknife ensemble average summed ratios
//   */
//   std::cout << "---> Writing ensemble of bias corrected jackknife ensemble average summed ratios" << std::endl;

//   tempKey3pt.key.npoint[1].t_slice = -60; // set the output ensem - w/ -60 indicating a summed ratio
//   std::string output = Hadron::ensemFileName(tempKey3pt.key);
//   output.erase(output.end()-4,output.end());
//   // std::string suffix = shortMom(threePt.normalize.keys[0].npoint[1].irrep.irrep_mom.mom,"");

  
//   std::ofstream outFile;
//   // outFile.open(output+".summedRatio_p"+suffix+".dat");
//   outFile.open(output+".summedRatio.dat");
//   outFile << func2pt.getCfgs() << " " << std::to_string(funcs3pt.size()) << " 1 0 1\n";
//   for ( int g = 0; g < func2pt.getCfgs(); ++g )
//     {
//       for ( std::vector<correlators>::iterator it = catAntiJkSummedRatio.begin();
// 	    it != catAntiJkSummedRatio.end(); ++it )
// 	{
// 	  // Catch the iterator value
// 	  auto c = it - catAntiJkSummedRatio.begin();
// 	  outFile << std::setprecision(10) << std::to_string(c) << " "
// 		  << (*it).ncor[0].real[g][0] << " " << (*it).ncor[0].imag[g][0] << "\n";
// 	}
//     }
//   outFile.close();
//   std::cout << "<--- Completed writing of ensemble of bias corrected jackknife ensemble average summed ratios"
// 	    << std::endl;


// // #if 0
// //   /*
// //     Create jackknife samples from bias corrected ensemble for each tsep
// //     & write to file for XMBF fitting
// //   */
// //   std::vector<correlators> jkRatioSamples(funcs3pt.size());
// //   for ( auto it = jkRatioSamples.begin(); it != jkRatioSamples.end(); ++it )
// //     {
// //       // Catch the iterator value
// //       auto c = it - jkRatioSamples.begin();
// //       std::string output = "summedRatio_tsnk"+std::to_string(funcs3pt[c].getT());
// //       it->ncor.resize(func2pt.getCfgs());
// //       for ( int j = 0; j < (*it).ncor.size(); ++j )
// // 	{
// // 	  (*it).ncor[j].real.resize(func2pt.getCfgs()-1);
// // 	  (*it).ncor[j].imag.resize(func2pt.getCfgs()-1);
// // 	  for ( int g = 0; g < (*it).ncor[0].real.size(); ++g )
// // 	    {
// // 	      // We are concatenating summedRatio jackknife samples here,
// // 	      // so insertion time dependence is removed - i.e. only 1 entry for each real/imag
// // 	      (*it).ncor[j].real[g].resize(1);
// // 	      (*it).ncor[j].imag[g].resize(1);
// // 	    }
// // 	}
      
// //       fill_jk_arrays(*it,catAntiJkSummedRatio[c]);

// //       // Write out the newly created jackknife samples here
// //       std::cout << "---> Writing jackknife ratio files" << std::endl;
// //       for ( int j = 0; j < func2pt.getCfgs(); ++j )
// // 	{
// // 	  outFile.open(output+"_jack"+std::to_string(j)+".dat");
// // 	  outFile << func2pt.getCfgs()-1 << " 1" << " 1 0 1\n";
// // 	  for ( int g = 0; g < (*it).ncor[0].real.size(); ++g )
// // 	    {
// // 	      outFile << std::setprecision(10) << "0 " << (*it).ncor[j].real[g][0]
// // 		      << " " << (*it).ncor[j].imag[g][0] << "\n";
// // 	    }
// // 	  outFile.close();
// // 	}
// //     }
// //   std::cout << "<--- All jackknife ratio files written..." << std::endl;  
// // #endif

  return 0;
}
