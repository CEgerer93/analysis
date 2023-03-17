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

// #include "summation.h"
#include "operators.h"
#include "pseudo_utils.h"
// #include "pseudo_structs.h"

#include "corr_utils.h"
#include "fit_util.h"
#include "threept_tr.h"

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

#include "hadron_npt_op_factory.h"
#include "single_hadron_factory.h"
#include "cont_hadron_op_factory.h"
#include "operators/spin_subduced_obj_factory.h"

// #include "hadron/irreps_cubic_factory.h"
// #include "hadron/irreps_cubic_oct_factory.h"
// #include "hadron/irreps_cubic_helicity_factory.h"
// #include "hadron/subduce_tables_oct_factory.h"
// #include "hadron/subduce_tables_lg_factory.h"
// #include "hadron/subduce_tables.h"
// #include "hadron/subduce_tables_factory.h"
// #include "hadron/single_hadron_coeffs.h"


using namespace ADATXML;
using namespace Pseudo;
using namespace NCOR;


typedef XMLArray::Array<XMLArray::Array<int> > XML2D;
typedef Hadron::KeyHadronSUNNPartNPtCorr_t       K;
// typedef ENSEM::VectorComplex                     V;
typedef std::vector<std::complex<double> > V;


// Define a struct for global props
global_t global;

// Define some structs to help determine what edbs to read
domain_t temporal3pt, temporal2ptFin, temporal2ptIni, temporal2ptRest;
info3pt db3ptInfo;
info2pt db2ptFinInfo, db2ptIniInfo, db2ptRestInfo;

// Constants
const std::complex<double> redFact(sqrt(2),0);


/*
  SOME CALLS TO HELP ADAT READ
*/
// Read a correlator edb
// std::vector<NCOR::corrEquivalence>
std::unordered_map<std::string, std::vector<NCOR::corrEquivalence> >
getCorrs(const std::vector<std::string>& dbases, std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& fetch)
{
  // Initialize map returned to main associating keys with correlators datatypes
  // n.b. cannot use a std::map here as I haven't created a comparator
  // ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> retMap;

  // vector of corrEquivalences forall 3pt tseps to return
  // std::vector<NCOR::corrEquivalence> ret( temporal3pt.numT() );




  // Order types of keys based on insertion
  std::unordered_map<std::string, std::vector<NCOR::corrEquivalence> > insKeysMap;


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

	  // If inserted operator name is 'gamma_z' then need another rho_rhoxDA__J1_T1mP
	  if ( k->key().npoint[2].irrep.op.ops[1].name == "gamma_z" )
	    {
	      ADATIO::SerialDBKey<K> rhoOpRow = *k;
	      std::vector<ADATIO::SerialDBData<V> > rhoDatRow2Fetch;

	      rhoOpRow.key().npoint[2].irrep.op.ops[1].name = "rho_rhoxDA__J1_T1mP";
	      rhoOpRow.key().npoint[2].irrep.irrep_mom.row  = 2;

	      ret = database.get(rhoOpRow,rhoDatRow2Fetch);

	      for ( int g = 0; g < rhoDatRow2Fetch.size(); ++g )
		{
		  std::vector<std::complex<double> > foo;
		  foo = rhoDatRow2Fetch[g].data();
		  // foo *= _I_;
		  foo *= -1.0;

		  ADATIO::SerialDBData<V> fooDBData;
		  fooDBData.data() = foo;

		  dataFetch[cnt].push_back(fooDBData);
		}
	    }
	  // If inserted operator name is 'gamma_x' or 'gamma_y' then linear combinations
	  // ..of rho_rhoxDA__J1_T1mP are needed
	  else if ( k->key().npoint[2].irrep.op.ops[1].name == "gamma_x" ||
		    k->key().npoint[2].irrep.op.ops[1].name == "gamma_y" )
	    {
	      
	      ADATIO::SerialDBKey<K> rhoOpRow = *k;
	      std::vector<ADATIO::SerialDBData<V> > rhoDatRow1Fetch;
	      std::vector<ADATIO::SerialDBData<V> > rhoDatRow3Fetch;
	      // std::vector< std::vector< std::complex<double> > > rhoDatRow1Fetch, rhoDatRow3Fetch;
	      
	      rhoOpRow.key().npoint[2].irrep.op.ops[1].name = "rho_rhoxDA__J1_T1mP";
	      rhoOpRow.key().npoint[2].irrep.irrep_mom.row  = 1;

	      // rhoDatRow1Fetch = database[rhoOpRow];
	      ret = database.get(rhoOpRow,rhoDatRow1Fetch);
	      
	      rhoOpRow.key().npoint[2].irrep.irrep_mom.row  = 3;
	      ret = database.get(rhoOpRow,rhoDatRow3Fetch);
#if 1	      
#warning "Trying arithmetic"
	      // Try some arithmetic
	      for ( int g = 0; g < rhoDatRow1Fetch.size(); ++g )
		{
		  std::vector<std::complex<double> > foo;
		  if ( k->key().npoint[2].irrep.op.ops[1].name == "gamma_x" )
		    {
		      foo = rhoDatRow1Fetch[g].data() - rhoDatRow3Fetch[g].data();
		      foo *= (1.0/sqrt(2));
		    }
		  if ( k->key().npoint[2].irrep.op.ops[1].name == "gamma_y" )
		    {
		      foo = rhoDatRow1Fetch[g].data() + rhoDatRow3Fetch[g].data();
		      foo *= (_mI_/sqrt(2));
		    }

		  ADATIO::SerialDBData<V> fooDBData;
		  fooDBData.data() = foo;
		  
		  dataFetch[cnt].push_back(fooDBData);
		}
#endif
	    }
	  else
	    {
	      ret = database.get(*k,dataFetch[cnt]);
	    }


	  // Attempt to access the key in the current dbase
	  if ( ret == 0 )
	    {
	      std::cout << "Successfully fetched key = " << k->key() << std::endl;
	      database.close();
	      break;
	    }
	  else
	    continue;
	  database.close();
	} // end dbase

      std::cout << "DataFetch[cnt].size = " << dataFetch[cnt].size() << std::endl;


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
      NCOR::correlator ensAdat;

      ensAdat.ensemble.ens.resize(enscfgs);
      // ensAdat.ncor[0].imag.resize(enscfgs);
      // Parse data for all configs for this key
      for ( int g = 0; g < dataFetch[k].size(); ++g )
	{
	  // Loop over the data elems for the config of this key
	  for ( int dd = 0; dd < dataFetch[k][g].data().size(); ++dd )
	    {
	  // // Loop over the data elems for the config of this key
	  // for ( int dd = 0; dd < dataFetch[k][g].data().numElem(); ++dd )
	  //   {
	      double r_, i_;

	      std::complex<double> dc = dataFetch[k][g].data()[dd];

	      // sIn << ENSEM::peekObs(ENSEM::real(dataFetch[k][g].data()),dd);
	      // sIn >> r_; sIn.clear(); sIn.str(std::string());
	      // sIn << ENSEM::peekObs(ENSEM::imag(dataFetch[k][g].data()),dd);
	      // sIn >> i_; sIn.clear(); sIn.str(std::string());
	      // std::complex<double> dc(r_,i_);

	      ensAdat.ensemble.ens[g].push_back(dc);
		// .real[g].push_back(r_);
	      // ensAdat.ncor[0].imag[g].push_back(i_);
	    }
	}
      // // Insert this key-correlator combo
      // RET[ (keyFetch[k].key().npoint[1].t_slice - temporal3pt.min)/temporal3pt.step ].
      // 	keyCorrMap.insert(keyFetch[k].key(),ensAdat);


      // NEW - Something hacky
      Hadron::KeyHadronSUNNPartNPtCorr_t tmp = keyFetch[k].key(); //convenience
      
      std::string insertedOp = tmp.npoint[2].irrep.op.ops[1].name;

#if 1
      // Insert new key-value pair if insertedOp doesn't yet exist
      if ( insKeysMap.count(insertedOp) == 0 )
	insKeysMap.insert(std::pair<std::string, std::vector<NCOR::corrEquivalence> >
			  (insertedOp,std::vector<NCOR::corrEquivalence> (db3ptInfo.rows.size())));


      if ( tmp.npoint[1].irrep.irrep_mom.row == 1 && tmp.npoint[3].irrep.irrep_mom.row == 1 )
	insKeysMap[insertedOp][0].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      if ( tmp.npoint[1].irrep.irrep_mom.row == 1 && tmp.npoint[3].irrep.irrep_mom.row == 2 )
	insKeysMap[insertedOp][1].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      if ( tmp.npoint[1].irrep.irrep_mom.row == 2 && tmp.npoint[3].irrep.irrep_mom.row == 1 )
	insKeysMap[insertedOp][2].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      if ( tmp.npoint[1].irrep.irrep_mom.row == 2 && tmp.npoint[3].irrep.irrep_mom.row == 2 )
	insKeysMap[insertedOp][3].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
#else
#warning "Trying Robert's idea"
      // Insert new key-value pair if insertedOp doesn't yet exist
      if ( insKeysMap.count(insertedOp) == 0 )
	insKeysMap.insert(std::pair<std::string, std::vector<NCOR::corrEquivalence> >
			  (insertedOp,std::vector<NCOR::corrEquivalence> (2)));


      if ( tmp.npoint[1].irrep.irrep_mom.row == 1 && tmp.npoint[3].irrep.irrep_mom.row == 1 )
	insKeysMap[insertedOp][0].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      else if ( tmp.npoint[1].irrep.irrep_mom.row == 2 && tmp.npoint[3].irrep.irrep_mom.row == 2 )
	insKeysMap[insertedOp][1].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      else
	continue;

#endif
    } // keyFetch iterator

  return insKeysMap;
}

// Reader for npt_props
void read(XMLReader& xml, const std::string path, global_t::nptProps_t& p)
{
  try {
    read(xml, path+"/creation_op", p.create);
    read(xml, path+"/smearedP", p.smear);
    read(xml, path+"/flavor", p.cgc);
    read(xml, path+"/names", p.cont_names);
    read(xml, path+"/disp_list", p.disp_list);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse nptProps_t of global properties " << e << std::endl;
  }
}

// Reader for global properties
void read(XMLReader& xml, const std::string path, global_t& g)
{
  try {
    read(xml, path+"/ensem", g.ensem);
    read(xml, path+"/cfgs", g.cfgs);
    read(xml, path+"/t2ptRows", g.t2ptRows);
    read(xml, path+"/observable", g.observable);
    read(xml, path+"/state", g.state);
    read(xml, path+"/projection", g.projector);
    read(xml, path+"/basic_op_structure/Irrep", g.basic_op);
    read(xml, path+"/npt_props/snk", g.snk);
    read(xml, path+"/npt_props/ins", g.ins);
    read(xml, path+"/npt_props/src", g.src);
    // read(xml, path+"/insertions/props", g.ins);
    read(xml, path+"/nvec", g.nvec);
    read(xml, path+"/Lt", g.Lt);
    read(xml, path+"/Lx", g.Lx);
    read(xml, path+"/rest", g.rest);
    read(xml, path+"/pfpi/pf", g.pf);
    read(xml, path+"/pfpi/pi", g.pi);
    read(xml, path+"/disp_list", g.disp_list);
    read(xml, path+"/dispNegate", g.dispNegate);
    read(xml, path+"/momNegate", g.momNegate);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse global properties " << e << std::endl;
  }
}

// Reader for temporal domains of data & any fits
void read(XMLReader& xml, const std::string path, domain_t& t)
{
  try {
    read(xml, path+"/tmin", t.min);
    read(xml, path+"/tstep", t.step);
    read(xml, path+"/tmax", t.max);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse domain_t struct from ini xml " << e << std::endl;
  }
}

void setType(const std::string type, NCOR::fitInfo_t& I)
{
  try {
    std::cout << "Read type string = " << type << std::endl;
    auto it = table.find(type);
    if ( it != table.end() )
      I.type = it->second;
  } catch ( std::string &e ) {
    std::cerr << "Whoops" << e << std::endl;
  }
}

// Reader for fit functions
void read(XMLReader& xml, const std::string path, NCOR::fitInfo_t& I)
{
  try {
    std::string type;
    read(xml, path+"/funcType", type); // Access type
    setType(type, I);
    // read(xml, path+"/funcType", type, I);

    read(xml, path+"/range", I.range);
    read(xml, path+"/bayes", I.bayesianFit);
    read(xml, path+"/imposeNonLinParamHierarchy", I.imposeNonLinParamHierarchy);
    read(xml, path+"/priors/prior", I.strParamValMaps.strPriorMap);
    read(xml, path+"/priors/width", I.strParamValMaps.strWidthMap);
    read(xml, path+"/minimizerProps/maxIters", I.initFitParams.maxIters);
    read(xml, path+"/minimizerProps/tolerance", I.initFitParams.tolerance);
    read(xml, path+"/minimizerProps/start", I.strParamValMaps.strStartMap);
    read(xml, path+"/minimizerProps/step", I.strParamValMaps.strStepMap);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse fitInfo_t from ini xml " << e << std::endl;
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
    read(xml, path+"/rows", d.rows);
    // read(xml, path+"/rowinfo", d.rowWeights);
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
    read(xml, path+"/rows", d.rows);
    // read(xml, path+"/rowinfo", d.rowWeights);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse db2ptInfo struct from ini xml " << e << std::endl;
    exit(1);
  }
}


// /*
//   Modify the src/snk row structure of passed key to allow a second RowCombo (RC2) to be accessed
// */
// void makeKeyRC2(Hadron::KeyHadronSUNNPartNPtCorr_t &kRC2, XMLArray::Array<int> &newRC)
// {
//   // newRC2 = db2ptInfo.rows[1]; i.e. an XML::Array<int> object
//   kRC2.npoint[1].irrep.irrep_mom.row = newRC[0]; // new snk row
//   kRC2.npoint[2].irrep.irrep_mom.row = newRC[1]; // new src row
// }


/*
  Accept a std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> of keys, and
  row/sign modifiers and create new keys to append
*/
void expandRowCombinations(std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& v,
			   XMLArray::Array<XMLArray::Array<int> >& r)
{
  // Make r.size local copies of the vector of keys
  std::vector<std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> > modV(r.size(),v);

  for ( int i = 0; i < r.size(); i++ )
    {

      // modV[i][0].npoint[1].irrep.irrep_mom.row = r[i][0];
      // modV[i][0].npoint[3].irrep.irrep_mom.row = r[i][1];

      // }
      


  // exit(10);

  // // Iterate over the *NEW* (starting at idx1) different row combinations
  // auto ri = r.begin();
  // while ( ri != r.end() )
  //   {

      for ( int j = 0; j < modV[0].size(); j++ )
	{

  // Perform a quick check on first key in &v to see if this row combo (ri) exists,
  // if so skip inputting this row combo
	  if ( modV[0][j].npoint[1].irrep.irrep_mom.row == r[i][0]
	       && modV[0][j].npoint[3].irrep.irrep_mom.row == r[i][1] )
	    {
	      continue;
	    }
	  else {
	    modV[i][j].npoint[1].irrep.irrep_mom.row = r[i][0];
	    modV[i][j].npoint[3].irrep.irrep_mom.row = r[i][1];
	    
	    // int idx = std::distance(r.begin(),ri);
	    // for ( auto i = modV[idx].begin(); i != modV[idx].end(); ++i )
	    //   {
	    //     // Modify the 1st npoint struct
	    //     i->npoint[1].irrep.irrep_mom.row = ri->first[0];
	    //     // Modify the 2nd npoint struct, unless we are iterating through a 3pt key
	    //     if ( i->npoint.size() == 2 ) { i->npoint[2].irrep.irrep_mom.row = ri->first[1]; }
	    //     else { i->npoint[3].irrep.irrep_mom.row = ri->first[1]; }
	    //   }
	    
	  }
	  
	}

      if ( i != 0 )
	// Now tack these modified keys onto the end of v
	v.insert(v.end(), modV[i].begin(), modV[i].end());
	
      // } // end check for existing row combo in first entry of &v
      
      // ri++; // iterate
    }
}
  

/*
  Utility to make 2pt db list
*/
std::vector<std::string> makeDBList(global_t& g, info2pt& I, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<std::string> s;
  
  // Get the momentum of the template key
  Array<int> thisMom = kTemplate->npoint[1].irrep.irrep_mom.mom;

  // A temporary db to set
  std::string tmp_db_plus = I.base+"/"+g.ensem+"/"+I.t0Tag+"/momXYZ."+shortMom(thisMom,".") \
    +"/EDB/"+g.ensem+"."+g.state+"_p"+shortMom(thisMom,"")+".n"+g.nvec+"."+I.t0Tag+".edb";
      
  // Push to db list
  s.push_back(tmp_db_plus);

  // Make a db name for potentially negated momentum combinations
  if ( shortMom(thisMom,"") != "000" && global.momNegate )
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
std::vector<std::string> makeDBList(global_t& g, info3pt& I, domain_t& t)
{
  // Vectors of all dbs to search from
  std::vector<std::string> s;

  // Loop over all (pf,pi) combinations
  for ( int mom = 0; mom < g.pf.size(); ++mom )
    {
      // Convenience
      const XMLArray::Array<int> _pf = g.pf[mom];
      const XMLArray::Array<int> _pi = g.pi[mom];
      
      for ( int ti = t.min; ti <= t.max; ti+=t.step )
	{
	  // A temporary db to set
	  std::string tmp_db_plus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/snk"+shortMom(_pf,".") \
	    +"_src"+shortMom(_pi,".")+"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+shortMom(_pf,"") \
	    +"_pi"+shortMom(_pi,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
	    I.tsnkTag+ti+"."+I.zTag+".edb";
	  
	  // Push to db list
	  s.push_back(tmp_db_plus);
	  // Check for any non-trivial 3pt momenta and make new db file names
	  if ( global.momNegate && ( shortMom(_pf,"") != "000" || shortMom(_pi,"") != "000" ) )
	    {
	      // A temporary db to set
	      std::string tmp_db_minus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/snk"+shortMom(_pf*-1,".") \
		+"_src"+shortMom(_pi*-1,".")+"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+shortMom(_pf*-1,"") \
		+"_pi"+shortMom(_pi*-1,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
		I.tsnkTag+ti+"."+I.zTag+".edb";
	      // Push to db list
	      s.push_back(tmp_db_minus);
	    }
	} // ti
    } // mom
  std::cout << "  Made 3pt db list of size = " << s.size() << std::endl;
  return s;
}


/*
  Make template keys
*/
std::vector<std::string> subName(const ADAT::MapObject<std::string, int> cont_ops, const Array<int> &mom)
{
  std::vector<std::string> sub_ops;
  // Iterate continuum names
  for ( auto it = cont_ops.begin(); it != cont_ops.end(); ++it )
    {
      // Register the continuum op
      auto cont_reg = Redstar::ContHadronOpEnv::TheHadronContOpRegInfoFactory::Instance().at(it->first);
      int cont_dim = cont_reg.twoJ + 1;

      // Register all subduced versions of continuum operator
      auto allSubOps = Redstar::ContHadronOpEnv::TheHadronContOpSubducedOpFactory::Instance().at(it->first);

      // Iterate over all subduced operators of continuum operator
      for ( auto sub_op = allSubOps.begin(); sub_op != allSubOps.end(); ++sub_op )
	{
	  // Register this specfic op
	  auto reg = Redstar::TheSingleHadronOpsRegInfoFactory::Instance().at(sub_op->first);
	  
	  // Init this specific op
	  auto spin_reg = Redstar::SpinObject::SubducedSpinObjRegFactory::Instance().at(reg.spin_id);
	  auto spin_op  = Redstar::SpinObject::TheSubducedSpinObjFactory::Instance().createObject(reg.spin_id,reg.spin_id);
	  
	  std::string lg = Hadron::getIrrepLG(spin_reg.irrep);
	  int        dim = Hadron::getIrrepDim(Hadron::removeHelicity(spin_reg.irrep));


	  // If little group doesn't match with what's induced by momentum, continue
	  if ( lg != Hadron::generateLittleGroup(mom) ) continue;
	  std::cout << "     MATCHED! Working w/ LG = " << lg << " and momentum = "
		    << mom << std::endl;
	  std::cout << "Sub Op name = " << sub_op->first << std::endl;

	  
	  sub_ops.push_back(sub_op->first);

	} // sub_op
    } // it
  return sub_ops;
}

std::vector<K> templateKeys(const global_t &g, int npt)
{
  // Vector of template keys
  std::vector<K> tmpKeys;

  // Loop over all (pf,pi) combinations
  for ( int mom = 0; mom < g.pf.size(); ++mom )
    {
      // Convenience
      const XMLArray::Array<int> _pf = g.pf[mom];
      const XMLArray::Array<int> _pi = g.pi[mom];

      Hadron::KeyCGCIrrepMom_t pfIrrepMom(1,_pf);
      Hadron::KeyCGCIrrepMom_t piIrrepMom(1,_pi);
      Hadron::KeyCGCIrrepMom_t qIrrepMom(1,_pf-_pi);
  
      /*
	Inserted op names
      */
      // auto subOps = subName(g.ins.cont_names,_pf-_pi); // will blow up w/o update to adat (reverting to old DA operator names)
      std::vector<std::string> insOps;
      for ( auto it = g.ins.cont_names.begin(); it != g.ins.cont_names.end(); ++it )
	insOps.push_back(it->first);
      // Manually duplicate rho entry
      // insOps.push_back("rho_rhoxDA__J1_T1mP");

  
      // Number of template keys for this (pf,pi) set by unique subduced insertions
      std::vector<K> tmp(insOps.size());

      // Fill the keys
      bool rho_visit = false;
      for ( auto k = tmp.begin(); k != tmp.end(); ++k )
	{
	  int idx = std::distance(tmp.begin(),k);
	  // Common
	  for ( int n = 1; n <= npt; ++n )
	    {
	      k->npoint.resize(3);
	      k->npoint[n].irrep = g.basic_op;
	    }
      
	  // Snk
	  k->npoint[1].t_slice              = 4;
	  k->npoint[1].irrep.creation_op    = g.snk.create;
	  k->npoint[1].irrep.smearedP       = g.snk.smear;
	  k->npoint[1].irrep.flavor         = g.snk.cgc;
	  k->npoint[1].irrep.irrep_mom      = pfIrrepMom;
	  k->npoint[1].irrep.op.ops[1] = Hadron::KeyParticleOp_t(subName(g.snk.cont_names,_pf)[0],
								 "", Hadron::canonicalOrder(_pf),
								 g.snk.disp_list);
      
	  // Ins
	  k->npoint[2].t_slice = -3;
	  k->npoint[2].irrep.creation_op    = g.ins.create;
	  k->npoint[2].irrep.smearedP       = g.ins.smear;
	  k->npoint[2].irrep.flavor         = g.ins.cgc;
	  k->npoint[2].irrep.irrep_mom      = qIrrepMom;
	  k->npoint[2].irrep.op.ops[1] = Hadron::KeyParticleOp_t(insOps[idx],"",
								 Hadron::canonicalOrder(_pf-_pi),
								 g.ins.disp_list);


	  // If this is second time visiting rho entry, change the row to 3
	  // ..this will ensure both \gamma_x & \gamma_y are accessible
	  if ( rho_visit && insOps[idx] == "rho_rhoxDA__J1_T1mP" )
	    k->npoint[2].irrep.irrep_mom.row = 3;
      
	  // Src
	  k->npoint[3].irrep.creation_op    = g.src.create;
	  k->npoint[3].irrep.smearedP       = g.src.smear;
	  k->npoint[3].irrep.flavor         = g.src.cgc;
	  k->npoint[3].irrep.irrep_mom      = piIrrepMom;
	  k->npoint[3].irrep.op.ops[1] = Hadron::KeyParticleOp_t(subName(g.src.cont_names,_pi)[0],
								 "", Hadron::canonicalOrder(_pi),
								 g.src.disp_list);
      
	  // Note if rho has been seen
	  if ( insOps[idx] == "rho_rhoxDA__J1_T1mP" )
	    rho_visit = true;
	} // auto k

      // Push this collection of tmp keys for this (pf,pi) into tmpKeys
      tmpKeys.insert(tmpKeys.end(), tmp.begin(), tmp.end());
    } // mom
  return tmpKeys;
}


/*
  Utility to make 2pt key list
*/
std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> makeKeyList(Hadron::KeyHadronSUNNPartNPtCorr_t &kTemplate)
{
  std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> s;

  s.push_back(kTemplate);
  // Default to row1-row1 combo (from template key); add the row2-row2 combo if t2ptRows = 0
  if ( global.t2ptRows == 0 )
    {
      // Assuming only a single 2pt key populates s, and that it is a row1-row1 piece
      s.push_back(kTemplate);
      s[1].npoint[1].irrep.irrep_mom.row = 2;
      s[1].npoint[2].irrep.irrep_mom.row = 2;
    }

  // Potentially include negated momentum combination
  if ( global.momNegate )
    checkKeyMomNegate(s,global);
  return s;
}


/*
  Utility to make 3pt key list
*/
std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> makeKeyList(domain_t& t, std::vector<K> *kTemplates)
{
  std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> s;

  for ( auto kTemplate : *kTemplates )
    {

      // With template 3pt key passed as pointer, iterate over tseps and make remaining 3pt keys
      for ( int ti = t.min; ti <= t.max; ti+=t.step )
	{
	  // Make a new key from the template
	  Hadron::KeyHadronSUNNPartNPtCorr_t tmp_k = kTemplate;
	  
	  // Change sink t_slice
	  tmp_k.npoint[1].t_slice=ti;
	  // Push this tmp key
	  s.push_back(tmp_k);
	  
	  if ( global.dispNegate )
	    {
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
	    } // global.dispNegate
	}
      
      // Potentially include negated momentum combination
      if ( global.momNegate )
	checkKeyMomNegate(s,global);

    } // auto kTemplate

  std::cout << "  Set all 3pt keys from template key" << std::endl;
  return s;
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


void Reads2pt(std::ifstream &inFile, NCOR::correlator& dum, info2pt& I, Hadron::KeyHadronSUNNPartNPtCorr_t *h)
{
  // // Look and check for existence of potentially several files
  // for ( auto d = I.base.begin(); d != I.base.end(); ++d )
  //   {  --> *d
      std::string dumb = I.base+"/"+I.momTag+"."+shortMom(h->npoint[1].irrep.irrep_mom.mom,".")
	+"/ENS/"+Hadron::ensemFileName(*h);


      NCOR::read(inFile,dumb,dum.ensemble);
      // arr_print(dum.ncor[0]);
}
#endif


// /*
//   Utility to ensure data has been parsed correctly from dbs
// */
// // void parseCheck(ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> A)
// void parseCheck(std::vector<corrEquivalence> T)
// {
//   for ( auto tt = T.begin(); tt != T.end(); ++tt )
//     {
//       for ( auto kc = tt->keyCorrMap.begin(); kc != tt->keyCorrMap.end(); ++kc )
// 	{
// 	  std::cout << "Performing simple ensemble average check on key = \n"
// 		    << "     " << kc->first << std::endl;
	  
// 	  corrFunc funcDum;
// 	  funcDum.data = kc->second;
// 	  std::cout << "We have dimensions: " << funcDum.data.ncor.size() << " X "
// 		    << funcDum.data.ncor[0].real.size() << " X "
// 		    << funcDum.data.ncor[0].real[0].size() << std::endl;
// 	  // arr_print(funcDum.data.ncor[0]);
// 	  ensemble funcDumAvg = corrAvg(funcDum.data);
// 	  std::cout << " Now the average" << std::endl;
// 	  arr_print(funcDumAvg);
// 	}
//     }
// }

/*
  Register all the ops
*/
bool registerAll(void)
{
  bool foo = true;

  try
    {
      // Operators                                                                                                                                                                                              
      foo &= Redstar::HadronNptOpEnv::registerAll();
    }
  catch(const std::out_of_range& e)
    {
      std::cerr << __PRETTY_FUNCTION__ << ": Caught out-of-range exception: " << e.what() << std::endl;
      exit(1);
    }
  catch(...)
    {
      std::cerr << __PRETTY_FUNCTION__ << ": caught generic exception" << std::endl;
      exit(1);
    }

  return foo;
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


  // Register all the necessary factories
  bool foo = registerAll();
  // // Hadron::IrrepsCubicHelicityEnv::registerAll();
  // Hadron::SubduceTablesOctEnv::registerAll();
  // Hadron::SubduceTablesLgEnv::registerAll();


  /*
    Determine global properties
  */
  read(xmlSR, "/PGITD/global", global);
  /*
    Determine the db info - sets where/what edbs to search for keys
  */
  read(xmlSR, "/PGITD/dbInfo/threePt/tseries/range", temporal3pt);
  read(xmlSR, "/PGITD/dbInfo/threePt", db3ptInfo);
  read(xmlSR, "/PGITD/dbInfo/twoPtFin/tseries/range", temporal2ptFin);
  read(xmlSR, "/PGITD/dbInfo/twoPtFin", db2ptFinInfo);
  read(xmlSR, "/PGITD/dbInfo/twoPtIni/tseries/range", temporal2ptIni);
  read(xmlSR, "/PGITD/dbInfo/twoPtIni", db2ptIniInfo);
  read(xmlSR, "/PGITD/dbInfo/twoPtRest/tseries/range", temporal2ptRest);
  read(xmlSR, "/PGITD/dbInfo/twoPtRest", db2ptRestInfo);
  // dumpRowInfo(db3ptInfo.rows,db3ptInfo.signs,3);
  // dumpRowInfo(db2ptInfo.rows,db2ptInfo.signs,2);


  // Set the number of tseps once and for all
  const int nTSeps = temporal3pt.numT();

  // Structure to hold all desired 3pt/2pt motion/normalizing keys
  struct keydb_t
  {
    std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> keys;
    std::vector<std::string>                        dbs;
    // Map of keys and correlators structs returned
    // std::map<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> keycorrMap; --> need a comparator
    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, NCOR::correlator> keycorrMap;
  };

  struct nptKeysValues_t
  {
    keydb_t keydbs;
  } threePt, twoPi, twoPf, twoPtRest;


  std::vector<K> tempKey3pt = templateKeys(global,3);
  for ( auto a : tempKey3pt )
    std::cout << a << std::endl;



  /*   Make the 3pt database structures to search   */
  threePt.keydbs.dbs    = makeDBList(global, db3ptInfo, temporal3pt); // , &tempKey3pt.key);
  // Make all 3pt keys from templates
  threePt.keydbs.keys   = makeKeyList(temporal3pt, &tempKey3pt);
  // Expand 3pt keys to include all desired row combinations
  expandRowCombinations(threePt.keydbs.keys, db3ptInfo.rows);
  /* Print the 3pt dbs & keys to read */
  // dumpDBs(threePt.keydbs.dbs);
  dumpKeys(threePt.keydbs.keys, 3);

  /*
    Access and store all three point functions
  */
  // std::vector<NCOR::corrEquivalence> funcs3pt;
  std::unordered_map<std::string, std::vector<NCOR::corrEquivalence> > funcs3pt;
  try {
    funcs3pt = getCorrs(threePt.keydbs.dbs,threePt.keydbs.keys);
  }
  catch (std::string e) {
    std::cout << "Failed access of 3pt functions : " << e << std::endl;
    exit(1);
  }
  std::cout << "...........3pt Correlators   ---   SUCCESS!" << std::endl;


  /*
    Do some checks of the three-pt functions
  */
#if 1
  // for ( std::unordered_map<std::string, std::vector<NCOR::corrEquivalence> >::iterator m = funcs3pt.begin(); m != funcs3pt.end(); ++m )
  for ( auto m = funcs3pt.begin(); m != funcs3pt.end(); ++m )
    {
      // for ( std::vector<NCOR::corrEquivalence>::iterator i = m->second.begin(); i != m->second.end(); ++i )
      for ( auto i = m->second.begin(); i != m->second.end(); ++i )
	{
	  std::cout << " This {op,abs(z),abs(p)} corr type has " << i->keyCorrMap.size() << " keys" << std::endl;
	  for ( auto it = i->keyCorrMap.begin(); it != i->keyCorrMap.end(); ++it )
	    {
	      std::cout << it->first << std::endl;
	      std::cout << it->second << std::endl;
	      // std::cout << it->second.getSrc().second << std::endl;
	      // std::cout << it->second.ensemble.ens << std::endl;

	      // toMerge[i-1].push_back(thisCorr.ensemble.ens);
	      // NCOR::correlator mergeCorr = NCOR::mergeCorrs(toMerge[j+i*srcOp.irrep_dim]);

	    }
	}
    }


  std::cout << "Number of tseps 3pt funcs are stored for = " << funcs3pt.size() << std::endl;
#if 0
  parseCheck(funcs3pt); std::cout << "\n";
#endif
  exit(8);
#endif



  /*
    GLOBAL PROPERTIES FOR CONSTRAINED SYSTEM TO SOLVE VIA SVD
  */
  const int GPD_RANK               = 8;                            // # AMPLITUDES EXTRACTED
  const int NUM_PFPI               = global.pf.size();             // # (PF,PI) TO SIMULTANEOUSLY CONSIDER
  const int NUM_CURRENTS           = global.ins.cont_names.size(); // # CURRENTS PER (PF,PI)
#if ROWONE
  const int NUM_ROWS_PER_INSERTION = 1;
#elif DIAGMATS
  const int NUM_ROWS_PER_INSERTION = 2;
#else
  const int NUM_ROWS_PER_INSERTION = 4;
#endif
  const int NUM_MATS = NUM_PFPI*NUM_CURRENTS*NUM_ROWS_PER_INSERTION;
  // const int NUM_MATS = 8;
  // const int MATS_PER_INSERTION = 4;
  
  // std::vector<Eigen::Matrix<std::complex<double>, NUM_MATS, 1> > MAT(global.cfgs);
  // std::vector<Eigen::MatrixXcd(NUM_MATS,1)> MAT(global.cfgs);
  std::vector<Eigen::MatrixXcd> MAT(global.cfgs,Eigen::MatrixXcd(NUM_MATS,1));
  std::vector<Eigen::Matrix<std::complex<double>, GPD_RANK, 1> > AMP(global.cfgs);
  const std::vector<int> DISP  = shortZ(global.disp_list);
  
  /*
    SPECIFY HOW COLLECTION OF (PF,PI), GAMMAS, AND ROWS ARE ASSEMBLED
  */
  // How each current is organized for a fixed (pf,pi)
  std::unordered_map<std::string, int> currentInMATOrder =
    {
#if ROWONE
      { "b_b0xDA__J0_A1pP", 0 }, { "gamma_x", 1 }, { "gamma_y", 2 },
#elif DIAGMATS
      { "b_b0xDA__J0_A1pP", 0 }, { "gamma_x", 2 }, { "gamma_y", 4 },
#else
      { "b_b0xDA__J0_A1pP", 0 }, { "gamma_x", 4 }, { "gamma_y", 8 }
      // { "b_b0xDA__J0_A1pP", 0 }, { "gamma_x", 4 }, { "gamma_y", 8 }, { "gamma_z", 12 }
#endif
    };
  // How rows are organized for a fixed (pf,pi,gamma) 
  std::map<std::pair<int,int>, int> matIDX;
#if ROWONE
  matIDX[std::make_pair(1,1)] = 0;
#elif DIAGMATS
  matIDX[std::make_pair(1,1)] = 0; matIDX[std::make_pair(2,2)] = 1;
#else
  matIDX[std::make_pair(1,1)] = 0; matIDX[std::make_pair(1,2)] = 1;
  matIDX[std::make_pair(2,1)] = 2; matIDX[std::make_pair(2,2)] = 3;
#endif
  




  /*
    PROCESS FIT INFORMATION FOR ALL 2PTS & 3PT
  */
  // Get the fit info for ini/fin 2pts & 3pt
  NCOR::fitInfo_t twoPtRestFitInfo, twoPtFinFitInfo, twoPtIniFitInfo, threePtFitInfo;
  read(xmlSR, "/PGITD/fitting/twoPtFin", twoPtFinFitInfo);
  read(xmlSR, "/PGITD/fitting/twoPtIni", twoPtIniFitInfo);
  read(xmlSR, "/PGITD/fitting/twoPtRest", twoPtRestFitInfo);
  read(xmlSR, "/PGITD/fitting/threePt", threePtFitInfo);
  // Parse the strParamValMaps
  twoPtFinFitInfo.parseParamMaps();
  twoPtIniFitInfo.parseParamMaps();
  twoPtRestFitInfo.parseParamMaps();
  threePtFitInfo.parseParamMaps();
  //*****************************************************************************************


  /*
    STANDARD 2PT REST FRAME STUFF
  */
  // Assuming same general structure for a template 2pt rest key as tempKey2Pi
  K tempKeyRest; tempKeyRest.npoint.resize(2);
  tempKeyRest.npoint[1].irrep.op.ops[1].name = "NucleonMG1g1MxD0J0S_J1o2_G1g1";
  tempKeyRest.npoint[1].irrep.irrep_mom.mom = global.rest;
  tempKeyRest.npoint[1].irrep.op.ops[1].mom_type = global.rest;
  tempKeyRest.npoint[2] = tempKeyRest.npoint[1];
  tempKeyRest.npoint[1].t_slice = -2;
  tempKeyRest.npoint[2].t_slice = 0;
  std::cout << tempKeyRest << std::endl;

  twoPtRest.keydbs.dbs  = makeDBList(global, db2ptRestInfo, &tempKeyRest);
  twoPtRest.keydbs.keys = makeKeyList(tempKeyRest);

  prop_t * props = new prop_t(global.cfgs,temporal2ptRest,tempKeyRest);
  props->npt = 2;
  NCOR::correlator rest2pt(*props);
  delete props;

  /*
    ACCESS 2PT REST FRAME - IF FIT RESULTS EXIST, READ THEM AND SKIP READING DBS & FITTING
  */
  if ( ! readCorrFitResH5(&rest2pt,"corr2pt-FitRes.h5") )
  // if ( ! rest2pt.fitExists )
    {
      for ( auto k = twoPtRest.keydbs.keys.begin(); k != twoPtRest.keydbs.keys.end(); ++k )
	Reads2pt(inFile, rest2pt, db2ptRestInfo, &(*k) );

      rest2pt.jackknife();  rest2pt.ensAvg();  rest2pt.Cov();
      // Pass correlator data covariance to set up fit
      // submatrix of covariance grabbed internally & inverse computed
      rest2pt.fit = NCOR::fitFunc_t(twoPtRestFitInfo,rest2pt.cov.dat,temporal2ptRest);

      NFIT::driver(&rest2pt, "real", true);  fitResW(&rest2pt, "real");
      writeCorr(&rest2pt);
    }
  // else
  //   {
  //     readCorrFitResH5(&rest2pt,"corr2pt-FitRes.h5");
  //   }
  std::cout << "Whoop" << std::endl;
  exit(10);
  //************************************************************************************************
  

//   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   /// MASSIVE LOOP OVER EACH (PF,PI) COMBINATION - EACH SUCCESSIVELY PROCESSED AND STACKED INTO SVD PROCEDURE ///
//   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
//   for ( int _MOM = 0; _MOM < g.pf.size(); ++_MOM )
//     {
//       // Convenience
//       const XMLArray::Array<int> _PF = g.pf[_MOM];
//       const XMLArray::Array<int> _PI = g.pi[_MOM];

//       // Template 2pt stuff
//       K tempKey2Pf, tempKey2Pi;
//       tempKey2Pf.npoint.resize(2); tempKey2Pi.npoint.resize(2);
//       for ( int i = 1; i <= 2; ++i )
// 	{
// 	  tempKey2Pf.npoint[i] = tempKey3pt[0].npoint[1];
// 	  tempKey2Pi.npoint[i] = tempKey3pt[0].npoint[3];
// 	}
//       tempKey2Pf.npoint[1].t_slice = tempKey2Pi.npoint[1].t_slice = -2;
//       tempKey2Pf.npoint[2].t_slice = tempKey2Pi.npoint[2].t_slice = 0;
      
      
      
//       std::cout << tempKey2Pf << std::endl;
//       std::cout << tempKey2Pi << std::endl;
      
      
  
//   /*
//     NOW LOAD THE 2PT FUNCTIONS
//   */
//   twoPi.keydbs.dbs      = makeDBList(global, db2ptIniInfo, &tempKey2Pi);
//   twoPi.keydbs.keys     = makeKeyList(tempKey2Pi);
//   twoPf.keydbs.dbs      = makeDBList(global, db2ptFinInfo, &tempKey2Pf);
//   twoPf.keydbs.keys     = makeKeyList(tempKey2Pf);
  


//   /*
//     Access and store all 2pt functions
//   */
//   // std::vector<NCOR::corrEquivalence> twoPtIni(1), twoPtFin(1);
//   prop_t * propsIni = new prop_t(global.cfgs,temporal2ptIni,tempKey2Pi);
//   propsIni->npt = 2;  
//   NCOR::correlator twoPtIni(*propsIni);
//   delete propsIni;

//   prop_t * propsFin = new prop_t(global.cfgs,temporal2ptFin,tempKey2Pf);
//   propsFin->npt = 2;
//   NCOR::correlator twoPtFin(*propsFin);
//   delete propsFin;




// #ifdef HAVE_DISPLIST_2PT
// #warning "Using getCorrs to read 2pt correlators"
//   twoPtIni = getCorrs(twoPi.keydbs.dbs,twoPi.keydbs.keys);
//   twoPtFin = getCorrs(twoPf.keydbs.dbs,twoPf.keydbs.keys);
// #else
// #warning ">>>>>>>>>>>  2pt correlators constructed w/o disp_list  -->  reverting to ascii reader"
//   std::ifstream inFile;
  
//   for ( auto k = twoPi.keydbs.keys.begin(); k != twoPi.keydbs.keys.end(); ++k )
//     Reads2pt(inFile, twoPtIni, db2ptIniInfo, &(*k) );
//   for ( auto k = twoPf.keydbs.keys.begin(); k != twoPf.keydbs.keys.end(); ++k )
//     Reads2pt(inFile, twoPtFin, db2ptFinInfo, &(*k) );

// #endif



//   twoPtIni.jackknife(); twoPtFin.jackknife(); 
//   twoPtIni.ensAvg();    twoPtFin.ensAvg();    
//   twoPtIni.Cov();       twoPtFin.Cov();       



  

//   // To set up fit properly, pass the correlator's data covariance
//   // Submatrix of covariance is internally grabbed, and its inverse computed
//   twoPtFin.fit = NCOR::fitFunc_t(twoPtFinFitInfo,twoPtFin.cov.dat,temporal2ptFin);
//   twoPtIni.fit = NCOR::fitFunc_t(twoPtIniFitInfo,twoPtIni.cov.dat,temporal2ptIni);
  

//   // Fire up the fits
//   std::vector<std::string> components(2);
//   components[0] = "real"; components[1] = "imag";
//   NFIT::driver(&twoPtFin, components[0], true); fitResW(&twoPtFin, components[0]);
//   NFIT::driver(&twoPtIni, components[0], true); fitResW(&twoPtIni, components[0]);
  


//   writeCorr(&twoPtFin);
//   writeCorr(&twoPtIni);


//   /*
//     Do some checks of 2pt functions
//   */
// #if 0
//   std::cout << "2pt INI KEYS" << std::endl; dumpKeys(twoPi.keydbs.keys,2);
//   std::cout << "2pt FIN KEYS" << std::endl; dumpKeys(twoPf.keydbs.keys,2);
//   parseCheck(twoPtIni); std::cout << "\n";
//   parseCheck(twoPtFin); std::cout << "\n";
//   exit(9);
// #endif


//   /*
//     Do some checks of the 2pt fit results
//   */
// #if 0
//   std::cout << "PARAMS FOR 2PT INI:" << std::endl;
//   for ( auto i = twoPtIni.res.params.begin(); i != twoPtIni.res.params.end(); ++i )
//     {
//       std::cout << "    " << i->first << " : ";
//       for ( auto p = i->second.begin(); p != i->second.end(); ++p )
// 	{
// 	  std::cout << *p << " ";
// 	}
//       std::cout << "\n";
//     }

//   std::cout << "PARAMS FOR 2PT FIN:" << std::endl;
//   for ( auto i = twoPtFin.res.params.begin(); i != twoPtFin.res.params.end(); ++i )
//     {
//       std::cout << "    " << i->first << " : ";
//       for ( auto p = i->second.begin(); p != i->second.end(); ++p )
//         {
// 	  std::cout << *p << " ";
//         }
//       std::cout << "\n";
//     }

//   std::cout << "PARAMS FOR 2PT REST:" << std::endl;
//   for ( auto i = rest2pt.res.params.begin(); i != rest2pt.res.params.end(); ++i )
//     {
//       std::cout << "    " << i->first << " : ";
//       for ( auto p = i->second.begin(); p != i->second.end(); ++p )
//         {
// 	  std::cout << *p << " ";
//         }
//       std::cout << "\n";
//     }
//   std::cout << "REST chi2: " << std::endl;
//   for ( auto c = rest2pt.res.chi2.begin(); c != rest2pt.res.chi2.end(); ++c )
//     std::cout << *c << " ";
//   std::cout << "\n";
// #endif




//   // Collect AvgP & Delta per jackknife sample to form correct lin.comb. of amps that project onto H & E
//   std::vector<std::vector<double> > collectAvgP(global.cfgs,std::vector<double>(4,0.0));
//   std::vector<std::vector<double> > collectDelta(global.cfgs,std::vector<double>(4,0.0));


  
//   std::vector<kinMatGPD_t> KIN(global.cfgs,
//   			       kinMatGPD_t(global.ins.cont_names.size()*MATS_PER_INSERTION,
// 					   GPD_RANK, current::VECTOR, DISP));

//   for ( auto k = KIN.begin(); k != KIN.end(); ++k )
//     {
//       int j = std::distance(KIN.begin(),k);
      
//       /*
// 	03/08/2023: try passing energies dictated by dispersion relation
//       */
//       double dispersionEf = sqrt( pow(rest2pt.res.params["E0"][j],2) + pow(2*PI/global.Lx,2)*(global.pf*global.pf) );
//       double dispersionEi = sqrt( pow(rest2pt.res.params["E0"][j],2) + pow(2*PI/global.Lx,2)*(global.pi*global.pi) );

//       double Ef = twoPtFin.res.params["E0"][j];
//       double Ei = twoPtIni.res.params["E0"][j];
//       // double Ef = dispersionEf;
//       // double Ei = dispersionEi;

//       // Initialize the final/initial state spinors for this jackknife sample
//       Spinor finSpin(tempKey2Pf.npoint[1].irrep.op.ops[1].name,
// 		     global.pf,Ef,rest2pt.res.params["E0"][j],global.Lx);
//       Spinor iniSpin(tempKey2Pi.npoint[1].irrep.op.ops[1].name,
// 		     global.pi,Ei,rest2pt.res.params["E0"][j],global.Lx);
//       // Build the spinors
//       finSpin.buildSpinors();
//       iniSpin.buildSpinors();




//       // Constant MASS variable for convenience
//       const double MASS = rest2pt.res.params["E0"][j];


//       // Init kinematic matrices for each of \gamma_4, \gamma_1, \gamma_2
//       kinMatGPD_t GPD_4(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,4,DISP);
//       kinMatGPD_t GPD_1(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,1,DISP);
//       kinMatGPD_t GPD_2(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,2,DISP);
//       kinMatGPD_t GPD_3(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,3,DISP);
//       // Assemble the kinematic matrices
//       GPD_4.assemble(true,MASS,&finSpin,&iniSpin);
//       GPD_1.assemble(true,MASS,&finSpin,&iniSpin);
//       GPD_2.assemble(true,MASS,&finSpin,&iniSpin);
//       GPD_3.assemble(true,MASS,&finSpin,&iniSpin);


//       // Hold avgP and Delta per config - for building actual amplitudes that project onto H/E
//       collectAvgP[j][0] = 0.5*(finSpin.getPhysMom()[0] + iniSpin.getPhysMom()[0]);
//       collectAvgP[j][1] = 0.5*(finSpin.getPhysMom()[1] + iniSpin.getPhysMom()[1]);
//       collectAvgP[j][2] = 0.5*(finSpin.getPhysMom()[2] + iniSpin.getPhysMom()[2]);
//       collectAvgP[j][3] = 0.5*(finSpin.getE() + iniSpin.getE());
//       collectDelta[j][0] = (finSpin.getPhysMom()[0] - iniSpin.getPhysMom()[0]);
//       collectDelta[j][1] = (finSpin.getPhysMom()[1] - iniSpin.getPhysMom()[1]);
//       collectDelta[j][2] = (finSpin.getPhysMom()[2] - iniSpin.getPhysMom()[2]);
//       collectDelta[j][3] = (finSpin.getE() - iniSpin.getE());

      
//       // Concatenate GPD_4,1,2 matrices into one large one for SVD
//       Eigen::MatrixXcd GPD(global.ins.cont_names.size()*MATS_PER_INSERTION, GPD_RANK);


//       // if ( global.ins.cont_names.begin()->first == "b_b0xDA__J0_A1pP" )
//       // 	for ( int i = 0; i < GPD_4.mat.rows(); ++i ) GPD.row(i) << GPD_4.mat.row(i);
//       // if ( global.ins.cont_names.begin()->first == "gamma_x" )
//       // 	for ( int i = 0; i < GPD_1.mat.rows(); ++i ) GPD.row(i) << GPD_1.mat.row(i);
//       // if ( global.ins.cont_names.begin()->first == "gamma_y" )
//       // 	for ( int i = 0; i < GPD_2.mat.rows(); ++i ) GPD.row(i) << GPD_2.mat.row(i);
//       // if ( global.ins.cont_names.begin()->first == "gamma_z" )
//       // 	for ( int i = 0; i < GPD_3.mat.rows(); ++i ) GPD.row(i) << GPD_3.mat.row(i);


//       // for ( int i = 0; i < GPD_4.mat.rows(); ++i ) GPD.row(i) << GPD_4.mat.row(i);
//       // for ( int i = MATS_PER_INSERTION; i < 2*MATS_PER_INSERTION; ++i )
//       // 	GPD.row(i) << GPD_2.mat.row(i-MATS_PER_INSERTION);

      

//       // Push GPD_4.mat, GPD_1.mat, GPD_2.mat into GPD
//       for ( int i = 0; i < GPD_4.mat.rows(); ++i ) GPD.row(i) << GPD_4.mat.row(i);
//       for ( int i = MATS_PER_INSERTION; i < 2*MATS_PER_INSERTION; ++i )
//       	GPD.row(i) << GPD_1.mat.row(i-MATS_PER_INSERTION);
//       for ( int i = 2*MATS_PER_INSERTION; i < 3*MATS_PER_INSERTION; ++i )
//       	GPD.row(i) << GPD_2.mat.row(i-2*MATS_PER_INSERTION);
//       // for ( int i = 3*MATS_PER_INSERTION; i < 4*MATS_PER_INSERTION; ++i )
//       // 	GPD.row(i) << GPD_3.mat.row(i-3*MATS_PER_INSERTION);
      


//       std::cout << "FINAL GPD = " << GPD << std::endl;
//       getSVs(&GPD);
//       getSingVecs(&GPD);

      
//       // Fill KIN matrix for this jackknife sample
//       k->mat = GPD;
//     }
  
  
//    } // MOM
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////




  
//   /*
//     Loop over inserted operators
//         -> formed summed ratios
// 	-> fit summed ratios and pack into 'MAT'
//   */
//   for ( auto OP = funcs3pt.begin(); OP != funcs3pt.end(); ++OP )
//     {
//       /*
// 	For this inserted operator, loop over all equivalent 3pt functions
// 	    --> diff. rows, but kinematically the same
//       */
//       for ( std::vector<NCOR::corrEquivalence>::iterator tsepItr = OP->second.begin();
// 	    tsepItr != OP->second.end(); ++tsepItr )
// 	{
// 	  // Convenience
// 	  int rowf = tsepItr->keyCorrMap.begin()->first.npoint[1].irrep.irrep_mom.row;
// 	  int rowi = tsepItr->keyCorrMap.begin()->first.npoint[3].irrep.irrep_mom.row;
// #if ROWONE
// 	  if ( rowf != 1 || rowi != 1 ) continue;
// #elif DIAGMATS
// 	  if ( rowf != rowi ) continue;
// #endif

      
// 	  /*
// 	    For each correlator stored in each corrEquivalence,
// 	    divide by 2pt function at same tsep
// 	    Done per jackknife ensemble average
// 	    Bias corrected after making ratio
	    
// 	    Sum up operator insertion time slice
// 	  */
// 	  std::vector<NCOR::correlator> ratio(nTSeps); // nTseps correlators of same rows/kinematics
// 	  for ( auto it = tsepItr->keyCorrMap.begin(); it != tsepItr->keyCorrMap.end(); ++it )
// 	    {
// 	      // tsepItr->keyCorrMap is not ordered by tsep yet --> form an index to store in ascending order
// 	      int idx = (it->first.npoint[1].t_slice - temporal3pt.min)/temporal3pt.step;
// 	      // Local copy of this tsep
// 	      const int TSEP = it->first.npoint[1].t_slice;
	      
	      
// 	      // Init this ratio --> REMEMBER, 3PTS HAVE TSLICES\IN[0,TSEP)
// 	      Pseudo::domain_t * d = new Pseudo::domain_t(0,1,it->first.npoint[1].t_slice-1);
// 	      prop_t * props = new prop_t(global.cfgs, *d, it->first);
// 	      props->npt     = 3;
// // #warning "Fix old chromaGamma member!"
// // 	      props->gamma   = -1; //global.chromaGamma;
	      
	      
// 	      // Now construct a new ratio
// 	      ratio[idx] = NCOR::correlator(*props, it->second.ensemble);
	      
// 	      // Jackknife this ratio so 3pt/2pt ratio can be formed
// 	      ratio[idx].jackknife();
// 	      // ratio[idx].ensAvg();
	      
	      
// 	      // Loop over insertion times for this 3pt TSEP
// 	      for ( auto tau = ratio[idx].ensemble.T.begin(); tau != ratio[idx].ensemble.T.end(); ++tau )
// 		{
// 		  // Loop over the jackknife ensemble averages
// 		  // ratio[idx] is for fixed T with tau varying
// 		  for ( int j = 0; j < global.cfgs; ++j )
// 		    {
// 		      // ratio[idx].ensemble.ens[j][*tau] =
// 		      // 	( ratio[idx].jack[j].avg[*tau] / twoPtFin.jack[j].avg[TSEP].real() )
// 		      // 	* sqrt(( twoPtIni.jack[j].avg[TSEP-*tau].real() * twoPtFin.jack[j].avg[*tau].real()
// 		      // 		 * twoPtFin.jack[j].avg[TSEP].real() )/
// 		      // 	       ( twoPtFin.jack[j].avg[TSEP-*tau].real() * twoPtIni.jack[j].avg[*tau].real()
// 		      // 		 * twoPtIni.jack[j].avg[TSEP].real() ));

// 		      ratio[idx].ensemble.ens[j][*tau] =
// 			( ratio[idx].jack[j].avg[*tau] / twoPtFin.jack[j].avg[TSEP] )
// 			* sqrt(( twoPtIni.jack[j].avg[TSEP-*tau] * twoPtFin.jack[j].avg[*tau]
// 				 * twoPtFin.jack[j].avg[TSEP] )/
// 			       ( twoPtFin.jack[j].avg[TSEP-*tau] * twoPtIni.jack[j].avg[*tau]
// 				 * twoPtIni.jack[j].avg[TSEP] ));
// 		    } // j
// 		} // tau
	      
	      
	      
// 	      // Include standard kinematic prefactors arising in forming optimized 3pt/2pt ratio
// 	      for ( int j = 0; j < global.cfgs; ++j )
// 		{
// 		  std::complex<double> commonKin(sqrt(4*twoPtFin.res.params["E0"][j]*twoPtIni.res.params["E0"][j]), 0.0);
		  
// 		  // Remove common kinematic factor & 1/\sqrt(2) from isovector current normalization
// 		  ratio[idx].ensemble.ens[j] *= (redFact*commonKin);
// 		} // j
	      
	      
// 	      // Get ensemble avg so bias removal can proceed
// 	      ratio[idx].ensAvg();
// #if 1
// 	      std::cout << "**************" << std::endl;
// 	      std::cout << "Ratio ens avg" << std::endl;
// 	      std::cout << ratio[idx] << std::endl;
// 	      std::cout << "**************" << std::endl;
// #endif
// 	      // Correct for bias in forming ratio
// 	      ratio[idx].removeBias();
	      
	      
// 	      // Summation of operator insertion
// 	      ratio[idx].summation();
	      
// 	      delete d;
// 	      delete props;
// 	    } // it
      

// 	  // Map the summed ratio data into correlator instance 'SR'
// 	  // --> so covariance/fitting members can be used
// 	  prop_t * xprops = new prop_t(global.cfgs, temporal3pt, ratio[0].key());
// 	  xprops->npt     = 3;
// // #warning "Fix old chromaGamma member!"
// // 	  xprops->gamma   = -1 ; //global.chromaGamma;


// 	  // Construct the single correlator instance 'SR'
// 	  // ---> this will be fit
// 	  NCOR::correlator SR(*xprops);
	  
// 	  for ( auto rptr = ratio.begin(); rptr != ratio.end(); ++rptr )
// 	    {
// 	      int ridx = std::distance(ratio.begin(), rptr);
// 	      for ( auto gg = rptr->ensemble.ens.begin(); gg != rptr->ensemble.ens.end(); ++gg )
// 		{
// 		  int gdx = std::distance(rptr->ensemble.ens.begin(), gg);
// 		  SR.ensemble.ens[gdx][ridx] = (*gg)[0];
// 		} // gg
// 	    } // rptr
// 	  /*
// 	    Now SR has been constructed
// 	  */

      
	  
// 	  /*
// 	    Perform linear fit to this summed ratio 'SR'
// 	    --> exposes matrix element that will be fed into SVD to extract amplitudes
// 	  */
// 	  SR.jackknife();
// 	  SR.ensAvg();
	  
// #if 1
// 	  std::cout << "With key = " << SR.key() << " ..." << std::endl;
// 	  std::cout << SR << std::endl;
// #endif
// #if 1
	  
// 	  // Make data covariance and initialize fit
// 	  SR.Cov();
// 	  SR.fit = NCOR::fitFunc_t(threePtFitInfo,SR.cov.dat,temporal3pt);
	  
// 	  writeCorr(&SR);
// 	  /*
// 	    Do the linear fits -- for both real/imag components
// 	  */
// 	  for ( auto f = components.begin(); f != components.end(); ++f )
// 	    {
// 	      NFIT::driver(&SR, *f, false);
	      
// 	      // Write out the fit results
// 	      fitResW(&SR, *f);
	      
	      
// 	      // Pipe fit results foreach jackknife sample into appropriate entry of MAT
// 	      for ( int g = 0; g < global.cfgs; ++g )
// 		{
// 		  std::pair<int,int> lookUp = std::make_pair(rowf,rowi);
// 		  if ( *f == "real" )
// 		    MAT[g](matIDX[lookUp] + currentInMATOrder[OP->first]).real(SR.res.params["b"][g]);
// 		  if ( *f == "imag" )
// 		    MAT[g](matIDX[lookUp] + currentInMATOrder[OP->first]).imag(SR.res.params["b"][g]);
// 		}
	      
// 	      // Destroy the stored fits values since they've been written
// 	      SR.res.chi2.clear(); SR.res.params.clear();
// 	    }
// #endif
// 	  delete xprops;
	  
// 	} // tsepItr

//     } // funcs3pt iterator


//       /*
// 	Now have summed ratios for inserted operators and src/snk row combinations
//       */


//   /*
//     Convenience for building two amplitudes that project onto H & E GPDs
//   */
  


//   // With MAT populated per jackknife ensemble avg
//   // Do the SVD per jackknife ensemble avg to extract amplitudes
//   extAmplitudes(&MAT,&KIN,&AMP);
//   std::cout << "What do these solutions look like?\n";
//   std::cout << "We have " << AMP.size() << " amplitudes" << std::endl;
//   // Put AMP results into a VectorXcd so writeAmplitudes can be reused - GPD_RANK + 2 to include two derived amplitudes that are correct combos we want
//   std::vector<Eigen::VectorXcd> finalAMP(global.cfgs,Eigen::VectorXcd(GPD_RANK+2));
//   for ( auto itr = AMP.begin(); itr != AMP.end(); ++itr )
//     {
//       int idx = std::distance(AMP.begin(), itr);

//       // Compute z\dot\Delta, z\dotP
//       double zDotDelta(0.0), zDotAvgP(0.0);
//       double zDotDeltaOverZDotAvgP = (-1.0*collectDelta[idx][2])/(-1.0*collectAvgP[idx][2]);
//       std::vector<int> disp(DISP); disp.push_back(0); // z^4 = 0 --> No time-like Wilson lines!
//       for ( int mu = 1; mu <=4; ++mu )
// 	{
// 	  for ( int nu = 1; nu <=4; ++nu )
// 	    {
// 	      zDotAvgP  += collectAvgP[idx][mu-1]*disp[nu-1]*metric(mu%4,nu%4);
// 	      zDotDelta += collectDelta[idx][mu-1]*disp[nu-1]*metric(mu%4,nu%4);
// 	    }
// 	}
//       // Done computing z\dot\Delta, z\dotP
//       // ------------------------------------------

//       for ( auto amp = 0; amp < itr->size(); ++amp )
// 	{
// 	  if ( amp < 8 )
// 	    {
// 	      std::cout << "A" << amp+1 << " = "  << (*itr)(amp) << "    ";
// 	      finalAMP[idx](amp) = (*itr)(amp);
// 	    }
// 	}

//       // G1 = A1 + (z\dot\Delta)/(2z\dot P)*A5
//       finalAMP[idx](8) = (*itr)(0) + 0.5*zDotDeltaOverZDotAvgP*(*itr)(4);
//       std::cout << "G1 = "  <<  finalAMP[idx](8) << "    ";

//       // G1 = A4 - (z\dot\Delta)/(2z\dot P)*A5 + (z\dot P)*A6 + (z\dot\Delta)*A7
//       finalAMP[idx](9) = (*itr)(3) - 0.5*zDotDeltaOverZDotAvgP*(*itr)(4) + zDotAvgP*(*itr)(5) + zDotDelta*(*itr)(6);
//       std::cout << "G2 = "  <<  finalAMP[idx](9) << "    ";
	    
//       std::cout << "\n";
//     }
//   std::cout << "\n";

//   std::cout << "This was the MAT[g=348] = ";
//   for ( int i = 0; i < NUM_MATS; ++i )
//     std::cout << MAT[348](i) << " ";
//   std::cout << "\n";



// #if 1
//   // Write the extracted amplitudes to h5
//   // -->clunky way to pass displacement
//   writeAmplitudes(&finalAMP,&global,&threePtFitInfo,&global.disp_list);
// #endif

  return 0;
}
