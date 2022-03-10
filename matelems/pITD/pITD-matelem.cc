/*
  ROUTINE TO READ IN 2PT/3PT CORRELATORS COMPUTED IN DISTILLATION,
  WHERE 3PT FUNCTION ARE MATELEMS OF PSEUDO-PDF FORMALISM

  SUMMATION OF OPERATOR INSERTION IS PERFORMED INTERNALLY,
  ALLOWING FOR SUMMATION METHOD TO EXTRACT MATRIX ELEMENT
*/
#include<vector>
#include<fstream>
#include<math.h>
#include<string>
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

#include "operators.h"
#include "pseudo_utils.h"
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

#include "hadron/irreps_cubic_factory.h"
#include "hadron/irreps_cubic_oct_factory.h"
#include "hadron/irreps_cubic_helicity_factory.h"
#include "hadron/subduce_tables_oct_factory.h"
#include "hadron/subduce_tables_lg_factory.h"
#include "hadron/subduce_tables.h"
#include "hadron/subduce_tables_factory.h"
#include "hadron/single_hadron_coeffs.h"

using namespace ADATXML;
using namespace Pseudo;
using namespace NCOR;

// Define a struct for global props
global_t global;

// Define some structs to help determine what edbs to read
domain_t temporal3pt, temporal2pt, temporal2ptRest;
info3pt db3ptInfo;
info2pt db2ptInfo, db2ptRestInfo;

// Constants
const std::complex<double> redFact(sqrt(2),0);

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

    opHelIrrepLG   = Hadron::getSingleHadronOpIrrep(name);
    opIrrepLG      = Hadron::removeHelicity(opHelIrrepLG);
    opIrrep        = Hadron::removeIrrepLG(opIrrepLG);
    opLG         = Hadron::getIrrepLG(opIrrepLG);
    // Don't know why getIrrepDim/getCubicRepNoParity fails when given irrep G1g1 - so do a hard code here
    if ( shortMom(m,"") == "000" )
      {
	opIrrepNoP   = Hadron::getCubicRepNoParity("G1g");
	irrep_dim    = Hadron::getIrrepDim("G1g");
      }
    else
      {
	opIrrepNoP   = Hadron::getCubicRepNoParity(opHelIrrepLG);
	irrep_dim    = Hadron::getIrrepDim(opIrrepLG);
      }

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
    std::cout << "  subduction table = " << std::endl;
    for ( int _i = 1; _i <= irrep_dim; ++_i )
      {
	std::cout << "                     ";
        for ( int _j = 1; _j <= irrep_dim; ++_j )
	  std::cout << (*H).operator()(_i,_j) << " ";
	std::cout << "\n";
      }
  }
};


/*
  Convert between lattice LG rows and helicity amplitudes
*/
std::complex<double> lgToHelicityAmps(const Hadron::KeyHadronSUNNPartNPtCorr_t *k)
{
  // Returned complex weight
  std::complex<double> weight;

  // Construct snk/src operator subduction info
  subduceInfo snkOp(k->npoint[1].irrep.op.ops[1].name,k->npoint[1].irrep.irrep_mom.mom);
  subduceInfo srcOp(k->npoint[3].irrep.op.ops[1].name,k->npoint[3].irrep.irrep_mom.mom);
  // Instantiate canonical rotation structs
  Hadron::CubicCanonicalRotation_t snkRot, srcRot;


  // Get the Euler angles, setting all to null for rest case
  if ( shortMom(k->npoint[1].irrep.irrep_mom.mom,"") != "000" )
    snkRot = Hadron::cubicCanonicalRotation(k->npoint[1].irrep.irrep_mom.mom);
  else
    {
      snkRot.alpha=0; snkRot.beta=0; snkRot.gamma=0;
    }
  if ( shortMom(k->npoint[3].irrep.irrep_mom.mom,"") != "000" )
    srcRot = Hadron::cubicCanonicalRotation(k->npoint[3].irrep.irrep_mom.mom);
  else
    {
      srcRot.alpha=0; srcRot.beta=0; srcRot.gamma=0;
    }

#if VERBOSITY>0
  // Irrep stuff
  std::cout << "SNK/SRC IRREP DIMS = " << snkOp.irrep_dim << "/" << srcOp.irrep_dim << std::endl;
  std::cout << "Snk Euler Angles:"
            << "        alpha = " << snkRot.alpha
            << "        beta  = " << snkRot.beta
            << "        gamma = " << snkRot.gamma << std::endl;
  std::cout << "Src Euler Angles:"
            << "        alpha = " << srcRot.alpha
            << "        beta  = " << srcRot.beta
            << "        gamma = " << srcRot.gamma << std::endl;
#endif

  // Inner subduction/Wigner-D matrix mults
  for ( int a = 1; a <= snkOp.irrep_dim; ++a )
    {
      int twoJZ_a = pow((-1),a-1);
      for ( int b = 1; b <= snkOp.irrep_dim; ++b )
        {
          int twoJZ_b = pow((-1),b-1);
          for ( int c = 1; c <= srcOp.irrep_dim; ++c )
            {
              int twoJZ_c = pow((-1),c-1);

#if VERBOSITY>2
	      std::cout << "NEW COMBO" << std::endl;
	      std::cout << "Snk Wigner-D: " << Hadron::Wigner_D(1,twoJZ_a,twoJZ_b,snkRot.alpha,snkRot.beta,snkRot.gamma) << std::endl;
	      std::cout << "Snk Subduction: " << (*snkOp.H).operator()(k->npoint[1].irrep.irrep_mom.row,a) << std::endl;
	      std::cout << "Src Wigner-D: " << Hadron::Wigner_D(1,twoJZ_b,twoJZ_c,srcRot.alpha,srcRot.beta,srcRot.gamma) << std::endl;
	      std::cout << "Src Subduction: " << (*srcOp.H).operator()(c,k->npoint[3].irrep.irrep_mom.row) << std::endl;
#endif

              // Build the weight
              weight += (*snkOp.H).operator()(k->npoint[1].irrep.irrep_mom.row,a)*
                Hadron::Wigner_D(1,twoJZ_a,twoJZ_b,
                                 snkRot.alpha,snkRot.beta,snkRot.gamma)*
                Hadron::Wigner_D(1,twoJZ_b,twoJZ_c,
                                 srcRot.alpha,srcRot.beta,srcRot.gamma)*
                (*srcOp.H).operator()(c,k->npoint[3].irrep.irrep_mom.row);


            } // c
        } // b
    } // a

#if VERBOSITY>0
  std::cout << "WEIGHT = " << weight << std::endl;
#endif

  return weight;                        
}

/*
  SOME CALLS TO HELP ADAT READ
*/
// Read a correlator edb
std::vector<NCOR::corrEquivalence>
getCorrs(const std::vector<std::string>& dbases, std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& fetch)
{
  typedef Hadron::KeyHadronSUNNPartNPtCorr_t       K;
  typedef ENSEM::VectorComplex                     V;

  // vector of corrEquivalences forall 3pt tseps to return
  std::vector<NCOR::corrEquivalence> ret( temporal3pt.numT() );


  // Make a vector large enough to hold just desired keys in fetch
  std::vector<ADATIO::SerialDBKey<K> > keyFetch(fetch.size());
  std::cout << "Will search for " << keyFetch.size() << " keys from databases" << std::endl;
  std::vector< std::vector<ADATIO::SerialDBData<V> > > dataFetch(keyFetch.size());

  // Now run over the keys
  for ( std::vector<ADATIO::SerialDBKey<K> >::iterator k = keyFetch.begin(); k != keyFetch.end(); k++ )
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
      // Initialize a prop_t struct for correlator initialization
#warning "Assuming foreach 3pt, trange = [0,tsep-1]"
      Pseudo::domain_t dom(0, 1, keyFetch[k].key().npoint[1].t_slice - 1);
      Pseudo::prop_t   p(global.cfgs, dom, keyFetch[k].key());
      // Make some array to associate with this key
      NCOR::correlator ensAdat(p);


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

              
	      std::complex<double> dc(r_,i_);

	      ensAdat.ensemble.ens[g][dd] = dc;
            }
        }
      // Insert this key-correlator combo
      ret[ (keyFetch[k].key().npoint[1].t_slice - temporal3pt.min)/temporal3pt.step ].
        keyCorrMap.insert(keyFetch[k].key(),ensAdat);
    }

  return ret;
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
    read(xml, path+"/insertion/gamma", g.chromaGamma);
    read(xml, path+"/nvec", g.nvec);
    read(xml, path+"/Lt", g.Lt);
    read(xml, path+"/Lx", g.Lx);
    read(xml, path+"/rest", g.rest);
    read(xml, path+"/pi", g.pi);
    read(xml, path+"/pf", g.pf);
    read(xml, path+"/momNegate", g.momNegate);
    read(xml, path+"/OpMap", g.opMomXML);
    // read(xml, path+"/t2ptNorm", g.t2ptNorm);
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

// Reader for fit functions
void read(XMLReader& xml, const std::string path, NCOR::fitInfo_t& I)
{
  try {
    read(xml, path+"/funcType", I.type);
    read(xml, path+"/range", I.range);
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
          }
        }

      if ( i != 0 )
        // Now tack these modified keys onto the end of v
        v.insert(v.end(), modV[i].begin(), modV[i].end());
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
  std::string tmp_db_plus = I.base+"/"+g.ensem+"/"+I.t0Tag+"/momXYZ."+shortMom(thisMom,".")
    +"/EDB/"+g.ensem+"."+g.state+"_p"+shortMom(thisMom,"")+".n"+g.nvec+"."+I.t0Tag+".edb";
      
  // Push to db list
  s.push_back(tmp_db_plus);
      
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
std::vector<std::string> makeDBList(global_t& g, info3pt& I, domain_t& t, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<std::string> s;
  for ( int ti = t.min; ti <= t.max; ti+=t.step )
    {
      // A temporary db to set
      std::string tmp_db_plus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/momXYZ."+shortMom(g.pf,".") \
        +"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+shortMom(g.pf,"") \
        +"_pi"+shortMom(g.pi,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
        I.tsnkTag+ti+"."+I.zTag+".edb";
          
      // Push to db list
      s.push_back(tmp_db_plus);
      // Check for any non-trivial 3pt momenta and make new db file names
      if ( global.momNegate && ( shortMom(g.pf,"") != "000" || shortMom(g.pi,"") != "000" ) )
        {
          // A temporary db to set
	  std::string tmp_db_minus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/momXYZ."+shortMom(g.pf*-1,".") \
            +"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+shortMom(g.pf*-1,"") \
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
std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> makeKeyList(domain_t& t, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> s;

  // With template 3pt key passed as pointer, iterate over tseps and make remaining 3pt keys
  for ( int ti = t.min; ti <= t.max; ti+=t.step )
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

  // Potentially include negated momentum combination
  if ( global.momNegate )
    checkKeyMomNegate(s,global);
  std::cout << "  Set all 3pt keys from template key" << std::endl;
  return s;
}


/*
  Try a new version of rowAvg
     - this will collapse keyCorrMap member of each corrEquivalence class to a single entry
     - single entry will have src/snk rows set to '0' (to denote average) and correlators averaged
     - this Adat MapObject will replace original
*/
void rowAvg(std::vector<NCOR::corrEquivalence>& v)
{
  for ( auto vi = v.begin(); vi != v.end(); ++vi )
    {
      // Convenience
      Hadron::KeyHadronSUNNPartNPtCorr_t tmp = vi->keyCorrMap.begin()->first;
      // Each corrEquivalence member must have same irrep momenta
      subduceInfo snkOp(tmp.npoint[1].irrep.op.ops[1].name,tmp.npoint[1].irrep.irrep_mom.mom);
      subduceInfo srcOp(tmp.npoint[3].irrep.op.ops[1].name,tmp.npoint[3].irrep.irrep_mom.mom);

      // Init canonical rotation structs
      Hadron::CubicCanonicalRotation_t snkRot, srcRot;
      // Set angles based on irrep momenta
      if ( shortMom(tmp.npoint[1].irrep.irrep_mom.mom,"") != "000" )
	snkRot = Hadron::cubicCanonicalRotation(tmp.npoint[1].irrep.irrep_mom.mom);
      else
	snkRot.alpha=0; snkRot.beta=0; snkRot.gamma=0;

      if ( shortMom(tmp.npoint[3].irrep.irrep_mom.mom,"") != "000" )
	srcRot = Hadron::cubicCanonicalRotation(tmp.npoint[3].irrep.irrep_mom.mom);
      else
	srcRot.alpha=0; srcRot.beta=0; srcRot.gamma=0;



      // Build the src subduction & Wigner-D matrices
      gsl_complex gc;
      std::complex<double> wig;
      gsl_matrix_complex * subduceSrc = gsl_matrix_complex_calloc(srcOp.irrep_dim,srcOp.irrep_dim);
      gsl_matrix_complex * wignerDSrc = gsl_matrix_complex_calloc(srcOp.irrep_dim,srcOp.irrep_dim);
      for ( int i = 1; i <= srcOp.irrep_dim; ++i )
	{
	  for ( int j = 1; j <= srcOp.irrep_dim; ++j )
	    {
	      gc = gsl_complex_rect((*srcOp.H).operator()(i,j).real(),
				    (*srcOp.H).operator()(i,j).imag());
	      gsl_matrix_complex_set(subduceSrc,i-1,j-1,gc);


	      wig = Hadron::Wigner_D(1,pow(-1,i+1),
				     pow(-1,j+1),srcRot.alpha,srcRot.beta,srcRot.gamma);
	      gc = gsl_complex_rect(wig.real(), wig.imag());
	      gsl_matrix_complex_set(wignerDSrc,i-1,j-1,gc);
	    }
	}
      
      // Build the snk subduction & Wigner-D matrices
      gsl_matrix_complex * subduceSnk = gsl_matrix_complex_calloc(snkOp.irrep_dim,snkOp.irrep_dim);
      gsl_matrix_complex * wignerDSnk = gsl_matrix_complex_calloc(snkOp.irrep_dim,snkOp.irrep_dim);
      for ( int i = 1; i <= snkOp.irrep_dim; ++i )
	{
	  for ( int j = 1; j <= snkOp.irrep_dim; ++j )
	    {
	      gc = gsl_complex_rect((*snkOp.H).operator()(i,j).real(),
				    (*snkOp.H).operator()(i,j).imag());
	      gsl_matrix_complex_set(subduceSnk,i-1,j-1,gc);


	      wig = Hadron::Wigner_D(1,pow(-1,i+1),
				     pow(-1,j+1),snkRot.alpha,snkRot.beta,snkRot.gamma);
	      gc = gsl_complex_rect(wig.real(), wig.imag());
	      gsl_matrix_complex_set(wignerDSnk,i-1,j-1,gc);
	    }
	}
      
      // Alloc & set tensor product of subduction & Wigner-D matrices
      gsl_matrix_complex *subduceTensProd, *wignerDTensProd;
      // res = gsl_matrix_complex_alloc
      subduceTensProd = LinAlg::tensorProd(subduceSnk,subduceSrc);
      wignerDTensProd = LinAlg::tensorProd(wignerDSnk,wignerDSrc);


      gsl_matrix_complex * subWigD = gsl_matrix_complex_calloc(subduceTensProd->size1,
							       subduceTensProd->size2);

      gsl_complex dum = gsl_complex_rect(1.0,0.0);
      gsl_complex dum2 = gsl_complex_rect(0.0,0.0);
      

      // std::cout << "Tensor product : ";
      // LinAlg::printMat(subduceTensProd);
      // std::cout << "Wigner Tensor Product : ";
      // LinAlg::printMat(wignerDTensProd);

      int stat = gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,dum,subduceTensProd,
				wignerDTensProd,dum2,subWigD);

      std::cout << "Result : ";
      LinAlg::printMat(subWigD);

      // Invert the subduction/Wigner-D products  -- maybe save this for later for now




      // Determine if +/- zsep need to be stored
      int disp2Store(1);
      std::vector<std::vector<int> > dispsContainer(1,tmp.npoint[2].irrep.op.ops[1].disp_list);
      if ( ! tmp.npoint[2].irrep.op.ops[1].disp_list.empty() )
	{
	  disp2Store = 2;
	  std::vector<int> dum = dispsContainer[0]; dum *= -1;
	  dispsContainer.push_back(dum);
	}


      // Collect the correlators to average
      std::vector<std::vector<NCOR::VVC> > toMerge(snkOp.irrep_dim * srcOp.irrep_dim);


      std::map<int, std::pair<int,int> > tensorProdMap; std::pair<int,int> tp;
      tp = std::make_pair(1,1); tensorProdMap[1] = tp;
      tp = std::make_pair(1,2); tensorProdMap[2] = tp;
      tp = std::make_pair(2,1); tensorProdMap[3] = tp;
      tp = std::make_pair(2,2); tensorProdMap[4] = tp;


      /*
	With subduction/wigner-D product determined, run through row combinations and
	reorganize as linear combination of original rows - this would make canonical states then
      */
      for ( int i = 1; i <= snkOp.irrep_dim*srcOp.irrep_dim; ++i )
	{
	  for ( int j = 1; j <= snkOp.irrep_dim*srcOp.irrep_dim; ++j )
	    {
	      for ( int d = 0; d < disp2Store; ++d )
		{
		  tmp.npoint[1].irrep.irrep_mom.row = tensorProdMap[j].first;
		  tmp.npoint[3].irrep.irrep_mom.row = tensorProdMap[j].second;
		  tmp.npoint[2].irrep.op.ops[1].disp_list = dispsContainer[d];

#if VERBOSITY>2
		  std::cout << "TMP = " << tmp << std::endl;
#endif

		  // Access this correlator
		  NCOR::correlator thisCorr = vi->keyCorrMap[tmp];
		  // Reweight
		  gsl_complex elem = gsl_matrix_complex_get(subWigD,i-1,j-1);
		  std::complex<double> weight(elem.dat[0],elem.dat[1]);

#if VERBOSITY>3
#warning "Hey! Head's up! Returned keyCorrMap from rowAvg call will be wrong!"
		  std::cout << "(" << i << "," << j << ") weight = " << weight << std::endl;

		  
		  thisCorr.jackknife(); thisCorr.ensAvg();

		  std::cout << "B4 = " << thisCorr << std::endl;

		  NCOR::correlator AFTER = thisCorr;
		  AFTER.ensemble = thisCorr.ensemble * weight;

		  AFTER.jackknife(); AFTER.ensAvg();
		  std::cout << "AFTER = " << AFTER << std::endl;
#endif

		  thisCorr.ensemble = thisCorr.ensemble * weight;

		  
		  // Get Ioffe-time
		  int ioffe(0);
		  if ( ! tmp.npoint[2].irrep.op.ops[1].disp_list.empty() )
		    ioffe = tmp.npoint[2].irrep.op.ops[1].disp_list.back();

		  // Conjugate if Ioffe-time is negative
		  if ( ioffe < 0 )
		    {
		      std::cout << "Ioffe = " << ioffe << " for key = ";
		      std::cout << tmp << std::endl;
		      NCOR::conj(&thisCorr.ensemble.ens);
		    }

		  /*
		    Only append this reweighted *.ensemble if weight != 0
		    Otherwise, merging below will average over non-zero 
		         & zero entries
		  */
		  // Pack this correlator for averaging
		  // toAvg[(j-1)*2+d].push_back(dum);
		  std::complex<double> null(0.0,0.0);
		  if ( weight != null )
		    toMerge[i-1].push_back(thisCorr.ensemble.ens);

		} // d
	    } // j
	  std::cout << "SIZE CHECK" << std::endl;
	  std::cout << "    toMerge[i-1].size() = " << toMerge[i-1].size() << std::endl;
	  for ( auto k = toMerge[i-1].begin(); k != toMerge[i-1].end(); ++k )
	    std::cout << "    toMerge[i-1][*].size() = " << k->size() << std::endl;
	} // i



      // Erase the original keyCorrMap
      vi->keyCorrMap.clear();
      // Make a new keyCorrMap
      ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, NCOR::correlator> cavg;


      // Loop once more over dimensions of src & snk subductions to assign rows
      for ( int i = 0; i < snkOp.irrep_dim; ++i )
	{
	  for ( int j = 0; j < srcOp.irrep_dim; ++j )
	    {
	      tmp.npoint[1].irrep.irrep_mom.row = (-1*(i+1));
	      tmp.npoint[3].irrep.irrep_mom.row = (-1*(j+1));
	      // Force disp key to be positive
	      if ( !tmp.npoint[2].irrep.op.ops[1].disp_list.empty()
		   && tmp.npoint[2].irrep.op.ops[1].disp_list.back() < 0 )
		tmp.npoint[2].irrep.op.ops[1].disp_list *= -1;


	      // Merge along toMerge[*] dimension
	      std::cout << "DEBUG - here we have " << toMerge[j+i*srcOp.irrep_dim].size() << " corrs to merge " << std::endl;
	      // NCOR::correlator mergeCorr = NCOR::mergeCorrs(toMerge[j+i*srcOp.irrep_dim]);


	      // Try a direct add
	      Pseudo::domain_t dumD(0,1,tmp.npoint[1].t_slice-1);
	      Pseudo::prop_t dumP(global.cfgs, dumD);
	      NCOR::correlator mergeCorr(dumP);

	      mergeCorr.ensemble.ens = toMerge[j+i*srcOp.irrep_dim][0];
	      if ( toMerge[j+i*srcOp.irrep_dim].size() > 1 )
		{
		  for ( int s = 1; s < toMerge[j+i*srcOp.irrep_dim].size(); ++s )
		    {
		      NCOR::dat_t dumDat; dumDat.ens = toMerge[j+i*srcOp.irrep_dim][s];
		      mergeCorr.ensemble += dumDat;
		      // NB(1/24/22): Is this where the doubling for z != 0 is coming from?
		    }
		}
	      // NB(1/24/22): Let's rescale
	      mergeCorr.ensemble *= (1.0/disp2Store);
		
	      cavg.insert(tmp, mergeCorr);
	    }
	}


      // Assign cavg to vi->keyCorrMap and proceed to next iteration
      vi->keyCorrMap = cavg;
      
      // Clean up
      gsl_matrix_complex_free(subduceSrc);
      gsl_matrix_complex_free(subduceSnk);
      gsl_matrix_complex_free(subduceTensProd);
    } // vi
}


#if 0
void rowAvg(std::vector<NCOR::corrEquivalence>& v)
{
  for ( auto a = v.begin(); a != v.end(); ++a )
    {
      // Access the first key in map to store info
      Hadron::KeyHadronSUNNPartNPtCorr_t tmp = a->keyCorrMap.begin()->first;
      
      // Obtain the snk/src subduction info
      subduceInfo snkOp(tmp.npoint[1].irrep.op.ops[1].name,tmp.npoint[1].irrep.irrep_mom.mom);
      subduceInfo srcOp(tmp.npoint[3].irrep.op.ops[1].name,tmp.npoint[3].irrep.irrep_mom.mom);

      
      


      // Determine if +/-zsep need to be stored
      int disp2Store(1);
      if ( ! tmp.npoint[2].irrep.op.ops[1].disp_list.empty() )
        disp2Store = 2;

      
      // Collect the correlators to average
      std::vector<NCOR::VVC> toAvg(snkOp.irrep_dim * srcOp.irrep_dim * disp2Store);
      

      // For this corrEquivance instance, access each correlator and pack into 'toAvg'
      for ( auto m = a->keyCorrMap.begin(); m != a->keyCorrMap.end(); ++m )
        {
          int idx = std::distance(a->keyCorrMap.begin(),m);

          // Get Ioffe-time
	  int ioffe(0);
          if ( ! m->first.npoint[2].irrep.op.ops[1].disp_list.empty() )
            ioffe = m->first.npoint[2].irrep.op.ops[1].disp_list.back();

          // Conjugate if Ioffe-time is negative
	  NCOR::VVC dum = m->second.ensemble.ens;
#if 1
          if ( ioffe < 0 )
            {
	      std::cout << "Ioffe = " << ioffe << " for key = ";
	      std::cout << m->first << std::endl;
	      NCOR::conj(&dum);
            }
#endif
                  
          // Pack this correlator
          toAvg[idx] = dum;
        }


      // Modify the tmp key to represent final row averaged correlator
      tmp.npoint[1].irrep.irrep_mom.row = 0;
      tmp.npoint[3].irrep.irrep_mom.row = 0;
      

      if ( !tmp.npoint[2].irrep.op.ops[1].disp_list.empty()
           && tmp.npoint[2].irrep.op.ops[1].disp_list.back() < 0 )
        tmp.npoint[2].irrep.op.ops[1].disp_list *= -1;


      // Erase the original keyCorrMap
      a->keyCorrMap.clear();
      
      // Remake the keyCorrMap
      ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, NCOR::correlator> cavg;
      NCOR::correlator mergeCorr = NCOR::mergeCorrs(toAvg);
      cavg.insert(tmp, mergeCorr);

      a->keyCorrMap = cavg;
    }
}
#endif

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
  std::string dumb = I.base+"/"+I.momTag+"."+shortMom(h->npoint[1].irrep.irrep_mom.mom,".")
    +"/ENS/"+Hadron::ensemFileName(*h);
  
  NCOR::read(inFile,dumb,dum.ensemble);
}
#endif



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
  Hadron::IrrepsCubicEnv::registerAll();
  Hadron::IrrepsCubicOctEnv::registerAll();
  Hadron::IrrepsCubicHelicityEnv::registerAll();
  Hadron::SubduceTablesOctEnv::registerAll();
  Hadron::SubduceTablesLgEnv::registerAll();


  /*
    Determine global properties
  */
  read(xmlSR, "/PITD/global", global);
  /*
    Determine the db info - sets where/what edbs to search for keys
  */
  read(xmlSR, "/PITD/dbInfo/threePt/tseries/range", temporal3pt);
  read(xmlSR, "/PITD/dbInfo/threePt", db3ptInfo);
  read(xmlSR, "/PITD/dbInfo/twoPt/tseries/range", temporal2pt);
  read(xmlSR, "/PITD/dbInfo/twoPt", db2ptInfo);
  read(xmlSR, "/PITD/dbInfo/twoPtRest/tseries/range", temporal2ptRest);
  read(xmlSR, "/PITD/dbInfo/twoPtRest", db2ptRestInfo);

#if 0
  // diracMat_t d1(1,false);
  // diracMat_t d2(2,false);
  // diracMat_t d3(3,false);
  // diracMat_t d4(4,false);
  // diracMat_t d5(5,false);
  // LinAlg::printMat(d1.gamma);
  // LinAlg::printMat(d2.gamma);
  // LinAlg::printMat(d3.gamma);
  // LinAlg::printMat(d4.gamma);
  // LinAlg::printMat(d5.gamma);

  // std::cout << "*&*&*&*&" << std::endl;

  polVec_t pol1 = polVec_t(1,true);
  polVec_t pol2 = polVec_t(2,true);
  polVec_t pol3 = polVec_t(3,true);
  polVec_t pol4 = polVec_t(4,true); // 0.5699
  std::cout << "S^1 = " << pol1.eval(global.pf,0.56989309716099856, 0.535,1,global.Lx) << " " << std::endl;
  std::cout << "S^2 = " << pol2.eval(global.pf,0.56989309716099856, 0.535,1,global.Lx) << " " << std::endl;
  std::cout << "S^3 = " << pol3.eval(global.pf,0.56989309716099856, 0.535,1,global.Lx) << " " << std::endl;
  std::cout << "S^4 = " << pol4.eval(global.pf,0.56989309716099856, 0.535,1,global.Lx) << " " << std::endl;

  std::cout << "Now the Jz=-1/2 polarization vectors" << std::endl;
  std::cout << "S^1 = " << pol1.eval(global.pf,0.56989309716099856, 0.535,2,global.Lx) << " " << std::endl;
  std::cout << "S^2 = " << pol2.eval(global.pf,0.56989309716099856, 0.535,2,global.Lx) << " " << std::endl;
  std::cout << "S^3 = " << pol3.eval(global.pf,0.56989309716099856, 0.535,2,global.Lx) << " " << std::endl;
  std::cout << "S^4 = " << pol4.eval(global.pf,0.56989309716099856, 0.535,2,global.Lx) << " " << std::endl;
  


  std::cout << "2E  = " << 2*0.6 << std::endl;
  std::cout << "E/m = " << 0.5699/0.535 << std::endl;


  exit(10);
#endif

  // Set the number of tseps once and for all
  const int nTSeps = temporal3pt.numT();

  // Structure to hold all desired 3pt/2pt motion/normalizing keys
  struct keydb_t
  {
    std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> keys;
    std::vector<std::string>                        dbs;
    // Map of keys and correlators structs returned
    // std::map<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> keyscorrs; --> need a comparator
    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, NCOR::correlator> keycorrMap;
  };

  struct nptKeysValues_t
  {
    keydb_t keydbs;
  } threePt, twoPi, twoPf, twoPtRest;


  struct keyTemplate_t
  {
    Hadron::KeyHadronSUNNPartNPtCorr_t key;
  } tempKey3pt, tempKey2Pi, tempKey2Pf, tempKeyRest;


  // Read the template keys from the ini xml
  try {
    read(xmlSR, "/PITD/Ops/threePt/keyTemplate", tempKey3pt.key);
    read(xmlSR, "/PITD/Ops/twoPt/keyTemplate", tempKey2Pi.key);
    read(xmlSR, "/PITD/Ops/twoPt/keyTemplate", tempKey2Pf.key);
    read(xmlSR, "/PITD/Ops/twoPtRest/keyTemplate", tempKeyRest.key);
  } catch ( std::string &e ) {
    std::cerr << "Unable to access key templates from ini xml " << e << std::endl;
  }
  std::cout << "READ THE TEMPLATE KEYS TO FETCH" << std::endl;


  /*
    Substitute correct momenta and src/snk operator names (based on passed global& global)
    into template 3pt & 2pt keys
  */
  setOpsMoms(&tempKey3pt.key, tempKey2Pf.key, tempKey2Pi.key, tempKeyRest.key, global);




  /*   Make the 3pt database structures to search   */
  threePt.keydbs.dbs    = makeDBList(global, db3ptInfo, temporal3pt, &tempKey3pt.key);
  // Make all 3pt keys from templates
  threePt.keydbs.keys   = makeKeyList(temporal3pt, &tempKey3pt.key);
  // Expand 3pt keys to include all desired row combinations
  expandRowCombinations(threePt.keydbs.keys, db3ptInfo.rows);
  /* Ensure 3pt keys conserve 3-momentum */
  conserveMom3PtKey(threePt.keydbs.keys);
  /* Print the 3pt dbs & keys to read */
  // dumpDBs(threePt.keydbs.dbs);
  // dumpKeys(threePt.keydbs.keys, 3);



  /*
    Access and store all three point functions
  */
  std::vector<NCOR::corrEquivalence> funcs3pt;
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
#if 0
  for ( auto i = funcs3pt.begin(); i != funcs3pt.end(); ++i )
    std::cout << " This {op,abs(z),abs(p)} corr type has " << i->keyCorrMap.size() << " keys" << std::endl;

  std::cout << "Number of tseps 3pt funcs are stored for = " << funcs3pt.size() << std::endl;
#if 0
  parseCheck(funcs3pt); std::cout << "\n";
#endif
  exit(8);
#endif



  /*
    NOW LOAD THE 2PT FUNCTIONS
  */
  twoPi.keydbs.dbs      = makeDBList(global, db2ptInfo, &tempKey2Pi.key);
  twoPi.keydbs.keys     = makeKeyList(tempKey2Pi.key);
  twoPf.keydbs.dbs      = makeDBList(global, db2ptInfo, &tempKey2Pf.key);
  twoPf.keydbs.keys     = makeKeyList(tempKey2Pf.key);
  twoPtRest.keydbs.dbs  = makeDBList(global, db2ptRestInfo, &tempKeyRest.key);
  twoPtRest.keydbs.keys = makeKeyList(tempKeyRest.key);


  /*
    Access and store all 2pt functions
  */
  // std::vector<NCOR::corrEquivalence> twoPtIni(1), twoPtFin(1);
  prop_t * props = new prop_t(global.cfgs,temporal2pt, tempKey2Pi.key);
  props->npt = 2;
  NCOR::correlator twoPtIni(*props); props->key = tempKey2Pf.key;
  NCOR::correlator twoPtFin(*props);

  props = new prop_t(global.cfgs,temporal2ptRest, tempKeyRest.key);
  props->npt = 2;
  NCOR::correlator rest2pt(*props);
  delete props;


#ifdef HAVE_DISPLIST_2PT
#warning "Using getCorrs to read 2pt correlators"
  twoPtIni = getCorrs(twoPi.keydbs.dbs,twoPi.keydbs.keys);
  twoPtFin = getCorrs(twoPf.keydbs.dbs,twoPf.keydbs.keys);
#else
#warning ">>>>>>>>>>>  2pt correlators constructed w/o disp_list  -->  reverting to ascii reader"
  std::ifstream inFile;
  
  for ( auto k = twoPi.keydbs.keys.begin(); k != twoPi.keydbs.keys.end(); ++k )
    {
      Reads2pt(inFile, twoPtIni, db2ptInfo, &(*k) ); // &twoPt.keydbs.keys[0]);
      // twoPtIni[0].keyCorrMap.insert(*k,dumCorr);
    }
  for ( auto k = twoPf.keydbs.keys.begin(); k != twoPf.keydbs.keys.end(); ++k )
    {
      Reads2pt(inFile, twoPtFin, db2ptInfo, &(*k) );
      // twoPtFin[0].keyCorrMap.insert(*k,dumCorr);
    }
  for ( auto k = twoPtRest.keydbs.keys.begin(); k != twoPtRest.keydbs.keys.end(); ++k )
    {
      Reads2pt(inFile, rest2pt, db2ptRestInfo, &(*k) );
    }
#endif


  twoPtIni.jackknife(); twoPtFin.jackknife(); rest2pt.jackknife();
  twoPtIni.ensAvg();    twoPtFin.ensAvg();    rest2pt.ensAvg();
  twoPtIni.Cov();       twoPtFin.Cov();       rest2pt.Cov();

  // LinAlg::printMat(twoPtFin.cov.dat["real"]);
  // LinAlg::printMat(twoPtFin.cov.inv["real"]);


  // Get the fit info for ini/fin 2pts & 3pt
  NCOR::fitInfo_t twoPtFitInfo, threePtFitInfo;
  read(xmlSR, "/PITD/fitting/twoPt", twoPtFitInfo);
  read(xmlSR, "/PITD/fitting/threePt", threePtFitInfo);

  // To set up fit properly, pass the correlator's data covariance
  // Submatrix of covariance is internally grabbed, and its inverse computed
  twoPtFin.fit = NCOR::fitFunc_t(twoPtFitInfo,twoPtFin.cov.dat,temporal2pt);
  twoPtIni.fit = NCOR::fitFunc_t(twoPtFitInfo,twoPtIni.cov.dat,temporal2pt);

  // Get the fit info for rest 2pt
  read(xmlSR, "/PITD/fitting/twoPtRest", twoPtFitInfo);
  rest2pt.fit = NCOR::fitFunc_t(twoPtFitInfo,rest2pt.cov.dat,temporal2pt);



  // Fire up the fits
  std::vector<std::string> components(2);
  components[0] = "real"; components[1] = "imag";
  NFIT::driver(&twoPtFin, components[0], true); fitResW(&twoPtFin, components[0]);
  NFIT::driver(&twoPtIni, components[0], true); fitResW(&twoPtIni, components[0]);
  NFIT::driver(&rest2pt, components[0], true);  fitResW(&rest2pt, components[0]);

  writeCorr(&twoPtFin);
  writeCorr(&twoPtIni);
  writeCorr(&rest2pt);


  /*
    Do some checks of 2pt functions
  */
#if 0
  std::cout << "2pt INI KEYS" << std::endl; dumpKeys(twoPi.keydbs.keys,2);
  std::cout << "2pt FIN KEYS" << std::endl; dumpKeys(twoPf.keydbs.keys,2);
  parseCheck(twoPtIni); std::cout << "\n";
  parseCheck(twoPtFin); std::cout << "\n";
  exit(9);
#endif


  /*
    Do some checks of the 2pt fit results
  */
#if 0
  std::cout << "PARAMS FOR 2PT INI:" << std::endl;
  for ( auto i = twoPtIni.res.params.begin(); i != twoPtIni.res.params.end(); ++i )
    {
      std::cout << "    " << i->first << " : ";
      for ( auto p = i->second.begin(); p != i->second.end(); ++p )
        {
	  std::cout << *p << " ";
        }
      std::cout << "\n";
    }

  std::cout << "PARAMS FOR 2PT FIN:" << std::endl;
  for ( auto i = twoPtFin.res.params.begin(); i != twoPtFin.res.params.end(); ++i )
    {
      std::cout << "    " << i->first << " : ";
      for ( auto p = i->second.begin(); p != i->second.end(); ++p )
        {
	  std::cout << *p << " ";
        }
      std::cout << "\n";
    }

  std::cout << "PARAMS FOR 2PT REST:" << std::endl;
  for ( auto i = rest2pt.res.params.begin(); i != rest2pt.res.params.end(); ++i )
    {
      std::cout << "    " << i->first << " : ";
      for ( auto p = i->second.begin(); p != i->second.end(); ++p )
        {
	  std::cout << *p << " ";
        }
      std::cout << "\n";
    }
  std::cout << "REST chi2: " << std::endl;
  for ( auto c = rest2pt.res.chi2.begin(); c != rest2pt.res.chi2.end(); ++c )
    std::cout << *c << " ";
  std::cout << "\n";
  exit(90);
#endif




  // for ( auto a = funcs3pt.begin(); a != funcs3pt.end(); ++a )
  //   {
  //     for ( auto aa = a->keyCorrMap.begin(); aa != a->keyCorrMap.end(); ++aa )
  // 	{
  // 	  std::cout << "NEW ONE" << std::endl;
  // 	  NCOR::correlator foo = aa->second;
  // 	  foo.jackknife();
  // 	  foo.ensAvg();
  // 	  std::cout << aa->first << std::endl;
  // 	  std::cout << foo << "\n" << std::endl;
  // 	}
  //   }

  // exit(8);

  
  /*
    3pt functions at this stage are reweighted
    Now average all 4 row combinations for each p/z combination
  */
  rowAvg(funcs3pt);

  // for ( auto a = funcs3pt.begin(); a != funcs3pt.end(); ++a )
  //   {
  //     for ( auto aa = a->keyCorrMap.begin(); aa != a->keyCorrMap.end(); ++aa )
  // 	{
  // 	  std::cout << "NEW ONE" << std::endl;
  // 	  NCOR::correlator foo = aa->second;
  // 	  foo.jackknife();
  // 	  foo.ensAvg();
  // 	  std::cout << aa->first << std::endl;
  // 	  std::cout << foo << "\n" << std::endl;
  // 	}
  //   }
  // exit(8);

  for ( auto a = funcs3pt.begin(); a != funcs3pt.end(); ++a )
    {
      // Copy reference of first key in each entry of funcs3pt
      Hadron::KeyHadronSUNNPartNPtCorr_t align = a->keyCorrMap.begin()->first;
      Hadron::KeyHadronSUNNPartNPtCorr_t antialign = align;
      align.npoint[1].irrep.irrep_mom.row = -1;
      align.npoint[3].irrep.irrep_mom.row = -1;
      antialign.npoint[1].irrep.irrep_mom.row = -2;
      antialign.npoint[3].irrep.irrep_mom.row = -2;


      std::complex<double> spinAvg(0.5,0.0);
      // Now form the unpolarized or polarized combination based on global projection
      switch ( global.chromaGamma )
	{
	case 8:
	  a->keyCorrMap[align].ensemble += a->keyCorrMap[antialign].ensemble;
	  a->keyCorrMap[align].ensemble *= (spinAvg * redFact);
	  break;
	case 11:
#warning "SUPER HACKY WAY TO MANAGE HELICITY! - FIX ME!!!!!!!!!"
	  if ( global.pf[2] >= 0 )
	    {
	      a->keyCorrMap[align].ensemble -= a->keyCorrMap[antialign].ensemble;
	      a->keyCorrMap[align].ensemble *= (spinAvg * redFact);
	    }
	  else if ( global.pf[2] < 0 )
	    {
	      a->keyCorrMap[antialign].ensemble -= a->keyCorrMap[align].ensemble;
	      a->keyCorrMap[align].ensemble = a->keyCorrMap[antialign].ensemble;
	      a->keyCorrMap[align].ensemble *= (spinAvg * redFact);
	    }
	  break;
	}

      a->keyCorrMap[align].ensAvg(); a->keyCorrMap[antialign].ensAvg();
      std::cout << "Align after spin-avg!" << a->keyCorrMap[align] << std::endl;


      // Unpol/pol combo now sits in the <-1,-1> row entry
      // Dump the remaining from the keyCorrMap
      a->keyCorrMap.erase(antialign);
      antialign.npoint[1].irrep.irrep_mom.row = -1;   a->keyCorrMap.erase(antialign);
      antialign.npoint[1].irrep.irrep_mom.row = -2;
      antialign.npoint[3].irrep.irrep_mom.row = -1;   a->keyCorrMap.erase(antialign);
    }

  
  // for ( auto a = funcs3pt.begin(); a != funcs3pt.end(); ++a )
  //   {
  //     for ( auto aa = a->keyCorrMap.begin(); aa != a->keyCorrMap.end(); ++aa )
  // 	{
  // 	  std::cout << "NEW ONE" << std::endl;
  // 	  NCOR::correlator foo = aa->second;
  // 	  foo.jackknife();
  // 	  foo.ensAvg();
  // 	  std::cout << aa->first << std::endl;
  // 	  std::cout << foo << "\n" << std::endl;
  // 	}
  //   }
  // exit(8);




  /*
    For each 3pt function, divide by 2pt function at same tsep
    Done per jackknife ensemble average
    Bias corrected after making ratio
  */
  std::vector<NCOR::correlator> ratio(funcs3pt.size());
  for ( auto it = funcs3pt.begin(); it != funcs3pt.end(); ++it )
    {
      int idx = std::distance(funcs3pt.begin(), it);
      // Init this ratio
      Pseudo::domain_t * d = new
        Pseudo::domain_t(0,1,it->keyCorrMap.begin()->second.ensemble.T.back());
      prop_t * props = new prop_t(global.cfgs, *d); //, it->keyCorrMap.begin()->first);
      props->npt     = 3;
      props->gamma   = global.chromaGamma;
      props->key     = it->keyCorrMap.begin()->first;

      // Lazy here - make sure the displacement being passed is positive
      if ( props->key.npoint[2].irrep.op.ops[1].disp_list.size() != 0
           && props->key.npoint[2].irrep.op.ops[1].disp_list[0] < 0 )
        props->key.npoint[2].irrep.op.ops[1].disp_list *= -1;
      
      // Now construct a new ratio
      ratio[idx] = NCOR::correlator(*props);
      delete d; delete props;


      // Should be only one correlator now...so the follow should be okay
      NCOR::correlator this3pt = it->keyCorrMap.begin()->second;


      this3pt.jackknife();
      this3pt.ensAvg();


      // Local copy of this tsep
      const int T = this3pt.ensemble.T.size();

      // Loop over time for this 3pt
      for ( auto tau = this3pt.ensemble.T.begin(); tau != this3pt.ensemble.T.end(); ++tau )
        {
          // Loop over the jackknife ensemble averages
          // this3pt is for fixed T with tau varying
          
          for ( int j = 0; j < global.cfgs; ++j )
            {
              ratio[idx].ensemble.ens[j][*tau] =
	      	( this3pt.jack[j].avg[*tau] / twoPtFin.jack[j].avg[T].real() ) 
                * sqrt(( twoPtIni.jack[j].avg[T-*tau].real() * twoPtFin.jack[j].avg[*tau].real()
                         * twoPtFin.jack[j].avg[T].real() )/
                       ( twoPtFin.jack[j].avg[T-*tau].real() * twoPtIni.jack[j].avg[*tau].real()
                         * twoPtIni.jack[j].avg[T].real() ));
            }
        } //tau

      // Get ensemble avg so bias removal can proceed
      ratio[idx].ensAvg();
      // std::cout << "Remove bias" << std::endl;
      // Correct for bias in forming ratio - i.e. "anti-jackknife"
      ratio[idx].removeBias();
      // std::cout << "Summation" << std::endl;
      ratio[idx].summation();


      // std::cout << "Raw ratio" << std::endl;
      // ratio[idx].jackknife(); ratio[idx].ensAvg();
      // std::cout << ratio[idx] << std::endl;

    } // it
  // exit(8);


  // Per jackknife sample, fold in standard kinematic factors
  for ( auto it = ratio.begin(); it != ratio.end(); ++it )
    {
      for ( int j = 0; j < global.cfgs; ++j )
        {
	  std::complex<double> kin(0.0,0.0);

          kin.real( ( 1 / ( 4*twoPtIni.res.params["E0"][j] *
                            (twoPtFin.res.params["E0"][j] + rest2pt.res.params["E0"][j]) )) *
                    sqrt( (twoPtIni.res.params["E0"][j]*
                           (twoPtFin.res.params["E0"][j]+rest2pt.res.params["E0"][j]))
                          / (twoPtFin.res.params["E0"][j]*(twoPtIni.res.params["E0"][j]+
                                                           rest2pt.res.params["E0"][j])) ) );

          it->ensemble.ens[j] *= (1.0/kin);
        }

      // it->jackknife();
      // it->ensAvg();
      // std::cout << *it << "\n" << std::endl;
    }
  // exit(90);

  // ratio now rescaled and basic kinefactor applied
  /*
    Remaining factor: kinematic factors from 3pt function trace
    Initialize the trace
  */
  projections::projSelect pr;
  if ( global.projector == 1 )
    pr = projections::projSelect::UNPOL;
  if ( global.projector == 2 )
    pr = projections::projSelect::POL;

  projector_t proj(pr,global.chromaGamma);
  std::cout << "Projector selected" << std::endl;

  for ( auto it = ratio.begin(); it != ratio.end(); ++it )
    {
      for ( int j = 0; j < global.cfgs; ++j )
        {
#if 1
	  std::complex<double> trace(0.0,0.0);
          trace = proj.eval(global.pf,global.pi,twoPtFin.res.params["E0"][j],
                            twoPtIni.res.params["E0"][j],
                            rest2pt.res.params["E0"][j],global.Lx);

#ifdef AMPPREFACTORS
	  std::cout << "Including amplitude prefactors in matrix element extractions!" << std::endl;
	  // Toss in prefactors of pseudo-ITDs
	  if ( global.chromaGamma == 8 )
	    trace *= (2*twoPtFin.res.params["E0"][j]); // 2*p^\mu prefactor
	  if ( global.chromaGamma == 11 )
	    {
	      polVec_t S(3,true);
	      // 2mS^\mu prefactor
	      trace *= (2*rest2pt.res.params["E0"][j]*S.eval(global.pf,twoPtFin.res.params["E0"][j],rest2pt.res.params["E0"][j],1,global.Lx));
	    }
#endif

#else
	  std::complex<double> trace(1.0,0.0);
#endif

          it->ensemble.ens[j] *= (1.0/trace);
        }

      // it->jackknife();
      // it->ensAvg();
      // std::cout << *it << std::endl;
      // std::cout << "*********************" << std::endl;
    }
  // exit(8);


#ifndef AMPPREFACTORS
  std::cout << "Not applying kinematic prefactors of pseudo-ITDs" << std::endl;
  if ( global.chromaGamma == 8 )
    {
      std::vector<std::complex<double> > pmu(global.cfgs,std::complex<double> (0.0,0.0));
      for ( auto a = pmu.begin(); a != pmu.end(); ++a )
	a->real(2*twoPtFin.res.params["E0"][std::distance(pmu.begin(),a)]);
      std::cout << "Writing 2p_mu prefactor" << std::endl;
      // Write 2p_mu prefactor for each jackknife bin to h5
      writePrefactor(3,4,global.cfgs,global.pf,"/2p_mu",pmu); // (3,3...) => npt=3 & mu=4
      std::cout << "Wrote 2p_mu prefactor" << std::endl;
    }
  if ( global.chromaGamma == 11 )
    {
      // Went through the trouble of allowing S^\mu evaluation,
      // so if gamma = 11, then compute S^\mu and store it in corr3pt-FitRes.h5
      std::vector<std::complex<double> > polVec(global.cfgs,std::complex<double> (0.0,0.0));
      for ( int j = 0; j < global.cfgs; ++j )
	{
	  polVec_t S(3,true);
	  polVec[j] = S.eval(global.pf,twoPtFin.res.params["E0"][j],
			     rest2pt.res.params["E0"][j],1,global.Lx);
	}
      
      std::cout << "Writing polarization vector" << std::endl;
      // Write the polarization vector for each jackknife bin to h5
      writePrefactor(3,3,global.cfgs,global.pf,"/polVec",polVec); // (3,3...) => npt=3 & mu=3
      std::cout << "Wrote polVec" << std::endl;
    }
#endif

  

  /*
    Concerning std::vector<NCOR::correlator> ratio:
      -- ensemble.ens formed & kinematic factors applied
      -- jackknifed & ens/jk ens avgs formed
  */
  // Now optionally sum up the operator insertion time slice for summation method
  prop_t * xprops = new prop_t(global.cfgs, temporal3pt);
  xprops->npt     = 3;
  xprops->gamma   = global.chromaGamma;
  xprops->key     = ratio[0].key();
  // xprops->mom.fin = global.pf;
  // xprops->mom.ini = global.pi;
  // xprops->disp    = ratio[0].getDisp();
  NCOR::correlator SR(*xprops);
  
  for ( auto it = ratio.begin(); it != ratio.end(); ++it )
    {
      int idx = std::distance(ratio.begin(), it);

      // Map the summed ratio data into the SR instance so covariance/fitting members can be used
      for ( auto gg = it->ensemble.ens.begin(); gg != it->ensemble.ens.end(); ++gg )
        {
          int gdx = std::distance(it->ensemble.ens.begin(), gg);
          SR.ensemble.ens[gdx][idx] = (*gg)[0];
        }
    }
  std::cout << "Summations done" << std::endl;

  SR.jackknife();
  SR.ensAvg();
#if 0
  std::cout << SR << std::endl;
  exit(90);
#endif
  SR.Cov();

  LinAlg::printMat(SR.cov.dat["real"]);
  LinAlg::printMat(SR.cov.dat["imag"]);

  // Init a linear fit for SR
  SR.fit = NCOR::fitFunc_t(threePtFitInfo,SR.cov.dat, temporal3pt);


  /*
    Do the linear fits
  */
  for ( auto f = components.begin(); f != components.end(); ++f )
    {
      NFIT::driver(&SR, *f, false);
      // Write out the fit results
      fitResW(&SR, *f);
      // Destroy the stored fits values since they've been written
      SR.res.chi2.clear(); SR.res.params.clear();
    }
  writeCorr(&SR);
  delete xprops;



  return 0;
}
