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
#include "rotations.h"
#include "shortcuts_gsl.h"

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
domain_t temporal3pt, temporal2ptFin, temporal2ptIni, temporal2ptRest;
info3pt db3ptInfo;
info2pt db2ptFinInfo, db2ptIniInfo, db2ptRestInfo;

// Constants
const std::complex<double> redFact(sqrt(2),0);


/*
  SOME CALLS TO HELP ADAT READ
*/
// Read a correlator edb
std::vector<NCOR::corrEquivalence>
getCorrs(const std::vector<std::string>& dbases, std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& fetch)
{
  typedef Hadron::KeyHadronSUNNPartNPtCorr_t       K;
  typedef ENSEM::VectorComplex                     V;

  // // vector of corrEquivalences forall 3pt tseps to return
  // std::vector<NCOR::corrEquivalence> ret( temporal3pt.numT() );

  // NEW
  std::vector<NCOR::corrEquivalence> RET(4);


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
//       // Initialize a prop_t struct for correlator initialization
// #warning "Assuming foreach 3pt, trange = [0,tsep-1]"
//       Pseudo::domain_t dom(0, 1, keyFetch[k].key().npoint[1].t_slice - 1);
//       Pseudo::prop_t   p(global.cfgs, dom, keyFetch[k].key());
//       // Make some array to associate with this key
//       NCOR::correlator ensAdat(p);


      NCOR::correlator ensAdat;
      ensAdat.ensemble.ens.resize(enscfgs);


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

	      // ensAdat.ensemble.ens[g][dd] = dc;
	      ensAdat.ensemble.ens[g].push_back(dc);
            }
        }
      // // Insert this key-correlator combo
      // ret[ (keyFetch[k].key().npoint[1].t_slice - temporal3pt.min)/temporal3pt.step ].
      //   keyCorrMap.insert(keyFetch[k].key(),ensAdat);

      // NEW - Something hacky
      Hadron::KeyHadronSUNNPartNPtCorr_t tmp = keyFetch[k].key(); //convenience
      if ( tmp.npoint[1].irrep.irrep_mom.row == 1 && tmp.npoint[3].irrep.irrep_mom.row == 1 )
	RET[0].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      if ( tmp.npoint[1].irrep.irrep_mom.row == 1 && tmp.npoint[3].irrep.irrep_mom.row == 2 )
	RET[1].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      if ( tmp.npoint[1].irrep.irrep_mom.row == 2 && tmp.npoint[3].irrep.irrep_mom.row == 1 )
	RET[2].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
      if ( tmp.npoint[1].irrep.irrep_mom.row == 2 && tmp.npoint[3].irrep.irrep_mom.row == 2 )
	RET[3].keyCorrMap.insert(keyFetch[k].key(),ensAdat);
    }

  return RET;
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
    read(xml, path+"/bayes", I.bayesianFit);
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
#if 1
      std::string tmp_db_plus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/snk"+shortMom(g.pf,".")
        +"_src"+shortMom(g.pi,".")+"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+
	shortMom(g.pf,"")+"_pi"+shortMom(g.pi,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
        I.tsnkTag+ti+"."+I.zTag+".edb";
#else
#warning "Lazy way to get q^2 FF data"
      std::string tmp_db_plus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/momXYZ"+shortMom(g.pf,".")
	+"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"+
	shortMom(g.pf,"")+"_pi"+shortMom(g.pi,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
	I.tsnkTag+ti+"."+I.zTag+".edb";
#endif


          
      // Push to db list
      s.push_back(tmp_db_plus);
      // Check for any non-trivial 3pt momenta and make new db file names
      if ( global.momNegate && ( shortMom(g.pf,"") != "000" || shortMom(g.pi,"") != "000" ) )
        {
          // A temporary db to set
	  std::string tmp_db_minus = I.base[0]+"/"+g.ensem+I.base[1]+ti+"/snk"+shortMom(g.pf*-1,".")
	    +"_src"+shortMom(g.pi*-1,".")+"/EDB/"+g.ensem+"."+g.state+"."+g.observable+"_pf"
	    +shortMom(g.pf*-1,"")+"_pi"+shortMom(g.pi*-1,"")+".n"+g.nvec+"."+I.t0Tag+"_"+
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
    }

  // Potentially include negated momentum combination
  if ( global.momNegate )
    checkKeyMomNegate(s,global);
  std::cout << "  Set all 3pt keys from template key" << std::endl;
  return s;
}


// /*
//   Try a new version of rowAvg
//      - this will collapse keyCorrMap member of each corrEquivalence class to a single entry
//      - single entry will have src/snk rows set to '0' (to denote average) and correlators averaged
//      - this Adat MapObject will replace original
// */
// void rowAvg(std::vector<NCOR::corrEquivalence>& v)
// {
//   for ( auto vi = v.begin(); vi != v.end(); ++vi )
//     {
//       // Convenience
//       Hadron::KeyHadronSUNNPartNPtCorr_t tmp = vi->keyCorrMap.begin()->first;
//       // Each corrEquivalence member must have same irrep momenta
//       subduceInfo snkOp(tmp.npoint[1].irrep.op.ops[1].name,tmp.npoint[1].irrep.irrep_mom.mom);
//       subduceInfo srcOp(tmp.npoint[3].irrep.op.ops[1].name,tmp.npoint[3].irrep.irrep_mom.mom);

//       // Init canonical rotation structs
//       Hadron::CubicCanonicalRotation_t snkRot, srcRot;
//       // Set angles based on irrep momenta
//       if ( shortMom(tmp.npoint[1].irrep.irrep_mom.mom,"") != "000" )
// 	snkRot = Hadron::cubicCanonicalRotation(tmp.npoint[1].irrep.irrep_mom.mom);
//       else
// 	snkRot.alpha=0; snkRot.beta=0; snkRot.gamma=0;

//       if ( shortMom(tmp.npoint[3].irrep.irrep_mom.mom,"") != "000" )
// 	srcRot = Hadron::cubicCanonicalRotation(tmp.npoint[3].irrep.irrep_mom.mom);
//       else
// 	srcRot.alpha=0; srcRot.beta=0; srcRot.gamma=0;



//       // Build the src subduction & Wigner-D matrices
//       gsl_complex gc;
//       std::complex<double> wig;
//       gsl_matrix_complex * subduceSrc = gsl_matrix_complex_calloc(srcOp.irrep_dim,srcOp.irrep_dim);
//       gsl_matrix_complex * wignerDSrc = gsl_matrix_complex_calloc(srcOp.irrep_dim,srcOp.irrep_dim);
//       for ( int i = 1; i <= srcOp.irrep_dim; ++i )
// 	{
// 	  for ( int j = 1; j <= srcOp.irrep_dim; ++j )
// 	    {
// 	      gc = gsl_complex_rect((*srcOp.H).operator()(i,j).real(),
// 				    (*srcOp.H).operator()(i,j).imag());
// 	      gsl_matrix_complex_set(subduceSrc,i-1,j-1,gc);


// 	      wig = Hadron::Wigner_D(1,pow(-1,i+1),
// 				     pow(-1,j+1),srcRot.alpha,srcRot.beta,srcRot.gamma);
// 	      gc = gsl_complex_rect(wig.real(), wig.imag());
// 	      gsl_matrix_complex_set(wignerDSrc,i-1,j-1,gc);
// 	    }
// 	}
      
//       // Build the snk subduction & Wigner-D matrices
//       gsl_matrix_complex * subduceSnk = gsl_matrix_complex_calloc(snkOp.irrep_dim,snkOp.irrep_dim);
//       gsl_matrix_complex * wignerDSnk = gsl_matrix_complex_calloc(snkOp.irrep_dim,snkOp.irrep_dim);
//       for ( int i = 1; i <= snkOp.irrep_dim; ++i )
// 	{
// 	  for ( int j = 1; j <= snkOp.irrep_dim; ++j )
// 	    {
// 	      gc = gsl_complex_rect((*snkOp.H).operator()(i,j).real(),
// 				    (*snkOp.H).operator()(i,j).imag());
// 	      gsl_matrix_complex_set(subduceSnk,i-1,j-1,gc);


// 	      wig = Hadron::Wigner_D(1,pow(-1,i+1),
// 				     pow(-1,j+1),snkRot.alpha,snkRot.beta,snkRot.gamma);
// 	      gc = gsl_complex_rect(wig.real(), wig.imag());
// 	      gsl_matrix_complex_set(wignerDSnk,i-1,j-1,gc);
// 	    }
// 	}
      
//       // Alloc & set tensor product of subduction & Wigner-D matrices
//       gsl_matrix_complex *subduceTensProd, *wignerDTensProd;
//       // res = gsl_matrix_complex_alloc
//       subduceTensProd = LinAlg::tensorProd(subduceSnk,subduceSrc);
//       wignerDTensProd = LinAlg::tensorProd(wignerDSnk,wignerDSrc);


//       gsl_matrix_complex * subWigD = gsl_matrix_complex_calloc(subduceTensProd->size1,
// 							       subduceTensProd->size2);

//       gsl_complex dum = gsl_complex_rect(1.0,0.0);
//       gsl_complex dum2 = gsl_complex_rect(0.0,0.0);
      

//       // std::cout << "Tensor product : ";
//       // LinAlg::printMat(subduceTensProd);
//       // std::cout << "Wigner Tensor Product : ";
//       // LinAlg::printMat(wignerDTensProd);

//       int stat = gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,dum,subduceTensProd,
// 				wignerDTensProd,dum2,subWigD);

//       std::cout << "Result : ";
//       LinAlg::printMat(subWigD);

//       // Invert the subduction/Wigner-D products  -- maybe save this for later for now




//       // Determine if +/- zsep need to be stored
//       int disp2Store(1);
//       std::vector<std::vector<int> > dispsContainer(1,tmp.npoint[2].irrep.op.ops[1].disp_list);
//       if ( ! tmp.npoint[2].irrep.op.ops[1].disp_list.empty() )
// 	{
// 	  disp2Store = 2;
// 	  std::vector<int> dum = dispsContainer[0]; dum *= -1;
// 	  dispsContainer.push_back(dum);
// 	}


//       // Collect the correlators to average
//       std::vector<std::vector<NCOR::VVC> > toMerge(snkOp.irrep_dim * srcOp.irrep_dim);


//       std::map<int, std::pair<int,int> > tensorProdMap; std::pair<int,int> tp;
//       tp = std::make_pair(1,1); tensorProdMap[1] = tp;
//       tp = std::make_pair(1,2); tensorProdMap[2] = tp;
//       tp = std::make_pair(2,1); tensorProdMap[3] = tp;
//       tp = std::make_pair(2,2); tensorProdMap[4] = tp;


//       /*
// 	With subduction/wigner-D product determined, run through row combinations and
// 	reorganize as linear combination of original rows - this would make canonical states then
//       */
//       for ( int i = 1; i <= snkOp.irrep_dim*srcOp.irrep_dim; ++i )
// 	{
// 	  for ( int j = 1; j <= snkOp.irrep_dim*srcOp.irrep_dim; ++j )
// 	    {
// 	      for ( int d = 0; d < disp2Store; ++d )
// 		{
// 		  tmp.npoint[1].irrep.irrep_mom.row = tensorProdMap[j].first;
// 		  tmp.npoint[3].irrep.irrep_mom.row = tensorProdMap[j].second;
// 		  tmp.npoint[2].irrep.op.ops[1].disp_list = dispsContainer[d];

// #if VERBOSITY>2
// 		  std::cout << "TMP = " << tmp << std::endl;
// #endif

// 		  // Access this correlator
// 		  NCOR::correlator thisCorr = vi->keyCorrMap[tmp];
// 		  // Reweight
// 		  gsl_complex elem = gsl_matrix_complex_get(subWigD,i-1,j-1);
// 		  std::complex<double> weight(elem.dat[0],elem.dat[1]);

// #if VERBOSITY>3
// #warning "Hey! Head's up! Returned keyCorrMap from rowAvg call will be wrong!"
// 		  std::cout << "(" << i << "," << j << ") weight = " << weight << std::endl;

		  
// 		  thisCorr.jackknife(); thisCorr.ensAvg();

// 		  std::cout << "B4 = " << thisCorr << std::endl;

// 		  NCOR::correlator AFTER = thisCorr;
// 		  AFTER.ensemble = thisCorr.ensemble * weight;

// 		  AFTER.jackknife(); AFTER.ensAvg();
// 		  std::cout << "AFTER = " << AFTER << std::endl;
// #endif

// 		  thisCorr.ensemble = thisCorr.ensemble * weight;

		  
// 		  // Get Ioffe-time
// 		  int ioffe(0);
// 		  if ( ! tmp.npoint[2].irrep.op.ops[1].disp_list.empty() )
// 		    ioffe = tmp.npoint[2].irrep.op.ops[1].disp_list.back();

// 		  // Conjugate if Ioffe-time is negative
// 		  if ( ioffe < 0 )
// 		    {
// 		      std::cout << "Ioffe = " << ioffe << " for key = ";
// 		      std::cout << tmp << std::endl;
// 		      NCOR::conj(&thisCorr.ensemble.ens);
// 		    }

// 		  /*
// 		    Only append this reweighted *.ensemble if weight != 0
// 		    Otherwise, merging below will average over non-zero 
// 		         & zero entries
// 		  */
// 		  // Pack this correlator for averaging
// 		  // toAvg[(j-1)*2+d].push_back(dum);
// 		  std::complex<double> null(0.0,0.0);
// 		  if ( weight != null )
// 		    toMerge[i-1].push_back(thisCorr.ensemble.ens);

// 		} // d
// 	    } // j
// 	  std::cout << "SIZE CHECK" << std::endl;
// 	  std::cout << "    toMerge[i-1].size() = " << toMerge[i-1].size() << std::endl;
// 	  for ( auto k = toMerge[i-1].begin(); k != toMerge[i-1].end(); ++k )
// 	    std::cout << "    toMerge[i-1][*].size() = " << k->size() << std::endl;
// 	} // i



//       // Erase the original keyCorrMap
//       vi->keyCorrMap.clear();
//       // Make a new keyCorrMap
//       ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, NCOR::correlator> cavg;


//       // Loop once more over dimensions of src & snk subductions to assign rows
//       for ( int i = 0; i < snkOp.irrep_dim; ++i )
// 	{
// 	  for ( int j = 0; j < srcOp.irrep_dim; ++j )
// 	    {
// 	      tmp.npoint[1].irrep.irrep_mom.row = (-1*(i+1));
// 	      tmp.npoint[3].irrep.irrep_mom.row = (-1*(j+1));
// 	      // Force disp key to be positive
// 	      if ( !tmp.npoint[2].irrep.op.ops[1].disp_list.empty()
// 		   && tmp.npoint[2].irrep.op.ops[1].disp_list.back() < 0 )
// 		tmp.npoint[2].irrep.op.ops[1].disp_list *= -1;


// 	      // Merge along toMerge[*] dimension
// 	      std::cout << "DEBUG - here we have " << toMerge[j+i*srcOp.irrep_dim].size() << " corrs to merge " << std::endl;
// 	      // NCOR::correlator mergeCorr = NCOR::mergeCorrs(toMerge[j+i*srcOp.irrep_dim]);


// 	      // Try a direct add
// 	      Pseudo::domain_t dumD(0,1,tmp.npoint[1].t_slice-1);
// 	      Pseudo::prop_t dumP(global.cfgs, dumD);
// 	      NCOR::correlator mergeCorr(dumP);

// 	      mergeCorr.ensemble.ens = toMerge[j+i*srcOp.irrep_dim][0];
// 	      if ( toMerge[j+i*srcOp.irrep_dim].size() > 1 )
// 		{
// 		  for ( int s = 1; s < toMerge[j+i*srcOp.irrep_dim].size(); ++s )
// 		    {
// 		      NCOR::dat_t dumDat; dumDat.ens = toMerge[j+i*srcOp.irrep_dim][s];
// 		      mergeCorr.ensemble += dumDat;
// 		      // NB(1/24/22): Is this where the doubling for z != 0 is coming from?
// 		    }
// 		}
// 	      // NB(1/24/22): Let's rescale
// 	      mergeCorr.ensemble *= (1.0/disp2Store);
		
// 	      cavg.insert(tmp, mergeCorr);
// 	    }
// 	}


//       // Assign cavg to vi->keyCorrMap and proceed to next iteration
//       vi->keyCorrMap = cavg;
      
//       // Clean up
//       gsl_matrix_complex_free(subduceSrc);
//       gsl_matrix_complex_free(subduceSnk);
//       gsl_matrix_complex_free(subduceTensProd);
//     } // vi
// }



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
  read(xmlSR, "/FormFact/global", global);
  /*
    Determine the db info - sets where/what edbs to search for keys
  */
  read(xmlSR, "/FormFact/dbInfo/threePt/tseries/range", temporal3pt);
  read(xmlSR, "/FormFact/dbInfo/threePt", db3ptInfo);
  read(xmlSR, "/FormFact/dbInfo/twoPtFin/tseries/range", temporal2ptFin);
  read(xmlSR, "/FormFact/dbInfo/twoPtFin", db2ptFinInfo);
  read(xmlSR, "/FormFact/dbInfo/twoPtIni/tseries/range", temporal2ptIni);
  read(xmlSR, "/FormFact/dbInfo/twoPtIni", db2ptIniInfo);
  read(xmlSR, "/FormFact/dbInfo/twoPtRest/tseries/range", temporal2ptRest);
  read(xmlSR, "/FormFact/dbInfo/twoPtRest", db2ptRestInfo);


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
    read(xmlSR, "/FormFact/Ops/threePt/keyTemplate", tempKey3pt.key);
    read(xmlSR, "/FormFact/Ops/twoPt/keyTemplate", tempKey2Pi.key);
    read(xmlSR, "/FormFact/Ops/twoPtRest/keyTemplate", tempKeyRest.key);
    tempKey2Pf = tempKey2Pi;
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
  twoPi.keydbs.dbs      = makeDBList(global, db2ptIniInfo, &tempKey2Pi.key);
  twoPi.keydbs.keys     = makeKeyList(tempKey2Pi.key);
  twoPf.keydbs.dbs      = makeDBList(global, db2ptFinInfo, &tempKey2Pf.key);
  twoPf.keydbs.keys     = makeKeyList(tempKey2Pf.key);
  twoPtRest.keydbs.dbs  = makeDBList(global, db2ptRestInfo, &tempKeyRest.key);
  twoPtRest.keydbs.keys = makeKeyList(tempKeyRest.key);


  /*
    Access and store all 2pt functions
  */
  // std::vector<NCOR::corrEquivalence> twoPtIni(1), twoPtFin(1);
  prop_t * propsIni = new prop_t(global.cfgs,temporal2ptIni,tempKey2Pi.key);
  propsIni->npt = 2;
  NCOR::correlator twoPtIni(*propsIni);
  delete propsIni;

  prop_t * propsFin = new prop_t(global.cfgs,temporal2ptFin,tempKey2Pf.key);
  propsFin->npt = 2;
  NCOR::correlator twoPtFin(*propsFin);
  delete propsFin;

  prop_t * props = new prop_t(global.cfgs,temporal2ptRest,tempKeyRest.key);
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
    Reads2pt(inFile, twoPtIni, db2ptIniInfo, &(*k) );
  for ( auto k = twoPf.keydbs.keys.begin(); k != twoPf.keydbs.keys.end(); ++k )
    Reads2pt(inFile, twoPtFin, db2ptFinInfo, &(*k) );
  for ( auto k = twoPtRest.keydbs.keys.begin(); k != twoPtRest.keydbs.keys.end(); ++k )
    Reads2pt(inFile, rest2pt, db2ptRestInfo, &(*k) );
#endif


  twoPtIni.jackknife(); twoPtFin.jackknife(); rest2pt.jackknife();
  twoPtIni.ensAvg();    twoPtFin.ensAvg();    rest2pt.ensAvg();
  twoPtIni.Cov();       twoPtFin.Cov();       rest2pt.Cov();

  // LinAlg::printMat(twoPtFin.cov.dat["real"]);
  // LinAlg::printMat(twoPtFin.cov.inv["real"]);

  // Get the fit info for ini/fin 2pts & 3pt
  NCOR::fitInfo_t twoPtRestFitInfo, twoPtFinFitInfo, twoPtIniFitInfo, threePtFitInfo;
  read(xmlSR, "/FormFact/fitting/twoPtFin", twoPtFinFitInfo);
  read(xmlSR, "/FormFact/fitting/twoPtIni", twoPtIniFitInfo);
  read(xmlSR, "/FormFact/fitting/twoPtRest", twoPtRestFitInfo);
  read(xmlSR, "/FormFact/fitting/threePt", threePtFitInfo);

  // Parse the strParamValMaps
  twoPtFinFitInfo.parseParamMaps();
  twoPtIniFitInfo.parseParamMaps();
  twoPtRestFitInfo.parseParamMaps();
  threePtFitInfo.parseParamMaps();

  std::cout << "PARSED" << std::endl;

  // To set up fit properly, pass the correlator's data covariance
  // Submatrix of covariance is internally grabbed, and its inverse computed
  twoPtFin.fit = NCOR::fitFunc_t(twoPtFinFitInfo,twoPtFin.cov.dat,temporal2ptFin);
  std::cout << "FIN FUNC SET" << std::endl;
  twoPtIni.fit = NCOR::fitFunc_t(twoPtIniFitInfo,twoPtIni.cov.dat,temporal2ptIni);
  std::cout << "INI FUNC SET" << std::endl;
  rest2pt.fit = NCOR::fitFunc_t(twoPtRestFitInfo,rest2pt.cov.dat,temporal2ptRest);
  std::cout << "REST FUNC SET" << std::endl;

  std::cout << "SET FIT FUNCS" << std::endl;

  // Fire up the fits
  std::vector<std::string> components(2);
  components[0] = "real"; components[1] = "imag";
  NFIT::driver(&twoPtFin, components[0], true); fitResW(&twoPtFin, components[0]);
  std::cout << "FIN FIT DONE" << std::endl;
  NFIT::driver(&twoPtIni, components[0], true); fitResW(&twoPtIni, components[0]);
  std::cout << "INI FIT DONE" << std::endl;
  NFIT::driver(&rest2pt, components[0], true);  fitResW(&rest2pt, components[0]);
  std::cout << "REST FIT DONE" << std::endl;

  std::cout << "DID THE FITS" << std::endl;

  writeCorr(&twoPtFin);
  writeCorr(&twoPtIni);
  writeCorr(&rest2pt);


  /*
    Do some checks of 2pt functions
  */
#if 0
  std::cout << "2pt BOOSTED KEYS" << std::endl; dumpKeys(twoPf.keydbs.keys,2);
  parseCheck(twoPt); std::cout << "\n";
  parseCheck(twoPt); std::cout << "\n";
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





  /*
    TRY SVD INSTEAD OF ALL THIS ROW AVERAGE BS
  */
  std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > MAT(global.cfgs);
  std::vector<ffMat_t> KIN(global.cfgs,ffMat_t(global.chromaGamma));
  for ( auto k = KIN.begin(); k != KIN.end(); ++k )
    {
      int j = std::distance(KIN.begin(),k);


      double inputFinE = 0.0;
      double inputIniE = 0.0;
#if 0
#warning "Energies passed to spinor constructors are from jackknife fits to C_2(pi) & C_2(pf)"
      inputFinE = twoPtFin.res.params["E0"][j];
      inputIniE = twoPtIni.res.params["E0"][j];
#else
#warning "Energies passed to spinor constructors will be determined from dispersion relation!"
      inputFinE = sqrt( pow(rest2pt.res.params["E0"][j],2) + pow(2*PI/global.Lx,2)*(global.pf*global.pf) );
      inputIniE = sqrt( pow(rest2pt.res.params["E0"][j],2) + pow(2*PI/global.Lx,2)*(global.pi*global.pi) );
#endif

      
      // Initialize the nucleon spinors for this jackknife sample
      Spinor finSpin(global.opMomXML[shortMom(global.pf,"")],
		     global.pf,inputFinE,rest2pt.res.params["E0"][j],global.Lx);
      Spinor iniSpin(global.opMomXML[shortMom(global.pi,"")],
		     global.pi,inputIniE,rest2pt.res.params["E0"][j],global.Lx);
      // Build the spinor
      finSpin.buildSpinors();
      iniSpin.buildSpinors();

      // Construct the kinematic matrix
      k->assemble(4,true,rest2pt.res.params["E0"][j],&finSpin,&iniSpin);

      std::cout << "KIN MAT!" << std::endl;
      std::cout << k->mat << std::endl;
    }
  std::vector<Eigen::Vector2cd> AMP(global.cfgs);



  /*
    Loop over all equivalent 3pt functions --> diff. rows, but kinematically the same
  */
  for ( std::vector<NCOR::corrEquivalence>::iterator tsepItr = funcs3pt.begin();
        tsepItr != funcs3pt.end(); ++tsepItr )
    {
      // int matIDX = std::distance(funcs3pt.begin(), tsepItr);

      // Convenience
      int rowf = tsepItr->keyCorrMap.begin()->first.npoint[1].irrep.irrep_mom.row;
      int rowi = tsepItr->keyCorrMap.begin()->first.npoint[3].irrep.irrep_mom.row;

      
      // Map the snk/src row combinations to a given element of correlator column vector
      std::map<std::pair<int,int>, int> matIDX;
      std::pair<int,int> rfri;

      rfri = std::make_pair(1,1); matIDX[rfri] = 0;
      rfri = std::make_pair(1,2); matIDX[rfri] = 1;
      rfri = std::make_pair(2,1); matIDX[rfri] = 2;
      rfri = std::make_pair(2,2); matIDX[rfri] = 3;

      /*
        For each correlator stored in each corrEquivalence,
        divide by 2pt function at same tsep
        Done per jackknife ensemble average
        Bias corrected after making ratio
        
        Sum up operator insertion time slice
      */
      std::vector<NCOR::correlator> ratio(nTSeps); // nTseps correlators of same rows/kinematics
      for ( auto it = tsepItr->keyCorrMap.begin(); it != tsepItr->keyCorrMap.end(); ++it )
        {
          // tsepItr->keyCorrMap is not ordered by tsep yet --> form an index to store in ascending order
          int idx = (it->first.npoint[1].t_slice - temporal3pt.min)/temporal3pt.step;
          // Local copy of this tsep
          const int TSEP = it->first.npoint[1].t_slice;


          // Init this ratio --> REMEMBER, 3PTS HAVE TSLICES\IN[0,TSEP)
	  Pseudo::domain_t * d = new Pseudo::domain_t(0,1,it->first.npoint[1].t_slice-1);
          prop_t * props = new prop_t(global.cfgs, *d, it->first);
          props->npt     = 3;
          props->gamma   = global.chromaGamma;

          
          // Now construct a new ratio
          ratio[idx] = NCOR::correlator(*props, it->second.ensemble);

          // Jackknife this ratio so 3pt/2pt ratio can be formed
          ratio[idx].jackknife();
          // ratio[idx].ensAvg();
	  // std::cout << ratio[idx] << std::endl;
          


          // Loop over insertion times for this 3pt TSEP
	  for ( auto tau = ratio[idx].ensemble.T.begin(); tau != ratio[idx].ensemble.T.end(); ++tau )
            {
              // Loop over the jackknife ensemble averages
              // ratio[idx] is for fixed T with tau varying
              for ( int j = 0; j < global.cfgs; ++j )
                {
		  // ratio[idx].ensemble.ens[j][*tau] = ratio[idx].jack[j].avg[*tau] / twoPtFin.jack[j].avg[TSEP].real();

		  ratio[idx].ensemble.ens[j][*tau] =
		    ( ratio[idx].jack[j].avg[*tau] / twoPtFin.jack[j].avg[TSEP].real() )
		    * sqrt(( twoPtIni.jack[j].avg[TSEP-*tau].real() * twoPtFin.jack[j].avg[*tau].real()
		  	     * twoPtFin.jack[j].avg[TSEP].real() )/
		  	   ( twoPtFin.jack[j].avg[TSEP-*tau].real() * twoPtIni.jack[j].avg[*tau].real()
		  	     * twoPtIni.jack[j].avg[TSEP].real() ));

		  

		  // ratio[idx].ensemble.ens[j][*tau] = ratio[idx].jack[j].avg[*tau] / twoPtFin.jack[j].avg[TSEP];

		  // ratio[idx].ensemble.ens[j][*tau] =
		  //   ( ratio[idx].jack[j].avg[*tau] / twoPtFin.jack[j].avg[TSEP] )
		  //   * sqrt(( twoPtIni.jack[j].avg[TSEP-*tau] * twoPtFin.jack[j].avg[*tau]
		  // 	     * twoPtFin.jack[j].avg[TSEP] )/
		  // 	   ( twoPtFin.jack[j].avg[TSEP-*tau] * twoPtIni.jack[j].avg[*tau]
		  // 	     * twoPtIni.jack[j].avg[TSEP] ));

                } // j
            } // tau
          


          // Include standard kinematic prefactors arising in forming optimized 3pt/2pt ratio
          for ( int j = 0; j < global.cfgs; ++j )
            {
	      // std::complex<double> commonKin(1.0,0.0);
	      std::complex<double> commonKin(sqrt(4*twoPtFin.res.params["E0"][j]*twoPtIni.res.params["E0"][j]), 0.0);

              // Remove common kinematic factor & 1/\sqrt(2) from isovector current normalization
	      ratio[idx].ensemble.ens[j] *= (redFact*commonKin);
            } // j


          // Get ensemble avg so bias removal can proceed
          ratio[idx].ensAvg();
#if 1
	  std::cout << "**************" << std::endl;
	  std::cout << "Ratio ens avg" << std::endl;
	  std::cout << ratio[idx] << std::endl;
	  std::cout << "**************" << std::endl;
#endif
          // Correct for bias in forming ratio
          ratio[idx].removeBias();


          // Summation of operator insertion
          ratio[idx].summation();

          delete d;
          delete props;
        } // it
      /*
        Now have summed ratio for a particular src/snk row combination
      */


      // Map the summed ratio data into correlator instance 'SR'
      // --> so covariance/fitting members can be used
      prop_t * xprops = new prop_t(global.cfgs, temporal3pt, ratio[0].key());
      xprops->npt     = 3;
      xprops->gamma   = global.chromaGamma;


      // Construct the single correlator instance 'SR'
      // ---> this will be fit
      NCOR::correlator SR(*xprops);
      
      for ( auto rptr = ratio.begin(); rptr != ratio.end(); ++rptr )
        {
          int ridx = std::distance(ratio.begin(), rptr);
          for ( auto gg = rptr->ensemble.ens.begin(); gg != rptr->ensemble.ens.end(); ++gg )
            {
              int gdx = std::distance(rptr->ensemble.ens.begin(), gg);
              SR.ensemble.ens[gdx][ridx] = (*gg)[0];
            } // gg
        } // rptr
      /*
        Now SR has been constructed
      */

      
      
      /*
        Perform linear fit to this summed ratio 'SR'
        --> exposes matrix element that will be fed into SVD to extract amplitudes
      */
      SR.jackknife();
      SR.ensAvg();

#if 1
      std::cout << "With key = " << SR.key() << " ..." << std::endl;
      std::cout << SR << std::endl;
#endif
#if 1

      // Make data covariance and initialize fit
      SR.Cov();
      SR.fit = NCOR::fitFunc_t(threePtFitInfo,SR.cov.dat,temporal3pt);

      writeCorr(&SR);
      /*
        Do the linear fits -- for both real/imag components
      */
      for ( auto f = components.begin(); f != components.end(); ++f )
        {
	  std::cout << "Before the fit" << std::endl;
	  // Do default fit
	  NFIT::driver(&SR, *f, false);

          // Write out the fit results
          fitResW(&SR, *f);


          // Pipe fit results foreach jackknife sample into appropriate entry of MAT
          for ( int g = 0; g < global.cfgs; ++g )
            {
	      std::pair<int,int> lookUp = std::make_pair(rowf,rowi);
              if ( *f == "real" )
                MAT[g](matIDX[lookUp]).real(SR.res.params["b"][g]);
              if ( *f == "imag" )
                MAT[g](matIDX[lookUp]).imag(SR.res.params["b"][g]);
            }

          // Destroy the stored fits values since they've been written
          SR.res.chi2.clear(); SR.res.params.clear();
        }
#endif
      delete xprops;

    } // funcs3pt iterator


  // With MAT populated per jackknife ensemble avg
  // Do the SVD per jackknife ensemble avg to extract amplitudes
  extAmplitudes(&MAT,&KIN,&AMP);
  std::cout << "What do these solutions look like?\n";

  for ( auto itr = AMP.begin(); itr != AMP.end(); ++itr )
    {
      int idx = std::distance(AMP.begin(),itr);
      
      double qSqrOver4MSqr = (1.0/(4*pow(rest2pt.res.params["E0"][idx],2)));
      qSqrOver4MSqr *= ( 2*(pow(rest2pt.res.params["E0"][idx],2)-
			    twoPtFin.res.params["E0"][idx]*twoPtIni.res.params["E0"][idx]
			    +pow((2*PI/global.Lx),2)*(global.pf*global.pi)) );

      std::complex<double> GE = (*itr)(0) + qSqrOver4MSqr*(*itr)(1);
      std::complex<double> GM = (*itr)(0) + (*itr)(1);

      qSqrOver4MSqr *= pow( hbarc/aLat, 2 ); // cast q^2/(4m^2) in physical units

      // Vector
      // std::cout << "F1 = " << (*itr)(0) << "     F2 = " << (*itr)(1) << std::endl;
      std::cout << "F1 = " << (*itr)(0) << "     F2 = " << (*itr)(1)
		<< "      GE = " << GE << "     GM = " << GM << "      q^2 (GeV^2) = " << qSqrOver4MSqr*(4*pow(rest2pt.res.params["E0"][idx],2)) << std::endl;
    }


#if 0  
  std::vector<Eigen::VectorXcd> finalAMP(global.cfgs,Eigen::VectorXcd(2));
  for ( auto itr = AMP.begin(); itr != AMP.end(); ++itr )
    {
      int idx = std::distance(AMP.begin(),itr);

      finalAMP[idx](0) = (*itr)(0);
      finalAMP[idx](1) = (*itr)(1);
    }

  writeAmplitudes(&finalAMP,&global,&threePtFitInfo,
		  &tempKey3pt.key.npoint[2].irrep.op.ops[1].disp_list);
#endif
  

  return 0;
}
