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

#include "summation.h"
#include "operators.h"
#include "fit_util.h"

#include <gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

// Bring in neccesary adat headers
#include "AllConfStoreDB.h"
#include "DBFunc.h"
#include "adat/map_obj.h"
#include "hadron/hadron_sun_npart_npt_corr.h"
#include "hadron/ensem_filenames.h"
#include "io/key_val_db.h"
#include "io/adat_xmlio.h"
#include "io/adat_xml_group_reader.h"

using namespace ADATXML;
using namespace Summation;
using namespace Pseudo;
using namespace FIT;

// Define a struct for global props
struct global
{
  int cfgs, t2pt, t2ptNorm;
} props;

// Define some structs to help determine what edbs to read
struct t3pt
{
  int tm3pt, ts3pt, tx3pt;
} temporal3pt;

struct info3pt
{
  // std::vector<std::string> base;
  std::string base;
  std::string momTag, tsnkTag, ensem, state, observable, nvec, t0Tag, zTag;
  // XMLArray::Array2d<int> rows;
  // XMLArray::Array2d<int> signs;

  XMLArray::Array<XMLArray::Array<int> > rows;
  XMLArray::Array<XMLArray::Array<int> > signs;
  
} db3ptInfo;

struct info2pt
{
  // std::vector<std::string> base;
  std::string base;
  std::string momTag, ensem, state, nvec, t0Tag;
  XMLArray::Array<XMLArray::Array<int> > rows;
  XMLArray::Array<XMLArray::Array<int> > signs;
} db2ptInfo;


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


// // Normalize a matelem of some non-trivial displacement by corresponding trivial displacement matelem
// void normalize(std::vector<correlators> &someNormedRatio, std::vector<correlators> &someRatio,
// 	       std::vector<correlators> &someRatio_noSep,int cfgs)
// {
//   std::vector<correlators>::iterator rSM;
//   //#pragma omp parallel for
// //   {
// // #pragma omp for
//   for ( rSM = someNormedRatio.begin(); rSM != someNormedRatio.end(); ++rSM )
//     {
//       auto c = rSM - someNormedRatio.begin();
      
//       rSM->ncor.resize(cfgs);
//       for ( int g = 0; g < cfgs; ++g )
// 	{
// 	  rSM->ncor[g].real.resize(1); rSM->ncor[g].imag.resize(1);
	  
// 	  double rd,id,rd0,id0;
// 	  rd  = someRatio[c].ncor[g].real[0][0];
// 	  id  = someRatio[c].ncor[g].imag[0][0];
// 	  rd0 = someRatio_noSep[c].ncor[g].real[0][0];
// 	  id0 = someRatio_noSep[c].ncor[g].imag[0][0];
	  
// 	  rSM->ncor[g].real[0].push_back( (rd*rd0+id*id0)/(pow(rd0,2)+pow(id0,2)) );
// 	  rSM->ncor[g].imag[0].push_back( (rd0*id-rd*id0)/(pow(rd0,2)+pow(id0,2)) );
// 	}
//     }
//   // } // #pragma omp parallel
// }



/*
  SOME CALLS TO HELP ADAT READ
*/
// Read a correlator edb
// std::map<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators>
ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators>
getCorrs(const std::vector<std::string>& dbases, std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t>& fetch)
{
  typedef Hadron::KeyHadronSUNNPartNPtCorr_t       K;
  typedef ENSEM::VectorComplex                     V;

  // Initialize map returned to main associating keys with correlators datatypes
  // n.b. cannot use a std::map here as I haven't created a comparator
  ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> retMap;



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
      retMap.insert(keyFetch[k].key(),ensAdat);
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

  return retMap;
}


// Reader for global properties
void read(XMLReader& xml, const std::string path, global& g)
{
  try {
    read(xml, path+"/cfgs", g.cfgs);
    read(xml, path+"/t2pt", g.t2pt);
    read(xml, path+"/t2ptNorm", g.t2ptNorm);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse global properties " << e << std::endl;
  }
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
    read(xml, path+"/ensem", d.ensem);
    read(xml, path+"/state", d.state);
    read(xml, path+"/observable", d.observable);
    read(xml, path+"/nvec", d.nvec);
    read(xml, path+"/t0Tag", d.t0Tag);
    read(xml, path+"/zTag", d.zTag);
    read(xml, path+"/rowinfo/rows", d.rows);
    read(xml, path+"/rowinfo/signs", d.signs);
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
    read(xml, path+"/ensem", d.ensem);
    read(xml, path+"/state", d.state);
    read(xml, path+"/nvec", d.nvec);
    read(xml, path+"/t0Tag", d.t0Tag);
    read(xml, path+"/rowinfo/rows", d.rows);
    read(xml, path+"/rowinfo/signs", d.signs);
  } catch ( std::string &e ) {
    std::cerr << "Unable to parse db2ptInfo struct from ini xml " << e << std::endl;
    exit(1);
  }
}



// Build a short version of momentum - inspired by adat/lib/hadron/irrep_util.cc
std::string shortMom(const Array<int>& mom, const std::string substr)
{
  std::ostringstream os;
  os << substr << mom[0] << substr << mom[1] << substr << mom[2];
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
  Utility to make 2pt db list
*/
std::vector<std::string> makeDBList(info2pt& I, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<std::string> s;
  
  // Get the momentum of the template key
  Array<int> momPlus = kTemplate->npoint[1].irrep.irrep_mom.mom;

  // // Loop over possible directory trees
  // for ( auto d = I.base.begin(); d != I.base.end(); ++d )
  //   {
      // A temporary db to set  -- *d
      std::string tmp_db_plus = I.base+"/"+I.momTag+shortMom(momPlus,".")+
	"/EDB/"+I.ensem+"."+I.state+"_p"+shortMom(momPlus,"")+".n"+I.nvec+"."+I.t0Tag+".edb";
      
      // Push to db list
      s.push_back(tmp_db_plus);
      
      if ( shortMom(momPlus,"") != "000" )
	{
	  Array<int> momMinus = momPlus*-1;
	  // A temporary db to set
	  std::string tmp_db_minus = I.base+"/"+I.momTag+shortMom(momMinus,".")+
	    "/EDB/"+I.ensem+"."+I.state+"_p"+shortMom(momMinus,"")+".n"+I.nvec+"."+I.t0Tag+".edb";
	  // Push to db list
	  s.push_back(tmp_db_minus);
	}
    // }
  return s;
}

/*
  Utility to make 3pt db list
*/
std::vector<std::string> makeDBList(info3pt& I, t3pt& t, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<std::string> s;
  for ( int ti = t.tm3pt; ti <= t.tx3pt; ti+=t.ts3pt )
    {
      // Get the momentum of the template key
      Array<int> momPlus = kTemplate->npoint[1].irrep.irrep_mom.mom;

      // // Loop over possible directory trees
      // for ( auto d = I.base.begin(); d != I.base.end(); ++d )
      // 	{
	  // A temporary db to set
	  std::string tmp_db_plus = I.base+"/"+I.tsnkTag+std::to_string(ti)+"/"+I.momTag+
	    shortMom(momPlus,".")+"/EDB/"+I.ensem+"."+I.state+"."+
	    I.observable+"_pf"+shortMom(momPlus,"")+"_pi"+
	    shortMom(momPlus,"")+".n"+I.nvec+"."+I.t0Tag+"_"+
	    I.tsnkTag+std::to_string(ti)+"."+I.zTag+".edb";
	  
	  // Push to db list
	  s.push_back(tmp_db_plus);
	  
	  if ( shortMom(momPlus,"") != "000" )
	    {
	      Array<int> momMinus = momPlus*-1;
	      // A temporary db to set
	      std::string tmp_db_minus = I.base+"/"+I.tsnkTag+std::to_string(ti)+"/"+I.momTag+
		shortMom(momMinus,".")+"/EDB/"+I.ensem+"."+I.state+"."+
		I.observable+"_pf"+shortMom(momMinus,"")+"_pi"+
		shortMom(momMinus,"")+".n"+I.nvec+"."+I.t0Tag+"_"+
		I.tsnkTag+std::to_string(ti)+"."+I.zTag+".edb";
	      // Push to db list
	      s.push_back(tmp_db_minus);
	    }
	// }
    }
  return s;
}


/*
  Utility to make 2pt key list
*/
std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> makeKeyList(Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> s;

  // Catch the momentum from the temporary key
  Array<int> mom = kTemplate->npoint[1].irrep.irrep_mom.mom;

  s.push_back(*kTemplate);

  // Copy and create the negative momentum keys, provided shortMom != 000
  if ( shortMom(mom,"") != "000" )
    {
      Hadron::KeyHadronSUNNPartNPtCorr_t h = *kTemplate;
      h.npoint[1].irrep.irrep_mom.mom *= -1;
      h.npoint[2].irrep.irrep_mom.mom *= -1;
      // Now add this negative momentum 2pt key
      s.push_back(h);
    }
  return s;
}


/*
  Utility to make 3pt key list
*/
std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> makeKeyList(t3pt& t, Hadron::KeyHadronSUNNPartNPtCorr_t *kTemplate)
{
  std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> s;
  // Catch the momentum from the temporary key
  Array<int> mom = kTemplate->npoint[1].irrep.irrep_mom.mom;

  // With template 3pt motion key stored, iterate over tseps and make remaining motion keys
  for ( int ti = t.tm3pt; ti <= t.tx3pt; ti+=t.ts3pt )
    {
      // Make a new key from the template
      Hadron::KeyHadronSUNNPartNPtCorr_t tmp_k_plus = *kTemplate;
      
      // Change sink t_slice
      tmp_k_plus.npoint[1].t_slice=ti;
      // Push this tmp key
      s.push_back(tmp_k_plus);

      // Make a key for opposing displacement, if disp_list != ''
      if ( tmp_k_plus.npoint[2].irrep.op.ops[1].disp_list.size() > 0 )
	{
	  for ( int d = 0; d < tmp_k_plus.npoint[2].irrep.op.ops[1].disp_list.size(); ++d )
	    {
	      tmp_k_plus.npoint[2].irrep.op.ops[1].disp_list[d]*=-1;
	    }
	  // insert this key with opposite displacement
	  s.push_back(tmp_k_plus);
	}
      // // Include the null displacement insertion as well
      // Hadron::KeyParticleOp_t * opPtr = &tmp_k_plus.npoint[2].irrep.op.ops[1]; // pointer to relevant stuff
      // tmp_k_plus.npoint[2].irrep.op.ops[1] = Hadron::KeyParticleOp_t(opPtr->name,opPtr->smear,opPtr->mom_type);

      // // Insert the null displacement key
      // s.push_back(tmp_k_plus);
    }

  // Copy and create the negative momentum keys, provided shortMom != 000
  if ( shortMom(mom,"") != "000" )
    {
      // Copy current keys to a new vector and negate src/snk momentum
      // Thereby creating all negative momentum keys
      std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> minusMomKeys(s);
      for ( auto a = minusMomKeys.begin(); a != minusMomKeys.end(); ++a )
	{
	  a->npoint[1].irrep.irrep_mom.mom *= -1;
	  a->npoint[3].irrep.irrep_mom.mom *= -1;
	  // Now add this negative momentum key to the std::vector
	  s.push_back(*a);
	}
    }

  return s;
}

/*
  Utility to retrive & order correlators stored in keyscorrs MapObjects
*/
void pop3ptFuncs(std::vector<corrFunc> &F, Hadron::KeyHadronSUNNPartNPtCorr_t &h,
		 ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators>& kc)
{
  for ( std::vector<corrFunc>::iterator f = F.begin(); f != F.end(); ++f )
    {
      // Local copy of this tsep & snk mom
      int t = f->getT();
      std::string mom = shortMom(h.npoint[1].irrep.irrep_mom.mom,"");
      
      // Set a key to the template key
      Hadron::KeyHadronSUNNPartNPtCorr_t sortKey = h;
      sortKey.npoint[1].t_slice=t;

      // Momentum combinations to lookup
      std::vector<Array<int> > momLookup(1,sortKey.npoint[1].irrep.irrep_mom.mom);
      if ( mom != "000" )
	{
	  Array<int> tmp = momLookup[0]; tmp *= -1;	  
	  momLookup.push_back(tmp);
	}

      // Displacement combinations to lookup
      std::vector<std::vector<int> > dispLookup(1,sortKey.npoint[2].irrep.op.ops[1].disp_list);
      if ( dispLookup[0].size() > 0 )
	{
	  std::vector<int> tmp;
	  for ( auto dd = dispLookup[0].begin(); dd != dispLookup[0].end(); ++dd )
	    {
	      tmp.push_back( (*dd)*-1 );
	    }
	  dispLookup.push_back(tmp);
	}
      

      // Now we know how many correlators need to be merged
      int num2Merge = momLookup.size()*dispLookup.size();
      correlators toMerge(num2Merge);

      std::cout << "3pt NUM2MERGE = " << num2Merge << std::endl;

      int cnt = 0; // talley of correlator included in toMerge total
      for ( int p = 0; p < momLookup.size(); ++p )
        {
          for ( int z = 0; z < dispLookup.size(); ++z )
            {
              // Set remainder of sortKey
              sortKey.npoint[1].irrep.irrep_mom.mom = momLookup[p];
              sortKey.npoint[3].irrep.irrep_mom.mom = momLookup[p];
              sortKey.npoint[2].irrep.op.ops[1].disp_list = dispLookup[z];

              // Retrieve the correlators datatype associated with this key
              ensemble retV = kc[sortKey].ncor[0];


	      // Sloppy capture of the sign of the ioffe time here!
	      int ioffe;
	      if ( dispLookup[z].empty() ) { ioffe = 0; }
	      else { ioffe = momLookup[p]*dispLookup[z]; }

	      // Conjugate if ioffe < 0
              if ( ioffe < 0 )
                {
                  cmplxConj(retV);
                }

              // Tack the retrieved correlators datatype onto toMerge object
              toMerge.ncor[cnt] = retV;
              cnt++;
            } // end z
        } // end p

      // Now we have all correlators at same tsep, and with same abs(pz) abs(z)
      // So now we merge them
      f->data = mergeCorrelators(toMerge);
    } // end F iterator
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
void parseCheck(ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> A)
{
  for ( auto kc = A.begin(); kc != A.end(); ++kc )
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


  /*
    Determine global properties
  */
  read(xmlSR, "/SummationPPDF/global", props);
  /*
    Determine the db info - sets where/what edbs to search for keys
  */
  read(xmlSR, "/SummationPPDF/dbInfo/threePt/tseries", temporal3pt);
  read(xmlSR, "/SummationPPDF/dbInfo/threePt", db3ptInfo);
  read(xmlSR, "/SummationPPDF/dbInfo/twoPt", db2ptInfo);


  // std::cout << "3PT ROW STUFF" << std::endl;
  // std::cout << db3ptInfo.rows[0][0] << " " << db3ptInfo.rows[0][1] << std::endl;
  // std::cout << db3ptInfo.rows[1][0] << " " << db3ptInfo.rows[1][1] << std::endl;
  // std::cout << "s: " << db3ptInfo.signs[0][0] << " " << db3ptInfo.signs[0][1] << std::endl;
  // std::cout << "s: " << db3ptInfo.signs[1][0] << " " << db3ptInfo.signs[1][1] << std::endl;
  // std::cout << "2PT ROW STUFF" << std::endl;
  // std::cout << db2ptInfo.rows[0][0] << " " << db2ptInfo.rows[0][1] << std::endl;
  // std::cout << db2ptInfo.rows[1][0] << " " << db2ptInfo.rows[1][1] << std::endl;
  // std::cout << "s: " << db2ptInfo.signs[0][0] << " " << db2ptInfo.signs[0][1] << std::endl;
  // std::cout << "s: " << db2ptInfo.signs[1][0] << " " << db2ptInfo.signs[1][1] << std::endl;
  // exit(8);


  // Structure to hold all desired 3pt/2pt motion/normalizing keys
  struct keydb
  {
    std::vector<Hadron::KeyHadronSUNNPartNPtCorr_t> keys;
    std::vector<std::string>                        dbs;
    // Map of keys and correlators structs returned
    // std::map<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> keyscorrs; --> need a comparator
    ADAT::MapObject<Hadron::KeyHadronSUNNPartNPtCorr_t, correlators> keyscorrs;
  };

  struct nptKVs
  {
    keydb keydbs;
  } threePt, twoPt, threePtRC2, twoPtRC2; // W/ RC2 = ROW COMBO 2


  struct keyTemplate
  {
    Hadron::KeyHadronSUNNPartNPtCorr_t key;
  } tempKey3pt, tempKey2pt, tempKey3ptRC2, tempKey2ptRC2;


  // Read the keys from the ini xml
  try {
    read(xmlSR, "/SummationPPDF/Ops/threePt/keyTemplate", tempKey3pt.key);
    read(xmlSR, "/SummationPPDF/Ops/twoPt/keyTemplate", tempKey2pt.key);
  } catch ( std::string &e ) {
    std::cerr << "Unable to access key templates from ini xml " << e << std::endl;
  }
  std::cout << "READ THE TEMPLATE KEYS TO FETCH" << std::endl;


  /*   Make the 3pt database structures to search   */
  threePt.keydbs.dbs    = makeDBList(db3ptInfo, temporal3pt, &tempKey3pt.key);
  threePtRC2.keydbs.dbs = threePt.keydbs.dbs;
  // tempKey3ptRC2 = tempKey3pt; // mod. template key struct to allow for 3pt RC2 to be accessed
  // Make all 3pt keys from templates
  threePt.keydbs.keys    = makeKeyList(temporal3pt, &tempKey3pt.key);

  /*   Make the 2pt database structures to search   */
  twoPt.keydbs.dbs      = makeDBList(db2ptInfo, &tempKey2pt.key);
  twoPt.keydbs.keys      = makeKeyList(&tempKey2pt.key); // make all 2pt keys from templates
  // Skip the row2-row2 2pt combo at rest
  if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
    {
      twoPtRC2.keydbs.dbs = twoPt.keydbs.dbs;
      tempKey2ptRC2 = tempKey2pt; // mod. template key struct to allow for 2pt RC2 to be accessed
      makeKeyRC2(tempKey2ptRC2.key, db2ptInfo.rows[1]); // pass second row combo and make new 2pt template key
      twoPtRC2.keydbs.keys   = makeKeyList(&tempKey2ptRC2.key); // make all RC2 2pt keys from templates
    }

# warning "-----------NO ROW AVERAGING IN 3PT FUNCTIONS (YET)!--------------"


  for ( int k = 0; k < threePt.keydbs.keys.size(); k++ )
    {
      std::cout << threePt.keydbs.keys[k] << std::endl;
    }


  // std::cout << "Number of RC1 2pt keys = " << twoPt.keydbs.keys.size() << std::endl;
  // for ( int k = 0; k < twoPt.keydbs.keys.size(); k++ )
  //   {
  //     std::cout << twoPt.keydbs.keys[k] << std::endl;
  //   }
  // std::cout << "Number of RC2 2pt keys = " << twoPtRC2.keydbs.keys.size() << std::endl;
  // for ( int k = 0; k < twoPt.keydbs.keys.size(); k++ )
  //   {
  //     std::cout << twoPtRC2.keydbs.keys[k] << std::endl;
  //   }
  // exit(9);



  /*
    Now go fetch all the correlators
  */
  threePt.keydbs.keyscorrs    = getCorrs(threePt.keydbs.dbs,threePt.keydbs.keys);
  // threePtRC2.keydbs.keyscorrs = getCorrs(threePtRC2.keydbs.dbs,threePtRC2.keydbs.keys);
#ifdef HAVE_DISPLIST_2PT
#warning "Using getCorrs to read 2pt correlators"
  twoPt.keydbs.keyscorrs      = getCorrs(twoPt.keydbs.dbs,twoPt.keydbs.keys);
  // Skip the row2-row2 2pt combo at rest
  if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
    {
      twoPtRC2.keydbs.keyscorrs   = getCorrs(twoPtRC2.keydbs.dbs,twoPtRC2.keydbs.keys);
    }
#endif


#if 0
  parseCheck(threePt.keydbs.keyscorrs); std::cout << "\n";
  parseCheck(twoPt.keydbs.keyscorrs); std::cout << "\n";
  exit(9);
#endif
  std::cout << "Successfully performed the reads!" << std::endl;




  /*
    Need to unpack the keys and push into corrFunc instances to reuse all code below

    funcs3pt:
          --> threePt.motion.keys  + momenta
	  --> threePt.motion.keys  - momenta
  */

  // A messy(ish) initialization of these std::vectors
  std::vector<corrFunc> funcs3pt;
  for ( int t = temporal3pt.tm3pt; t <= temporal3pt.tx3pt; t+=temporal3pt.ts3pt )
    {
      funcs3pt.push_back(corrFunc(props.cfgs,t));
    }
  pop3ptFuncs(funcs3pt, tempKey3pt.key, threePt.keydbs.keyscorrs);

  

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

 

  /*
    Initialize the 2pt functions
  */
  std::cout << "Initializing the 2pt functions" << std::endl;

  
  // Start by initializing single master 2pt correlator
  corrFunc func2pt;
  // Initialize a correlators database to hold all 2pt correlators to merge - defaults to ncor.size() = 1
  correlators toMerge2pts; // minimally 1-set to merge (e.g. rest w/ 1 RCs)



  // Create corrFunc objects for rest (or +p mom) RC 1 & 2
  corrFunc func2ptRC1(props.cfgs,props.t2pt), func2ptRC2(props.cfgs,props.t2pt);
  func2ptRC1.initialize(); func2ptRC2.initialize();
#ifndef HAVE_DISPLIST_2PT
#warning "   2pt correlators constructed w/o disp_list  -->  reverting to dumb 2pt reader"
  std::ifstream inFile;
  Reads2pt(inFile, func2ptRC1, db2ptInfo, &twoPt.keydbs.keys[0]);
  if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
    {
      Reads2pt(inFile, func2ptRC2, db2ptInfo, &twoPtRC2.keydbs.keys[0]);
    }
#endif

  // Push into toMerge2pts
  toMerge2pts.ncor[0] = func2ptRC1.data.ncor[0];

  
  // If shortMom != 000, then we must extend toMerge2pts once for the RC2 combo of 2pt w/ pz > 0
  // & initialize two more corrFunc objects for pz < 0 [ RC 1 & 2 ]
  if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
    {
      // Insert pz > 0 RC2
      toMerge2pts.ncor.push_back(func2ptRC2.data.ncor[0]);

      // Now prep stuff for the 2pts w/ pz < 0
      corrFunc func2ptMinusRC1(props.cfgs,props.t2pt),func2ptMinusRC2(props.cfgs,props.t2pt);
      func2ptMinusRC1.initialize(); func2ptMinusRC2.initialize();
#ifndef HAVE_DISPLIST_2PT
#warning "   2pt correlators constructed w/o disp_list  -->  reverting to dumb 2pt reader"
      std::ifstream inFile;
      Reads2pt(inFile, func2ptMinusRC1, db2ptInfo, &twoPt.keydbs.keys[1]);
      Reads2pt(inFile, func2ptMinusRC2, db2ptInfo, &twoPtRC2.keydbs.keys[1]);
#endif
      // Append these to the toMerge2pts object
      toMerge2pts.ncor.push_back(func2ptMinusRC1.data.ncor[0]);
      toMerge2pts.ncor.push_back(func2ptMinusRC2.data.ncor[0]);
      
#if 0
      std::cout << "2pt Minus Parse Check!" << std::endl;
      ensemble dum = corrAvg(func2ptMinusRC1.data);
      arr_print(dum);
      dum = corrAvg(func2ptMinusRC2.data);
      arr_print(dum);
#endif
    } // if shortMom(2ptkey"--") != 000

#if 0
  ensemble dum = corrAvg(func2ptRC1.data);
  std::cout << "2pt Parse Check!" << std::endl;
  arr_print(dum);
  if ( shortMom(tempKey2pt.key.npoint[1].irrep.irrep_mom.mom,"") != "000" )
    { dum = corrAvg(func2ptRC2.data); arr_print(dum); }
  // exit(9);
#endif

  // To get all info, reset func2pt equal to func2ptRC1
  func2pt = func2ptRC1;
  std::cout << "---> Averaging all rows/moms 2pt correlators..." << std::endl;
  std::cout << toMerge2pts.ncor.size() << std::endl;
  std::cout << toMerge2pts.ncor[0].real.size() << std::endl;
  std::cout << toMerge2pts.ncor[0].real[0].size() << std::endl;
  func2pt.data.ncor[0] = mergeCorrelators(toMerge2pts).ncor[0];
  std::cout << "<--- Finished averaging all rows/moms 2pt correlators..." << std::endl;


  /*
    ###########################################################################################
    BELOW HERE SHOULD BE ALL AVERAGED/CONJUGATED CORRELATORS...
    ###########################################################################################
  */

  // Print averages of 2pt & 3pt correlator
#if 0
  ensemble C2pt_avg = corrAvg(func2pt.data);
  std::vector<ensemble> C3pt_avgs(funcs3pt.size());
  std::cout << "C2pt_avg = " << std::endl;
  arr_print(C2pt_avg);
  std::cout << "******************" << std::endl;
  for ( int i = 0; i < funcs3pt.size(); ++i )
    {
      C3pt_avgs[i] = corrAvg(funcs3pt[i].data);
      arr_print(C3pt_avgs[i]);
      std::cout << "*************" << std::endl;
    }
  exit(8);
#endif



  /*
    Make 2pt jackknife samples for each correlator
  */
  std::cout << "---> Making 2pt jackknife samples..." << std::endl;
  func2pt.initializeJks();
  fill_jk_arrays(func2pt.dataJk, func2pt.data);
  std::cout << "---> Determining 2pt jackknife ensemble averages..." << std::endl;
  func2pt.initializeJkEnsAvgs();
  jkEnsAvg(func2pt.dataJk,func2pt.dataJkEnsAvg);

  /*
    Make 3pt jackknife samples for each correlator
  */
  std::cout << "---> Making 3pt jackknife samples for each correlator..." << std::endl;
  std::cout << "+++> & Determining ensemble averages per jackknife sample..." << std::endl;
  // Iterate over all 3pt data
  //#pragma omp parallel for
//   {
// #pragma omp for
  for ( std::vector<corrFunc>::iterator it = funcs3pt.begin(); it != funcs3pt.end(); ++it )
    {
      auto c = it - funcs3pt.begin();
      //*****************************
      // Initialize the 3pts
      (*it).initializeJks();
      fill_jk_arrays((*it).dataJk,(*it).data);
      (*it).initializeJkEnsAvgs();
      jkEnsAvg((*it).dataJk,(*it).dataJkEnsAvg);
    }
  // } // #pragma omp parallel

  std::cout << "<--- Completed jackknifing/jk averaging of 3pt data..." << std::endl;
  

  std::cout << "---> Performing the summation of current insertions per jackknife ensemble average"
	    << std::endl;
  // Now sum up the insertions for each jackknife ensemble average
  summation(funcs3pt);
  std::cout << "<--- Finished the summation of current insertions per jackknife ensemble average"
	    << std::endl;



  /*
    Form the summed ratio for each jackknife ensemble average
  */
  Ratios SR(props.cfgs, temporal3pt.tm3pt, temporal3pt.ts3pt, temporal3pt.tx3pt);
  SR.makeSR(funcs3pt, func2pt);
  /*
    Determine covariance + inverse of covariance matrices
  */
  SR.mean();
  SR.makeCovs();
  SR.makeInvCovs();


  /*
    Do the linear fit
    But first make a good name for the file
  */
  tempKey3pt.key.npoint[1].t_slice = -60; // set the output ensem - w/ -60 indicating a summed ratio
  std::string fitFileName = Hadron::ensemFileName(tempKey3pt.key);
  fitFileName.erase(fitFileName.end()-4,fitFileName.end());
  fitFileName += ".summedRatio";

  SR.fit("LINEAR", fitFileName);


#if 0
  /*
    Correct for bias in each jackknife ensemble average estimate of summed ratio
    & concatenate corrected jackknife ensemble average estimates per tsep
  */
  std::vector<correlators> antiJkSummedRatio(funcs3pt.size());
  std::vector<correlators> catAntiJkSummedRatio(funcs3pt.size());
  //#pragma omp parallel for
//   {
// #pragma omp for
  for ( std::vector<correlators>::iterator it = antiJkSummedRatio.begin(); it != antiJkSummedRatio.end(); ++it )
    {
      // Again catch the iterator value
      auto c = it-antiJkSummedRatio.begin();
      // Correct for the bias
      *it = antiJkRatio(ratioSumAvgEst[c],ratioSum[c]);
      // Concatenate bias corrected jackknife estimates
      catAntiJkSummedRatio[c] = catAntiJkEnsRatio(*it);
    }
  // } // #pragma omp parallel


  /*
    Write out ensemble of bias corrected jackknife ensemble average summed ratios
  */
  std::cout << "---> Writing ensemble of bias corrected jackknife ensemble average summed ratios" << std::endl;

  tempKey3pt.key.npoint[1].t_slice = -60; // set the output ensem - w/ -60 indicating a summed ratio
  std::string output = Hadron::ensemFileName(tempKey3pt.key);
  output.erase(output.end()-4,output.end());
  // std::string suffix = shortMom(threePt.normalize.keys[0].npoint[1].irrep.irrep_mom.mom,"");

  
  std::ofstream outFile;
  // outFile.open(output+".summedRatio_p"+suffix+".dat");
  outFile.open(output+".summedRatio.dat");
  outFile << func2pt.getCfgs() << " " << std::to_string(funcs3pt.size()) << " 1 0 1\n";
  for ( int g = 0; g < func2pt.getCfgs(); ++g )
    {
      for ( std::vector<correlators>::iterator it = catAntiJkSummedRatio.begin();
	    it != catAntiJkSummedRatio.end(); ++it )
	{
	  // Catch the iterator value
	  auto c = it - catAntiJkSummedRatio.begin();
	  outFile << std::setprecision(10) << std::to_string(c) << " "
		  << (*it).ncor[0].real[g][0] << " " << (*it).ncor[0].imag[g][0] << "\n";
	}
    }
  outFile.close();
  std::cout << "<--- Completed writing of ensemble of bias corrected jackknife ensemble average summed ratios"
	    << std::endl;



  /*
    Create jackknife samples from bias corrected ensemble for each tsep
    & write to file for XMBF fitting
  */
  std::vector<correlators> jkRatioSamples(funcs3pt.size());
  for ( auto it = jkRatioSamples.begin(); it != jkRatioSamples.end(); ++it )
    {
      // Catch the iterator value
      auto c = it - jkRatioSamples.begin();
      std::string output = "summedRatio_tsnk"+std::to_string(funcs3pt[c].getT());
      it->ncor.resize(func2pt.getCfgs());
      for ( int j = 0; j < (*it).ncor.size(); ++j )
	{
	  (*it).ncor[j].real.resize(func2pt.getCfgs()-1);
	  (*it).ncor[j].imag.resize(func2pt.getCfgs()-1);
	  for ( int g = 0; g < (*it).ncor[0].real.size(); ++g )
	    {
	      // We are concatenating summedRatio jackknife samples here,
	      // so insertion time dependence is removed - i.e. only 1 entry for each real/imag
	      (*it).ncor[j].real[g].resize(1);
	      (*it).ncor[j].imag[g].resize(1);
	    }
	}
      
      fill_jk_arrays(*it,catAntiJkSummedRatio[c]);

      // Write out the newly created jackknife samples here
      std::cout << "---> Writing jackknife ratio files" << std::endl;
      for ( int j = 0; j < func2pt.getCfgs(); ++j )
	{
	  outFile.open(output+"_jack"+std::to_string(j)+".dat");
	  outFile << func2pt.getCfgs()-1 << " 1" << " 1 0 1\n";
	  for ( int g = 0; g < (*it).ncor[0].real.size(); ++g )
	    {
	      outFile << std::setprecision(10) << "0 " << (*it).ncor[j].real[g][0]
		      << " " << (*it).ncor[j].imag[g][0] << "\n";
	    }
	  outFile.close();
	}
    }
  std::cout << "<--- All jackknife ratio files written..." << std::endl;  
#endif

  return 0;
}
