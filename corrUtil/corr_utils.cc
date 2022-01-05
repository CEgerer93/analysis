/*
  Member functions of correlator class
 */

#include<sstream>
#include "corr_utils.h"
#include "pseudo_utils.h"


// namespace NCOR
// {
/*
  OPERATORS
*/
// Define operator to easily print contents of a gsl_vector
std::ostream& operator<<(std::ostream& os, const gsl_vector *v)
{
  if ( v->size > 0 )
    {
      os << gsl_vector_get(v,0);
      for ( int i = 1; i < v->size; ++i )
	os << "," << gsl_vector_get(v,i);
    }
  return os;
}

namespace NCOR { 
// template<class C, typename T>
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T>& v)
{
  if ( v.size() > 0 )
    {
      for ( auto iv = v.begin(); iv != v.end(); ++iv )
	os << *iv << " ";
      os << "\n";
    }
  return os;
}

// template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<std::complex<double> >& v)
{
  if ( v.size() > 0 )
    {
      for ( auto iv = v.begin(); iv != v.end(); ++iv )
	os << *iv << " ";
      os << "\n";
    }
  return os;
}
}


namespace NCOR
{

  // // template<typename T>
  // std::vector<int>& operator*=(std::vector<int>& v, int i)
  // {
  //   for ( std::vector<int>::iterator it = v.begin(); it != v.end(); ++it )
  //     (*it)*i;
  //   return v;
  // }

  // template<typename T>
  std::vector<std::complex<double> > operator*=(std::vector<std::complex<double> >& v, double d)
  {
    for ( auto it = v.begin(); it != v.end(); ++it )
      (*it)*d;
    return v;
  }

  // template<typename T>
  std::vector<std::complex<double> > operator*=(std::vector<std::complex<double> >& v,
						std::complex<double> c)
  {
    for ( auto vi = v.begin(); vi != v.end(); ++vi )
      (*vi) *= c;
    return v;
  }

  std::vector<std::complex<double> > operator+=(std::vector<std::complex<double> >& v1,
						std::vector<std::complex<double> >& v2)
  {
    // Check dimensions of vectors are equal
    if ( v1.size() != v2.size() )
      {
	std::cout << "Vectors of unequal size!" << std::endl;
	exit(2);
      }

    for ( int n = 0; n < v1.size(); ++n )
      v1[n] += v2[n];

    return v1;
  }

  template<typename T>
  std::vector<std::complex<T> > operator-=(std::vector<std::complex<T> >& v1,
					   std::vector<std::complex<T> >& v2)
  {
    // Check dimensions of vectors are equal
    if ( v1.size() != v2.size() )
      {
	std::cout << "Vectors of unequal size!" << std::endl;
	exit(2);
      }

    for ( int n = 0; n < v1.size(); ++n )
      v1[n] -= v2[n];

    return v1;
  }

  std::ostream& operator<<(std::ostream& os, correlator &c)
  {
    os << "{cfgs, Nt} = {" << c.getCfgs() << ", " << c.getNt() << "}\n";
    os << "Ensemble dimensions:\n";
    os << "    - Temporal: size    = " << c.ensemble.T.size() << "\n";
    os << "                entries = " << c.ensemble.T << "\n";
    os << "    - ens:      sizes   = " << c.ensemble.ens.size() << "x"
       << c.ensemble.ens[0].size() << "\n";
    os << "    - avg:      size    = " << c.ensemble.avg.size() << "\n";
    if ( c.ensemble.avg.size() > 0 )
      {
	os << "    - avg:      entries = " << c.ensemble.avg;

#warning "Delete the below commented stuff once errors of SR data are comfirmed correct"	
	// for ( int i = 0; i < c.cov.dat["real"]->size1; ++i )
	//   {
	//     for ( int j = 0; j < c.cov.dat["real"]->size2; ++j )
	//       {
	// 	os << "    - err:      entries = ";
	// 	if ( i == j )
	// 	  os << "(" << sqrt(gsl_matrix_get(c.cov.dat["real"],i,j))
	// 	     << ", " << sqrt(gsl_matrix_get(c.cov.dat["imag"],i,j)) << ") ";
	//       }
	//   }
      }
    return os;
  }

  dat_t operator*(dat_t& d, std::complex<double> c)
  {
    for ( auto vi = d.ens.begin(); vi != d.ens.end(); ++vi )
      (*vi) *= c;
    d.avg *= c;

    return d;
  }

  dat_t operator*=(dat_t& d, std::complex<double> c)
  {
    for ( auto vi = d.ens.begin(); vi != d.ens.end(); ++vi )
      (*vi) *= c;
    d.avg *= c;

    return d;
  }

  // Direct add of ensembles
  dat_t operator+=(dat_t& d1, dat_t& d2)
  {
    // Check first dimensions of dat_t.ens are equal
    if ( d1.ens.size() != d2.ens.size() )
      {
	std::cout << "Ensembles not of same size!" << std::endl;
	exit(2);
      }
    
    for ( int n = 0; n < d1.ens.size(); ++n )
      d1.ens[n] += d2.ens[n];
    d1.avg += d2.avg;

    return d1;
  }

  // Direct subtract of ensembles
  dat_t operator-=(dat_t& d1, dat_t& d2)
  {
    // Check first dimensions of dat_t.ens are equal
    if ( d1.ens.size() != d2.ens.size() )
      {
	std::cout << "Ensembles not of same size!" << std::endl;
	exit(2);
      }
    
    for ( int n = 0; n < d1.ens.size(); ++n )
      d1.ens[n] -= d2.ens[n];
    d1.avg -= d2.avg;

    return d1;
  }

      

  ////////////////////////////////////////////////////////////////////////////////////////////

  void Avg(std::vector<std::complex<double> > &a, VVC &e, int _t, int _g)
  {
    a.resize(_t);
    for ( int t = 0; t < _t; ++t )
      {
	std::complex<double> _c(0.0,0.0);
	for ( int g = 0; g < _g; ++g )
	  {
	    _c += e[g][t];
	  }
	a[t].real(_c.real()/(1.0*_g));
	a[t].imag(_c.imag()/(1.0*_g));
      }
  }

  // Ensemble average
  void correlator::ensAvg()
  {
    Avg(ensemble.avg, ensemble.ens, getNt(), getCfgs() );
  }

  // Jackknife the data
  void correlator::jackknife()
  {
    jack.resize( getCfgs() );
    
    // Loop to exclude gth gauge config
    for ( int g = 0; g < getCfgs() ; ++g )
      {
	jack[g].T = ensemble.T; // same time domain
	jack[g].ens.resize( getCfgs()-1 );
	// Fill entries of jackknifed array with index lower than 'g'
	if ( g > 0 )
	  {
	    for ( int j = 0; j < g; ++j )
	      {
		for ( int t = 0; t < getNt(); ++t )
		  jack[g].ens[j].push_back(ensemble.ens[j][t]);
	      } // end j
	  }
	
	// Fill entries of jackknifed array starting at index 'g'
	for ( int j = g; j < getCfgs()-1; ++j )
	  {
	    for ( int t = 0; t < getNt(); ++t )
	      jack[g].ens[j].push_back(ensemble.ens[j+1][t]);
	  } // j

	// Form ensemble avg for this jackknife sample
	Avg(jack[g].avg, jack[g].ens, getNt(), jack[g].ens.size());
      } // g
  } // jackknife



  // Compute data covariance & its inverse
  void correlator::Cov()
  {
    // Convenience
    int Nt = getNt(); int cfgs = getCfgs();

    cov.dat["real"] = gsl_matrix_calloc(Nt,Nt); cov.inv["real"] = gsl_matrix_calloc(Nt,Nt);
    cov.dat["imag"] = gsl_matrix_calloc(Nt,Nt); cov.inv["imag"] = gsl_matrix_calloc(Nt,Nt);

    for ( int ti = 0; ti < Nt; ++ti )
      {
	for ( int tj = 0; tj < Nt; ++tj )
	  {
	    double _r(0.0), _i(0.0);
	    for ( int g = 0; g < cfgs; ++g )
	      {
		_r += ( jack[g].avg[ti].real() - ensemble.avg[ti].real() )*
		  ( jack[g].avg[tj].real() - ensemble.avg[tj].real() );
		_i += ( jack[g].avg[ti].imag() - ensemble.avg[ti].imag() )*
		  ( jack[g].avg[tj].imag() - ensemble.avg[tj].imag() );
	      } // g

	    gsl_matrix_set(cov.dat["real"], ti, tj, (( cfgs - 1 )/(1.0*cfgs))*_r);
	    gsl_matrix_set(cov.dat["imag"], ti, tj, (( cfgs - 1 )/(1.0*cfgs))*_i);

	  } // tj
      } // ti

    std::cout << "Attempt LinAlg::matrixInv call" << std::endl;
    // Now that data covariance is formed, compute the inverse
    cov.svs["real"] = LinAlg::matrixInv(cov.dat["real"],cov.inv["real"]);
    cov.svs["imag"] = LinAlg::matrixInv(cov.dat["imag"],cov.inv["imag"]);
#if 0
    std::cout << "C: "; LinAlg::printMat(cov.dat["real"]);
    std::cout << "IC: "; LinAlg::printMat(cov.inv["real"]);
    gsl_matrix * id = gsl_matrix_alloc(Nt,Nt);
    gsl_matrix_set_zero(id);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,cov.dat["real"],cov.inv["real"],0.0,id);
    std::cout << "ID?: "; LinAlg::printMat(id);
#endif
    
    std::cout << "Passed LinAlg::matrixInv call" << std::endl;

  } // Cov


  // Sum up temporal (e.g. summation method)
  void correlator::summation()
  {
    // Collapse VVC from N x T to and N x 1 vector
    // Maybe do this with a pointer, so VVC is left unaltered.
    for ( int g = 0; g < getCfgs(); ++g )
      {
	// Explicitly ignore t=0 contact term - add t\in[2,Nt-1] terms to t=1 term
	for ( int t = 2; t < getNt(); ++t )
	  ensemble.ens[g][1] += ensemble.ens[g][t];
	// Replace t=0 value with summed value
	ensemble.ens[g][0] = ensemble.ens[g][1];
	ensemble.ens[g].erase(ensemble.ens[g].begin()+1,ensemble.ens[g].end());
      }
    // Temporal range is now only a single element w/ value of correlator.Nt before summation
    ensemble.T.resize(1);
    ensemble.T[0] = getNt(); aspects.Nt = 1;
  }

  // Correct for bias in forming ratios of Npt/2pt funcs (per jk ensemble avg)
  void correlator::removeBias()
  {
    for ( auto ee = ensemble.ens.begin(); ee != ensemble.ens.end(); ++ee )
      {
	for ( auto tt = ee->begin(); tt != ee->end(); ++tt )
	  {
	    int tdx = std::distance(ee->begin(), tt);
	    // std::cout << "B4: " << *tt << std::endl;
	    // std::cout << "   ens avg: " << ensemble.avg[tdx] << std::endl;
	    tt->real(  ( getCfgs()*ensemble.avg[tdx].real() - ( getCfgs() - 1 )*tt->real()  )   );
	    tt->imag(  ( getCfgs()*ensemble.avg[tdx].imag() - ( getCfgs() - 1 )*tt->imag()  )   );
	    // std::cout << "After: " << *tt << std::endl;
	  }
      }
  }

  /*
    Other methods
  */
  /*
  Read a correlator from file and store in an 'ensemble' structure
  */
  // Read the correlator from file
  void read( std::ifstream& input, std::string& file, dat_t& _d )
  {
    // Read the file into memory
    try {
      std::cout << "Attempting to open file" << std::endl;
      input.open(file);
    }
    catch ( std::ifstream::failure e ) {
      std::cerr << "Caught exception in read..." << std::endl;
      exit(999);
    }

    std::cout << "Successfully read file..." << std::endl;

    // variables for reading in data
    std::string line;
    int t;
    double _r, _i;
    // Local copies
    int _g = _d.ens.size();
    int _n = _d.ens[0].size();


    // variables to track line location
    int cntr = 0; // counter to track which line in file is currently under consideration
    int config_num = 0; // counter of which gauge config is being looked at
    


    std::cout << "Reading correlator = " << file << std::endl;
    // Extract the correlator from file
    if ( input.is_open() )
      {
	while ( getline(input,line) && config_num < _g )
	  {
	    // Read in data line-by-line from file
	    input >> t >> _r >> _i;
	    _d.ens[config_num][t] = std::complex<double>(_r,_i);
	    // _d.real[config_num][t] = _r;
	    // _d.imag[config_num][t] = _i;

	    cntr++; // Iterate the counter tracking line location in file

	    if ( cntr%_n == 0 )
	      {
		config_num++; // increase the gauge_config tally by one
	      }
	  }
	input.close();
      }

  } // read


  correlator mergeCorrs(std::vector<VVC>& v)
  {
    int n(v.size());
    int g(v[0].size());
    int t(v[0][0].size());
    
    dat_t datMerge;
    datMerge.ens.resize(g);
    
    for ( int gg = 0; gg < g; ++gg )
      {
	datMerge.ens[gg].resize(t); // , std::complex<double> (0.0,0.0))
	for ( int tt = 0; tt < t; ++tt )
	  {
	    std::complex<double> dum(0.0,0.0);
	    for ( int nn = 0; nn < n; ++nn )
	      datMerge.ens[gg][tt] += v[nn][gg][tt];
	    datMerge.ens[gg][tt] /= n;
	  }
      }

#warning "Assuming each 3pt function has temporal range [0,t-1] in steps of one!"
    Pseudo::domain_t domMerge(0,1,t-1);
    Pseudo::prop_t dum(g, domMerge);
    correlator merge(dum);
    merge.ensemble.ens = datMerge.ens;
    return merge;
  } // mergeCorrs

  void conj(VVC* c)
  {
    for ( VVC::iterator g = c->begin(); g != c->end(); ++g )
      {
	for ( std::vector<std::complex<double> >::iterator it = g->begin(); it != g->end(); ++it )
	  (*it).imag( (*it).imag() * -1 );
      }
  }


  /*
    H5 Support
  */
  void H5Read(char *inH5, correlator *c)
  {
  }

  
  void writeCorr(correlator *c)
  {
    Exception::dontPrint(); // silence auto-printing of failures - try/catch to manage them
    
    /*
      Set up
    */
    const std::string& outH5 = "corr"+std::to_string(c->npt())+"pt-Data.h5";
    H5File h5;
    
    PredType DTYPE(PredType::IEEE_F64LE);
    std::string DATASET = "data";

    hsize_t dims[] = {c->getCfgs(), c->getNt()};
    DataSpace * Space = new DataSpace(2, dims, NULL);

    //---------------------------------------------------------------------------------
    // Groups to make/open
    std::vector<std::string> grps;
    std::string ops = c->getSnk().first+"-"+c->getSrc().first;
    std::string rows = "rows_"+std::to_string(c->getSnk().second)+std::to_string(c->getSrc().second);
    switch(c->npt())
      {
      case 2:
	grps.push_back("p"+Pseudo::shortMom(c->getPi(),""));
	grps.push_back("t0_avg");
	grps.push_back(ops);
	grps.push_back(rows);
	break;
      case 3:
	grps.push_back("pf" + Pseudo::shortMom(c->getPf(),"") + "_pi"
		       + Pseudo::shortMom(c->getPi(),""));
	grps.push_back("zsep" + std::to_string(Pseudo::shortZ(c->getDisp())[0])
		       + std::to_string(Pseudo::shortZ(c->getDisp())[1])
		       + std::to_string(Pseudo::shortZ(c->getDisp())[2]));
	grps.push_back("gamma-" + std::to_string(c->getGamma()));
	grps.push_back("t0_avg");
	grps.push_back(ops);
	grps.push_back(rows);
	break;
      }
    //----------------------------------------------------------------------------------


    try {
      h5.openFile(outH5.c_str(), H5F_ACC_RDWR);
      std::cout << "     Opened h5: " << outH5 << std::endl;
    } catch (H5::FileIException &file_dne) {
      h5 = H5File(outH5.c_str(), H5F_ACC_TRUNC);
      std::cout << "     Created h5: " << outH5 << std::endl;
    }


    // Try to open a first group
    Group root;
    try {
      root = h5.openGroup(&grps[0][0]);
      std::cout << "Root group exists...opening..." << std::endl;
    } catch(ReferenceException &groupDNE) {
      root = h5.createGroup(&grps[0][0]);
      std::cout << "Created new root dir..." << std::endl;
    } catch (FileIException &groupDNE) {
      root = h5.createGroup(&grps[0][0]);
      std::cout << "Created new root dir..." << std::endl;
    }

    auto x = grps.begin(); x++;
    while ( x != grps.end() )
      {
	// Create the remaining subgroups
	try {
	  root = root.openGroup(&(*x)[0]);
	  std::cout << "Opened group = " << *x << std::endl;
	} catch(GroupIException &groupDNE) {
	  root = root.createGroup(&(*x)[0]);
	  std::cout << "Created group = " << *x << std::endl;
	}
	++x;
      }


    // Components
    const char * const comps[] = {"real", "imag"};


    double realDATA[c->getCfgs()][c->getNt()];
    double imagDATA[c->getCfgs()][c->getNt()];
    
    for ( auto j = c->ensemble.ens.begin(); j != c->ensemble.ens.end(); ++j )
      {
	int jdx = std::distance(c->ensemble.ens.begin(), j);
	for ( auto ti = j->begin(); ti != j->end(); ++ti )
	  {
	    int tdx = std::distance(j->begin(), ti);
	    realDATA[jdx][tdx] = ti->real();
	    imagDATA[jdx][tdx] = ti->imag();
	  } // ti
      } // j



    //--------------------------
    // Onto the writing
    DataSet datR, datI;
    Group real, imag;
    try {
      real = root.openGroup(comps[0]);
      imag = root.openGroup(comps[1]);
    } catch (GroupIException &groupDNE) {
      real = root.createGroup(comps[0]);
      imag = root.createGroup(comps[1]);
    }

    try {
      datR = real.openDataSet(&DATASET[0]);
      datI = imag.openDataSet(&DATASET[0]);
      std::cout << "Open existing dataset" << std::endl;
    } catch(GroupIException &dsetDNE) {
      datR = real.createDataSet(&DATASET[0], DTYPE, *Space);
      datI = imag.createDataSet(&DATASET[0], DTYPE, *Space);
      std::cout << "Create new dataset" << std::endl;
    }

    try {
      datR.write(realDATA, DTYPE);
      datI.write(imagDATA, DTYPE);
    } catch (DataSetIException &d) {
      std::cout << "Failed to write dataset" << std::endl;
      exit(1);
    }
    //---------------------------------------------------------


    root.close();
    h5.close();
  }

  void fitResW(correlator *c, std::string& comp) // const char * comp)
  {
    Exception::dontPrint(); // silence auto-printing of failures - try/catch to manage them
    
    /*
      Set up
    */
    const std::string& outH5 = "corr"+std::to_string(c->npt())+"pt-FitRes.h5";
    H5File h5;
    
    // Other groups to make ----------------------------------------------------------------
    // const char *cmzg[] = {}; // [c->npt()dum[];
    std::vector<std::string> cmzg(2);
    cmzg[0] = comp;
    switch(c->npt())
      {
      case 2:
	cmzg[1] = "p"+Pseudo::shortMom(c->getPi(),"");
	break;
      case 3:
	std::string pfpi = "pf" + Pseudo::shortMom(c->getPf(),"") + "_pi" +
	  Pseudo::shortMom(c->getPi(),"");
	std::string z    = "zsep" + std::to_string(Pseudo::shortZ(c->getDisp())[0])
	  + std::to_string(Pseudo::shortZ(c->getDisp())[1])
	  + std::to_string(Pseudo::shortZ(c->getDisp())[2]);
	std::string gamma = "gamma-" + std::to_string(c->getGamma());

	cmzg[1] = pfpi; cmzg.push_back(z); cmzg.push_back(gamma);
	break;
      } // switch
    // -------------------------------------------------------------------------------------


    try {
      h5.openFile(outH5.c_str(), H5F_ACC_RDWR);
      std::cout << "     Opened h5: " << outH5 << std::endl;
    } catch (H5::FileIException &file_dne) {
      h5 = H5File(outH5.c_str(), H5F_ACC_TRUNC);
      std::cout << "     Created h5: " << outH5 << std::endl;
    }


    // Define the datatypes that will be written
    PredType DTYPE(PredType::IEEE_F64LE);

    // Some helpers
    const char * const collections[] = {"mean", "bins"};

    std::string fit = "/" + c->fit.theFit.verbose();
    std::string DATASET = "tfit_" + std::to_string(c->fit.theFit.range.min) + "-" +
      std::to_string(c->fit.theFit.range.max);


    // Make root group based on fit function
    Group root;
    try {
      // root.openGroup(&fit[0]);
      root = h5.openGroup(&fit[0]);
      std::cout << "Opened existing group" << std::endl;
    } catch (ReferenceException &groupDNE) {
      root = h5.createGroup(&fit[0]);
      std::cout << "Created new group" << std::endl;
    } catch (FileIException &groupDNE) {
      root = h5.createGroup(&fit[0]);
      std::cout << "Created new group" << std::endl;
    }


    // Save an extra set of write loops by making local params map that includes chi2
    std::map<std::string, std::vector<double> > _p(c->res.params);
    std::pair<std::string, std::vector<double> > chi = std::make_pair("chi2", c->res.chi2);
    _p.insert(chi);


    // Loop over named parameters
    for ( auto p = _p.begin(); p != _p.end(); ++p )
      {
	Group param, means, bins;

	try {
	  param = root.openGroup(&p->first[0]);
	  means = param.openGroup(collections[0]);
	  bins  = param.openGroup(collections[1]);
	  std::cout << "Opened existing param/means/bins/" << std::endl;
	} catch (GroupIException &groupDNE) {
	  param = root.createGroup(&p->first[0]);
	  means = param.createGroup(collections[0]);
	  bins  = param.createGroup(collections[1]);
	} catch (H5::ReferenceException &r) {
	  std::cout << "Caught reference exception:";
	}


	Group subMeans = means;
	Group subBins  = bins;
	// for ( int x = 0; x < sizeof(cmzg)/sizeof(cmzg[0]); ++x )
	for ( auto x = cmzg.begin(); x != cmzg.end(); ++x )
	  {
	    // Create more subgroups
	    try {
	      subMeans = subMeans.openGroup(&(*x)[0]); // cmzg[x]);
	      subBins  = subBins.openGroup(&(*x)[0]); // cmzg[x]);
	      std::cout << "Opened existing subMeans/Bins" << std::endl;
	    } catch (GroupIException &groupDNE) {
	      subMeans = subMeans.createGroup(&(*x)[0]); // cmzg[x]);
	      subBins  = subBins.createGroup(&(*x)[0]); // cmzg[x]);
	      std::cout << "Created new subMeans/Bins" << std::endl;
	    }
	  }

	// Spaces for mean/bins datasets
	hsize_t binDim[1]  = {c->getCfgs()};
	hsize_t meanDim[1] = {2};
	
	// hid_t
	DataSpace * binSpace  = new DataSpace(1,binDim,NULL);
	DataSpace * meanSpace = new DataSpace(1,meanDim,NULL);

	
	// Pull out and organize data into arrays
	double MEAN[2] = {0.0, 0.0};
	double BINS[c->getCfgs()];
	
	for ( int j = 0; j < c->getCfgs(); ++j )
	  {
	    BINS[j] = p->second[j];
	    MEAN[0] += ( p->second[j] / (1.0*c->getCfgs()) );
	  }

	for ( int j = 0; j < c->getCfgs(); ++j )
	  MEAN[1] += pow( p->second[j] - MEAN[0], 2);
	MEAN[1] = sqrt( (( c->getCfgs() - 1 )/(1.0*c->getCfgs())) * MEAN[1] );
	

	DataSet bin, mean;
	try {
	  bin = subBins.openDataSet(&DATASET[0]);
	  mean = subMeans.openDataSet(&DATASET[0]);
	  std::cout << "Open existing dataset" << std::endl;
	} catch(GroupIException &dsetDNE) {
	  bin = subBins.createDataSet(&DATASET[0], DTYPE, *binSpace);
	  mean = subMeans.createDataSet(&DATASET[0], DTYPE, *meanSpace);
	  std::cout << "Create new dataset" << std::endl;
	}

	try {
	  bin.write(BINS, DTYPE);
	  mean.write(MEAN, DTYPE);
	} catch (DataSetIException &d) {
	  std::cout << "Failed to write dataset" << std::endl;
	  exit(1);
	}

	delete binSpace;
	delete meanSpace;

	bins.close();
	means.close();
	param.close();
      } // param map
    
  root.close();
  h5.close();
  } // fitResW

}
