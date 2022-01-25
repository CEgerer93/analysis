/*
  Utilities to help manage three pt function traces
*/
#include "threept_tr.h"
#include "hdf5.h"
#include "H5Cpp.h"

#include "pseudo_utils.h"

using namespace H5;

void writePolVec(int npt, int mu, int cfgs, const XMLArray::Array<int> &mom,
		 std::vector<std::complex<double> > &polVec)
{
  Exception::dontPrint(); // silence auto-printing of failures - try/catch to manage them

  // Loop over the real/imaginary components of S^mu internally
  std::vector<std::string> comps = {"real", "imag"};

  for ( auto c = comps.begin(); c != comps.end(); ++c )
    {
      // Rip out desired component into a separate vector
      std::vector<double> S;
      for ( auto a = polVec.begin(); a != polVec.end(); ++a )
	{
	  if ( *c == "real" )
	    {
	      S.push_back(a->real());
	  }
	  else if ( *c == "imag" )
	    {
	      S.push_back(a->imag());
	    }
	} // accessing done

      /*
	Set up
      */
      const std::string& outH5 = "corr"+std::to_string(npt)+"pt-FitRes.h5";
      H5File h5;

      // Groups to make
      std::string strRoot = "/polVec"; ///idx";
      std::string idx  = std::to_string(mu);
      std::string DATASET = "pf" + Pseudo::shortMom(mom,"") + "_pi" + Pseudo::shortMom(mom,"");

      // Define the datatypes that will be written
      PredType DTYPE(PredType::IEEE_F64LE);

      // Some helpers
      const char * const collections[] = {"mean", "bins"};


      // Open/create
      try {
	h5.openFile(outH5.c_str(), H5F_ACC_RDWR);
	std::cout << "     Opened h5: " << outH5 << std::endl;
      } catch (H5::FileIException &file_dne) {
	h5 = H5File(outH5.c_str(), H5F_ACC_TRUNC);
	std::cout << "     Created h5: " << outH5 << std::endl;
      }

      // Make groups
      Group root;
      try {
	root = h5.openGroup(&strRoot[0]);
	std::cout << "Opened existing group" << std::endl;
      } catch (ReferenceException &groupDNE) {
	root = h5.createGroup(&strRoot[0]);
	std::cout << "Created new group" << std::endl;
      } catch (FileIException &groupDNE) {
	root = h5.createGroup(&strRoot[0]);
	std::cout << "Created new group" << std::endl;
      }

      
      Group index;
      try {
	index = root.openGroup(&idx[0]);
	std::cout << "Opened existing idx" << std::endl;
      } catch (GroupIException &groupDNE) {
	index = root.createGroup(&idx[0]);
      } catch (H5::ReferenceException &r) {
	std::cout << "Caught reference exception:";
      }

      Group means, bins;
      try {
	means = index.openGroup(collections[0]);
	bins  = index.openGroup(collections[1]);
	std::cout << "Opened existing means/bins" << std::endl;
      } catch ( GroupIException &groupDNE) {
	means = index.createGroup(collections[0]);
	bins  = index.createGroup(collections[1]);
      } catch (H5::ReferenceException &r ) {
	std::cout << "Caught reference exception" << std::endl;
      }

      // Now make the groups for real/imag
      Group subMeans, subBins;
      try {
	subMeans = means.openGroup(&(*c)[0]); //c[0]);
	subBins  = bins.openGroup(&(*c)[0]);
	std::cout << "Opened existing subMeans/Bins" << std::endl;
      } catch (GroupIException &groupDNE) {
	subMeans = means.createGroup(&(*c)[0]);
	subBins  = bins.createGroup(&(*c)[0]);
	std::cout << "Created new subMeans/Bins" << std::endl;
      }
      
  
      // Spaces for mean/bins datasets
      hsize_t binDim[1] = {cfgs};
      hsize_t meanDim[1] = {2};

      // hid_t
      DataSpace * binSpace = new DataSpace(1,binDim,NULL);
      DataSpace * meanSpace = new DataSpace(1,meanDim,NULL);

      // Pull out and organize data into arrays
      double MEAN[2] = {0.0, 0.0};
      double BINS[cfgs];

      for ( int j = 0; j < cfgs; ++j )
	{
	  BINS[j] = S[j];
	  MEAN[0] += ( S[j] / (1.0*cfgs) );
	}

      for ( int j = 0; j < cfgs; ++j )
	MEAN[1] += pow( S[j] - MEAN[0], 2);
      MEAN[1] = sqrt( (( cfgs - 1 )/(1.0*cfgs)) * MEAN[1] );


      // Attempt to make datasets
      DataSet bin, mean;
      try {
	bin  = subBins.openDataSet(&DATASET[0]);
	mean = subMeans.openDataSet(&DATASET[0]);
	std::cout << "Open existing dataset" << std::endl;
      } catch (GroupIException &dsetDNE) {
	bin  = subBins.createDataSet(&DATASET[0], DTYPE, *binSpace);
	mean = subMeans.createDataSet(&DATASET[0], DTYPE, *meanSpace);
	std::cout << "Create new dataset" << std::endl;
      }
      
      // Write the data
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
      index.close();
      root.close();
      h5.close();
    } // comps iter
} // writePolVec
