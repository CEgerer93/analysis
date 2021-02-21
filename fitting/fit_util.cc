/*
  Define classes/structs/methods needed to handle matrix element extraction
  From adat based correlators; fitting with gsl
*/
#include "fit_util.h"
#include "summation.h"

#include<iostream>
#include<gsl/gsl_permutation.h> // permutation header for matrix inversions
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h> // linear algebra

namespace FIT
{
  // Generic gsl_matrix viewer
  void printMat(gsl_matrix *g)
  {
    std::cout << "{";
    for ( size_t i = 0; i < g->size1; i++ ) {
      std::cout << "{";
      for ( size_t j = 0; j < g->size1; j++ ) {
	std::cout << gsl_matrix_get(g,i,j) << ",";
      }
      std::cout << "},\n";
    }    
  }


  /*
    PERFORM INVERSION OF PASSED MATRIX - RETURN # OF SVs REMOVED
  */
  int matrixInv(gsl_matrix * M, gsl_matrix * MInv)
  {
    size_t dataDim = M->size1;
    gsl_matrix * toInvert = gsl_matrix_alloc(dataDim,dataDim);
    gsl_matrix_memcpy(toInvert,M); // make a copy of M, as toInvert is modified below
    gsl_matrix * V = gsl_matrix_alloc(dataDim,dataDim);
    gsl_vector *S = gsl_vector_alloc(dataDim);
    gsl_vector *work = gsl_vector_alloc(dataDim);

    /*
      PERFORM THE SINGULAR VALUE DECOMPOSITION OF DATA COVARIANCE MATRIX (A)
      
      A = USV^T
          ---> A is an MxN matrix
    	  ---> S is the singular value matrix (diagonal w/ singular values along diag - descending)
    	  ---> V is returned in an untransposed form
    */
    gsl_linalg_SV_decomp(toInvert,V,S,work); // SVD decomp;  'toInvert' replaced w/ U on return

    // Define an svd cut
    double svdCut = 1e-16;
    // Initialize the inverse of the S diag
    gsl_vector *pseudoSInv = gsl_vector_alloc(dataDim);
    gsl_vector_set_all(pseudoSInv,0.0);

    // Vector of singular values that are larger than specified cut
    std::vector<double> aboveCutVals;
    
    std::cout << "The singular values above SVD Cut = " << svdCut << " are..." << std::endl;
    for ( int s = 0; s < dataDim; s++ )
      {
    	double dum = gsl_vector_get(S,s);
    	if ( dum >= svdCut ) { aboveCutVals.push_back(dum); }
      }
    
    // Assign the inverse of aboveCutVals to the pseudoSInv vector
    for ( std::vector<double>::iterator it = aboveCutVals.begin(); it != aboveCutVals.end(); ++it )
      {
    	gsl_vector_set(pseudoSInv,it-aboveCutVals.begin(),1.0/(*it));
      }

    // Promote this pseudoSInv vector to a matrix, where entries are placed along diagonal
    gsl_matrix * pseudoSInvMat = gsl_matrix_alloc(dataDim,dataDim);
    gsl_matrix_set_zero(pseudoSInvMat);
    for ( int m = 0; m < dataDim; m++ )
      {
    	gsl_matrix_set(pseudoSInvMat,m,m,gsl_vector_get(pseudoSInv,m));
    	std::cout << gsl_vector_get(pseudoSInv,m) << std::endl;
      }

    /*
      With singular values that are zero, best we can do is construct a pseudo-inverse

      In general, the inverse we are after is
                               VS^(-1)U^T
		     with S the diagonal matrix of singular values
    */
    // In place construct the transpose of U ('toInvert' was modified in place to U above)
    gsl_matrix_transpose(toInvert);

    gsl_matrix * SinvUT = gsl_matrix_alloc(dataDim,dataDim); gsl_matrix_set_zero(SinvUT);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,pseudoSInvMat,toInvert,0.0,SinvUT);

    // Now make the inverse of 'toInvert'
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,SinvUT,0.0,MInv);


    return pseudoSInv->size - aboveCutVals.size();
  }




  double chi2Linear(const gsl_vector * x, void *data)
  {
    
    // Get a pointer to the void data structure
    linFit_t * jfitCpy = (linFit_t *)data;
    // The current fit params
    FitRes_t fit( "LINEAR", gsl_vector_get(x,0), gsl_vector_get(x,1) );


    // Evaluate the linear fit for these parameters at each time
    std::vector<double> predict;
    for ( auto t = jfitCpy->T.begin(); t != jfitCpy->T.end(); ++t )
      {
	predict.push_back( fit.func(*t) );
      }


    // Begin chi2 computation
    double chi2(0.0);
    gsl_vector *iDiffVec = gsl_vector_alloc(jfitCpy->T.size());
    gsl_vector *jDiffVec = gsl_vector_alloc(jfitCpy->T.size());

    for ( auto l = jfitCpy->m.begin(); l != jfitCpy->m.end(); ++l )
      {
	int idx = std::distance(jfitCpy->m.begin(),l);
	gsl_vector_set(iDiffVec, idx, predict[idx] - *l );
      }
    // Copy the difference vector
    gsl_vector_memcpy(jDiffVec, iDiffVec);

    // Initialize cov^-1 right multiplying jDiffVec
    gsl_vector *invCovRightMult = gsl_vector_alloc(jfitCpy->m.size());
    gsl_blas_dgemv(CblasNoTrans,1.0,jfitCpy->covinv,jDiffVec,0.0,invCovRightMult);

    // Form the scalar dot product of iDiffVec & result of invCov x jDiffVec
    gsl_blas_ddot(iDiffVec,invCovRightMult,&chi2);

    // Free some memory
    gsl_vector_free(iDiffVec);
    gsl_vector_free(jDiffVec);
    gsl_vector_free(invCovRightMult);

    return chi2;

  } // chi2Linear
		       
    


  
// #if 0  
//   /*
//     READER FOR PASSED H5 FILES
//   */
//   void H5Read(char *inH5, reducedPITD *dat, int gauge_configs, int zmin, int zmax, int pmin,
// 	      int pmax, std::string dTypeName)
//   {
//     /*
//     OPEN THE H5 FILE CONTAINING ENSEM/JACK RESULTS OF pPDF
//     */
//     std::cout << "     READING H5 FILE = " << inH5 << std::endl;
//     hid_t h5File = H5Fopen(inH5,H5F_ACC_RDONLY,H5P_DEFAULT);
//     /*
//     Some h5 handles for accessing data
//     */
//     hid_t space, h5Pitd;
//     herr_t h5Status;
//     // The name of data actually stored in h5 file
//     const char * DATASET = &dTypeName[0];

//     // Access the first group entry within root - i.e. the current
//     hid_t h5Current = H5Gopen(h5File, "/b_b0xDA__J0_A1pP", H5P_DEFAULT);

//     // Other group headings w/in h5 file
//     const char * const momenta[] = {"pz0","pz1","pz2","pz3","pz4","pz5","pz6"};
//     const char * const comp[] = {"Re", "Im"}; // {"1","2"};
//     const char * const zsep[] = {"zsep0","zsep1","zsep2","zsep3","zsep4","zsep5","zsep6","zsep7","zsep8",
// 				 "zsep9","zsep10","zsep11","zsep12","zsep13","zsep14","zsep15","zsep16"};


//     // Iterator through stored displacements
//     for ( int z = zmin; z <= zmax; z++ )
//       {
// 	zvals dumZ;

// 	hid_t h5Zsep = H5Gopen(h5Current, zsep[z], H5P_DEFAULT);

// 	// Iterate through each momentum stored
// 	for ( int m = pmin; m <= pmax; m++ )
// 	  {
// 	    hid_t h5Mom = H5Gopen(h5Zsep, momenta[m], H5P_DEFAULT);


// 	    // // Start by grabbing handle to ensemble group
// 	    // hid_t h5Ensem = H5Gopen(h5Mom, "ensemble", H5P_DEFAULT);
	    
// 	    // for ( int c = 0; c < ReIm; c++ )
// 	    //   {
// 	    //     // Get the component handle
// 	    //     hid_t h5Comp = H5Gopen(h5Ensem,comp[c], H5P_DEFAULT);
	    
// 	    //     for ( int z = zmin; z <= zmax; z++ )
// 	    //       {
// 	    // 	// Get the zsep handle
// 	    // 	hid_t h5Zsep = H5Gopen(h5Comp, zsep[z], H5P_DEFAULT);
	    
// 	    // 	double ensemPitd[3];
// 	    // 	// Grab the dataset handle
// 	    // 	h5Pitd  = H5Dopen(h5Zsep, DATASET, H5P_DEFAULT);
// 	    // 	// Grab the data!
// 	    // 	h5Status = H5Dread(h5Pitd, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ensemPitd);
	    
	    
// 	    // 	if ( c == 0 )
// 	    // 	  {
// 	    // 	    dat->data.real.disps[z].ensem.IT[m-pmin]=ensemPitd[0];
// 	    // 	    ens->real.disps[z].ensem.avgM[m-pmin]=ensemPitd[1];
// 	    // 	    ens->real.disps[z].ensem.errM[m-pmin]=ensemPitd[2];
// 	    // 	  }
// 	    // 	if ( c == 1 )
// 	    // 	  {
// 	    // 	    ens->imag.disps[z].ensem.IT[m-pmin]=ensemPitd[0];
// 	    // 	    ens->imag.disps[z].ensem.avgM[m-pmin]=ensemPitd[1];
// 	    // 	    ens->imag.disps[z].ensem.errM[m-pmin]=ensemPitd[2];
// 	    // 	  }
// 	    //       } // end z
// 	    //     H5Gclose(h5Comp);
// 	    //   } // end c
// 	    // H5Gclose(h5Ensem); // Finished parsing the ensemble groups
	    
// 	    // Now grab a handle to the jack group
// 	    hid_t h5Jack = H5Gopen(h5Mom, "jack", H5P_DEFAULT);


// 	    // A container for this {p,z} data
// 	    momVals dumMomVal(gauge_configs);

	
// 	    // for ( int c = 1; c > -1; c-- )
// 	    for ( int c = 0; c < 2; c++ )
// 	      {
// 		// Get the component handle
// 		hid_t h5Comp = H5Gopen(h5Jack,comp[c], H5P_DEFAULT);
		

// 		std::cout << "/b_b0xDA__J0_A1pP/" << zsep[z] << "/" << momenta[m]
// 			  << "jack/" << comp[c] << "/" << dTypeName << std::endl;
		
// 		// // Get the zsep handle
// 		// hid_t h5Zsep = H5Gopen(h5Comp, zsep[z], H5P_DEFAULT);
		
// 		// Initialize a buffer to read jackknife dataset
// 		// double read_buf[gauge_configs*3];
// 		double read_buf[gauge_configs*2];
// 		/* double *read_buf = NULL; */
// 		/* read_buf = (double*) malloc(sizeof(double)*gauge_configs*3); */
		
		
// 		hid_t dset_id = H5Dopen(h5Comp, DATASET, H5P_DEFAULT);
// 		herr_t status = H5Dread(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, read_buf);
		
		
		
// 		// // A container for this {p,z} data
// 		// momVals dumMomVal(gauge_configs);
// 		dumMomVal.IT = read_buf[0];    // Same Ioffe-time
// 		double _avgR(0.0), _avgI(0.0); // Talley the averages
		
// 		// Fill the dumMomVal Ioffe-time data for each jackknife sample
// 		for ( int J = 0; J < gauge_configs; J++ )
// 		  {
// 		    if ( c == 0 )
// 		      {
// 			dumMomVal.mat[J].real(read_buf[1+2*J]);
// 			_avgR += read_buf[1+2*J];
// 		      }
// 		    if ( c == 1 )
// 		      {
// 			dumMomVal.mat[J].imag(read_buf[1+2*J]);
// 			_avgI += read_buf[1+2*J];
// 		      }
// 		  }
// 		// Compute the averge Mat for the {p, z}
// 		if ( c == 0 )
// 		  dumMomVal.matAvg.real( _avgR / gauge_configs );
// 		if ( c == 1 )
// 		  dumMomVal.matAvg.imag( _avgI / gauge_configs );


// 		H5Gclose(h5Comp);
// 	      } // end c jack

// 		// Associate a std::string with this collection of IT data
// 		std::pair<std::string, momVals> amom ( momenta[m], dumMomVal );
		
// 		dumZ.moms.insert(amom);
		
// 		// // Associate an int w/ zvals object
// 		// std::pair<int, zvals> az ( z, dumZ );
// 		// // Insert this zvals map into the pitd.disps map
// 		// dat->data.disps.insert( az );

// 		// std::cout << "XXX = " << _avgI / gauge_configs << std::endl;
// 		// std::cout << "XXX = " << dumZ.moms["pz1"].mat[0].real() << std::endl;

		
// 		// End of parsing jack data from read_buf


		
// 	      // 	H5Gclose(h5Comp);
// 	      // } // end c jack



// 	    H5Gclose(h5Jack);
// 	    H5Gclose(h5Mom);

// 	  } // end mom

// 	// Associate an int w/ zvals object
// 	std::pair<int, zvals> az ( z, dumZ );
// 	// Insert this zvals map into the pitd.disps map
// 	dat->data.disps.insert( az );
	
//       } // end z
//     H5Gclose(h5Current);
//     H5Fclose(h5File);

//     std::cout << "READ H5 FILE SUCCESS" << std::endl;
//   } // end H5Read

  
//   /*
//     WRITER FOR MAKING NEW H5 FILES - E.G. EVOLVED/MATCHED DATASETS
//   */
//   void H5Write(char *outH5, reducedPITD *dat, int gauge_configs, int zmin, int zmax, int pmin,
// 	       int pmax, std::string dTypeName)
//   {
//     /*
//       OPEN THE H5 FILE FOR WRITING
//     */
//     std::cout << "    WRITING H5 FILE = " << outH5 << std::endl;
//     hid_t h5File = H5Fcreate(outH5,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
//     /*
//       Some h5 handles for writing data
//     */
//     hid_t space, h5Pitd;
//     herr_t h5Status;
//     // The name of data actually stored in h5 file
//     const char * DATASET = &dTypeName[0];

//     // Make the first group entry within root - i.e. the current
//     hid_t h5Current = H5Gcreate(h5File, "/b_b0xDA__J0_A1pP", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
//     // Other group headings w/in h5 file
//     const char * const momenta[] = {"pz0","pz1","pz2","pz3","pz4","pz5","pz6"};
//     const char * const comp[] = {"Re", "Im"}; // {"1","2"};
//     const char * const zsep[] = {"zsep0","zsep1","zsep2","zsep3","zsep4","zsep5","zsep6","zsep7","zsep8",
//   				 "zsep9","zsep10","zsep11","zsep12","zsep13","zsep14","zsep15","zsep16"};
//     // const char * const zsep[] = {"0","1","2","3","4","5","6","7","8",
//     // 				 "9","10","11","12","13","14","15","16"};


//     // Iterator through displacements to store
//     for ( int z = zmin; z <= zmax; z++ )
//       {
// 	zvals dumZ;

// 	hid_t h5Zsep = H5Gcreate(h5Current, zsep[z], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

// 	// Iterate through each momentum to store
// 	for ( int m = pmin; m <= pmax; m++ )
// 	  {

// 	    hid_t h5Mom = H5Gcreate(h5Zsep, momenta[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// 	    // if ( H5Aexists(h5Current, momenta[m] ) )
// 	    //   {
// 	    // 	std::cout << "Group exists, opening" << std::endl;
// 	    // 	h5Mom = H5Gopen1(h5Current, momenta[m]);
// 	    //   }
// 	    // else
// 	    //   {
// 	    // 	std::cout << "Group doesn't exist, creating ---- > " << momenta[m] << std::endl;
// 	    // 	h5Mom = H5Gcreate(h5Current, momenta[m], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// 	    //   }
	

// 	    // Now grab a handle to the jack group
// 	    hid_t h5Jack = H5Gcreate(h5Mom, "jack", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// 	    // if ( H5Aexists(h5Mom, "jack") )
// 	    //   {
// 	    // 	std::cout << "jack exists, opening" << std::endl;
// 	    // 	h5Jack = H5Gopen1(h5Mom, "jack");
// 	    //   }
// 	    // else
// 	    //   {
// 	    // 	std::cout << "jack doesn't exist, creating" << std::endl;
// 	    // 	h5Jack = H5Gcreate(h5Mom, "jack", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// 	    //   }



// 	    // This set of {p,z} data
// 	    momVals dumMomVal = dat->data.disps[z].moms[momenta[m]];
// 	    std::cout << "Got to dumMomVal" << std::endl;

// 	    for ( int c = 0; c != 2; c++ )
// 	    // for ( int c = 1; c > -1; c-- )
// 	      {
// 		// Get the component handle
// 		hid_t h5Comp = H5Gcreate(h5Jack,comp[c], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// 		// if ( H5Aexists(h5Jack, comp[c]) )
// 		//   {
// 		//     std::cout << "Comp exists, opening" << std::endl;
// 		//     h5Comp = H5Gopen1(h5Jack, comp[c]);
// 		//   }
// 		// else
// 		//   {
// 		//     std::cout << "Comp doesn't exist, creating" << std::endl;
// 		//     H5Gcreate(h5Jack,comp[c], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// 		//   }


//   		// // Get the zsep handle
//   		// hid_t h5Zsep;
// 		// if ( H5Aexists(h5Comp, zsep[z]) )
// 		//   {
// 		//     std::cout << "zsep exists, opening" << std::endl;
// 		//     h5Zsep = H5Gopen1(h5Comp, zsep[z]);
// 		//   }
// 		// else
// 		//   {
// 		//     std::cout << "zsep doesn't exist, creating" << std::endl;
// 		//     H5Gcreate(h5Comp, zsep[z], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// 		//   }


//   		// Prepare buffer for writing
//   		double jackBuff[2*gauge_configs];
// 		// jackBuff[0] = dumMomVal.IT;
//   		for ( int J = 0; J < gauge_configs; J++ )
//   		  {
// 		    jackBuff[2*J] = dumMomVal.IT;
//   		    if ( c == 0 )
// 		      jackBuff[2*J+1]=dumMomVal.mat[J].real();
//   		    if ( c == 1 )
// 		      jackBuff[2*J+1]=dumMomVal.mat[J].imag();
//   		  } // end J

//   		// Grab the dataset handle
//   		hsize_t dims[2] = {gauge_configs,2};
// 		std::cout << "BEFORE DATASPACE SET" << std::endl;
//   		hid_t DATASPACE = H5Screate_simple(2,dims,NULL);
//   		hid_t dset_id = H5Dcreate(h5Comp, DATASET, H5T_IEEE_F64LE, DATASPACE,
//   					  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   		// Push data to file
//   		herr_t status = H5Dwrite(dset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, jackBuff);


// 		H5Gclose(h5Comp);
// 	      } // end c jack

// 	    H5Gclose(h5Jack);
// 	    H5Gclose(h5Mom);
// 	  } // end mom

// 	H5Gclose(h5Zsep);
//       } // end z
//     H5Gclose(h5Current);
//     H5Fclose(h5File);

//     std::cout << "WRITE H5 FILE SUCCESS" << std::endl;
//   } // end H5Write
// #endif
}

