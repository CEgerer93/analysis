/*
  Utilities to help manage three pt function traces
*/
#include "threept_tr.h"
#include "hdf5.h"
#include "H5Cpp.h"
#include "cov_utils.h"

// #include "pseudo_utils.h"
// #include "operators.h"

using namespace H5;
using namespace Pseudo;


#define gc_rect(r,i) gsl_complex_rect(r,i)

/*
  Common
*/
gc zero = gsl_complex_rect(0.0,0.0);
gc one  = gsl_complex_rect(1.0,0.0);
gc mone = gsl_complex_rect(-1.0,0.0);
gc I    = gsl_complex_rect(0.0,1.0);
gc mI   = gsl_complex_rect(0.0,-1.0);

//! Zero out bits of a complex number - n.b. not in an adat header
std::complex<double> zeroComplex(const std::complex<double>& w)
{
  const double fuzz = 1.0e-11;
  // Pretty print the complex number. Replace noise with 0.
  double op_r = real(w);
  double op_i = imag(w);

  if (fabs(op_r) < fuzz)
    op_r = 0;
  if (fabs(op_i) < fuzz)
    op_i = 0;

  return std::complex<double>(op_r,op_i);
}


/*
  pauli_t constructor
*/
pauli_t::pauli_t(int n)
{
  gc c;
  switch(n)
    {
    case 1:
      c = gc_rect(1.0,0);
      gsl_matrix_complex_set(m,0,1,c); gsl_matrix_complex_set(m,1,0,c); break;
    case 2:
      c = gsl_complex_rect(0,1.0);
      gsl_matrix_complex_set(m,1,0,c);
      c = gsl_complex_rect(0,-1.0);
      gsl_matrix_complex_set(m,0,1,c); break;
    case 3:
      c = gsl_complex_rect(1.0,0);
      gsl_matrix_complex_set(m,0,0,c);
      c = gsl_complex_rect(-1.0,0);
      gsl_matrix_complex_set(m,1,1,c); break;
    default:
      std::cout << "Assuming you won't need a pauli matrix - unsupported index passed" << std::endl;
    }
}

/*
  diracMat_t constructor
*/
diracMat_t::diracMat_t(int n, bool MINK)
{
  pauli_t s(n);
  switch(n)
    {
    case 1:
      for ( size_t i = 0; i < s.m->size1; ++i )
	{
	  for ( size_t j = 0; j < s.m->size2; ++j )
	    {
	      if ( MINK )
		{
		  gsl_matrix_complex_set(gamma,i,2+j,gsl_matrix_complex_get(s.m,i,j));
		  gsl_matrix_complex_set(gamma,i+2,j,
					 gsl_complex_mul(mone,gsl_matrix_complex_get(s.m,i,j)));
		}
	      else if ( !MINK )
		{
		  gsl_matrix_complex_set(gamma,i,2+j,
					 gsl_complex_mul(mI,gsl_matrix_complex_get(s.m,i,j)));
		  gsl_matrix_complex_set(gamma,i+2,j,
					 gsl_complex_mul(I,gsl_matrix_complex_get(s.m,i,j)));
		}
	    }
	}
      break;
    case 2:
      for ( size_t i = 0; i < s.m->size1; ++i )
	{
	  for ( size_t j = 0; j < s.m->size2; ++j )
	    {
	      if ( MINK )
		{
		  gsl_matrix_complex_set(gamma,i,2+j,gsl_matrix_complex_get(s.m,i,j));
		  gsl_matrix_complex_set(gamma,i+2,j,
					 gsl_complex_mul(mone,gsl_matrix_complex_get(s.m,i,j)));
		}
	      else if ( !MINK )
		{
		  gsl_matrix_complex_set(gamma,i,2+j,
					 gsl_complex_mul(mI,gsl_matrix_complex_get(s.m,i,j)));
		  gsl_matrix_complex_set(gamma,i+2,j,
					 gsl_complex_mul(I,gsl_matrix_complex_get(s.m,i,j)));
		}
	    }
	}
      break;
    case 3:
      for ( size_t i = 0; i < s.m->size1; ++i )
	{
	  for ( size_t j = 0; j < s.m->size2; ++j )
	    {
	      if ( MINK )
		{
		  gsl_matrix_complex_set(gamma,i,2+j,gsl_matrix_complex_get(s.m,i,j));
		  gsl_matrix_complex_set(gamma,i+2,j,
					 gsl_complex_mul(mone,gsl_matrix_complex_get(s.m,i,j)));
		}
	      else if ( !MINK )
		{
		  gsl_matrix_complex_set(gamma,i,2+j,
					 gsl_complex_mul(mI,gsl_matrix_complex_get(s.m,i,j)));
		  gsl_matrix_complex_set(gamma,i+2,j,
					 gsl_complex_mul(I,gsl_matrix_complex_get(s.m,i,j)));
		}
	    }
	}
      break;
    case 4:
      for ( size_t i = 0; i < 2; ++i )
	gsl_matrix_complex_set(gamma,i,i,one);
      for ( size_t i = 2; i < 4; ++i )
	gsl_matrix_complex_set(gamma,i,i,mone);
      break;
    case 5:
      for ( size_t i = 0; i < 4; ++i )
	{
	  if ( MINK )
	    gsl_matrix_complex_set(gamma,i,(i+2)%4,one);
	  else if ( !MINK )
	    gsl_matrix_complex_set(gamma,i,(i+2)%4,mone);
	}     
      break;
    default:
      std::cerr << n << " is not a valid Dirac matrix!" << std::endl;
      exit(2);
    }
}


//=============================================================================================
/*
  spinor_t methods
*/
template<typename T>
void spinor_t::build(XMLArray::Array<T>& p, double E, double m, int L)
{
  pauli_t sx(1), sy(2), sz(3);

  gmc * ID2d = gsl_matrix_complex_alloc(2,2); gsl_matrix_complex_set_identity(ID2d);
  gvc * basis1 = gsl_vector_complex_calloc(2); gsl_vector_complex_set_basis(basis1,0);
  gvc * basis2 = gsl_vector_complex_calloc(2); gsl_vector_complex_set_basis(basis2,1);

  // Split \vec{p} into std::vector of gsl_complexes
  std::vector<gsl_complex> pc(3);
  for ( auto a = pc.begin(); a != pc.end(); ++a )
    *a = gsl_complex_rect(((2*M_PI/L)*p[std::distance(pc.begin(),a)])/(E+m),0);

  //----------------------------------------------------------------------
  // Form \sigma\cdot\vec{p}
  gmc * sDotP = gsl_matrix_complex_calloc(2,2);
  // px \dot sigma_x OVER (E+m)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc[0],sx.m,ID2d,zero,sDotP);
  // PLUS py \dot sigma_y OVER (E+m)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc[1],sy.m,ID2d,one,sDotP);
  // PLUS pz \dot sigma_z OVER (E+m)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,pc[2],sz.m,ID2d,one,sDotP);
  //----------------------------------------------------------------------


  // Make the spinor for twoJz = +/-1
  for ( int jz = 1; jz > -2; jz-=2 )
    {
      gvc * components = gsl_vector_complex_calloc(4);
      /* gsl_vector_complex_set_zero(components); */
      gvc * tmp = gsl_vector_complex_calloc(2);

      switch(jz)
	{
	case 1:
	  // Mult 2-comp tmp by sDotP  --  this shall be lower 2 components
	  gsl_blas_zgemv(CblasNoTrans,one,sDotP,basis1,zero,tmp);
	  for ( int s = 0; s < 2; ++s )
	    {
	      gsl_vector_complex_set(components,s,gsl_vector_complex_get(basis1,s));
	      gsl_vector_complex_set(components,s+2,gsl_vector_complex_get(tmp,s));
	    }
	  break;
	case -1:
	  // Mult 2-comp tmp by sDotP  --  this shall be lower 2 components
	  gsl_blas_zgemv(CblasNoTrans,one,sDotP,basis2,zero,tmp);
	  for ( int s = 0; s < 2; ++s )
	    {
	      gsl_vector_complex_set(components,s,gsl_vector_complex_get(basis2,s));
	      gsl_vector_complex_set(components,s+2,gsl_vector_complex_get(tmp,s));
	    }
	  break;
	}
      gsl_vector_complex_free(tmp);
        
      // Rescale spinor components by \sqrt( (E+m)/(2m) )
      gsl_blas_zdscal(sqrt((E+m)/(2*m)),components);

      // Insert into map
      std::pair<int, gvc> jzSpinor = std::make_pair(jz,*components);
      twoJz.insert(jzSpinor);

      /* gsl_vector_complex_free(components); */
    }
  //---------------- Now we have spinor ---------------------------------
} // build
//==============================================================================================



//==============================================================================================
/*
  ugu_t methods
*/
std::complex<double> ugu_t::eval(gsl_vector_complex * left, gsl_vector_complex * right)
{
  std::complex<double> val(0.0,0.0);

  gc res        = gc_rect(0.0,0.0);               // to collect result from gsl
  gvc * matVec  = gsl_vector_complex_calloc(4);
  gmc * matProd = gsl_matrix_complex_calloc(4,4); // To hold \gamma^4 \times \gamma^\mu

  // Multiply \gamma^4\gamma^\mu
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,g4.gamma,d.gamma,zero,matProd);
  // Matrix vector product (gamma^\mu on right-spinor)
  gsl_blas_zgemv(CblasNoTrans,one,matProd,right,zero,matVec);
  // Inner product of left^\dagger & (gamma^\mu \times right)
  gsl_blas_zdotc(left,matVec,&res);

  // Push real/imag components of gsl result into std::complex val
  val.real(GSL_REAL(res)); val.imag(GSL_IMAG(res));
  
  gsl_vector_complex_free(matVec);
  gsl_matrix_complex_free(matProd);

  return val;
}
//==============================================================================================



//**********************************************************************************************
/*
  Spinor class methods
*/
void Spinor::buildSpinors()
{
  // Canonical momentum is |\vec{p}| in +z-direction
  XMLArray::Array<double> canonMom(3);
  canonMom[0]=0; canonMom[1]=0; canonMom[2]=absMom();

  // Start by building absolute & canonical spinors: u(\vec{p},s) u(|\vec{p}|\hat{z},s)
  absolute.build(mom,E,m,L);
  canon.build(canonMom,E,m,L);

#if 1
  LinAlg::printMat(wig);
  LinAlg::printMat(coeffS);

  // Build subduced spinor from canonical spinor
  gmc * prod = gsl_matrix_complex_calloc(subduce.irrep_dim,twoJ+1);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,coeffS,wig,zero,prod);
  LinAlg::printMat(prod);
  for ( int i = twoJ; i >= -twoJ; i -= subduce.irrep_dim )
    {

      // Placeholder for up/down vector views to be operated on
      gvc * tmpUp = gsl_vector_complex_calloc(2);
      gvc * fooUp = gsl_vector_complex_alloc(2);
      gvc * tmpDown = gsl_vector_complex_calloc(2);
      gvc * fooDown = gsl_vector_complex_alloc(2);
      
      // Access upper/lower components of canonical spinor
      gsl_vector_complex_view up = gsl_vector_complex_subvector(&canon.twoJz[i],0,2);
      gsl_vector_complex_view down = gsl_vector_complex_subvector(&canon.twoJz[i],2,2);


      // Hold up.vector & down.vector
      gsl_vector_complex_memcpy(tmpUp,&up.vector);
      gsl_vector_complex_memcpy(fooUp,tmpUp);
      gsl_vector_complex_memcpy(tmpDown,&down.vector);
      gsl_vector_complex_memcpy(fooDown,tmpDown);


      gsl_blas_zgemv(CblasNoTrans,one,prod,fooUp,zero,tmpUp);
      gsl_blas_zgemv(CblasNoTrans,one,prod,fooDown,zero,tmpDown);

      // // gsl_vector_complex_memcpy(dum,&up.vector);

      // // L-mult up.vector by Subduction/Wigner-D product
      // gsl_blas_zgemv(CblasNoTrans,one,prod,dum,zero,&up.vector);

      // // Hold down.vector
      // gsl_vector_complex_memcpy(dum,&down.vector);
      // // L-muly down.vector by Subduction/Wigner-D product
      // gsl_blas_zgemv(CblasNoTrans,one,prod,dum,zero,&down.vector);
      
#if 0
      std::cout << "\nup view check: size = " << up.vector.size << std::endl;
      for ( auto j = 0; j < up.vector.size; ++j )
	std::cout << gsl_vector_complex_get(&up.vector,j).dat[0] << " " << gsl_vector_complex_get(&up.vector,j).dat[1] << " ";
      std::cout << "\ndown view check" << std::endl;
      for ( auto j = 0; j < down.vector.size; ++j )
      	std::cout << gsl_vector_complex_get(&down.vector,j).dat[0] << " " << gsl_vector_complex_get(&down.vector,j).dat[1] << " ";
      std::cout << "\n";
      exit(1000);
#endif

#if 1
      gsl_vector_complex * foo = gsl_vector_complex_calloc(4);
      gsl_vector_complex_set(foo,0,gsl_vector_complex_get(tmpUp,0));
      gsl_vector_complex_set(foo,1,gsl_vector_complex_get(tmpUp,1));
      gsl_vector_complex_set(foo,2,gsl_vector_complex_get(tmpDown,0));
      gsl_vector_complex_set(foo,3,gsl_vector_complex_get(tmpDown,1));


      gsl_vector_complex * toIns = gsl_vector_complex_alloc(4);
      gsl_vector_complex_memcpy(toIns,foo);

      std::pair<int, gvc> aSubducedSpinor(i,*toIns);
      subduced.twoJz.insert(aSubducedSpinor); // [i] = *foo;
      // gsl_vector_complex_memcpy(&subduced.twoJz[i],foo);

      gsl_vector_complex_free(foo);
#endif
    }
#endif
} // buildSpinors

void Spinor::initSubduce(const std::string& opName)
{
  // Gather subduction info for this operator
  subduce = subduceInfo(opName,mom);
  // Set the ranks of subduction coefficient
  coeffS = gsl_matrix_complex_calloc(subduce.irrep_dim,twoJ+1);
  

  // Build the Wigner-D
  for ( int i = 0; i < subduce.irrep_dim; ++i )
    {
      int twoJz_i = pow((-1),i);
      for ( int j = 0; j < subduce.irrep_dim; ++j )
	{
	  int twoJz_j = pow((-1),j);
	  gsl_matrix_complex_set(wig,i,j,gc_rect(Hadron::Wigner_D(twoJ,twoJz_i,twoJz_j,rot.alpha,rot.beta,rot.gamma).real(),
						 Hadron::Wigner_D(twoJ,twoJz_i,twoJz_j,rot.alpha,rot.beta,rot.gamma).imag()));
	}
    }
  // for ( int i = twoJ; i >= -twoJ; i -= subduce.irrep_dim )
  //   {
  //     for ( int j = twoJ; j >= -twoJ; j -= subduce.irrep_dim )
  // 	gsl_matrix_complex_set(wig,i,j,gc_rect(Hadron::Wigner_D(twoJ,i,j,rot.alpha,
  // 								rot.beta,rot.gamma).real(),
  // 					       Hadron::Wigner_D(twoJ,i,j,rot.alpha,
  // 								rot.beta,rot.gamma).imag()));
  //   }

  // Build the subduction matrix
  for ( int row = 1; row <= subduce.irrep_dim; ++row )
    {
      for ( int h = 1; h <= twoJ+1; ++h )
	{
	  gsl_matrix_complex_set(coeffS,row-1,h-1,
				 gc_rect((*subduce.H).operator()(row,h).real(),
					 (*subduce.H).operator()(row,h).imag()));
	}
    }
}
// End of Spinor class methods
//************************************************************************************************

void writePrefactor(int npt, int mu, int cfgs, const XMLArray::Array<int> &mom, std::string strRoot,
		    std::vector<std::complex<double> > &P)
{
  Exception::dontPrint(); // silence auto-printing of failures - try/catch to manage them

  // Loop over the real/imaginary components of S^mu internally
  std::vector<std::string> comps = {"real", "imag"};

  for ( auto c = comps.begin(); c != comps.end(); ++c )
    {
      // Rip out desired component into a separate vector
      std::vector<double> v;
      for ( auto a = P.begin(); a != P.end(); ++a )
	{
	  if ( *c == "real" )
	    {
	      v.push_back(a->real());
	    }
	  else if ( *c == "imag" )
	    {
	      v.push_back(a->imag());
	    }
	} // accessing done

      /*
	Set up
      */
      const std::string& outH5 = "corr"+std::to_string(npt)+"pt-FitRes.h5";
      H5File h5;

      // Groups to make
      // std::string strRoot = "/polVec"; ///idx"; CMT
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
	  BINS[j] = v[j];
	  MEAN[0] += ( v[j] / (1.0*cfgs) );
	}

      for ( int j = 0; j < cfgs; ++j )
	MEAN[1] += pow( v[j] - MEAN[0], 2);
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
} // writePrefactor

void writeSpinorContract(int npt, int mu, int cfgs, const XMLArray::Array<int>& pf,
			 const XMLArray::Array<int>& pi, int twoJz_f, int twoJz_i,
			 std::vector<std::complex<double> >& C)
{
  Exception::dontPrint(); // silence auto-printing of failures - try/catch to manage them

  // Loop over the real/imaginary components of S^mu internally
  std::vector<std::string> comps = {"real", "imag"};

  for ( auto c = comps.begin(); c != comps.end(); ++c )
    {
      // Rip out desired component into a separate vector
      std::vector<double> uGu;
      for ( auto a = C.begin(); a != C.end(); ++a )
        {
          if ( *c == "real" )
            {
              uGu.push_back(a->real());
	    }
          else if ( *c == "imag" )
            {
              uGu.push_back(a->imag());
            }
        } // accessing done

      /*
        Set up
      */
      const std::string& outH5 = "corr"+std::to_string(npt)+"pt-FitRes.h5";
      H5File h5;

      // Groups to make
      std::string strRoot = "/ufbar*gamma_mu*ui";
      std::vector<std::string> idxMoms(2);
      idxMoms[0] = std::to_string(mu);
      idxMoms[1] = "pf" + Pseudo::shortMom(pf,"") + "_pi" + Pseudo::shortMom(pi,"");
      std::string DATASET = "twoJfz_" + std::to_string(twoJz_f)+"__twoJiz_"
	+std::to_string(twoJz_i);

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

      // Make root
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


      Group g = root;
      // Make remaining groups
      for ( auto it = idxMoms.begin(); it != idxMoms.end(); ++it )
	{
	  try {
	    g = g.openGroup(&(*it)[0]);
	    std::cout << "Opened existing idx/mom" << std::endl;
	  } catch (GroupIException &groupDNE) {
	    g = g.createGroup(&(*it)[0]);
	  } catch (H5::ReferenceException &r) {
	    std::cout << "Caught reference exception:";
	  }
	}

      Group means, bins;
      try {
        means = g.openGroup(collections[0]);
        bins  = g.openGroup(collections[1]);
	std::cout << "Opened existing means/bins" << std::endl;
      } catch ( GroupIException &groupDNE) {
        means = g.createGroup(collections[0]);
        bins  = g.createGroup(collections[1]);
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
          BINS[j] = uGu[j];
          MEAN[0] += ( uGu[j] / (1.0*cfgs) );
        }

      for ( int j = 0; j < cfgs; ++j )
        MEAN[1] += pow( uGu[j] - MEAN[0], 2);
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
      g.close();
      root.close();
      h5.close();
    } // comps iter
} // writeSpinorContract
