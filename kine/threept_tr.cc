/*
  Utilities to help manage three pt function traces
*/
#include "threept_tr.h"
#include "hdf5.h"
#include "H5Cpp.h"
#include "cov_utils.h"

// #include "shortcuts_gsl.h"

using namespace H5;
using namespace Pseudo;
using namespace LinAlg;

// Print elements of real-valued gsl_vector
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

// Print elements of complex-valued gsl_vector
std::ostream& operator<<(std::ostream& os, const gsl_vector_complex *v)
{
  if ( v->size > 0 )
    {
      for ( int i = 0; i < v->size; ++i )
        {
          for ( int c = 0; c < 2; ++c )
            os << gsl_vector_complex_get(v,i).dat[c] << " ";
          if ( i != v->size-1 )
            os << ", ";
        }
    }
  return os;
}


// Print elements of real-valued gsl_matrix
std::ostream& operator<<(std::ostream& os, const gsl_matrix * m)
{
  int r, c;
  for ( r = 0; r < m->size1; ++r )
    {
      os << "\n";
      for ( c = 0; c < m->size2; ++c )
	os << gsl_matrix_get(m,r,c) << ", ";
    }
  return os;
}

// Print elements of complex-valued gsl_matrix
std::ostream& operator<<(std::ostream& os, const gsl_matrix_complex * m)
{
  int r, c;
  for ( r = 0; r < m->size1; ++r )
    {
      os << "\n";
      for ( c = 0; c < m->size2; ++c )
	os << "[" << GSL_REAL(gsl_matrix_complex_get(m,r,c)) << ", "
	   << GSL_IMAG(gsl_matrix_complex_get(m,r,c)) << "] ";
    }
	
  return os;
}

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

  // Reset any "-0"'s to zero
  for ( int row = 0; row < gamma->size1; ++row )
    {
      for ( int col = 0; col < gamma->size2; ++col )
	{
	  gc holder = gsl_matrix_complex_get(gamma,row,col);
	  if ( holder.dat[0] == -0 )
	    holder = gsl_complex_rect( 0, holder.dat[1] );
	  if ( holder.dat[1] == -0 )
	    holder = gsl_complex_rect( holder.dat[0], 0 );

	  gsl_matrix_complex_set(gamma,row,col,holder);
	}
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
    *a = gsl_complex_rect(((2*PI/L)*p[std::distance(pc.begin(),a)])/(E+m),0);

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
#if 0
      // Rescale spinor components by \sqrt( (E+m)/(2m) )
      gsl_blas_zdscal(sqrt((E+m)/(2*m)),components);
#else
      // Rescale spinor components by \sqrt( (E+m) )
      gsl_blas_zdscal(sqrt(E+m),components);
#endif

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
  polVec_t methods
*/
std::complex<double> polVec_t::eval(gsl_vector_complex * left, gsl_vector_complex * right, double mass)
// std::complex<double> polVec_t::eval(gsl_vector_complex * left, gsl_vector_complex * right)
{
  std::complex<double> val(0.0,0.0);
  
  gc res = gc_rect(0.0,0.0); // to collect result from gsl
  gvc * matVec = gsl_vector_complex_calloc(4);
  gmc * matProd1 = gsl_matrix_complex_calloc(4,4);
  gmc * matProd2 = gsl_matrix_complex_calloc(4,4);

  // Multiply \gamma^4\gamma^\mu
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,g4.gamma,d.gamma,zero,matProd1);
  // gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,g4.gamma,g5.gamma,zero,matProd1);
  // Right multiply by \gamma^5
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,matProd1,g5.gamma,zero,matProd2);
  // gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,matProd1,d.gamma,zero,matProd2);
  // Matrix vector product (gamma's right mult'd by right spinor)
  gsl_blas_zgemv(CblasNoTrans,one,matProd2,right,zero,matVec);
  // Inner product of left^\dagger & (gammas * right product)
  gsl_blas_zdotc(left,matVec,&res);

  // Push real/imag components of gsl result into std::complex val
  val.real(GSL_REAL(res)); val.imag(GSL_IMAG(res));

#if 1
// #warning "This won't compile because 'm' is no longer passed"
  // Normalize val by 1/(2m)
  val *= (1.0/(2*mass));
#endif

  gsl_vector_complex_free(matVec);
  gsl_matrix_complex_free(matProd1);
  gsl_matrix_complex_free(matProd2);

  return val;
}
//==============================================================================================

//==============================================================================================
/*
  u1u_t methods
*/
std::complex<double> u1u_t::eval(gsl_vector_complex * left, gsl_vector_complex * right)
{
  std::complex<double> val(0.0,0.0);

  gc res        = gc_rect(0.0,0.0);               // to collect result from gsl
  gvc * matVec  = gsl_vector_complex_calloc(4);

  // Matrix vector product (gamma^4 on right-spinor)
  gsl_blas_zgemv(CblasNoTrans,one,g4.gamma,right,zero,matVec);
  // Inner product of left^\dagger & (gamma^4 \times right)
  gsl_blas_zdotc(left,matVec,&res);

  // Push real/imag components of gsl result into std::complex val
  val.real(GSL_REAL(res)); val.imag(GSL_IMAG(res));

  gsl_vector_complex_free(matVec);

  return val;
}
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


//==============================================================================================
/*
  utu_t methods
*/
std::complex<double> utu_t::eval(gsl_vector_complex * left, gsl_vector_complex * right)
{
  std::complex<double> val(0.0,0.0);

  gc res = gc_rect(0.0,0.0); // to collect result from gsl
  gc two = gc_rect(2.0,0.0); // for rescaling
  gc Ihalf = gc_rect(0.0,0.5);

  gvc * matVec = gsl_vector_complex_calloc(4);
  gmc * matMuNu = gsl_matrix_complex_calloc(4,4); // To hold \gamma^\mu x \gamma^\nu
  gmc * matNuMu = gsl_matrix_complex_calloc(4,4); // To hold \gamma^\nu x \gamma^\nu
  // gmc * matProd4Mu = gsl_matrix_complex_calloc(4,4); // To hold \gamma^4 \times \gamma^mu
  gmc * matProd = gsl_matrix_complex_calloc(4,4); // To hold (\gamma^4\gamma^\mu) \times\gamma^\nu


  // Multiply \gamma^\mu\gamma^\nu & \gamma^\nu\gamma^\mu
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,dl.gamma,dr.gamma,zero,matMuNu);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,dr.gamma,dl.gamma,zero,matNuMu);

#if 0
  LinAlg::printMat(matMuNu);
  LinAlg::printMat(matNuMu);
#endif
  // Subtract \gamma^\nu\gamma^\mu from \gamma^\mu\gamma^\nu ( result stored in matMuNu )
  gsl_matrix_complex_sub(matMuNu,matNuMu);

#if 0
  std::cout << "Commutator" << std::endl;
  LinAlg::printMat(matMuNu);
#endif
  
  // Left multiply resulting commutator of gamma matrices by \gamma^4 from \bar{u}
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,g4.gamma,matMuNu,zero,matProd);

  // // sigma^{\mu\nu} = (I/2) * [ \gamma^\mu, \gamma^\nu ]
  // //                = I * [ \gamma_\mu \gamma_\nu - delta_{\mu\nu} ]

  // // Possibly subtract \gamma_4 (coming from {\bar u})
  // if ( mu == nu )
  //   gsl_matrix_complex_sub(matProdFin,g4.gamma);


  // rescale by I/2
  gsl_matrix_complex_scale(matProd,Ihalf);

#if 0
  std::cout << "After L-mult gamma_4 & mult by I/2" << std::endl;
  LinAlg::printMat(matProd);
#endif

  // Now we have formed \gamma_4 * \sigma_{\mu\nu}
  // Proceed to form inner product with \left^\dagger & \right
  

  // Matrix vector product (on right-spinor)
  gsl_blas_zgemv(CblasNoTrans,one,matProd,right,zero,matVec);

#if 0
  for ( int i = 0; i < matVec->size; ++i )
    std::cout << gsl_vector_complex_get(matVec,i).dat[0] << " " << gsl_vector_complex_get(matVec,i).dat[1] << std::endl;
#endif

#if 0
  std::cout << "Res before: " << GSL_REAL(res) << " " << GSL_IMAG(res) << std::endl;
  std::cout << "Left spinor:" << std::endl;
  for ( int i = 0; i < left->size; ++i )
    std::cout << gsl_vector_complex_get(left,i).dat[0] << " " << gsl_vector_complex_get(left,i).dat[1] << std::endl;
#endif

  // Inner product of left^\dagger & (gamma^4\sigma^{\mu\nu} \times right)
  gsl_blas_zdotc(left,matVec,&res);

#if 0
  std::cout << "In sigma eval, we have res = ";
  std::cout << GSL_REAL(res) << " " << GSL_IMAG(res) << std::endl;
#endif

  // // And lastly, mult. result by (-1/2) if \mu\nu=4j or (-I/2) if \mu\nu=ij
  // if (( mu == 4 && nu != 4 ) || ( mu != 4 && nu == 4 ))
  //   res = gsl_complex_mul(res,mhalf);
  // else
  //   res = gsl_complex_mul(res,gsl_complex_mul(mhalf,I));

  // Push real/imag components of gsl result into std::complex val
  val.real(GSL_REAL(res)); val.imag(GSL_IMAG(res));

  gsl_vector_complex_free(matVec);
  gsl_matrix_complex_free(matMuNu);
  gsl_matrix_complex_free(matNuMu);
  gsl_matrix_complex_free(matProd);

  return val;
}
//==============================================================================================

//==============================================================================================
/*
  kinMatPDF_t methods
*/
void kinMatPDF_t::assemble(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini)
{
  if ( mu != 4 )
    {
      std::cout << "Sigma^{mu*} for mu = " << mu << " is unsupported" << std::endl;
      exit(3);
    }

  if ( fin->getL() != ini->getL() )
    {
      std::cerr << "Big problem: initial/final state spinors have different L's!" << std::endl;
      exit(4);
    }

  if ( fin->getMom() != ini->getMom() )
    {
      std::cerr << "kinMatPDF_t expects initial/final state spinors of equal momenta!" << std::endl;
      exit(4);
    }

  // The spinor contraction w/ \gamma_\mu
  ugu_t ugu(mu,MINK);

  // Map the rows of spinor at snk/src to rows of gmc matrix
  std::map<int, std::pair<int,int> > rowMap;

  // Forward case we only need two row combinations - diagonals!
  rowMap[0] = std::make_pair(fin->getTwoJ(),ini->getTwoJ());
  rowMap[1] = std::make_pair(-fin->getTwoJ(),-ini->getTwoJ());
  
  for ( auto r = rowMap.begin(); r != rowMap.end(); ++r )
    {
      std::complex<double> foo(0.0,0.0);
      foo = ugu.eval(&(fin->subduced.twoJz[r->second.first]),
		     &(ini->subduced.twoJz[r->second.second]));

      mat(r->first,0) = zeroComplex(foo);
      foo = std::complex<double>(0.0,0.0);

      // Sigma is in fact not needed. Should evaluate to zero from qrho = 0 below
      for ( int nu = 1; nu <= 4; ++nu )
	{
	  utu_t sig(mu,nu,MINK);

	  for ( int rho = 1; rho <= 4; ++rho )
	    {
	      double qrho(0.0);

	      if ( rho == 4 )
		qrho = ( fin->getE() - ini->getE() );
	      else
		qrho = ( (2*PI/fin->getL())*( fin->getMom()[rho-1] - ini->getMom()[rho-1] ) );

	      foo += sig.eval(&(fin->subduced.twoJz[r->second.first]),
			      &(ini->subduced.twoJz[r->second.second]))*metric(nu%4,rho%4)*qrho;
	    } // rho
	} // nu

      // Scale by (I/2m)
      foo *= std::complex<double>(0.0,1.0/(2*mass));

      mat(r->first,1) = zeroComplex(foo);      
    } // auto r
}
//==============================================================================================


//==============================================================================================
/*
  kinMat_t methods
*/
void kinMat_t::assemble(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini,
			const std::vector<int> &disp)
{
  if ( mu != 4 )
    {
      std::cout << "Sigma^{mu*} for mu = " << mu << " is unsupported" << std::endl;
      exit(3);
    }

  if ( fin->getL() != ini->getL() )
    {
      std::cerr << "Big problem: initial/final state spinors have different L's!" << std::endl;
      exit(4);
    }

  // The spinor contraction w/ \gamma_\mu
  ugu_t ugu(mu,MINK);

  // Map the rows of spinor at snk/src to rows of gmc matrix
  std::map<int, std::pair<int,int> > rowMap;

  // Off-forward case we need all four row combinations
  rowMap[0] = std::make_pair(fin->getTwoJ(),ini->getTwoJ());
  rowMap[1] = std::make_pair(fin->getTwoJ(),-ini->getTwoJ());
  rowMap[2] = std::make_pair(-fin->getTwoJ(),ini->getTwoJ());
  rowMap[3] = std::make_pair(-fin->getTwoJ(),-ini->getTwoJ());

  for ( auto r = rowMap.begin(); r != rowMap.end(); ++r )
    {
      std::complex<double> foo(0.0,0.0);
      foo = ugu.eval(&(fin->subduced.twoJz[r->second.first]),
		     &(ini->subduced.twoJz[r->second.second]));

      mat(r->first,0) = zeroComplex(foo);
      foo = std::complex<double>(0.0,0.0);

      for ( int nu = 1; nu <= 4; ++nu )
	{
	  utu_t sig(mu,nu,MINK);

	  for ( int rho = 1; rho <= 4; ++rho )
	    {
	      double qrho(0.0);

	      if ( rho == 4 )
		qrho = ( fin->getE() - ini->getE() );
	      else
		qrho = ( (2*PI/fin->getL())*( fin->getMom()[rho-1] - ini->getMom()[rho-1] ) );

	      foo += sig.eval(&(fin->subduced.twoJz[r->second.first]),
			      &(ini->subduced.twoJz[r->second.second]))*metric(nu%4,rho%4)*qrho;
	    } // rho
	} // nu

      // Scale by (I/2m)
      foo *= std::complex<double>(0.0,1.0/(2*mass));

      mat(r->first,1) = zeroComplex(foo);


      /*
	Mon. Oct. 24, 2022 - Explore adding in a q^\mu term
      */
      foo = std::complex<double>(0.0,0.0);

      // The spinor contraction w/ identity
      u1u_t u1u(MINK);

      foo = u1u.eval(&(fin->subduced.twoJz[r->second.first]),
		     &(ini->subduced.twoJz[r->second.second]));
      
      if ( mu == 4 )
	foo *= ( fin->getE() - ini->getE() );
      else
	foo *= ( (2*PI/fin->getL())*( fin->getMom()[mu] - ini->getMom()[mu] ) );
      // Rescale by 1/2m so energy dimensions agree
      foo *= std::complex<double>(-1.0/mass,0.0); // "-" so we have (p1 - p2)

      mat(r->first,2) = zeroComplex(foo);


      /*
	Tues. Oct. 25, 2022 - Explore adding in a im\sigma^{\mu\nu}z_\nu term
      */
      // Dead in the water --> Need eight matrix elements to constrain!
      foo = std::complex<double>(0.0,0.0);

      for ( int nu = 1; nu <= 4; ++nu )
	{
	  utu_t sig(mu,nu,MINK);
	  for ( int rho = 1; rho <= 4; ++rho )
	    {
	      double zrho(0.0);

	      if ( rho == 4 )
		zrho = 0.0;
	      else
		zrho = disp[rho-1];

	      foo += sig.eval(&(fin->subduced.twoJz[r->second.first]),
			      &(ini->subduced.twoJz[r->second.second]))*metric(rho%4,nu%4)*zrho;
	    } // rho
	} //nu
      
      // Scale by (I*m)
      foo *= std::complex<double>(0.0,1.0*mass);

      mat(r->first,3) = zeroComplex(foo);
    } // auto r
}
//==============================================================================================

/*
  kinMat3_t methods
*/
void kinMat3_t::assemble(int mu, bool MINK, double mass, Spinor *s, const std::vector<int> &disp)
{
  // The polarization vector associated with mu
  polVec_t S(mu,MINK);

  // Map the rows of spinor at snk/src to rows of gmc matrix
  std::map<int, std::pair<int,int> > rowMap;
  rowMap[0] = std::make_pair(s->getTwoJ(),s->getTwoJ());
  rowMap[1] = std::make_pair(s->getTwoJ(),-s->getTwoJ());
  rowMap[2] = std::make_pair(-s->getTwoJ(),s->getTwoJ());
  rowMap[3] = std::make_pair(-s->getTwoJ(),-s->getTwoJ());

#if 0
  // For vector insertions
  rowMap[0] = std::make_pair(s->getTwoJ(),s->getTwoJ());
  rowMap[1] = std::make_pair(-s->getTwoJ(),-s->getTwoJ());
#endif
#if 0
  // For axial insertions
  if ( shortMom(s->getMom(),"") == "000" )
    {
      rowMap[0] = std::make_pair(s->getTwoJ(),s->getTwoJ());
      rowMap[1] = std::make_pair(-s->getTwoJ(),-s->getTwoJ());
    }
  else
    {
      rowMap[0] = std::make_pair(s->getTwoJ(),-s->getTwoJ());
      rowMap[1] = std::make_pair(-s->getTwoJ(),s->getTwoJ());
    }
#endif


  double pref = 1.0;

  for ( auto r = rowMap.begin(); r != rowMap.end(); ++r )
    {
#if 1
      // Use foo to collect spinor contractions
      std::complex<double> foo(pref*mass*S.eval(&(s->subduced.twoJz[r->second.first]),
						&(s->subduced.twoJz[r->second.second]),mass));

      // Mat[*,0] = -1*<<gamma^mu gamma^5>>
      mat(r->first,0) = foo*std::complex<double>(0.0,-1.0); // _I_;
      // mat(r->first,0) = foo;

      foo = 0.0;
      // Loop over Lorentz indices to compute z_nu \cdot <<gamma^nu\gamma^5>>
      // Skip 4th component, since \vec{z}^4 is always zero
      for ( int nu = 1; nu < 4; ++nu )
	{
	  polVec_t tmpS(nu,MINK);

	  foo += tmpS.eval(&(s->subduced.twoJz[r->second.first]),
			   &(s->subduced.twoJz[r->second.second]),mass)*disp[nu-1]; //*_I_;
	}

      // Mat[*,1] = 2 m^3 z^mu (z_nu \cdot <<gamma^nu\gamma^5>>)
      // mat(r->first,1) = pref*pow(mass,3)*disp[mu-1]*foo;
      mat(r->first,1) = pref*pow(mass,3)*disp[mu-1]*foo*std::complex<double>(0.0,-1.0);
#endif
    } // auto r
}
//=========


void kinMat3_t::assembleBig(int mu, bool MINK, double mass, Spinor *s, const std::vector<int> &disp)
{
  // The polarization vector associated with mu
  polVec_t S(mu,MINK);

  // Map the rows of spinor at snk/src to rows of gmc matrix
  std::map<int, std::pair<int,int> > rowMap;
  rowMap[0] = std::make_pair(s->getTwoJ(),s->getTwoJ());
  rowMap[1] = std::make_pair(s->getTwoJ(),-s->getTwoJ());
  rowMap[2] = std::make_pair(-s->getTwoJ(),s->getTwoJ());
  rowMap[3] = std::make_pair(-s->getTwoJ(),-s->getTwoJ());

  std::complex<double> prefactor(0.0,0.0);
  std::complex<double> polVecEval(0.0,0.0);

  for ( auto r = rowMap.begin(); r != rowMap.end(); ++r )
    {
      polVecEval = S.eval(&(s->subduced.twoJz[r->second.first]),
			  &(s->subduced.twoJz[r->second.second]),mass);
      
      mat(r->first,0) = zeroComplex(polVecEval*std::complex<double>(-2*mass,0.0));
      

      // Reset
      polVecEval = std::complex<double>(0.0,0.0);
      // Loop over Lorentz indices to compute z_\nu \cdot <<gamma^nu\gamma^5>>
      for ( int nu = 1; nu < 4; ++nu )
	{
	  std::complex<double> foo;
	  polVec_t tmpS(nu,MINK);
	  
	  foo = tmpS.eval(&(s->subduced.twoJz[r->second.first]),
			  &(s->subduced.twoJz[r->second.second]),mass);

	  for ( int rho = 1; rho <= 4; ++rho )
	    polVecEval += metric(nu%4,rho%4)*disp[rho-1]*foo;
	}

      double pMu = ( mu < 4 ) ? s->getMom()[mu-1] : s->getE();

      mat(r->first,1) = zeroComplex(polVecEval*std::complex<double>(0.0,-2*mass*pMu));
      mat(r->first,2) = zeroComplex(polVecEval*std::complex<double>(2*pow(mass,3)*disp[mu-1]));
    } // auto r
}


/*
  ffMat_t methods
*/
void ffMat_t::assemble(int mu, bool MINK, double mass, Spinor *fin, Spinor *ini)
{
  // Prefactor
  double pref = 1.0;
  // Generic
  ugu_t ff1 = ugu_t(mu,MINK);

  // // Scalar current selected
  // if ( gamma == 0 )
  //   ff1 = u1u_t(MINK);
  // // Vector current selected
  // else if ( gamma == 1 || gamma == 2 || gamma == 4 || gamma == 8 )
  //   ff1 = ugu_t(mu,MINK);
  // // Axial current selected
  // else if ( gamma == 7 || gamma == 11 || gamma == 13 || gamma == 14 )
  //   ff1 = polVec_t(mu,MINK);


  // Map the rows of spinor at snk/src to rows of gmc matrix
  std::map<int, std::pair<int,int> > rowMap;
  rowMap[0] = std::make_pair(fin->getTwoJ(),ini->getTwoJ());
  rowMap[1] = std::make_pair(fin->getTwoJ(),-ini->getTwoJ());
  rowMap[2] = std::make_pair(-fin->getTwoJ(),ini->getTwoJ());
  rowMap[3] = std::make_pair(-fin->getTwoJ(),-ini->getTwoJ());

  // Iterate over rowMap pairs and assemble the kinematic matrix
  for ( auto r = rowMap.begin(); r != rowMap.end(); ++r )
    {
      std::complex<double> foo(0.0,0.0);
      foo = pref*ff1.eval(&(fin->subduced.twoJz[r->second.first]),
			  &(ini->subduced.twoJz[r->second.second]));

      mat(r->first,0) = zeroComplex(foo);
      foo = std::complex<double>(0.0,0.0);

      for ( int nu = 1; nu <= 4; ++nu )
	{
	  utu_t ff2(mu,nu,MINK);

	  for ( int rho = 1; rho <= 4; ++rho )
	    {
	      double qrho(0.0);
	      
	      if ( rho == 4 )
		qrho = ( fin->getE() - ini->getE() );
	      else
		qrho = ( (2*PI/fin->getL())*( fin->getMom()[rho-1] - ini->getMom()[rho-1] ) );


	      foo += ff2.eval(&(fin->subduced.twoJz[r->second.first]),
			      &(ini->subduced.twoJz[r->second.second]))*metric(nu%4,rho%4)*qrho;
	    } // rho
	} //nu

      // Scale by (I/2m)
      foo *= std::complex<double>(0.0,1.0/(2*mass));

      mat(r->first,1) = zeroComplex(foo);
    } // auto r	  
}

// void ffMat_t::assemble(int mu, int nu, bool MINK, double mass, Spinor *s)
// {
//   auto ff1 = utu_t(mu,nu,MINK);
// }

//**********************************************************************************************
/*
  Spinor class methods
*/
void Spinor::projector()
{
  gmc * extProd = gsl_matrix_complex_calloc(4,4);
  gmc * projector = gsl_matrix_complex_calloc(4,4);
  diracMat_t g4 = diracMat_t(4,false);
  
  // Compute u(p,s1)u^\dagger(p,s1) + u(p,s2)u^\dagger(p,s2)
  for ( int s = 1; s > -2; s-=2 )
    gsl_blas_zgeru(one,&canon.twoJz[s],&canon.twoJz[s],extProd);

  // Right multiply by \gamma_4 to complete \sum_s u(p,s)\bar{u}(p,s)
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,extProd,g4.gamma,one,projector);

  std::cout << "PROJ!" << std::endl;
  std::cout << projector << std::endl;

  gsl_matrix_complex_free(extProd);
  gsl_matrix_complex_free(projector);
  gsl_matrix_complex_free(g4.gamma);
}

void Spinor::buildSpinors()
{
  // Canonical momentum is |\vec{p}| in +z-direction
  XMLArray::Array<double> canonMom(3);
  canonMom[0]=0.0; canonMom[1]=0.0; canonMom[2]=absMom();

  // Start by building absolute & canonical spinors: u(\vec{p},s) u(|\vec{p}|\hat{z},s)
  absolute.build(mom,E,m,L);
  canon.build(canonMom,E,m,L);

  // // Build helicity spinors by applying Euler rotations to upper/lower comps. of canon spinors
  // helicity = canon;


  std::cout << "Canonical spinors(1): " << &canon.twoJz[1] << std::endl;
  std::cout << "Canonical spinors(2): " << &canon.twoJz[-1] << std::endl;
  // Interate over both twoJz components of helicity spinor
  for ( int i = twoJ; i >= -twoJ; i -= subductInfo.irrep_dim )
    {
      gsl_vector_complex_view sliceUp   = gsl_vector_complex_subvector(&canon.twoJz[i],0,2);
      gsl_vector_complex_view sliceDown = gsl_vector_complex_subvector(&canon.twoJz[i],2,2);
      gvc * up   = gsl_vector_complex_calloc(2);
      gvc * down = gsl_vector_complex_calloc(2);
      // gsl_vector_complex_memcpy(up,&sliceUp.vector);
      // gsl_vector_complex_memcpy(down,&sliceDown.vector);
      
      // Left-mult upper/lower components of this twoJz spinor by Euler rotation matrix
      gsl_blas_zgemv(CblasNoTrans,one,eulerRot,&sliceUp.vector,one,up);
      gsl_blas_zgemv(CblasNoTrans,one,eulerRot,&sliceDown.vector,one,down);

      // Concatenate up/down pieces to form 4-component helicity spinor
      gvc * hel = gsl_vector_complex_calloc(4);
      gsl_vector_complex_set(hel,0,gsl_vector_complex_get(up,0));
      gsl_vector_complex_set(hel,1,gsl_vector_complex_get(up,1));
      gsl_vector_complex_set(hel,2,gsl_vector_complex_get(down,0));
      gsl_vector_complex_set(hel,3,gsl_vector_complex_get(down,1));

      // Pack and insert
      // helicity.twoJz[i] = *hel;
      std::pair<int, gvc> helSpinor(i,*hel);
      helicity.twoJz.insert(helSpinor);
    } // i
  // Now have helicity spinors


  std::cout << "Helicity spinors(1): " << &helicity.twoJz[1] << std::endl;
  std::cout << "Helicity spinors(2): " << &helicity.twoJz[-1] << std::endl;

#if 1
#warning "Packing full spinor"
#else
#warning "Packing non-relativistic spinor"

  gsl_vector_complex * goo = gsl_vector_complex_calloc(4);
  // gsl_vector_complex * toIns = gsl_vector_complex_alloc(4);
  // gsl_vector_complex_memcpy(toIns,&helicity.twoJz[1]);
  gsl_vector_complex_memcpy(goo,&helicity.twoJz[1]);

  gmc * ID4d = gsl_matrix_complex_calloc(4,4); gsl_matrix_complex_set_identity(ID4d);
  gmc * nrelProj = gsl_matrix_complex_calloc(4,4);
  gsl_matrix_complex_set_identity(nrelProj);

  diracMat_t dd(4,true);

  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,half,dd.gamma,ID4d,half,nrelProj);

  // Now left-mult the non-relativistic projector onto spinor
  gsl_blas_zgemv(CblasNoTrans,one,nrelProj,goo,zero,&helicity.twoJz[1]);

  
  gsl_vector_complex_memcpy(goo,&helicity.twoJz[-1]);
  gsl_blas_zgemv(CblasNoTrans,one,nrelProj,goo,zero,&helicity.twoJz[-1]);
#endif




  gvc * subSpinor = gsl_vector_complex_calloc(4);
  gsl_blas_zaxpy(gsl_matrix_complex_get(coeffS,0,1),&helicity.twoJz[-1],subSpinor);
  gsl_blas_zaxpy(gsl_matrix_complex_get(coeffS,0,0),&helicity.twoJz[1],subSpinor);

  std::pair<int, gvc> subIns(1,*subSpinor);
  subduced.twoJz.insert(subIns);


  gvc * subSpinor2 = gsl_vector_complex_calloc(4);
  gsl_blas_zaxpy(gsl_matrix_complex_get(coeffS,1,1),&helicity.twoJz[-1],subSpinor2);
  gsl_blas_zaxpy(gsl_matrix_complex_get(coeffS,1,0),&helicity.twoJz[1],subSpinor2);
  subIns = std::make_pair(-1,*subSpinor2);
  subduced.twoJz.insert(subIns);

  std::cout << "Subduced spinors(1): " << &subduced.twoJz[1] << std::endl;
  std::cout << "Subduced spinors(2): " << &subduced.twoJz[-1] << std::endl;

  // gsl_vector_complex_free(subSpinor);
  // gsl_vector_complex_free(subSpinor2);


#if 0
  LinAlg::printMat(coeffS);

  // Build subduced spinor from canonical spinor
  gmc * prod = gsl_matrix_complex_calloc(subductInfo.irrep_dim,twoJ+1);
  gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,one,coeffS,eulerRot,zero,prod);
  LinAlg::printMat(prod);
  for ( int i = twoJ; i >= -twoJ; i -= subductInfo.irrep_dim )
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


      // Left mult. product of Subduction/spin=1/2 Euler rotation matrix onto upper/lower components
      gsl_blas_zgemv(CblasNoTrans,one,prod,fooUp,zero,tmpUp);
      gsl_blas_zgemv(CblasNoTrans,one,prod,fooDown,zero,tmpDown);

     
#if 0
      std::cout << "\nup view check: size = " << up.vector.size << std::endl;
      std::cout << &up.vector << std::endl;
      std::cout << "\ndown view check: size = " << down.vector.size << std::endl;
      std::cout << &down.vector << std::endl;
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

      std::cout << "Inserted spinor = " << toIns << std::endl;

      std::pair<int, gvc> aSubducedSpinor(i,*toIns);
      subduced.twoJz.insert(aSubducedSpinor); // [i] = *foo;
      // gsl_vector_complex_memcpy(&subduced.twoJz[i],foo);

      gsl_vector_complex_free(foo);
#endif
    } // end i
#endif
} // buildSpinors

void Spinor::initSubduce(const std::string& opName)
{
  // Gather subduction info for this operator
  subductInfo = subduceInfo(opName,mom);
  // Set the ranks of subduction coefficient
  coeffS = gsl_matrix_complex_calloc(subductInfo.irrep_dim,twoJ+1);
  

  // // Build the Wigner-D
  // for ( int i = 0; i < subductInfo.irrep_dim; ++i )
  //   {
  //     int twoJz_i = pow((-1),i);
  //     for ( int j = 0; j < subductInfo.irrep_dim; ++j )
  // 	{
  // 	  int twoJz_j = pow((-1),j);
  // 	  gsl_matrix_complex_set(wig,i,j,gc_rect(Hadron::Wigner_D(twoJ,twoJz_i,twoJz_j,rot.alpha,rot.beta,rot.gamma).real(),
  // 						 Hadron::Wigner_D(twoJ,twoJz_i,twoJz_j,rot.alpha,rot.beta,rot.gamma).imag()));
  // 	}
  //   }


  // Build the subduction matrix
  for ( int row = 1; row <= subductInfo.irrep_dim; ++row )
    {
      for ( int h = 1; h <= twoJ+1; ++h )
	{
	  gsl_matrix_complex_set(coeffS,row-1,h-1,
				 gc_rect((*subductInfo.H).operator()(row,h).real(),
					 (*subductInfo.H).operator()(row,h).imag()));
	}
    }
}
// End of Spinor class methods
//************************************************************************************************



/*
  EXTRACT INVARIANT AMPLITUDES USING (IN GENERAL) AN SVD DECOMPOSITION
*/
// //--------- Unpol. PDF */
// void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * MAT,
// 		   std::vector<kinMat_t> * KIN,
// 		   std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * AMP)
//--------- Unpol. GPD
void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * MAT,
		   std::vector<kinMat_t> * KIN,
		   std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * AMP)
{
  for ( auto a = AMP->begin(); a != AMP->end(); ++a )
    {
      int idx = std::distance(AMP->begin(),a);
      Eigen::JacobiSVD<Eigen::MatrixXcd> SVD((*KIN)[idx].mat,
					     Eigen::ComputeFullU | Eigen::ComputeFullV);
#ifdef PRINT_SVS
      std::cout << "      --> singular values = " << SVD.singularValues() << std::endl;
#endif
      // SVD.matrixU();
      // SVD.matrixV();

      // Solution of SVD to AMP
      *a = SVD.solve((*MAT)[idx]);
    } // a
}

void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * MAT,
		   std::vector<kinMat3_t> * KIN,
		   std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * AMP)
// void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * MAT,
// 		   std::vector<kinMat3_t> * KIN,
// 		   std::vector<Eigen::Vector2cd> * AMP)
{
  for ( auto a = AMP->begin(); a != AMP->end(); ++a )
    {
      int idx = std::distance(AMP->begin(),a);
      Eigen::JacobiSVD<Eigen::MatrixXcd> SVD((*KIN)[idx].mat,
					     Eigen::ComputeFullU | Eigen::ComputeFullV);
#ifdef PRINT_SVS
      std::cout << "      --> singular values = " << SVD.singularValues() << std::endl;
#endif

      // Solution of SVD to AMP
      *a = SVD.solve((*MAT)[idx]);
    } // a
}

void extAmplitudes(std::vector<Eigen::Matrix<std::complex<double>, 4, 1> > * MAT,
		   std::vector<ffMat_t> * KIN,
		   std::vector<Eigen::Matrix<std::complex<double>, 2, 1> > * AMP)
{
  for ( auto a = AMP->begin(); a != AMP->end(); ++a )
    {
      int idx = std::distance(AMP->begin(),a);
      Eigen::JacobiSVD<Eigen::MatrixXcd> SVD((*KIN)[idx].mat,
      					     Eigen::ComputeFullU | Eigen::ComputeFullV);
#ifdef PRINT_SVS
      std::cout << "      --> singular values = " << SVD.singularValues() << std::endl;
#endif

      // Solution of SVD to AMP
      *a = SVD.solve((*MAT)[idx]);
    } // a
}


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
