/*
  Test validity of Gordon Identity w/in Code
*/
#include "three_pt.h"

int main (int argc, char *argv[])
{
  double M(0.535);
  double Ef(0.63394355156848614);
  double Ei(0.66365470597820764);
  u1u_t id(true);
  ugu_t gam(4,true);

  std::complex<double> rhs(0.0,0.0);

  Spinor fin(global.opMomXML[shortMom(global.pf,"")],
	     global.pf,Ef,M,global.Lx);
  Spinor ini(global.opMomXML[shortMom(global.pi,"")],
	     global.pi,Ei,M,global.Lx);
  ini.buildSpinors();
  fin.buildSpinors();


  /*
    \bar{u}(p,s')\gamma^\mu u(p,s) = \bar{u}(p,s')(p^\mu/M) u(p,s)

    std::cout << "<<gamma_mu>> = " << gam.eval(&fin.absolute.twoJz[1],&fin.absolute.twoJz[1]) << std::endl;
    std::cout << "rhs = " << ((Ef/M)*id.eval(&fin.absolute.twoJz[1],&fin.absolute.twoJz[1])) << std::endl;
    
    PASSED
  */

  /*
    <<\sigma_{43}>> = 0.18897240910531166
    w/ Ef(1,1,1)=0.63394355156848614, Ei(0,0,2)=0.66365470597820764, M=0.535

    utu_t sig(4,3,true);
    std::cout << sig.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1]) << std::endl;

    PASSED
  */
  utu_t sig_x(4,1,true); utu_t sig_y(4,2,true); utu_t sig_z(4,3,true);
  std::cout << sig_x.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1]) << std::endl;
  std::cout << sig_y.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1]) << std::endl;
  std::cout << sig_z.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1]) << std::endl;

  /*
    \bar{u}(p,s')u(p,s) = 2m\delta_{s's}

    std::cout << "2m = " << 2*M << std::endl;
    for ( int i = 1; i < 3; ++i )
    {
    for ( int j = 1; j < 3; ++j )
    {
    std::cout << "ID CHECK[ij] = " << id.eval(&fin.absolute.twoJz[pow(-1,i+1)],&fin.absolute.twoJz[pow(-1,j+1)])
    << std::endl;
    
    PASSED
  */
  
  for ( int i = 1; i <=4; ++i )
    {
      diracMat_t D(i,true);
      LinAlg::printMat(D.gamma);
    }
  

  std::cout << "2m*<<gamma_mu>> = "
  	    << (2*M*gam.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1]))
  	    << std::endl;

  rhs += (id.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1])*(Ef+Ei));
  // std::cout << "RHS b4 sigma = " << rhs << std::endl;

  std::complex<double> sigTot(0.0,0.0);

  for ( int j = 1; j < 4; ++j )
    {
      utu_t sig(4,j,true);

      // std::cout << "This sigma eval = " << sig.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1])
      // 		<< std::endl;
     
      // Minus from metric tensor!
      sigTot -= (sig.eval(&fin.absolute.twoJz[1],&ini.absolute.twoJz[1])*(2*M_PIl/fin.getL())*( fin.getMom()[j-1] - ini.getMom()[j-1] ));

    }
  // std::cout << "sigTot = " << sigTot << std::endl;


  sigTot *= std::complex<double>(0,1.0);
  // std::cout << "sigTot = " << sigTot << std::endl;

  rhs += sigTot;
  std::cout << "RHS = " << rhs << std::endl;

  /*
    GORDAN IDENTITY IS SATISFIED! JUST DONT FORGET THE METRIC IN THE SIGMA * (PF-PI) PIECE
  */

  return 0;
}
