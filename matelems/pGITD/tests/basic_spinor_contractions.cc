/*
  Debug:
     -> construction of spinors
     -> contraction of fin/ini spinors with \gamma_3\gamma_5, \gamma_4, \sigma_{\mu\nu}
*/
#include "three_pt.h"

int main (int argc, char *argv[])
{
  Spinor dumF(global.opMomXML[shortMom(global.pf,"")],global.pf,0.66365470597820764,0.535,global.Lx); // (011)
  Spinor dumI(global.opMomXML[shortMom(global.pi,"")],global.pi,0.56989309716099856,0.535,global.Lx); // (001)
  dumF.buildSpinors();
  dumI.buildSpinors();

  std::cout << &(dumF.absolute.twoJz[1]) << std::endl;
  std::cout << &(dumF.absolute.twoJz[-1]) << std::endl;
  std::cout << &(dumF.canon.twoJz[1]) << std::endl;
  std::cout << &(dumF.canon.twoJz[-1]) << std::endl;
  std::cout << &(dumF.subduced.twoJz[1]) << std::endl;
  std::cout << &(dumF.subduced.twoJz[-1]) << std::endl;

  
  polVec_t P(3,true);
  std::cout << "PolVec[11] = " << P.eval(&(dumI.subduced.twoJz[1]),&(dumI.subduced.twoJz[1])) << std::endl;
  std::cout << "PolVec[12] = " << P.eval(&(dumI.subduced.twoJz[1]),&(dumI.subduced.twoJz[-1])) << std::endl;
  std::cout << "PolVec[21] = " << P.eval(&(dumI.subduced.twoJz[-1]),&(dumI.subduced.twoJz[1])) << std::endl;
  std::cout << "PolVec[22] = " << P.eval(&(dumI.subduced.twoJz[-1]),&(dumI.subduced.twoJz[-1])) << std::endl;
  exit(80);


  ugu_t foo(4,false);
  std::cout << "[a1]Inner prod.= " << foo.eval(&(dumF.absolute.twoJz[1]),&(dumF.absolute.twoJz[1])) << std::endl;
  std::cout << "[a-1]Inner prod.= " << foo.eval(&(dumF.absolute.twoJz[-1]),&(dumF.absolute.twoJz[-1])) << std::endl;
  std::cout << "[c1]Inner prod.= " << foo.eval(&(dumF.canon.twoJz[1]),&(dumF.canon.twoJz[1])) << std::endl;
  std::cout << "[c-1]Inner prod.= " << foo.eval(&(dumF.canon.twoJz[-1]),&(dumF.canon.twoJz[-1])) << std::endl;
  std::cout << "[s11]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[1]),&(dumF.subduced.twoJz[1])) << std::endl;
  std::cout << "[s1-1]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[1]),&(dumF.subduced.twoJz[-1])) << std::endl;
  std::cout << "[s-11]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[-1]),&(dumF.subduced.twoJz[1])) << std::endl;
  std::cout << "[s-1-1]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[-1]),&(dumF.subduced.twoJz[-1])) << std::endl;

  std::cout << "vvvvvvvvvv Off-forward contractions vvvvvvvvv" << std::endl;
  std::cout << "[a1]Inner prod.= " << foo.eval(&(dumF.absolute.twoJz[1]),&(dumI.absolute.twoJz[1])) << std::endl;
  std::cout << "[a-1]Inner prod.= " << foo.eval(&(dumF.absolute.twoJz[-1]),&(dumI.absolute.twoJz[-1])) << std::endl;
  std::cout << "[c1]Inner prod.= " << foo.eval(&(dumF.canon.twoJz[1]),&(dumI.canon.twoJz[1])) << std::endl;
  std::cout << "[c-1]Inner prod.= " << foo.eval(&(dumF.canon.twoJz[-1]),&(dumI.canon.twoJz[-1])) << std::endl;
  std::cout << "[s11]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[1]),&(dumI.subduced.twoJz[1])) << std::endl;
  std::cout << "[s1-1]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[1]),&(dumI.subduced.twoJz[-1])) << std::endl;
  std::cout << "[s-11]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[-1]),&(dumI.subduced.twoJz[1])) << std::endl;
  std::cout << "[s-1-1]Inner prod.= " << foo.eval(&(dumF.subduced.twoJz[-1]),&(dumI.subduced.twoJz[-1])) << std::endl;



  std::cout << "Canonical test sigma" << std::endl;
  utu_t sigTest(4,2,false);
  std::cout << sigTest.eval(&(dumF.absolute.twoJz[1]),&(dumI.absolute.twoJz[1])) << std::endl;

  std::cout << "Contract w/ sigma" << std::endl;
  for ( int i = 1; i < 4; ++i )
    {
      utu_t sig(4,i,false);
      std::cout << "(4" << i << "): s[11] = " << sig.eval(&(dumF.subduced.twoJz[1]),&(dumI.subduced.twoJz[1])) << std::endl;
      std::cout << "(4" << i << "): s[1-1] = " << sig.eval(&(dumF.subduced.twoJz[1]),&(dumI.subduced.twoJz[-1])) << std::endl;
      std::cout << "(4" << i << "): s[-11] = " << sig.eval(&(dumF.subduced.twoJz[-1]),&(dumI.subduced.twoJz[1])) << std::endl;
      std::cout << "(4" << i << "): s[-1-1] = " << sig.eval(&(dumF.subduced.twoJz[-1]),&(dumI.subduced.twoJz[-1])) << std::endl;
      std::cout << "---------------------------" << std::endl;
    }


  return 0;
