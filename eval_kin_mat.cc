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

#include "hadron_npt_op_factory.h"
#include "single_hadron_factory.h"
#include "cont_hadron_op_factory.h"
#include "operators/spin_subduced_obj_factory.h"


using namespace ADATXML;
using namespace Pseudo;
using namespace NCOR;

typedef std::complex<double> cd;

/*
  Register all the ops
*/
bool registerAll(void)
{
  bool foo = true;

  try
    {
      // Operators                                                                                                                                                                                              
      foo &= Redstar::HadronNptOpEnv::registerAll();
    }
  catch(const std::out_of_range& e)
    {
      std::cerr << __PRETTY_FUNCTION__ << ": Caught out-of-range exception: " << e.what() << std::endl;
      exit(1);
    }
  catch(...)
    {
      std::cerr << __PRETTY_FUNCTION__ << ": caught generic exception" << std::endl;
      exit(1);
    }

  return foo;
}


void getSVs(Eigen::MatrixXcd *h)
{
  Eigen::JacobiSVD<Eigen::MatrixXcd> s(*h,Eigen::ComputeFullU | Eigen::ComputeFullV);
  std::cout << s.singularValues() << std::endl;
}

std::vector<std::string> subName(const ADAT::MapObject<std::string, int> cont_ops, const Array<int> &mom)
{
  std::vector<std::string> sub_ops;
  // Iterate continuum names
  for ( auto it = cont_ops.begin(); it != cont_ops.end(); ++it )
    {
      // Register the continuum op
      auto cont_reg = Redstar::ContHadronOpEnv::TheHadronContOpRegInfoFactory::Instance().at(it->first);
      int cont_dim = cont_reg.twoJ + 1;

      // Register all subduced versions of continuum operator
      auto allSubOps = Redstar::ContHadronOpEnv::TheHadronContOpSubducedOpFactory::Instance().at(it->first);

      // Iterate over all subduced operators of continuum operator
      for ( auto sub_op = allSubOps.begin(); sub_op != allSubOps.end(); ++sub_op )
	{
	  // Register this specfic op
	  auto reg = Redstar::TheSingleHadronOpsRegInfoFactory::Instance().at(sub_op->first);
	  
	  // Init this specific op
	  auto spin_reg = Redstar::SpinObject::SubducedSpinObjRegFactory::Instance().at(reg.spin_id);
	  auto spin_op  = Redstar::SpinObject::TheSubducedSpinObjFactory::Instance().createObject(reg.spin_id,reg.spin_id);
	  
	  std::string lg = Hadron::getIrrepLG(spin_reg.irrep);
	  int        dim = Hadron::getIrrepDim(Hadron::removeHelicity(spin_reg.irrep));


	  // If little group doesn't match with what's induced by momentum, continue
	  if ( lg != Hadron::generateLittleGroup(mom) ) continue;
	  std::cout << "     MATCHED! Working w/ LG = " << lg << " and momentum = "
		    << mom << std::endl;
	  std::cout << "Sub Op name = " << sub_op->first << std::endl;

	  
	  sub_ops.push_back(sub_op->first);

	} // sub_op
    } // it
  return sub_ops;
}

int main(int argc, char *argv[])
{
  if ( argc != 11 )
    {
      std::cout << "\nUsage: $0 <pfx> <pfy> <pfz> <pix> <piy> <piz> <dx> <dy> <dz> <mu>\n" << std::endl;
      exit(1);
    }
  // Register all the necessary factories for relating lattice matrix elements and helicity amplitudes
  // Register all the necessary factories
  bool foo = registerAll();

  int pfx, pfy, pfz, pix, piy, piz, dx, dy, dz, mu;

  std::stringstream ss;
  // Get ints from CL
  ss << argv[1]; ss >> pfx; ss.clear();
  ss << argv[2]; ss >> pfy; ss.clear();
  ss << argv[3]; ss >> pfz; ss.clear();
  ss << argv[4]; ss >> pix; ss.clear();
  ss << argv[5]; ss >> piy; ss.clear();
  ss << argv[6]; ss >> piz; ss.clear();
  ss << argv[7]; ss >> dx; ss.clear();
  ss << argv[8]; ss >> dy; ss.clear();
  ss << argv[9]; ss >> dz; ss.clear();
  ss << argv[10]; ss >> mu; ss.clear();
  

  double mass = 0.535;
  int Lx      = 32;

  XMLArray::Array<int> mom_i;
  mom_i.resize(3);
  mom_i[0] = pix; mom_i[1] = piy; mom_i[2] = piz;
  XMLArray::Array<int> mom_f;
  mom_f.resize(3);
  mom_f[0] = pfx; mom_f[1] = pfy; mom_f[2] = pfz;

  double Ei    = sqrt( pow(mass,2) + pow(2*PI/Lx,2)*(mom_i*mom_i) );
  double Ef    = sqrt( pow(mass,2) + pow(2*PI/Lx,2)*(mom_f*mom_f) );

  std::cout << "Ef = " << Ef << std::endl;
  std::cout << "Ei = " << Ei << std::endl;
  
  // // std::string opName = "NucleonMG1g1MxD0J0S_J1o2_G1g1";
  // std::string opName = "NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1";
  // std::string opNameD2 = "NucleonMG1g1MxD0J0S_J1o2_H1o2D2E";



  // std::vector<K> tempKey3pt = templateKeys(global,3);
  // K tempKey2Pf; tempKey2Pf.npoint.resize(2);
  // K tempKey2Pi; tempKey2Pi.npoint.resize(2);
  // for ( int i = 1; i <= 2; ++i )
  //   {
  //     tempKey2Pf.npoint[i] = tempKey3pt[0].npoint[1];
  //     tempKey2Pi.npoint[i] = tempKey3pt[0].npoint[3];
  //   }

  // tempKey2Pi.npoint[1].irrep.op.ops[1].name;
  // tempKey2Pf.npoint[1].irrep.op.ops[1].name;



  std::vector<int> disp_list = { dx, dy, dz };

  // Hadron::KeyCGCIrrepMom_t pfIrrepMom(1,g.pf);
  // Hadron::KeyCGCIrrepMom_t piIrrepMom(1,g.pi);

  // Hadron::KeyParticleOp_t(subName("NucleonMG1g1MxD0J0S_J1o2",mom_i)[0],
  // 			  "", Hadron::canonicalOrder(mom_i),disp_list);
  


  ADAT::MapObject<std::string, int> dumMap;
  dumMap.insert("NucleonMG1g1MxD0J0S_J1o2", 1);

  std::string opPi = subName(dumMap,mom_i)[0];
  std::string opPf = subName(dumMap,mom_f)[0];
  
  std::vector<int> disp(disp_list); disp.push_back(0);


  Spinor si(opPi,mom_i,Ei,mass,Lx);
  si.buildSpinors();
  Spinor sf(opPf,mom_f,Ef,mass,Lx);
  sf.buildSpinors();
 

  // diracMat_t g4(4,true); diracMat_t gx(1,true); diracMat_t gy(2,true); diracMat_t gz(3,true);
  diracMat_t g4(4,false); diracMat_t gx(1,false); diracMat_t gy(2,false); diracMat_t gz(3,false);
  diracMat_t g5(5,true);
  
  std::cout << "G4 = " << g4.gamma << std::endl;
  std::cout << "GX = " << gx.gamma << std::endl;
  std::cout << "GY = " << gy.gamma << std::endl;
  std::cout << "GZ = " << gz.gamma << std::endl;
  std::cout << "G5 = " << g5.gamma << std::endl;



  std::cout << "Si ref Angles = " << si.getRefRot().alpha << " " << si.getRefRot().beta
	    << " " << si.getRefRot().gamma << std::endl;
  std::cout << "Si lat Angles = " << si.getLatRot().alpha << " " << si.getLatRot().beta
	    << " " << si.getLatRot().gamma << std::endl;
  std::cout << "Sf ref Angles = " << sf.getRefRot().alpha << " " << sf.getRefRot().beta
	    << " " << sf.getRefRot().gamma << std::endl;
  std::cout << "Sf lat Angles = " << sf.getLatRot().alpha << " " << sf.getLatRot().beta
	    << " " << sf.getLatRot().gamma << std::endl;

  for ( int mm = 1; mm >= -1; mm-=2 )
    {
      for ( int m = 1; m >= -1; m-=2 )
	{
	  std::complex<double> w_siRef = Hadron::Wigner_D(1,m,mm,si.getRefRot().alpha,
							  si.getRefRot().beta,si.getRefRot().gamma);
	  std::complex<double> w_siLat = Hadron::Wigner_D(1,m,mm,si.getLatRot().alpha,
							  si.getLatRot().beta,si.getLatRot().gamma);
	  std::complex<double> w_sfRef = Hadron::Wigner_D(1,m,mm,sf.getRefRot().alpha,
							  sf.getRefRot().beta,sf.getRefRot().gamma);
	  std::complex<double> w_sfLat = Hadron::Wigner_D(1,m,mm,sf.getLatRot().alpha,
							  sf.getLatRot().beta,sf.getLatRot().gamma);

	  
	  std::cout << "--------------------------" << std::endl;
	  std::cout << " (m,mm) = " << m << " , " << mm << std::endl;
	  std::cout << "     w_siRef = " << w_siRef << std::endl;
	  std::cout << "     w_siLat = " << w_siLat << std::endl;
	  std::cout << "     w_sfRef = " << w_sfRef << std::endl;
	  std::cout << "     w_sfLat = " << w_sfLat << std::endl;
	}
    }



  u1u_t idCHECK(true);
  // ABSOLUTE SPINORS ARE BRANDON KRIESTEN's 'HELICITY' SPINORS
  std::cout << "ID[11]" << idCHECK.eval(&sf.helicity.twoJz[1],
					&si.helicity.twoJz[1]) << std::endl;
  std::cout << "ID[12]" << idCHECK.eval(&sf.helicity.twoJz[1],
					&si.helicity.twoJz[-1]) << std::endl;
  std::cout << "ID[21]" << idCHECK.eval(&sf.helicity.twoJz[-1],
					&si.helicity.twoJz[1]) << std::endl;
  std::cout << "ID[22]" << idCHECK.eval(&sf.helicity.twoJz[-1],
					&si.helicity.twoJz[-1]) << std::endl;



  ugu_t g4CHECK(4,true);
  std::cout << "G4[11]" << g4CHECK.eval(&sf.helicity.twoJz[1],
					&si.helicity.twoJz[1]) << std::endl;
  std::cout << "G4[12]" << g4CHECK.eval(&sf.helicity.twoJz[1],
					&si.helicity.twoJz[-1]) << std::endl;
  std::cout << "G4[21]" << g4CHECK.eval(&sf.helicity.twoJz[-1],
					&si.helicity.twoJz[1]) << std::endl;
  std::cout << "G4[22]" << g4CHECK.eval(&sf.helicity.twoJz[-1],
					&si.helicity.twoJz[-1]) << std::endl;

  
  ugu_t gXCHECK(1,true);
  std::cout << "GX[11]" << gXCHECK.eval(&sf.helicity.twoJz[1],
					&si.helicity.twoJz[1]) << std::endl;
  std::cout << "GX[12]" << gXCHECK.eval(&sf.helicity.twoJz[1],
					&si.helicity.twoJz[-1]) << std::endl;
  std::cout << "GX[21]" << gXCHECK.eval(&sf.helicity.twoJz[-1],
					&si.helicity.twoJz[1]) << std::endl;
  std::cout << "GX[22]" << gXCHECK.eval(&sf.helicity.twoJz[-1],
					&si.helicity.twoJz[-1]) << std::endl;






  const int NUM_MATS           = 12;
  const int MATS_PER_INSERTION = 4;
  const int GPD_RANK           = 8;

  kinMatGPD_t GPD_4(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,4,disp);
  kinMatGPD_t GPD_1(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,1,disp);
  kinMatGPD_t GPD_2(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,2,disp);
  kinMatGPD_t GPD_3(MATS_PER_INSERTION,GPD_RANK,current::VECTOR,3,disp);
  // Assemble the kinematic matrices
  GPD_4.assemble(true,mass,&sf,&si);
  GPD_1.assemble(true,mass,&sf,&si);
  GPD_2.assemble(true,mass,&sf,&si);
  GPD_3.assemble(true,mass,&sf,&si);
  // Concatenate GPD_4,1,2 matrices into one large one for SVD
  Eigen::MatrixXcd GPD(3*MATS_PER_INSERTION, GPD_RANK);

  // Push GPD_4.mat, GPD_1.mat, GPD_2.mat into GPD
  for ( int i = 0; i < GPD_4.mat.rows(); ++i ) GPD.row(i) << GPD_4.mat.row(i);
  for ( int i = MATS_PER_INSERTION; i < 2*MATS_PER_INSERTION; ++i )
    GPD.row(i) << GPD_1.mat.row(i-MATS_PER_INSERTION);
  for ( int i = 2*MATS_PER_INSERTION; i < 3*MATS_PER_INSERTION; ++i )
    GPD.row(i) << GPD_2.mat.row(i-2*MATS_PER_INSERTION);
      

  std::cout << "FINAL GPD = " << GPD << std::endl;

  
  std::cout << "Ward ID Okay?" << std::endl;
  for ( int s = 0; s < 4; ++s )
    {
      std::complex<double> DIFF = GPD(s,0)*(Ef-Ei)-GPD(s+4,0)*(2*PI/Lx)*(mom_f[0]-mom_i[0])-
	GPD(s+8,0)*(2*PI/Lx)*(mom_f[1]-mom_i[1])-GPD_3.mat(s,0)*(2*PI/Lx)*(mom_f[2]-mom_i[2]);
      std::cout << "Ward (sdx = " << s << " ) = " << DIFF << std::endl;
    }



  ugu_t foo_x(1,true); ugu_t foo_y(2,true); ugu_t foo_z(3,true); ugu_t foo_4(4,true);

  std::cout << "Simple \gamma_x checks" << std::endl;
  std::cout << foo_x.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_x.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[-1]) << std::endl;
  std::cout << foo_x.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_x.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[-1]) << std::endl;
  std::cout << "Simple \gamma_y checks" << std::endl;
  std::cout << foo_y.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_y.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[-1]) << std::endl;
  std::cout << foo_y.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_y.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[-1]) << std::endl;
  std::cout << "Simple \gamma_z checks" << std::endl;
  std::cout << foo_z.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_z.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[-1]) << std::endl;
  std::cout << foo_z.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_z.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[-1]) << std::endl;
  std::cout << "Simple \gamma_4 checks" << std::endl;
  std::cout << foo_4.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_4.eval(&sf.subduced.twoJz[1],&si.subduced.twoJz[-1]) << std::endl;
  std::cout << foo_4.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[1]) << std::endl;
  std::cout << foo_4.eval(&sf.subduced.twoJz[-1],&si.subduced.twoJz[-1]) << std::endl;
  
  
  


  exit(3);


#if 0

  // K.assembleBig(mu,true,mass,&s,disp);

  // std::cout << "Kinematic matrix evaluates to: " << std::endl;
  // std::cout << K.mat << std::endl;

  // // Singular values
  // getSVs(&K.mat);


  // Hard code some svd check
  Eigen::MatrixXcd h(12,3);
  // Eigen::MatrixXcd h(8,3);
  // Eigen::MatrixXcd h(6,3);

  // p = (0,0,1)
  h(0,0) = cd(0,0);        h(0,1) = cd(0,0);       h(0,2) = cd(0,0);
  h(1,0) = cd(-1.13979,0); h(1,1) = cd(0,1.13979); h(1,2) = cd(-0.326235,0);
  h(2,0) = cd(-1.13979,0); h(2,1) = cd(0,1.13979); h(2,2) = cd(-0.326235,0);
  h(3,0) = cd(0,0);        h(3,1) = cd(0,0);       h(3,2) = cd(0,0);

  // p = (0,1,1)
  h(4,0) = cd(0,0);        h(4,1) = cd(0,0);       h(4,2) = cd(0,0);
  h(5,0) = cd(-0.852445,0.756604);  h(5,1)=cd(0.756604,0.852445);  h(5,2)=cd(-0.243991,0.216559);
  h(6,0) = cd(-0.852445,-0.756604); h(6,1)=cd(-0.756604,0.852445); h(6,2)=cd(-0.243991,-0.216559);
  h(7,0) = cd(0,0);        h(7,1) = cd(0,0);       h(7,2) = cd(0,0);

  // // p = (1,1,1)
  h(4,0) = cd(0.873651,0);   h(4,1) = cd(0,-0.873651); h(4,2) = cd(0.250061,0);
  h(5,0) = cd(-0.732015,0);  h(5,1) = cd(0,0.732015);  h(5,2) = cd(-0.209521,0);
  h(6,0) = cd(-0.732015,0);  h(6,1) = cd(0,0.732015);  h(6,2) = cd(-0.209521,0);
  h(7,0) = cd(-0.873651,0);  h(7,1) = cd(0,0.873651);  h(7,2) = cd(-0.250061,0);

  std::cout << "For big matrix, SVs = " << std::endl;
  getSVs(&h);


  // Hard code another
  Eigen::MatrixXcd hsmall(6,3);

  h(0,0) = cd(-1.13979,0); h(0,1) = cd(0,1.13979); h(0,2) = cd(-0.326235,0);
  h(1,0) = cd(-1.13979,0); h(1,1) = cd(0,1.13979); h(1,2) = cd(-0.326235,0);
  h(2,0) = cd(-0.852445,0.756604); h(2,1)=cd(0.756604,0.852445); h(2,2)=cd(-0.243991,0.216559);
  h(3,0) = cd(-0.852445,-0.756604); h(3,1)=cd(-0.756604,0.852445); h(3,2)=cd(-0.243991,-0.216559);
  h(4,0) = cd(-0.732015,0); h(4,1)=cd(0,0.732015); h(4,2)=cd(-0.209521,0);
  h(5,0) = cd(-0.732015,0);  h(5,1) = cd(0,0.732015); h(5,2) = cd(-0.209521,0);

  std::cout << "For small matrix, SVs = " << std::endl;
  getSVs(&hsmall);



  kinMatPDF_t V(2,8,current::VECTOR,mu,disp);
  kinMatPDF_t A(2,2,current::AXIAL,mu,disp);
  
  // kinMatPDF_t A(2,3,current::AXIAL,mu,disp);
  // kinMatPDF_t T(2,3,current::TENSOR,mu,4,disp);
  

  V.assemble(true,mass,&si,&si);
  std::cout << "Vector = " << V.mat << std::endl;

  A.assemble(true,mass,&si,&si);
  std::cout << "Axial = " << A.mat << std::endl;

  // T.assemble(true,mass,&s,&s);
  // std::cout << "Tensor = " << T.mat << std::endl;
  
  


  int gpdRank = 6; // 8

  kinMatGPD_t GPD(4,gpdRank,current::VECTOR,mu,disp);
  GPD.assemble(true,mass,&sf,&si);
  std::cout << "GPD = " << GPD.mat << std::endl;
  getSVs(&GPD.mat);


  kinMatGPD_t GPD_4(4,gpdRank,current::VECTOR,4,disp);
  kinMatGPD_t GPD_1(4,gpdRank,current::VECTOR,1,disp);
  kinMatGPD_t GPD_2(4,gpdRank,current::VECTOR,2,disp);
  GPD_4.assemble(true,mass,&sf,&si);
  GPD_1.assemble(true,mass,&sf,&si);
  GPD_2.assemble(true,mass,&sf,&si);

  // Try to concatenate GPD kinematic matrices
  Eigen::MatrixXcd bigGPD(12,gpdRank);
  std::cout << GPD_4.mat.row(0) << std::endl;
  std::cout << bigGPD.row(0) << std::endl;
  bigGPD.row(0) << GPD_4.mat.row(0);
  std::cout << bigGPD.row(0) << std::endl;
  for ( int i = 0; i < 4; ++i )  bigGPD.row(i) << GPD_4.mat.row(i);
  for ( int i = 4; i < 8; ++i )  bigGPD.row(i) << GPD_1.mat.row(i-4);
  for ( int i = 8; i < 12; ++i ) bigGPD.row(i) << GPD_2.mat.row(i-8);
  std::cout << "Big GPD = " << bigGPD << std::endl;
  getSVs(&bigGPD);

#endif

  return 0;
}
