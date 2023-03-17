CXX=icpc

EIGEN=/u/home/cegerer/src/eigen

# -DAMPPREFACTORS
# -DHAVE_DISPLIST
# -DPRINT_SVS

# -DMYDECOMP
# -DMYDECOMP2
# -DARDECOMP
# -DMARTHADECOMP
# -DDIAGMATS
# -DROWONE

CXXFLAGS=-g -qopenmp -fma -falias -std=c++14 -DMYDECOMP2 -DPRINT_SVS -DHAVE_DISPLIST -DVERBOSITY=1 -I/dist/gmp-6.0.0/include -I/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/adat/include -I/usr/include/libxml2 -I/home/cegerer/analysis/pPDF -I/home/cegerer/analysis/fitting -I/home/cegerer/analysis/linAlg -I/home/cegerer/analysis/kine -I/home/cegerer/analysis/corrUtil -I/u/home/cegerer/gsl-2.6/include -I/u/home/cegerer/hdf5-1.12.1/include -I$(EIGEN) -I/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/redstar/include -I/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/colorvec/include -I/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/hadron/include -I/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/tensor/include

LDFLAGS=-L/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/adat/lib -L/u/home/cegerer/gsl-2.6/lib -L/home/cegerer/analysis/pPDF -L/u/home/cegerer/hdf5-1.12.1/lib -L/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/hadron/lib -L/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/colorvec/lib -L/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/redstar/lib -L/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/tensor/lib

LDLIBS=-lredstar -ladat -lcolorvec -lhadron -ltensor -lfiledb -lfilehash -lXPathReader -lxmlWriter -lMinuit2Base -lxml2 -lz -lm -ldl /dist/gmp-6.0.0/lib/libgmpxx.a /dist/gmp-6.0.0/lib/libgmp.a -lgsl -lgslcblas -lhdf5 -lhdf5_cpp -lhdf5_hl -lhdf5_hl_cpp


# CXXFLAGS=-g -qopenmp -fma -falias -std=c++14 -DPRINT_SVS -DHAVE_DISPLIST -DVERBOSITY=1 -I/dist/gmp-6.0.0/include -I/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/adat/include -I/usr/include/libxml2 -I/home/cegerer/analysis/pPDF -I/home/cegerer/analysis/fitting -I/home/cegerer/analysis/linAlg -I/home/cegerer/analysis/kine -I/home/cegerer/analysis/corrUtil -I/u/home/cegerer/gsl-2.6/include -I/u/home/cegerer/hdf5-1.12.1/include -I$(EIGEN)

# LDFLAGS=-L/home/cegerer/arch/jlab/scalar-knl_PDFs-no-AVX512/install/adat/lib -L/u/home/cegerer/gsl-2.6/lib -L/home/cegerer/analysis/pPDF -L/u/home/cegerer/hdf5-1.12.1/lib

# LDLIBS=-ladat -lfiledb -lfilehash -lXPathReader -lxmlWriter -lMinuit2Base -lxml2 -lz -lm -ldl /dist/gmp-6.0.0/lib/libgmpxx.a /dist/gmp-6.0.0/lib/libgmp.a -lgsl -lgslcblas -lhdf5 -lhdf5_cpp -lhdf5_hl -lhdf5_hl_cpp




pGITD-matelem: pGITD-matelem.o operators.o pseudo_utils.o fit_util.o corr_utils.o cov_utils.o varpro.o threept_tr.o rotations.o old_irrep_util_angles.o
	$(CXX) $(CXXFLAGS) -o pGITD-matelem pGITD-matelem.o operators.o pseudo_utils.o fit_util.o corr_utils.o cov_utils.o varpro.o threept_tr.o rotations.o old_irrep_util_angles.o $(LDFLAGS) $(LDLIBS)

pGITD-matelem.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/matelems/pGITD/pGITD-matelem.cc

operators.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/pPDF/operators.cc

pseudo_utils.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/pPDF/pseudo_utils.cc

fit_util.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/fitting/fit_util.cc

varpro.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/fitting/varpro.cc

corr_utils.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/corrUtil/corr_utils.cc

cov_utils.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/linAlg/cov_utils.cc

threept_tr.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/kine/threept_tr.cc

rotations.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/kine/rotations.cc

old_irrep_util_angles.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/kine/old_irrep_util_angles.cc
clean:
	rm -f pGITD-matelem pGITD-matelem.o threept_tr.o
clean-all:
	rm -f pGITD-matelem pGITD-matelem.o personal_fns.o operators.o pseudo_utils.o fit_util.o corr_utils.o cov_utils.o varpro.o threept_tr.o rotations.o old_irrep_util_angles.o
