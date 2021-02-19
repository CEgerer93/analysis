CXX=icpc

CXXFLAGS=-Wall -Wextra -g -fsanitize=undefined -qopenmp -D_REENTRANT -O3 -finline-functions -fma -falias -std=c++14 -xAVX -qopt-report=0 -qopt-report-phase=vec -restrict -I/dist/gmp-6.0.0/include -I/volatile/JLabLQCD/cegerer/qcd_software/jlab/scalar-knl_PDFs/install/adat/include -I/usr/include/libxml2 -I/dist/gsl-1.15/include/gsl -I/home/cegerer/src -I/home/cegerer/analysis/fitting


LDFLAGS=-L/volatile/JLabLQCD/cegerer/qcd_software/jlab/scalar-knl_PDFs/install/adat/lib -L/dist/gsl-1.15/lib

LDLIBS=-ladat -lfiledb -lfilehash -lXPathReader -lxmlWriter -lMinuit2Base -lxml2 -lz -lm -ldl /dist/gmp-6.0.0/lib/libgmpxx.a /dist/gmp-6.0.0/lib/libgmp.a -lgsl -lgslcblas


summation-matelem: summation-matelem.o summation.o personal_fns.o operators.o fit_util.o
	$(CXX) $(CXXFLAGS) -o summation-matelem summation-matelem.o summation.o personal_fns.o operators.o fit_util.o $(LDFLAGS) $(LDLIBS)

summation-matelem.o:
	$(CXX) $(CXXFLAGS) -c summation-matelem.cc

summation.o:
	$(CXX) $(CXXFLAGS) -c summation.cc

operators.o:
	$(CXX) $(CXXFLAGS) -c operators.cc

personal_fns.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/src/personal_fns.cc

fit_util.o:
	$(CXX) $(CXXFLAGS) -c /home/cegerer/analysis/fitting/fit_util.cc

clean:
	rm -f summation-matelem summation-matelem.o summation.o personal_fns.o fit_util.o
