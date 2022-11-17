#!/bin/bash

source ~/arch/jlab/chromaform/env.sh
export LD_LIBRARY_PATH=/u/home/cegerer/gsl-2.6/lib:${LD_LIBRARY_PATH}

#eval `modulecmd bash load hdf5-1.8.12`
export PATH=/u/home/cegerer/hdf5-1.12.1/bin:${PATH}
export LD_LIBRARY_PATH=/u/home/cegerer/hdf5-1.12.1/lib:${LD_LIBRARY_PATH}
