#!/dist/anaconda/bin/python

#####
# CONVERT CHRISTOS HDF5 FILE FORMAT INTO MY OWN
#####

import numpy as np
import h5py
import matplotlib.pyplot as plt
import pylab # to save figures to file
import sys,optparse
from collections import OrderedDict
#from eff_charge_util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-c", "--ensemble_cfgs", type="int", default=1,
                  help='Number of configurations (default = 1)')


parser.add_option("-H", "--christoFile", type="string", default="",
                  help='Christo H5 file (default = "")')
parser.add_option("-C", "--christoInsertion", type="string", default="",
                  help='Insertion in Christo H5 file (default = "")')
parser.add_option("-I", "--insertion", type="string", default="",
                  help='Gamma matrix to consider (default = "")')



# Parse the input arguments
(options, args) = parser.parse_args()
gauge_configs=options.ensemble_cfgs

####################################
# ACCESS FILE HANDLE FOR CHRISTO'S H5 FILE
####################################
christoH5 = h5py.File(options.christoFile,'r')
cIns = options.christoInsertion



##############################
# INITIALIZE HDF5 STORAGE
##############################
#h5File = h5py.File('%s.h5'%options.insertion, 'w')
h5File = h5py.File('%s_%s.h5'%(options.insertion,cIns), 'a')

h5Dirac = h5File.get(options.insertion)
# If the insertion group doesn't exist, make it
if h5Dirac == None:
    h5Dirac = h5File.create_group(options.insertion)
# Now have "/insertion" in h5 file...


for p in range(1,7):
    h5Mom = h5Dirac.get("pz%s"%str(p))
    # If this momentum group doesn't exist, make it
    if h5Mom == None:
        h5Mom = h5Dirac.create_group("pz%s"%str(p))
    # Now have "/insertion/momz" in h5 file...


    for comp in [1, 2]:
        compToStr=''
        if comp == 1:
            compToStr="Re"
        if comp == 2:
            compToStr="Im"

        for zsep in range(0,9):
            ioffeTime=((2*np.pi)/(32.0))*p*zsep

            # Check to see if ensemble group exists, and make it none
            h5Ensem = h5Mom.get("ensemble/%s/%s"%(str(comp),str(zsep)))
            if h5Ensem == None:
                h5Ensem = h5Mom.create_group("ensemble/%s/%s"%(str(comp),str(zsep)))


            # Check to see if jack group exists, and make it none
            h5Jack = h5Mom.get("jack/%s/%s"%(str(comp),str(zsep)))
            if h5Jack == None:
                h5Jack = h5Mom.create_group("jack/%s/%s"%(str(comp),str(zsep)))



            # Now have "/insertion/ensemble/comp/zsep" in h5 file
            # Now have "/insertion/jack/comp/zsep" in h5 file
            
            # Create the ensem data structure
            #   -- holding 1 config; ioffeTime, avgM, errM
            h5EnsemData = h5Ensem.create_dataset("pitd", (1,3), 'd')
            
            # Create the jackknife data structure
            #   -- holding N configs; ioffetime, avgM, errM
            h5JackData = h5Jack.create_dataset("pitd", (gauge_configs,3), 'd')


            h5EnsemData[0,0]=ioffeTime
            h5EnsemData[0,1]=christoH5['/%s/Pz_%d/z_%d/%s/%s'\
                                       %(cIns,p,zsep,compToStr,"ave")][0]
            h5EnsemData[0,2]=christoH5['/%s/Pz_%d/z_%d/%s/%s'\
                                       %(cIns,p,zsep,compToStr,"ave")][1]


            for g in range(0,gauge_configs):
                h5JackData[g,0] = ioffeTime
                h5JackData[g,1] = christoH5['/%s/Pz_%d/z_%d/%s/%s'\
                                             %(cIns,p,zsep,compToStr,"jack")][g]
                # Full data covariance that will be determined from jackknife samples
                # provides best estimate of data covariance, and will thus be used
                # as covariance within each jackknife bin
                # So, set jackknife error to EnsemData error
                h5JackData[g,2] = christoH5['/%s/Pz_%d/z_%d/%s/%s'\
                                            %(cIns,p,zsep,compToStr,"ave")][1]
############################################################################


