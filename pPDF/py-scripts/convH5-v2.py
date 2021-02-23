#!/dist/anaconda/bin/python

#####
# CONVERT CHRISTOS HDF5 FILE FORMAT INTO MY OWN
#
# where CK format is: "/bins/<L-summ/Plat>/mom_0_0_$p/disp_z$z/insertion_gt/<Re/Im>"
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
parser.add_option("-x", "--Nx", type="int", default=0,
                  help='Spatial extent of lattice (default = 0)')

parser.add_option("-H", "--christoFile", type="string", default="",
                  help='Christo H5 file (default = "")')
parser.add_option("-f", "--fittype", type="string", default="",
                  help='Fit scheme used to obtain rpITD <eg. Plat/L-summ > (default = '')')
parser.add_option("-I", "--insertion", type="int", default="-1",
                  help='Gamma matrix to consider (default = -1)')
parser.add_option("-z", "--zmax", type="int", default=8,
                  help='Max zsep in h5 (default = 8)')
parser.add_option("-p", "--momRange", type="string", default="",
                  help='Min/Max Momenta in h5 <pmin>.<pmax> (default = '')')


insertions = { 8: { 'Redstar': 'b_b0xDA__J0_A1pP', 'CK' : 'insertion_gt' } }



# Parse the input arguments
(options, args) = parser.parse_args()
gauge_configs=options.ensemble_cfgs

####################################
# ACCESS FILE HANDLE FOR CHRISTO'S H5 FILE
####################################
christoH5 = h5py.File(options.christoFile,'r')
# Grab just the jcakknife bins
cCollection = "bins"


##################################
# INITIALIZE HDF5 STORAGE TO MAKE
##################################
h5File = h5py.File('%s.CK_%s.all_data.h5'%(insertions[options.insertion]['Redstar'],
                                            options.fittype),'a')
# h5File = h5py.File('%s.CK_%s.mom_%s-%s.h5'%(insertions[options.insertion]['Redstar'],
#                                             options.fittype,options.momRange.split('.')[0],
#                                             options.momRange.split('.')[1]), 'a')

h5Dirac = h5File.get(insertions[options.insertion]['Redstar'])
# If the insertion group doesn't exist, make it
if h5Dirac == None:
    h5Dirac = h5File.create_group(insertions[options.insertion]['Redstar'])
# Now have "/insertion" in h5 file...


# for p in range(int(options.momRange.split('.')[0]),int(options.momRange.split('.')[1])):
#     pTag=str(p)
#     if p > 0:
#         pTag="+"+pTag


#     h5Mom = h5Dirac.get("pz%s"%str(p))
#     # If this momentum group doesn't exist, make it
#     if h5Mom == None:
#         h5Mom = h5Dirac.create_group("pz%s"%str(p))
#     # Now have "/insertion/momz" in h5 file...


#     for comp in [1, 2]:
#         compToStr=''
#         if comp == 1:
#             compToStr="Re"
#         if comp == 2:
#             compToStr="Im"

#         for zsep in range(0,options.zmax+1):
#             zTag=str(zsep)
#             if zsep > 0:
#                 zTag="z+%s"%str(zsep)

#             ioffeTime=((2*np.pi)/(1.0*options.Nx))*p*zsep


#             # Check to see if jack group exists, and make it none
#             h5Jack = h5Mom.get("jack/%s/zsep%s"%(compToStr,str(zsep)))
#             if h5Jack == None:
#                 h5Jack = h5Mom.create_group("jack/%s/zsep%s"%(compToStr,str(zsep)))
#             # Now have "/insertion/pz$p/jack/comp/zsep" in h5 file

            
#             # Create the jackknife data structure
#             #   -- holding N configs; ioffetime, jkM
#             h5JackData = h5Jack.create_dataset("pitd", (gauge_configs,2), 'd')

#             print "/%s/%s/mom_0_0_%s/disp_%s/%s/%s"%(cCollection,options.fittype,pTag,
#                                                      zTag,insertions[options.insertion]['CK'],
#                                                      compToStr)

#             for g in range(0,gauge_configs):
#                 h5JackData[g,0] = ioffeTime
#                 # h5JackData[g,1] = christoH5['/%s/%s/mom_0_0_%s/disp_%s/%s/%s'\
#                 #                              %(cCollection,options.fittype,pTag,
#                 #                                zTag,insertions[options.insertion]['CK'],
#                 #                                compToStr)][g]
#                 h5JackData[g,1] = christoH5['/%s/Pz_%d/z_%d/%s/jack'\
#                                              %(options.fittype,p,
#                                                zsep,compToStr)][g]
############################################################################


for zsep in range(0,options.zmax+1):
    zTag=str(zsep)
    if zsep > 0:
        zTag="z+%s"%str(zsep)

    h5Z = h5Dirac.get("zsep%s"%str(zsep))
    if h5Z == None:
        h5Z = h5Dirac.create_group("zsep%s"%str(zsep)) 


    for p in range(int(options.momRange.split('.')[0]),int(options.momRange.split('.')[1])+1):
        pTag=str(p)
        if p > 0:
            pTag="+"+pTag

        h5Mom = h5Z.get("pz%s"%str(p))
        # If this momentum group doesn't exist, make it
        if h5Mom == None:
            h5Mom = h5Z.create_group("pz%s"%str(p))
            

        ioffeTime=((2*np.pi)/(1.0*options.Nx))*p*zsep 

            
        for comp in [1, 2]:
            compToStr=''
            if comp == 1:
                compToStr="Re"
            if comp == 2:
                compToStr="Im"


            # Check to see if jack group exists, and make it none
            h5Jack = h5Mom.get("jack/%s"%compToStr)
            if h5Jack == None:
                h5Jack = h5Mom.create_group("jack/%s"%compToStr)
            # Now have "/insertion/pz$p/jack/comp/zsep" in h5 file

            
            # Create the jackknife data structure
            #   -- holding N configs; ioffetime, jkM
            h5JackData = h5Jack.create_dataset("pitd", (gauge_configs,2), 'd')

            # print "/%s/%s/mom_0_0_%s/disp_%s/%s/%s"%(cCollection,options.fittype,pTag,
            #                                          zTag,insertions[options.insertion]['CK'],
            #                                          compToStr)

            for g in range(0,gauge_configs):
                h5JackData[g,0] = ioffeTime
                h5JackData[g,1] = christoH5['/%s/%s/mom_0_0_%s/disp_%s/%s/%s'\
                                             %(cCollection,options.fittype,pTag,
                                               zTag,insertions[options.insertion]['CK'],
                                               compToStr)][g]
                # h5JackData[g,1] = christoH5['/%s/Pz_%d/z_%d/%s/jack'\
                #                              %(options.fittype,p,
                #                                zsep,compToStr)][g]
