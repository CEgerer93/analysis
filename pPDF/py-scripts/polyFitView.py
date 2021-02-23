#!/dist/anaconda/bin/python

#####
# PLOT JACK pITD DATA STORED IN AN H5 FILE
# TOGETHER WITH POLYNOMIAL FITS
#####

import numpy as np
import h5py
import matplotlib.pyplot as plt
import pylab # to save figures to file
import sys,optparse
from collections import OrderedDict
from pitd_util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-H", "--h5File", type="string", default="",
                  help='H5 file to plot (default = "")')
parser.add_option("-P", "--polyFitFile", type="string", default="",
                  help='Polynomial fits file (default = "")')
parser.add_option("-z", "--zRange", type="string", default="",
                  help='Min/Max zsep in h5 <zmin>.<zmax> (default = '')')
parser.add_option("-p", "--momRange", type="string", default="",
                  help='Min/Max Momenta in h5 <pmin>.<pmax> (default = '')')
parser.add_option("-I", "--insertion", type="int", default="-1",
                  help='Gamma matrix to consider (default = -1)')
parser.add_option("-c", "--cfgs", type="int", default=1,
                  help='Number of configurations (default = 1)')


insertions = { 8: { 'Redstar': 'b_b0xDA__J0_A1pP', 'CK' : 'insertion_gt' } }
labels = { 'evoKernel': { 'Re': r'$\mathfrak{Re} B\otimes\mathfrak{M}\left(\nu,z^2\right)$',
                          'Im': r'$\mathfrak{Im} B\otimes\mathfrak{M}\left(\nu,z^2\right)$' },
           'matchingKernel' : { 'Re': r'$\mathfrak{Re} L\otimes\mathfrak{M}\left(\nu,z^2\right)$',
                                'Im': r'$\mathfrak{Im} L\otimes\mathfrak{M}\left(\nu,z^2\right)$' },
           'pitd': {'Re': '', 'Im': ''} }


# Parse the input arguments
(options, args) = parser.parse_args()
cfgs=options.cfgs
zmin=int(options.zRange.split('.')[0])
zmax=int(options.zRange.split('.')[1])
pmin=int(options.momRange.split('.')[0])
pmax=int(options.momRange.split('.')[1])

sizeFont=16

# Instantiate some figures
plt.rcParams["mathtext.fontset"]="stix"
fig_real = plt.figure()
fig_imag = plt.figure()
axesR = fig_real.gca()
axesI = fig_imag.gca()
for n, f in enumerate([axesR, axesI]):
    f.set_xlabel(r'$\nu$',fontsize=sizeFont)
axesR.set_ylabel(r'$\mathfrak{Re} \mathfrak{M}\left(\nu,z^2\right)$',fontsize=sizeFont)
axesI.set_ylabel(r'$\mathfrak{Im} \mathfrak{M}\left(\nu,z^2\right)$',fontsize=sizeFont)
axesR.set_ylim([-0.5,1.1])
axesI.set_ylim([-0.5,1.1])
axesR.axhline(y=0.0,ls=':',color='gray')
axesI.axhline(y=0.0,ls=':',color='gray')

####################################
# ACCESS FILE HANDLE(S)
####################################
for h5file in [options.h5File]:
    h5In = h5py.File(h5file,'r')
    dtypeName = "pitd"


    for z in range(zmin,zmax+1):
        ztag="zsep%d"%z
        for m in range(pmin,pmax+1):
            ptag="pz%d"%m

            for comp in ["Re", "Im"]:

                ioffeTime = -1
                avgMat = 0.0
                avgMatErr = 0.0
                for g in range(0,cfgs):
                    ioffeTime, mat = h5In['/%s/%s/%s/jack/%s/%s'%(insertions[options.insertion]['Redstar'],\
                                                                  ztag,ptag,comp,dtypeName)][g]
                    avgMat += mat
                    
                avgMat *= (1.0/cfgs)
                    
                for g in range(0,cfgs):
                    ioffeTime, mat = h5In['/%s/%s/%s/jack/%s/%s'%(insertions[options.insertion]['Redstar'],\
                                                                  ztag,ptag,comp,dtypeName)][g]
                    avgMatErr += np.power( mat - avgMat, 2)

                avgMatErr = np.sqrt( ((1.0*(cfgs-1))/cfgs)*avgMatErr )


                if comp == "Re":
                    axesR.errorbar(ioffeTime, avgMat, yerr=avgMatErr, fmt='o',\
                                   color=mainColors[z],mec=mainColors[z],mfc=mainColors[z],label="z=%s"%z)
                if comp == "Im":
                    axesI.errorbar(ioffeTime, avgMat, yerr=avgMatErr, fmt='o',\
                                   color=mainColors[z],mec=mainColors[z],mfc=mainColors[z],label="z=%s"%z)



########################################
# ACCESS THE POLY FIT RESULTS
#
# <zsep> <jk> <real/imag> <a> <b> <c> <rchi2>
#
########################################
# Big dictionary for all fit params
poly={z: {0 : {'jks': []}, 1: {'jks': []}} for z in range(1,17)}


with open(options.polyFitFile) as ptr:
    for cnt, line in enumerate(ptr):
        # Capture lines; remove spaces
        z=int(line.split(' ')[0])
        jk=int(line.split(' ')[1])
        comp=int(line.split(' ')[2])
        A=float(line.split(' ')[3])
        B=float(line.split(' ')[4])
        C=float(line.split(' ')[5])
        chi2=float(line.split(' ')[6].rstrip())

        params=[p for p in [A, B, C, chi2]]

        # Pack the params
        poly[z][comp]['jks'].append(params)

        
# Make dictionary of polyFit classes
thePolyFits={ z: polyFit(poly[z][0], poly[z][1], z) for z in range(1,17)}

# Determine the average and covariances of all fit parameters
for z in range(1,17):
    thePolyFits[z].getAvgParams()
    thePolyFits[z].getParamCov()

# # A verbose check
# for n,v in thePolyFits[2].dumPoly.paramOrder.items():
#     print "%.6f +/- %.6f"%(thePolyFits[2].avgParams[0][v],np.sqrt(thePolyFits[2].cov[0][(v,v)]))



########################
# Make the fit bands
########################
for z in range(zmin,zmax+1):
    thePolyFits[z].plotFit(axesR,axesI)



# Now the real plot labels
handles, labels = axesR.get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
fig_real.legend(by_label.values(), by_label.keys(),fontsize=12,loc='right')
# Now the imag plot labels
handles, labels = axesI.get_legend_handles_labels() 
by_label = OrderedDict(zip(labels, handles))
fig_imag.legend(by_label.values(), by_label.keys(),fontsize=12,loc='right')

fig_real.savefig("polyFitView_real.png", dpi=400,bbox_includes='tight',pad_inches=0)
fig_imag.savefig("polyFitView_imag.png", dpi=400,bbox_includes='tight',pad_inches=0)
plt.show()
