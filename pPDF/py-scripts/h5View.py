#!/dist/anaconda/bin/python

#####
# PLOT JACK DATA STORED IN AN H5 FILE
#####

import numpy as np
import h5py
import matplotlib.pyplot as plt
import pylab # to save figures to file
import sys,optparse
from collections import OrderedDict
from pitd_util import mainColors

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-H", "--h5File", type="string", default="",
                  help='H5 file(s) to plot <x>:<x> (default = "x:x")')
parser.add_option("-d", "--dtypeName", type="string", default="",
                  help='Datatype name(s) to access <d>:<d> (default = "d:d")')
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


# Prefactor of convolution
# prefactor = (0.303*4/3)/(2*np.pi)
prefactor = 1.0


# Parse the input arguments
(options, args) = parser.parse_args()
cfgs=options.cfgs
zmin=int(options.zRange.split('.')[0])
zmax=int(options.zRange.split('.')[1])
pmin=int(options.momRange.split('.')[0])
pmax=int(options.momRange.split('.')[1])

# Instantiate some figures
plt.rcParams["mathtext.fontset"]="stix"
fig_real, axesR = plt.subplots(len(options.h5File.split(':')),1,sharex=True) # ,sharey=True)
fig_imag, axesI = plt.subplots(len(options.h5File.split(':')),1,sharex=True) # ,sharey=True)
fig_real.subplots_adjust(hspace=0.05)
fig_imag.subplots_adjust(hspace=0.05)
# fig_real.set_xlabel(r'$\nu$',fontsize=16)
# fig_imag.xlabel(r'$\nu$',fontsize=16)
if len(fig_real.get_axes()) > 1:
    axesR[-1].set_xlabel(r'$\nu$',fontsize=16)
    axesI[-1].set_xlabel(r'$\nu$',fontsize=16)
else:
    axesR.set_xlabel(r'$\nu$',fontsize=16)
    axesI.set_xlabel(r'$\nu$',fontsize=16)

    



####################################
# ACCESS FILE HANDLE(S)
####################################
for nH5, h5file in enumerate(options.h5File.split(':')):
    h5In = h5py.File(h5file,'r')
    dtypeName = options.dtypeName.split(':')[nH5]


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
                    avgMat += mat*prefactor
                    
                avgMat *= (1.0/cfgs)
                    
                for g in range(0,cfgs):
                    ioffeTime, mat = h5In['/%s/%s/%s/jack/%s/%s'%(insertions[options.insertion]['Redstar'],\
                                                                  ztag,ptag,comp,dtypeName)][g]
                    avgMatErr += np.power( mat*prefactor - avgMat, 2)

                avgMatErr = np.sqrt( ((1.0*(cfgs-1))/cfgs)*avgMatErr )


                if comp == "Re":
                    if len(fig_real.get_axes()) > 1:
                        axesR[nH5].errorbar(ioffeTime, avgMat, yerr=avgMatErr, fmt='o',\
                                            color=mainColors[z],mec=mainColors[z],mfc=mainColors[z],label="z=%s"%z)
                    else:
                        axesR.errorbar(ioffeTime, avgMat, yerr=avgMatErr, fmt='o',\
                                       color=mainColors[z],mec=mainColors[z],mfc=mainColors[z],label="z=%s"%z)
                if comp == "Im":
                    if len(fig_real.get_axes()) > 1:
                        axesI[nH5].errorbar(ioffeTime, avgMat, yerr=avgMatErr, fmt='o',\
                                            color=mainColors[z],mec=mainColors[z],mfc=mainColors[z],label="z=%s"%z)
                    else:
                        axesI.errorbar(ioffeTime, avgMat, yerr=avgMatErr, fmt='o',\
                                       color=mainColors[z],mec=mainColors[z],mfc=mainColors[z],label="z=%s"%z)


for n, d in enumerate(options.dtypeName.split(':')):
    if len(fig_real.get_axes()) > 1:
        axesR[n].set_ylabel(labels[d]['Re'],fontsize=16)
        axesI[n].set_ylabel(labels[d]['Im'],fontsize=16)

        axesR[n].axhline(y=0.0,ls=':',color='gray')
        axesI[n].axhline(y=0.0,ls=':',color='gray')

        yrangeR = axesR[n].get_ylim()[1] - axesR[n].get_ylim()[0]
        yrangeI = axesI[n].get_ylim()[1] - axesI[n].get_ylim()[0]
    
        axesR[n].set_ylim([axesR[n].get_ylim()[0]-yrangeR*0.05, axesR[n].get_ylim()[1]+yrangeR*0.05])
        axesI[n].set_ylim([axesI[n].get_ylim()[0]-yrangeI*0.05, axesI[n].get_ylim()[1]+yrangeI*0.05])

    else:
        axesR.set_ylabel(labels[d]['Re'],fontsize=16)
        axesI.set_ylabel(labels[d]['Im'],fontsize=16)

        axesR.axhline(y=0.0,ls=':',color='gray')
        axesI.axhline(y=0.0,ls=':',color='gray')

        yrangeR = axesR.get_ylim()[1] - axesR.get_ylim()[0]
        yrangeI = axesI.get_ylim()[1] - axesI.get_ylim()[0]

        axesR.set_ylim([axesR.get_ylim()[0]-yrangeR*0.05, axesR.get_ylim()[1]+yrangeR*0.05])
        axesI.set_ylim([axesI.get_ylim()[0]-yrangeI*0.05, axesI.get_ylim()[1]+yrangeI*0.05])

        axesR.set_ylim([-0.5,1.1])
        axesI.set_ylim([-0.5,1.1])


# axesR[-1].set_ylim([axesR[0].get_ylim()[0], axesR[0].get_ylim()[1]])
# axesI[-1].set_ylim([axesI[0].get_ylim()[0], axesI[0].get_ylim()[1]])


if len(fig_real.get_axes()) > 1:                
    handles, labels = axesR[-1].get_legend_handles_labels()
else:
    handles, labels = axesR.get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
fig_real.legend(by_label.values(), by_label.keys(),fontsize=12,loc='right')
if len(fig_real.get_axes()) > 1:
    handles, labels = axesI[-1].get_legend_handles_labels()
else:
    handles, labels = axesI.get_legend_handles_labels() 
by_label = OrderedDict(zip(labels, handles))
fig_imag.legend(by_label.values(), by_label.keys(),fontsize=12,loc='right')

fig_real.savefig("viewh5_real.png", dpi=400,bbox_includes='tight',pad_inches=0)
fig_imag.savefig("viewh5_imag.png", dpi=400,bbox_includes='tight',pad_inches=0)
plt.show()
