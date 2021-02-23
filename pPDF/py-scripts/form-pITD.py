#!/dist/anaconda/bin/python

import sys,optparse
sys.path.append('/u/home/cegerer/src')

import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
# import scipy as scipy
import pylab # to save figures to file
from collections import OrderedDict

from pitd_util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-m", "--boostMom", type="str", default='',
                  help='Momentum of rpITD numerator X.X.X (default = '')')
parser.add_option("-n", "--normMom", type="str", default='',
                  help='Momentum of rpITD denominator X.X.X (default = '')')
parser.add_option("-z", "--zsep", type="int", default=-99,
                  help='Displacement of rpITD numerator (default = -99)')
parser.add_option("-I", "--insertion", type="int", default="-1",
                  help='Gamma matrix to consider (default = -1)')
parser.add_option("-c", "--cfgs", type="int", default=1,
                  help='Number of configurations (default = 1)')

# Parse the input arguments
(options, args) = parser.parse_args()

insertions = { 8: { 'Redstar': 'b_b0xDA__J0_A1pP', 'CK' : 'insertion_gt' } }


bx=int(options.boostMom.split('.')[0])
by=int(options.boostMom.split('.')[1])
bz=int(options.boostMom.split('.')[2])
nx=int(options.normMom.split('.')[0])
ny=int(options.normMom.split('.')[1])
nz=int(options.normMom.split('.')[2])
zsep=makeZSep(options.zsep)

nucBoostOp = 'NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1'
nucNormOp = 'NucleonMG1g1MxD0J0S_J1o2_G1g1'
if nx != 0 or ny != 0 or nz != 0:
    nucNormOp = nucBoostOp
if bx == 0 and by == 0 and bz == 0:
    nucBoostOp = 'NucleonMG1g1MxD0J0S_J1o2_G1g1'

current = 'b_b0xDA__J0_A1pP'

insOp = nptOp('tm3','fI2Y0i0',1,mom2Str(0,0,0),current,zsep); insOp.setName()
snkOpB = nptOp('tm60','fI1Y3i1',1,mom2Str(bx,by,bz),nucBoostOp); snkOpB.setName()
srcOpB = nptOp('tm3','fI1Y3i1',1,mom2Str(bx,by,bz),nucBoostOp); srcOpB.setName()
snkOpN = nptOp('tm60','fI1Y3i1',1,mom2Str(nx,ny,nz),nucNormOp); snkOpN.setName()
srcOpN = nptOp('tm3','fI1Y3i1',1,mom2Str(nx,ny,nz),nucNormOp); srcOpN.setName()
# Boost src op could change (i.e diff. tsrc stamp), so make a another for the no disp boosted correlator
srcOpBnd = nptOp('tm3','fI1Y3i1',1,mom2Str(bx,by,bz),nucBoostOp); srcOpBnd.setName()


# (12/03/2020) Apparently long-z pz = 4,5,6 data has tsrc time label as 'tm3' rather than 't0'
if bz > 3 and options.zsep > 8:
    srcOpB.t = 'tm3'; srcOpB.setName()
if nx > 0:
    srcOpN.t = 'tm3'; srcOpN.setName()

srcOpB.t = 'tm3'
srcOpN.t = 'tm3'

boost=correlator(snkOpB,insOp,srcOpB,'summedRatio','Summed-Matelems/%s/RES/'%current)
norm=correlator(snkOpN,insOp,srcOpN,'summedRatio','Summed-Matelems/%s/RES/'%current)
# Alter insertion to no displacement and make boosted/norm no disp corrs
insOp.disp = ''; insOp.setName()
boostNoDisp=correlator(snkOpB,insOp,srcOpBnd,'summedRatio','Summed-Matelems/%s/RES/'%current)
normNoDisp=correlator(snkOpN,insOp,srcOpN,'summedRatio','Summed-Matelems/%s/RES/'%current)

print "\n"
print boost.name
print boostNoDisp.name
print norm.name
print normNoDisp.name
print "\n"

# Big parse of all res files
for A in [boost, boostNoDisp, norm, normNoDisp]:
    A.resRead()
    for a, m in A.amplitude.items():
        m.parseParams()
        m.getParams()
        m.jk(); m.avgjk()

# Form the rpITD by combining fit results per jackknife sample of each amplitude
# NO JACKKNIFING HERE!!!
rpITD=pitd(boost,boostNoDisp,norm,normNoDisp)

# Take the jk ensemble averages and form the pITD
rpITD.fillJkAvgs()
rpITD.avg()
rpITD.err()

print "\nReal = %.7f +/- %.7f"%(float(rpITD.pITD[1]['avg']),float(rpITD.pITD[1]['err']))
print "Imag = %.7f +/- %.7f\n"%(float(rpITD.pITD[2]['avg']),float(rpITD.pITD[2]['err']))


# Define the ioffe time for this momentum/zsep combination
ioffeTime=ioffeTime(options.boostMom,options.zsep,32.0)
print "IoffeTime = %.4f"%ioffeTime


# Append asymptotic value of matelem + error + |zsep|
# W=open("pITD/NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1."+options.insertion+".pITD.dat",'a')
# W=open("TEST-NEW-NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1."+current+".pITD.dat",'a')
# W=open("TEST-NEW-XNORM-NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1."+current+".pITD.dat",'a')
# W=open("TEST-NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1."+current+".pITD.p001_norm.dat",'a')
# W=open("TEST-NEW_bugcorrected-NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1."+current+".pITD.dat",'a')
# for c in [1, 2]:
#     W.write(str(ioffeTime)+" "+str(c)+" "\
#             +str(float(rpITD.pITD[c]['avg']))+" "\
#             +str(float(rpITD.pITD[c]['err']))+" "+str(options.zsep)+"\n")
# W.close()


##############################
# INITIALIZE HDF5 STORAGE
##############################
#h5File = h5py.File('%s.h5'%options.insertion, 'w')
h5File = h5py.File('%s.CPE_L-summ_rpITD.h5'%insertions[options.insertion]['Redstar'], 'a')
 
h5Dirac = h5File.get(insertions[options.insertion]['Redstar'])
# If the insertion group doesn't exist, make it
if h5Dirac == None:
    h5Dirac = h5File.create_group(insertions[options.insertion]['Redstar'])
# Now have "/insertion" in h5 file...

h5Zsep = h5Dirac.get("zsep%s"%str(options.zsep))
# If this zsep group doesn't exist, make it
if h5Zsep == None:
    h5Zsep = h5Dirac.create_group("zsep%s"%str(options.zsep))

h5Mom = h5Zsep.get("pz%s"%str(bz))
# If this momentum group doesn't exist, make it
if h5Mom == None:
    h5Mom = h5Zsep.create_group("pz%s"%str(bz))
# Now have "/insertion/momz" in h5 file...

for n, component in enumerate([ "Re", "Im" ]):
    # Check to see if ensemble group exists, and make if none
    h5Ensem = h5Mom.get("ensemble/%s"%component)
    if h5Ensem == None:
        h5Ensem = h5Mom.create_group("ensemble/%s"%component)

    # Check to see if jack group exists, and make if none
    h5Jack = h5Mom.get("jack/%s"%component)
    if h5Jack == None:
        h5Jack = h5Mom.create_group("jack/%s"%component)

    # Now have "/insertion/ensemble/comp/zsep" in h5 file
    # Now have "/insertion/jack/comp/zsep" in h5 file


    # Create the ensem data structure
    #   -- holding 1 config; ioffeTime, avgM, errM
    h5EnsemMatelem = h5Ensem.create_dataset("pitd", (1,3), 'd')
    h5EnsemMatelem[0,0]=ioffeTime
    h5EnsemMatelem[0,1]=float(rpITD.pITD[n+1]['avg'])
    h5EnsemMatelem[0,2]=float(rpITD.pITD[n+1]['err'])

    # Create the jackknife data structure
    #   -- holding N configs; ioffetime, jackM
    h5JackMatelem = h5Jack.create_dataset("pitd", (options.cfgs,2), 'd')
    for g in range(0,options.cfgs):
        h5JackMatelem[g,0] = ioffeTime
        h5JackMatelem[g,1] = float(rpITD.pITD[n+1]['jks'][g])
    ############################################################################
