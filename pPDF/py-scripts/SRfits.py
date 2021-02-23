#!/dist/anaconda/bin/python

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
parser.add_option("-x", "--complexity", type="int", default=0,
                  help='Component of Npt function to analyze [Real=1 Imag=2] (default = 1)')
parser.add_option("-t", "--numTSeps", type="int", default=0,
                  help='Number of tseps present in SR file (default = 0)')
parser.add_option("-z", "--zsep", type="string", default='0',
                  help='Z Displacement of Wilson line (default = "0")')
parser.add_option("-p", "--mom", type="str", default="X.X.X",
                  help='Momentum (default = X.X.X)')
parser.add_option("-n", "--SRCorrFile", type="string", default="",
                  help='File containing full path/name for SR correlator (default="")')
parser.add_option("-I", "--insertion", type="string", default="",
                  help='Gamma matrix to consider (default = "")')
parser.add_option("-r", "--resFile", type="string", default="",
                  help='Name of fit results file (default = "")')
parser.add_option("-w", "--dataTMin", type="int", default="1",
                  help='Min tsep in SR file (default = 0)')
parser.add_option("-y", "--dataTMax", type="int", default="1",
                  help='Max tsep in SR file (default = 1)')
parser.add_option("-E", "--excludedTFits", type="string", default="",
                  help='One-based Times excluded from fit tseries (default="")')
parser.add_option("-H", "--makeH5", type="int", default=0,
                  help='Make hdf5 of matelem results (default = 0)')
parser.add_option("-l", "--lightBkgd", type="int", default=1,
                  help='Format figs for light (1) or dark (0) background (default = 1)')
parser.add_option("-R", "--showFitRes", type="int", default=1,
                  help='Show fit results on plot (default = 1)')

# Parse the input arguments
(options, args) = parser.parse_args()

################
# INITIALIZE GLOBAL PROPERTIES OF FIGURES
################
# Finalize the figures
plt.rc('text', usetex=True, fontsize=18)
plt.rcParams["mathtext.fontset"]="stix"
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rc('xtick.major',size=5)
plt.rc('ytick.major',size=5)
plt.rc('xtick',labelsize=18)
plt.rc('ytick',labelsize=18)
plt.rc('axes',labelsize=18)
truthTransparent=False
FrameAlpha=0.8
# Optionally swap default black labels for white
if options.lightBkgd == 0:
    truthTransparent=True
    plt.rcParams['text.color'] = 'white'
    plt.rcParams['axes.edgecolor'] = 'white'
    plt.rc('axes',edgecolor='white')
    plt.rc('axes',labelcolor='white')
    plt.rc('xtick',color='white')
    plt.rc('ytick',color='white')
    plt.rc('text',color='white')
    FrameAlpha=0



gauge_configs=options.ensemble_cfgs
complexity=options.complexity
# gamma=options.chromaGamma
paramsFile=options.resFile

component=""
if options.complexity == 1:
    component="Re"
if options.complexity == 2:
    component="Im"

momx=int(options.mom.split('.')[0])
momy=int(options.mom.split('.')[1])
momz=int(options.mom.split('.')[2])


# Define the ioffe time for this momentum/zsep combination
ioffeTime=((2*np.pi)/(32.0))*(momz)*int(options.zsep)


# Access the summed ratio file
dataSR=np.loadtxt(options.SRCorrFile, delimiter=' ', skiprows=1)
# Set the generated summed ratio plot
output_name="%s_comp%i"%(options.SRCorrFile.replace('.dat',''),options.complexity)



if options.makeH5:
    ##############################
    # INITIALIZE HDF5 STORAGE
    ##############################
    #h5File = h5py.File('%s.h5'%options.insertion, 'w')
    h5File = h5py.File('%s.CPE_L-summ_fits.h5'%options.insertion, 'a')

    h5Dirac = h5File.get(options.insertion)
    # If the insertion group doesn't exist, make it
    if h5Dirac == None:
        h5Dirac = h5File.create_group(options.insertion)
    # Now have "/insertion" in h5 file...

    h5Zsep = h5Dirac.get("zsep%s"%str(options.zsep))
    # If this zsep group doesn't exist, make it
    if h5Zsep == None:
        h5Zsep = h5Dirac.create_group("zsep%s"%str(options.zsep))

    h5Mom = h5Zsep.get("pz%s"%str(momz))
    # If this momentum group doesn't exist, make it
    if h5Mom == None:
        h5Mom = h5Zsep.create_group("pz%s"%str(momz))
    # Now have "/insertion/momz" in h5 file...

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
    # Create the ensem intercept data structure
    h5EnsemIntercept = h5Ensem.create_dataset("intercept", (1,3), 'd')
    h5EnsemIntercept[0,0] = ioffeTime

    # Create the jackknife data structure
    #   -- holding N configs; ioffetime, jackM
    h5JackMatelem = h5Jack.create_dataset("pitd", (gauge_configs,2), 'd')
    h5JackIntercept = h5Jack.create_dataset("intercept", (gauge_configs,2), 'd')
    for g in range(0,gauge_configs):
        h5JackMatelem[g,0] = ioffeTime
        h5JackIntercept[g,0] = ioffeTime
    ############################################################################



# The Fit function to SR data
# def SR(A,B,C,dE,T):
#     # return A+B*np.exp(-dE*T)+C*T*np.exp(-dE*T)
#     return B+A*T+C*T*np.exp(-dE*T)

def SR(A,B,T):
    return A+B*T

####################
# DEFINE A DICTIONARY TO HOLD ORDERING OF KEY/VALUES
####################
# paramOrder={0: 'A', 1: 'B', 2: 'C', 3: 'dE'}
paramOrder={0: 'A', 1: 'B', 2: 'rChi2'}


#####################
# GET THE FITTED PARAMETERS AND ERRORS (TO JACKKNIFE ENSEMBLES)
#####################
numFitParams=len(paramOrder)
print "Parsing %i distinct fit parameters from file"%numFitParams
params={}
PARAMS=np.loadtxt(paramsFile,dtype="string")
for k,v in paramOrder.items():
    avgFitParam=0.0
    for g in range(0,gauge_configs):
        avgFitParam+=float(PARAMS[k+numFitParams*g][1])

        if options.makeH5:
            # Only add the intercept & asymptotic value of matrix element (i.e. fit params A & B)
            if k == 0:
                h5JackIntercept[g,1]=float(PARAMS[numFitParams*g][1])
            if k == 1:
                h5JackMatelem[g,1]=float(PARAMS[k+numFitParams*g][1])

    avgFitParam*=(1.0/(1.0*gauge_configs))
    # Add the value to the dictionary
    params.update({v: avgFitParam})

    if options.makeH5:
        # Only add the intercept & asymptotic value of matrix element (i.e. fit params A & B)
        if k == 0:
            h5EnsemIntercept[0,1]=avgFitParam
        if k == 1:
            h5EnsemMatelem[0,1]=avgFitParam



##########################
# NOW ESTIMATE THE FIT PARAMETER COVARIANCE DIRECTLY FROM ALL JACKKNIFE FIT RESULTS
##########################

COV={}
for pk,pv in paramOrder.items():
    for qk,qv in paramOrder.items():
        # Make the tuple key for the covariance
        key=(pv,qv)
        dum=0.0
        for g in range(0,gauge_configs):

            dum+=(float(PARAMS[pk+numFitParams*g][1])-params.get(pv))*(float(PARAMS[qk+numFitParams*g][1])-params.get(qv))

        dum*=((1.0*(gauge_configs-1))/gauge_configs)

        COV.update({key : dum})


############################
# SAVE THE AVG FIT PARAMS BASED ON EACH JACKKNIFE FIT RESULT
############################
A=params['A']
B=params['B']
# C=params['C']
# dE=params['dE']
print "A = %.7f"%params['A']
print "B = %.7f"%params['B']
# print "C = %.7f"%params['C'] 
# print "dE = %.7f"%params['dE'] 
#######################################################

# # Dump out the fit params and errors
# for Q in 'A','B','C','E0','E1','a','b':
#     print "%s +/- %s"%(params[Q]['avg'],params[Q]['err'])
# sys.exit()

# Set the fit parameter covariance for easier access
covAA=COV[('A', 'A')]
covAB=COV[('A', 'B')]
# covAC=COV[('A', 'C')]
# covAdE=COV[('A', 'dE')]
covBA=COV[('B', 'A')]
covBB=COV[('B', 'B')]
# covBC=COV[('B', 'C')]
# covBdE=COV[('B', 'dE')]
# covCA=COV[('C', 'A')]
# covCB=COV[('C', 'B')]
# covCC=COV[('C', 'C')]
# covCdE=COV[('C', 'dE')]
# covdEA=COV[('dE', 'A')]
# covdEB=COV[('dE', 'B')]
# covdEC=COV[('dE', 'C')]
# covdEdE=COV[('dE', 'dE')]

dA=np.sqrt(covAA)
dB=np.sqrt(covBB)
# dC=np.sqrt(covCC)
# ddE=np.sqrt(covdEdE)

print "A = %.7f +/- %.7f"%(params['A'],dA)
print "B = %.7f +/- %.7f"%(params['B'],dB)
# print "C = %.7f +/- %.7f"%(params['C'],dC)
# print "dE = %.7f +/- %.7f"%(params['dE'],ddE)



if options.makeH5:
    #########################################################
    # Push the dA value to the ensemble & jack h5 files
    #########################################################
    h5EnsemIntercept[0,2]=dA
    h5EnsemMatelem[0,2]=dB
    #########################################################



# Fill and hold the full fit parameter covariance matrix
covariance=[]
for l,v in paramOrder.items():
    tmp=[]
    for k,w in paramOrder.items():
        tmp.append(COV[(v,w)])
    covariance.append(tmp)


class DERIV_SR:
    def __init__(self,params):
        self.PX=params
        self.DFDA=None
        self.DFDB=None
        self.DFDC=None
        self.DFDdE=None
        self.DFDT=None

    def SR(self,T):
        # return A+B*np.exp(-self.PX['dE']*T)+C*T*np.exp(-self.PX['dE']*T)
        # return self.PX['B']+self.PX['A']*T+self.PX['C']*T*np.exp(-self.PX['dE']*T)
        return self.PX['A']+self.PX['B']*T

    def populatePartials(self,T):
        # self.DFDA=1
        # self.DFDB=np.exp(-self.PX['dE']*T)
        # self.DFDC=T*np.exp(-self.PX['dE']*T)
        # self.DFDdE=-self.PX['B']*T*np.exp(-self.PX['dE']*T)-self.PX['C']*(T**2)*np.exp(-self.PX['dE']*T)
        # self.DFDT=self.PX['C']*np.exp(-self.PX['dE']*T)-self.PX['B']*self.PX['dE']*np.exp(-self.PX['dE']*T)-self.PX['C']*self.PX['dE']*np.exp(-self.PX['dE']*T)*T

        self.DFDA=1
        self.DFDB=T
        self.DFDT=self.PX['B']
        # self.DFDC=T*np.exp(-self.PX['dE']*T)
        # self.DFDdE=-T*self.PX['C']*T*np.exp(-self.PX['dE']*T)
        # self.DFDT=self.PX['A']+self.PX['C']*(np.exp(-self.PX['dE']*T)+-self.PX['dE']*T*np.exp(-self.PX['dE']*T))



# DETERMINE THE ERROR OF FIT TO SR DATA
def DerivSR(T,covariance):
    # Make a new instance of DERIV_SR class
    derivTracker=DERIV_SR(params)
    derivTracker.populatePartials(T)

    partials={ 'A': derivTracker.DFDA, 'B': derivTracker.DFDB}
               # 'C' : derivTracker.DFDC, 'dE': derivTracker.DFDdE }  
    
    error=0.0
    for i in ['A','B']: #,'C','dE']:
        for j in ['A','B']: #,'C','dE']:
            error+=partials[i]*COV[(i, j)]*partials[j]

    return np.sqrt(error)



def jackknifeSamples(data,complexity,numT):
    dataJks=np.zeros((gauge_configs,gauge_configs-1,numT))
    # Make jackknife samples for an input array
    for n in range(0,gauge_configs):
        if n > 0:
            for l in range(0,n):
                for t in range(0,numT):
                    dataJks[n,l,t]=data[t+l*numT,complexity]
                    
        for l in range(n,gauge_configs-1):
            for t in range(0,numT):
                dataJks[n,l,t]=data[t+(l+1)*numT,complexity]

    return dataJks


def avgJackknifeSamples(data,numT):
    dataJkAvg=np.zeros((gauge_configs,numT))
    # Form an average per jackknife sample
    for n in range(0,gauge_configs):
        for t in range(0,numT):
            dum=0.0
            for g in range(0,gauge_configs-1):
                dum+=data[n,g,t]

            dataJkAvg[n,t]=(1.0/(1.0*(gauge_configs-1)))*dum

    return dataJkAvg


def avgData(data,comp,numT):
    dataAvg=np.zeros(numT)
    for t in range(0,numT):
        dum=0.0
        for g in range(0,gauge_configs):
            dum+=data[g*numT+t,comp]

        dataAvg[t]=(1.0/(1.0*gauge_configs))*dum
    return dataAvg


def dataErrEst(dataAvg,dataJkAvg,numT):
    dataAvgErr=np.zeros(numT)
    for t in range(0,numT):
        dum=0.0
        for g in range(0,gauge_configs):
            dum+=np.power(dataJkAvg[g,t]-dataAvg[t],2)
        dataAvgErr[t]=np.sqrt(((1.0*(gauge_configs-1))/(1.0*gauge_configs))*dum)
    return dataAvgErr


# GET AVERAGE OF SR DATA FOR EACH TSEP
print "Determining average summed ratio data for each tsep"
dataSRAvg=avgData(dataSR,complexity,options.numTSeps)

# MAKE THE JACKKNIFE SAMPLES OF SR DATA
print "Making jackknife samples"
dataSRJks=jackknifeSamples(dataSR,complexity,options.numTSeps)

# MAKE AVERAGES PER JACKKNIFE SAMPLE OF SR DATA
print "Calculating averages per jackknife sample"
data2JkAvg=avgJackknifeSamples(dataSRJks,options.numTSeps)

# MAKE ERROR ESTIMATE OF SR DATA WITH JACKKNIFE SAMPLE AVERAGES
print "Determining error estimates for summed ratio data for each tsep"
dataSRAvgErr=dataErrEst(dataSRAvg,data2JkAvg,options.numTSeps)



# Initialize and populate the extracted matrix element and its error
Matelem=np.zeros(3000)
MatelemError=np.zeros(3000)
Matelem[:]=B
MatelemError[:]=np.sqrt(COV[('B', 'B')])
print "Matelem = %.7f"%Matelem[0]
print "Matelem error = %.7f"%MatelemError[0]

timeSRFit=np.linspace(2,16,3000)
# Plot the asymptotic value/error for the matrix element
plt.plot(timeSRFit,Matelem,'black')
# plt.fill_between(timeSRFit,Matelem+MatelemError,Matelem-MatelemError,color='#d8dcd6')
plt.fill_between(timeSRFit,Matelem+MatelemError,Matelem-MatelemError,color='grey',alpha=0.6)




# # Append asymptotic value of matelem + error + |zsep|
# W=open("pITD/NucleonMG1g1MxD0J0S_J1o2_H1o2D4E1."+options.insertion+".pITD.dat",'a')
# W.write(str(ioffeTime)+" "+str(options.complexity)+" "+str(A)+" "+str(MatelemError[0])+" "+options.zsep+"\n")
# W.close()


# Determine the central fit
fitSR=np.linspace(2,16,3000)
for i in range(0,len(timeSRFit)-1):
    fitSR[i]=SR(params['A'],params['B'],timeSRFit[i])
    # fitSR[i]=SR(params['A'],params['B'],params['C'],params['dE'],timeSRFit[i])
# Determine error of fit    
fitSRError=np.linspace(2,26,3000)
for i in range(0,len(timeSRFit)-1):
    fitSRError[i]=DerivSR(timeSRFit[i],covariance)
# Add fit to plot    
plt.plot(timeSRFit,fitSR[:],'c',label=r'Fit')
# plt.fill_between(timeSRFit,fitSR+fitSRError,fitSR-fitSRError,color='#B2E5FF',alpha=0.6)
plt.fill_between(timeSRFit,fitSR+fitSRError,fitSR-fitSRError,color='c',alpha=0.25)

# Adjust vertical range
# if options.complexity==1:
#     plt.ylim([0.75,1.25])
# if options.complexity==2:
#     plt.ylim([-0.25,0.25])


# if MatelemError[0] == 0:
#     plt.ylim([0.99998,1.00002])
# else:
#     plt.ylim([Matelem[0]-4*MatelemError[0],Matelem[0]+4*MatelemError[0]])


# For better vertical range specification, find the min and max errorbar
minErr=min(dataSRAvg-dataSRAvgErr)
maxErr=max(dataSRAvg+dataSRAvgErr)
# Set a 10% buffer
Ymin=minErr-(maxErr-minErr)
Ymax=maxErr+(maxErr-minErr)
plt.ylim([Ymin,Ymax])
# plt.ylim([0.998,0.9998])
# plt.ylim([-0.04,0.06])
# plt.ylim([0.8,7.2])

# plt.ylim([0.6,3.8])

ax=plt.gca()
xext=abs(ax.get_xlim()[1]-ax.get_xlim()[0]) # xrange
yext=abs(ax.get_ylim()[1]-ax.get_ylim()[0]) # yrange

# sizefont=18

# Plot the summed ratio data - color according to which data pts were fit
Time=np.linspace(options.dataTMin,options.dataTMax,options.numTSeps)
excludeTime=options.excludedTFits.split(',')
excludedTimes=[int(q) for q in excludeTime] # Cast as ints
for tt in range(0,len(Time)):
    if tt in excludedTimes:
        plt.errorbar(Time[tt],dataSRAvg[tt],yerr=dataSRAvgErr[tt],fmt='ro',
                     label="",mfc='r',mec='r',alpha=0.3)
    else:
         plt.errorbar(Time[tt],dataSRAvg[tt],yerr=dataSRAvgErr[tt],fmt='yo',
                      label=r'$R_{sum}$',mfc='y',mec='y')

# Make the figure
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
# plt.legend(by_label.values(), by_label.keys(), fontsize=sizefont)
plt.xlabel(r'$T/a$') #,fontsize=sizefont)
# plt.ylabel(r'$\mathcal{M}_{\rm{eff}}\left(\nu,z^2\right)$', fontsize=sizefont)
plt.ylabel(r'$R(T)$') # fontsize=sizefont)

# Add the momentum and displacement of this matelem
plt.text(ax.get_xlim()[1]-0.7*xext,ax.get_ylim()[0]+0.1*yext,\
         r'$\vec{p}=\left(%d,%d,%d\right)\quad\vec{z}=\left(0,0,%s\right)$'\
         %(momx,momy,momz,options.zsep)) #,fontsize=sizefont)

# Add fit function and fit results
if options.showFitRes:
    # plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.1*yext,\
        #          r'$f(T)=B+MT+CTe^{-\Delta ET}$', fontsize=14)
    plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.1*yext,\
             r'$f(T)=A+BT$') #, fontsize=14)
    plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.15*yext,\
             r'$B = %.4f(%s)$'%(params['B'],\
                                str(round(dB,4)).lstrip('0').lstrip('.').lstrip('0')))#,fontsize=14)
    plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.2*yext,\
             r'$A = %.4f(%s)$'%(params['A'],\
                                str(round(dA,5)).lstrip('0').lstrip('.').lstrip('0')))#,fontsize=14)
    # plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.25*yext,\
        #          r'$C = %.4f(%s)$'%(params['C'],\
        #                             str(round(dC,5)).lstrip('0').lstrip('.').lstrip('0')),\
        #          fontsize=14)
    # plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.3*yext,\
        #          r'$\Delta E = %.4f(%s)$'%(params['dE'],\
        #                             str(round(ddE,5)).lstrip('0').lstrip('.')),\
        #          fontsize=14)

# plt.savefig(output_name+'.png',transparent=truthTransparent,dpi=600)
# plt.show()

