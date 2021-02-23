#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pylab # to save figures to file
import sys,optparse
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
parser.add_option("-p", "--mom", type="string", default="X.X.X",
                  help='Momentum (default = X.X.X)')
# parser.add_option("-q", "--Q", type="string", default="0.0.0",
#                   help='Momentum transfer <qx>.<qy>.<qz> (default = 0.0.0)')
# parser.add_option("-p", "--pf", type="string", default="0.0.0",
#                   help='Sink momentum <pfx>.<pfy>.<pfz> (default = 0.0.0)')
# parser.add_option("-N", "--nptCorrsFile", type="string", default="",
#                   help='List file containing full path/name for Npt correlator(s) (default="")')
parser.add_option("-n", "--SRCorrFile", type="string", default="",
                  help='File containing full path/name for 2pt correlator (default="")')
parser.add_option("-I", "--insertion", type="string", default="",
                  help='Gamma matrix to consider (default = "")')
# parser.add_option("-o", "--output", type="string", default="",
#                   help='Override name of generated pITD file (default = "p<momx><momy><momz>_g<gamma>_rows<snkRow><srcRow>")')
parser.add_option("-r", "--resFile", type="string", default="",
                  help='Name of fit results file (default = "")')
parser.add_option("-v", "--cov", type="string", default="",
                  help='Fit parameter covariance matrix (default = "")')
parser.add_option("-w", "--srcRow", type="int", default="1",
                  help='Source interpolator row (default = 1)')
parser.add_option("-y", "--snkRow", type="int", default="1",
                  help='Sink interpolator row (default = 1)')
# parser.add_option("-s", "--cutNoisy", type="int", default="0",
#                   help='Flag to cut noisy data (tsep>10) default = 0')
# parser.add_option("-f", "--show2ptFitToData", type="int", default="0",
#                   help='Plot 2pt effective energy & fit following simultaneous fit - 0 or 1 (default = 0)')
# parser.add_option("-F", "--fitRange2pt", type="string", default="",
#                   help='The fit range used to fit 2pt function <tmin>.<tstep>.<tmax> - only relevant if "-f" option is true (default = "")')

# Parse the input arguments
(options, args) = parser.parse_args()
gauge_configs=options.ensemble_cfgs
complexity=options.complexity
# gamma=options.chromaGamma
paramsFile=options.resFile
covFile=options.cov

momx=int(options.mom.split('.')[0])
momy=int(options.mom.split('.')[1])
momz=int(options.mom.split('.')[2])
# snkMom=options.pf
# insMom=options.Q
# Grab and read the SR correlator
# with open(options.SRCorrFile) as pt:
#     SRCorr=pt.read()
# dataSR=np.loadtxt(SRCorr, delimiter=' ', skiprows=1)


dataSR=np.loadtxt(options.SRCorrFile, delimiter=' ', skiprows=1)


# # Method to generate distinct momentum dictionaries
# def makeMom(mom_tag):
#     momVal = []
#     for p in mom_tag.split('.'):
#         momVal.append(int(p))
#     momType = [abs(m) for m in momVal]
#     momType.sort(reverse=True)

#     return {"val": momVal, "type": momType}

# srcMomInfo=makeMom(srcMom)
# insMomInfo=makeMom(insMom)
# snkMomInfo=makeMom(snkMom)
# # A nested momentum dictionary
# momenta = {"Sink" : snkMomInfo, "Ins" : insMomInfo, "Src" : srcMomInfo}


# Set the generated summed ratio plot
output_name="%s_comp%i"%(options.SRCorrFile.replace('.dat',''),options.complexity)

    # output_name="p%s%s%s_g%i_rows%i%i"%(momenta["Src"]["val"][0],momenta["Src"]["val"][1],momenta["Src"]["val"][2],gamma,options.snkRow,options.srcRow)





# The Fit function to SR data
# def SR(A,B,C,dE,T):
#     return A+B*np.exp(-dE*T)+C*T*np.exp(-dE*T)
def SR(A,B,T):
    return B+A*T

####################
# DEFINE A DICTIONARY TO HOLD ORDERING OF KEY/VALUES
####################
# paramOrder={0: 'A', 1: 'B', 2: 'C', 3: 'dE'}
paramOrder={0: 'A', 1: 'B'}



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

    avgFitParam*=(1.0/(1.0*gauge_configs))
    # Add the value to the dictionary
    params.update({v: avgFitParam})


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


# Set the fit params and errors for easier access
A=params['A']
B=params['B']
# C=params['C']
# dE=params['dE']

print "A = %.7f"%params['A']
print "B = %.7f"%params['B']
# print "C = %.7f"%params['C'] 
# print "dE = %.7f"%params['dE'] 



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
        return self.PX['B']+self.PX['A']*T

    def populatePartials(self,T):
        # self.DFDA=1
        # self.DFDB=np.exp(-self.PX['dE']*T)
        # self.DFDC=T*np.exp(-self.PX['dE']*T)
        # self.DFDdE=-self.PX['B']*T*np.exp(-self.PX['dE']*T)-self.PX['C']*(T**2)*np.exp(-self.PX['dE']*T)
        # self.DFDT=self.PX['C']*np.exp(-self.PX['dE']*T)-self.PX['B']*self.PX['dE']*np.exp(-self.PX['dE']*T)-self.PX['C']*self.PX['dE']*np.exp(-self.PX['dE']*T)*T
        self.DFDA=T
        self.DFDB=1
        self.DFDT=self.PX['A']


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
Matelem[:]=A
MatelemError[:]=np.sqrt(COV[('A', 'A')])
print "Matelem   = %.7f +/- %.7f"%(Matelem[0],MatelemError[0])
print "Intercept = %.7f +/- %.7f"%(B,dB)

timeSRFit=np.linspace(2,16,3000)
# Plot the asymptotic value/error for the matrix element
plt.plot(timeSRFit,Matelem,'black')
# plt.fill_between(timeSRFit,Matelem+MatelemError,Matelem-MatelemError,color='#d8dcd6',alpha=0.6)
plt.fill_between(timeSRFit,Matelem+MatelemError,Matelem-MatelemError,color='grey',alpha=0.6)



ioffeTime=((2*np.pi)/(32.0))*(momz)*int(options.zsep)


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
fitSRError=np.linspace(2,16,3000)
for i in range(0,len(timeSRFit)-1):
    fitSRError[i]=DerivSR(timeSRFit[i],covariance)
# Add fit to plot    
plt.plot(timeSRFit,fitSR[:],'b',label=r'Fit')
# plt.fill_between(timeSRFit,fitSR+fitSRError,fitSR-fitSRError,color='#B2E5FF',alpha=0.6)
plt.fill_between(timeSRFit,fitSR+fitSRError,fitSR-fitSRError,color='b',alpha=0.25)

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

# plt.ylim([0.9988,1])

# Plot the summed ratio data
Time=np.linspace(4,14,6)
plt.errorbar(Time,dataSRAvg,yerr=dataSRAvgErr,fmt='bo',label=r'$\mathcal{M}$',mfc='b',mec='b')
# Make the figure
plt.legend(fontsize=10)
plt.rc('text', usetex=True)
plt.xlabel(r'T_{\rm{sep}}',fontsize=14)
plt.ylabel(r'$\mathcal{M}_{\rm{eff}}\left(\nu,z^2\right)$', fontsize=14)


ax=plt.gca()
xext=abs(ax.get_xlim()[1]-ax.get_xlim()[0]) # xrange
yext=abs(ax.get_ylim()[1]-ax.get_ylim()[0]) # yrange

# Add the momentum and displacement of this matelem
plt.text(ax.get_xlim()[1]-0.7*xext,ax.get_ylim()[0]+0.1*yext,\
         r'$\vec{p}=\left(%d,%d,%d\right)\quad\vec{z}=\left(0,0,%s\right)$'\
         %(momx,momy,momz,options.zsep))

plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.1*yext,\
         r'$f(T)=B+MT$') #, fontsize=14)
plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.15*yext,\
         r'$M = %.4f(%s)$'%(params['A'],\
                            str(round(dA,4)).lstrip('0').lstrip('.').lstrip('0')))
plt.text(ax.get_xlim()[0]+0.1*xext,ax.get_ylim()[1]-0.2*yext,\
         r'$B = %.4f(%s)$'%(params['B'],\
                            str(round(dB,5)).lstrip('0').lstrip('.').lstrip('0')))

#plt.savefig(output_name+'.png',format='png')#,dpi=1200)
plt.show()

