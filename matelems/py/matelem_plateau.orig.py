#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import pylab # to save figures to file
import sys,optparse
#from eff_charge_util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-c", "--ensemble_cfgs", type="int", default=1,
                  help='Number of configurations (default = 1)')
parser.add_option("-x", "--complexity", type="int", default=0,
                  help='Component of Npt function to analyze [Real=1 Imag=2] (default = 1)')
parser.add_option("-k", "--pi", type="string", default="0.0.0",
                  help='Source momentum <pix>.<piy>.<piz> (default = 0.0.0)')
parser.add_option("-q", "--Q", type="string", default="0.0.0",
                  help='Momentum transfer <qx>.<qy>.<qz> (default = 0.0.0)')
parser.add_option("-p", "--pf", type="string", default="0.0.0",
                  help='Sink momentum <pfx>.<pfy>.<pfz> (default = 0.0.0)')
parser.add_option("-N", "--nptCorrsFile", type="string", default="",
                  help='List file containing full path/name for Npt correlator(s) (default="")')
parser.add_option("-n", "--twoptCorrsFile", type="string", default="",
                  help='File containing full path/name for 2pt correlator (default="")')
parser.add_option("-g", "--chromaGamma", type="int", default="-1",
                  help='Gamma matrix to consider (default = -1)')
parser.add_option("-o", "--output", type="string", default="",
                  help='Override name of generated effective charge plot (default = "p<momx><momy><momz>_g<gamma>_rows<snkRow><srcRow>")')
parser.add_option("-r", "--resFile", type="string", default="",
                  help='Name of fit results file (default = "")')
parser.add_option("-v", "--cov", type="string", default="",
                  help='Fit parameter covariance matrix (default = "")')
parser.add_option("-w", "--srcRow", type="string", default="X",
                  help='Source interpolator row (default = X)')
parser.add_option("-y", "--snkRow", type="string", default="X",
                  help='Sink interpolator row (default = X)')
parser.add_option("-s", "--cutNoisy", type="int", default="0",
                  help='Cut noisy data tsep>10 (1) -or- Cut noisy data tsep=16 (2)  (default = 0)')
parser.add_option("-f", "--show2ptFitToData", type="int", default="0",
                  help='Plot 2pt effective energy & fit following simultaneous fit - 0 or 1 (default = 0)')
parser.add_option("-F", "--fitRange2pt", type="string", default="",
                  help='The fit range used to fit 2pt function <tmin>.<tstep>.<tmax> - only relevant if "-f" option is true (default = "")')
parser.add_option("-R", "--rawRatio", type="int", default=0,
                  help='Plot raw 3pt/2pt lattice ratios -or- effective charges [default] (default = 0)')

# Parse the input arguments
(options, args) = parser.parse_args()
gauge_configs=options.ensemble_cfgs
complexity=options.complexity
gamma=options.chromaGamma
paramsFile=options.resFile
covFile=options.cov
srcMom=options.pi
snkMom=options.pf
insMom=options.Q
# Grab and read the 2pt correlator
with open(options.twoptCorrsFile) as pt:
    twoptCorr=pt.read()
data2pt=np.loadtxt(twoptCorr, delimiter=' ', skiprows=1)
# Grab and read Npt correlators
with open(options.nptCorrsFile) as tpt:
    nptCorrs=tpt.read().splitlines()
data3pt=[]
for f in nptCorrs:
    data3pt.append(np.loadtxt(f, delimiter=' ', skiprows=1))
# Set the momenta throughout


# Method to generate distinct momentum dictionaries
def makeMom(mom_tag):
    momVal = []
    for p in mom_tag.split('.'):
        momVal.append(int(p))
    momType = [abs(m) for m in momVal]
    momType.sort(reverse=True)

    return {"val": momVal, "type": momType}

srcMomInfo=makeMom(srcMom)
insMomInfo=makeMom(insMom)
snkMomInfo=makeMom(snkMom)
# A nested momentum dictionary
momenta = {"Sink" : snkMomInfo, "Ins" : insMomInfo, "Src" : srcMomInfo}


# Set the generated effective charge plot with fits
if options.output == "":
    output_name="p%i%i%i_g%i_rows%s%s"%(momenta["Src"]["val"][0],momenta["Src"]["val"][1],momenta["Src"]["val"][2],gamma,options.snkRow,options.srcRow)
else:
    output_name=options.output

# Just some temp stuff
dataT6=data3pt[0]
dataT8=data3pt[1]
dataT10=data3pt[2]
if options.cutNoisy == 0:
    dataT12=data3pt[3]
    dataT14=data3pt[4]
    dataT16=data3pt[5]
elif options.cutNoisy == 2:
    dataT12=data3pt[3]
    dataT14=data3pt[4]



# SET RENORMALIZATION CONSTANTS
renormalizations = { 0 : 0.800, 1 : 0.834, 2 : 0.834, 4 : 0.834, 8 : 0.834, 7 : 0.887, 11 : 0.887,
                     13 : 0.887, 14 : 0.887, 3 : 0.925, 12 : 0.925, 15: 1.000 }

print "Gamma  = %d"%gamma
print "Renorm = %.7f"%renormalizations[gamma]



def twoState(a,b,E0,E1,T):
    return np.exp(-E0*T)*(a+b*np.exp(-(E1-E0)*T))

# def fit3pt2ptRatio(A,B,C,E0,E1,a,b,t,T):
#     return (np.exp(-E0*T)*(A+B*np.exp(-(E1-E0)*T)+C*np.exp(-(E1-E0)*(T/2.0))*np.cosh((E1-E0)*(t-(T/2.0)))))/(np.exp(-E0*T)*(a+b*np.exp(-(E1-E0)*T)))

def fit3pt2ptRatio(paramsDict,t,T):
    A=paramsDict['A']
    B=paramsDict['B']
    C=paramsDict['C']
    E0=paramsDict['E0']
    E1=paramsDict['E1']
    a=paramsDict['a']
    b=paramsDict['b']

    if options.srcRow==options.snkRow:
        return (np.exp(-E0*T)*(A+B*np.exp(-(E1-E0)*T)+C*np.exp(-(E1-E0)*(T/2.0))*np.cosh((E1-E0)*(t-(T/2.0)))))/(np.exp(-E0*T)*(a+b*np.exp(-(E1-E0)*T)))
    if options.srcRow!=options.snkRow:
        D=paramsDict['D']
        return (np.exp(-E0*T)*(A+B*np.exp(-(E1-E0)*T)+C*np.exp(-(E1-E0)*t)+D*np.exp(-(E1-E0)*T)*np.exp((E1-E0)*t)))/(np.exp(-E0*T)*(a+b*np.exp(-(E1-E0)*T)))



####################
# DEFINE A DICTIONARY TO HOLD ORDERING OF KEY/VALUES
####################
paramOrder={}
if options.srcRow==options.snkRow:
    paramOrder={0: 'A', 1: 'B', 2: 'C', 3: 'E0', 4: 'E1', 5: 'a', 6: 'b'}
if options.srcRow!=options.snkRow:
    paramOrder={0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E0', 5: 'E1', 6: 'a', 7: 'b'}


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

# print params['a']
# print np.sqrt(COV[('a', 'a')])
# print params['A']
# print np.sqrt(COV[('A', 'A')])
# print COV[('A', 'a')]
# print "***************"

# Set the fit params and errors for easier access
A=params['A']
B=params['B']
C=params['C']
E0=params['E0']
E1=params['E1']
a=params['a']
b=params['b']
if options.srcRow!=options.snkRow:
    D=params['D']

print "A = %.7f"%params['A']
print "B = %.7f"%params['B']
print "C = %.7f"%params['C'] 
print "E0 = %.7f"%params['E0'] 
print "E1 = %.7f"%params['E1'] 
print "a = %.7f"%params['a'] 
print "b = %.7f"%params['b'] 
if options.srcRow!=options.snkRow:
    print "D = %.7f"%params['D'] 




# A=params['A']['avg']
# dA=params['A']['err']
# B=params['B']['avg']
# dB=params['B']['err']
# C=params['C']['avg']
# dC=params['C']['err']
# E0=params['E0']['avg']
# dE0=params['E0']['err']
# E1=params['E1']['avg']
# dE1=params['E1']['err']
# a=params['a']['avg']
# da=params['a']['err']
# b=params['b']['avg']
# db=params['b']['err']


# # Dump out the fit params and errors
# for Q in 'A','B','C','E0','E1','a','b':
#     print "%s +/- %s"%(params[Q]['avg'],params[Q]['err'])
# sys.exit()

# Set the fit parameter covariance for easier access
covAA=COV[('A', 'A')]
covAB=COV[('A', 'B')]
covAC=COV[('A', 'C')]
covAE0=COV[('A', 'E0')]
covAE1=COV[('A', 'E1')]
covAa=COV[('A', 'a')]
covAb=COV[('A', 'b')]
covBA=COV[('B', 'A')]
covBB=COV[('B', 'B')]
covBC=COV[('B', 'C')]
covBE0=COV[('B', 'E0')]
covBE1=COV[('B', 'E1')]
covBa=COV[('B', 'a')]
covBb=COV[('B', 'b')]
covCA=COV[('C', 'A')]
covCB=COV[('C', 'B')]
covCC=COV[('C', 'C')]
covCE0=COV[('C', 'E0')]
covCE1=COV[('C', 'E1')]
covCa=COV[('C', 'a')]
covCb=COV[('C', 'b')]
covE0A=COV[('E0', 'A')]
covE0B=COV[('E0', 'B')]
covE0C=COV[('E0', 'C')]
covE0E0=COV[('E0', 'E0')]
covE0E1=COV[('E0', 'E1')]
covE0a=COV[('E0', 'a')]
covE0b=COV[('E0', 'b')]
covE1A=COV[('E1', 'A')]
covE1B=COV[('E1', 'B')]
covE1C=COV[('E1', 'C')]
covE1E0=COV[('E1', 'E0')]
covE1E1=COV[('E1', 'E1')]
covE1a=COV[('E1', 'a')]
covE1b=COV[('E1', 'b')]
covaA=COV[('a', 'A')]
covaB=COV[('a', 'B')]
covaC=COV[('a', 'C')]
covaE0=COV[('a', 'E0')]
covaE1=COV[('a', 'E1')]
covaa=COV[('a', 'a')]
covab=COV[('a', 'b')]
covbA=COV[('b', 'A')]
covbB=COV[('b', 'B')]
covbC=COV[('b', 'C')]
covbE0=COV[('b', 'E0')]
covbE1=COV[('b', 'E1')]
covba=COV[('b', 'a')]
covbb=COV[('b', 'b')]
if options.srcRow != options.snkRow:
    covDD=COV[('D', 'D')]

dA=np.sqrt(covAA)
dB=np.sqrt(covBB)
dC=np.sqrt(covCC)
dE0=np.sqrt(covE0E0)
dE1=np.sqrt(covE1E1)
da=np.sqrt(covaa)
print da
db=np.sqrt(covbb)
if options.srcRow!=options.srcRow:
    dD=np.sqrt(covDD)

# Fill and hold the full fit parameter covariance matrix
covariance=[]
for l,v in paramOrder.items():
    tmp=[]
    for k,w in paramOrder.items():
        tmp.append(COV[(v,w)])
    covariance.append(tmp)


# print covariance[0][0]
# print covAA
# print covariance[1][1]
# print covBB
# print covariance[2][2]
# print covCC
# print covariance[3][3]
# print covE0E0
# print covariance[4][4]
# print covE1E1
# print covariance[5][5]
# print covaa
# print covariance[6][6]
# print covbb

# print covariance[0][0]
# print covAA
# print covariance[1][1]
# print covBB
# print covariance[2][2]
# print covCC
# print covariance[3][3]
# print covDD
# print covariance[4][4]
# print covE0E0
# print covariance[5][5]
# print covE1E1
# print covariance[6][6]
# print covaa
# print covariance[7][7]
# print covbb

# sys.exit()

class DERIV2pt:
    def __init__(self,params):
        self.PX=params
        self.dFda=None
        self.dFdb=None
        self.dFdE0=None
        self.dFdE1=None

    def C2(self,T):
        return np.exp(-self.PX['E0']*T)*(self.PX['a']+self.PX['b']*np.exp(-(self.PX['E1']-self.PX['E0'])*T))

    def populatePartials(self,T,tdiff):

        # self.dFda=(2*self.PX['a']+self.PX['b']*np.exp(self.PX['E0']*T-self.PX['E1']*(1+T))*(np.exp(self.PX['E0'])+np.exp(self.PX['E1'])))/((self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(1+T)))*tdiff)

        # self.dFdb=(np.exp(self.PX['E0']*T-self.PX['E1']*(1+T))*(2*self.PX['b']*np.exp(self.PX['E0']+self.PX['E0']*T-self.PX['E1']*T)+self.PX['a']*(np.exp(self.PX['E0'])+np.exp(self.PX['E1']))))/((self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(1+T)))*tdiff)

        # self.dFdE0=(self.PX['a']**2+2*(self.PX['b']**2)*np.exp((self.PX['E0']-self.PX['E1'])*(1+2*T))*(1+T)+self.PX['a']*self.PX['b']*np.exp(self.PX['E0']*T-self.PX['E1']*(1+T))*(np.exp(self.PX['E1'])*(1+T)+np.exp(self.PX['E0'])*(2+T)))/((self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(1+T)))*tdiff)

        # self.dFdE1=-((self.PX['b']*np.exp(self.PX['E0']*T)*(self.PX['b']*np.exp(self.PX['E0']+self.PX['E0']*T)*(1+2*T)+self.PX['a']*np.exp(self.PX['E1']*T)*(np.exp(self.PX['E1'])*T+np.exp(self.PX['E0'])*(1+T))))/((self.PX['b']*np.exp(self.PX['E0']*T)+self.PX['a']*np.exp(self.PX['E1']*T))*(self.PX['a']*np.exp(self.PX['E1']*(1+T))+self.PX['b']*np.exp(self.PX['E0']+self.PX['E0']*T))*tdiff))


	self.dFda=(1.0/tdiff)*(1/(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))-1/(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(tdiff+T))))

	

	self.dFdb=(self.PX['a']/(self.PX['b']*tdiff))*(-(1/(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T)))+1/(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(tdiff+T))))


	self.dFdE0=(self.PX['a']*(self.PX['a']*tdiff+self.PX['b']*np.exp(self.PX['E0']*T-self.PX['E1']*(tdiff+T))*(-np.exp(tdiff*self.PX['E0'])*T+np.exp(tdiff*self.PX['E1'])*(tdiff+T))))/(tdiff*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(tdiff+T))))



	self.dFdE1=(self.PX['b']*(self.PX['b']*tdiff*np.exp((self.PX['E0']-self.PX['E1'])*(tdiff+2*T))-self.PX['a']*np.exp((self.PX['E0']-self.PX['E1'])*T)*T+self.PX['a']*np.exp((self.PX['E0']-self.PX['E1'])*(tdiff+T))*(tdiff+T)))/(tdiff*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(tdiff+T))))





        # self.dFda=(self.PX['b']*np.exp((self.PX['E0']+self.PX['E1'])*T)*(np.exp(self.PX['E0'])-np.exp(self.PX['E1'])))/((self.PX['b']*np.exp(self.PX['E0']*T)+self.PX['a']*np.exp(self.PX['E1']*T))*(self.PX['a']*np.exp(self.PX['E1']*(1+T))+self.PX['b']*np.exp(self.PX['E0']+self.PX['E0']*T))*tdiff)


        # self.dFdb=(self.PX['a']*np.exp(self.PX['E0']*T-self.PX['E1']*(1+T))*(-np.exp(self.PX['E0'])+np.exp(self.PX['E1'])))/((self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(1+T)))*tdiff)


        # self.dFdE0=(self.PX['a']*(self.PX['a']+self.PX['b']*np.exp(self.PX['E0']*T-self.PX['E1']*(1+T))*(-np.exp(self.PX['E0'])*T+np.exp(self.PX['E1'])*(1+T))))/((self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(1+T)))*tdiff)


        # self.dFdE1=(self.PX['b']*np.exp(self.PX['E0']*T-self.PX['E1']*(1+T))*(self.PX['b']*np.exp(self.PX['E0']+self.PX['E0']*T-self.PX['E1']*T)-self.PX['a']*np.exp(self.PX['E1'])*T+self.PX['a']*np.exp(self.PX['E0'])*(1+T)))/((self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T))*(self.PX['a']+self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*(1+T)))*tdiff)





class DERIV:
    def __init__(self,snkrow,srcrow,params):
        self.snkrow=snkrow
        self.srcrow=srcrow
        self.PX=params
        self.dRdA=None
        self.dRdB=None
        self.dRdC=None
        self.dRdD=None
        self.dRdE0=None
        self.dRdE1=None
        self.dRda=None
        self.dRdb=None

    def C2(self,T):
        return self.PX['a']+self.PX['b']*np.exp(-(self.PX['E1']-self.PX['E0'])*T)

    def C4(self,T,t):
        if self.srcrow==self.snkrow:
            return self.PX['A']+self.PX['B']*np.exp(-(self.PX['E1']-self.PX['E0'])*T)+self.PX['C']*np.exp(-(self.PX['E1']-self.PX['E0'])*(T/2.0))*np.cosh((self.PX['E1']-self.PX['E0'])*(t-(T/2.0)))
        if self.srcrow!=self.snkrow:
            return self.PX['A']+self.PX['B']*np.exp(-(self.PX['E1']-self.PX['E0'])*T)+self.PX['C']*np.exp(-(self.PX['E1']-self.PX['E0'])*t)+self.PX['D']*np.exp(-(self.PX['E1']-self.PX['E0'])*T)*np.exp((self.PX['E1']-self.PX['E0'])*t)

    def populatePartials(self,T,t):
        self.dRdA=1/self.C2(T)
        self.dRdB=(np.exp(-(self.PX['E1']-self.PX['E0'])*T))/(self.C2(T))
        if self.srcrow==self.snkrow:
            self.dRdC=(np.exp(-(self.PX['E1']-self.PX['E0'])*(T/2))*np.cosh((self.PX['E1']-self.PX['E0'])*(t-T/2)))/(self.C2(T))
            self.dRdE0=(self.C2(T)*(T*self.PX['B']*np.exp(-(self.PX['E1']-self.PX['E0'])*T)+(T/2)*self.PX['C']*np.exp(-(self.PX['E1']-self.PX['E0'])*(T/2))*np.cosh((self.PX['E1']-self.PX['E0'])*(t-T/2))-(t-T/2)*self.PX['C']*np.exp(-(self.PX['E1']-self.PX['E0'])*(T/2))*np.sinh((self.PX['E1']-self.PX['E0'])*(t-T/2)))-self.C4(T,t)*self.PX['b']*np.exp(-(self.PX['E1']-self.PX['E0'])*T))/(self.C2(T)**2)
            self.dRdE1=(self.C2(T)*(-T*self.PX['B']*np.exp(-(self.PX['E1']-self.PX['E0'])*T)-(T/2)*self.PX['C']*np.exp(-(self.PX['E1']-self.PX['E0'])*(T/2))*np.cosh((self.PX['E1']-self.PX['E0'])*(t-T/2))+(t-T/2)*self.PX['C']*np.exp(-(self.PX['E1']-self.PX['E0'])*(T/2))*np.sinh((self.PX['E1']-self.PX['E0'])*(t-T/2)))+self.C4(T,t)*self.PX['b']*np.exp(-(self.PX['E1']-self.PX['E0'])*T))/(self.C2(T)**2)
        if self.srcrow!=self.snkrow:
            self.dRdC=np.exp(-(self.PX['E1']-self.PX['E0'])*t)/self.C2(T)
            self.dRdD=np.exp((self.PX['E1']-self.PX['E0'])*t-(self.PX['E1']-self.PX['E0'])*T)/self.C2(T)
            self.dRdE0=-((self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T)*(self.PX['A']+self.PX['C']*np.exp((self.PX['E0']-self.PX['E1'])*t)+self.PX['B']*np.exp((self.PX['E0']-self.PX['E1'])*T)+self.PX['D']*np.exp((-self.PX['E0']+self.PX['E1'])*t+(self.PX['E0']-self.PX['E1'])*T))*T)/(self.C2(T))**2)+(self.PX['C']*np.exp((self.PX['E0']-self.PX['E1'])*t)*t+self.PX['B']*np.exp((self.PX['E0']-self.PX['E1'])*T)*T+self.PX['D']*np.exp((-self.PX['E0']+self.PX['E1'])*t+(self.PX['E0']-self.PX['E1'])*T)*(-t + T))/(self.C2(T))

            self.dRdE1=(self.PX['b']*np.exp((self.PX['E0']-self.PX['E1'])*T)*(self.PX['A']+self.PX['C']*np.exp((self.PX['E0']-self.PX['E1'])*t)+self.PX['B']*np.exp((self.PX['E0']-self.PX['E1'])*T)+self.PX['D']*np.exp((-self.PX['E0']+self.PX['E1'])*t+(self.PX['E0']-self.PX['E1'])*T))*T)/(self.C2(T))**2+(-self.PX['C']*np.exp((self.PX['E0']-self.PX['E1'])*t)*t+self.PX['D']*np.exp((-self.PX['E0']+self.PX['E1'])*t+(self.PX['E0']-self.PX['E1'])*T)*(t-T)-self.PX['B']*np.exp((self.PX['E0']-self.PX['E1'])*T)*T)/(self.C2(T))

        self.dRda=-self.C4(T,t)/(self.C2(T)**2)
        self.dRdb=(-self.C4(T,t)*np.exp(-(self.PX['E1']-self.PX['E0'])*T))/(self.C2(T)**2)


# DETERMINE THE ERROR OF FIT TO 2PT FN
def Deriv2pt(T,tdiff,covariance):
    # Make a new instance of Deriv2pt class
    derivTracker=DERIV2pt(params)
    derivTracker.populatePartials(T,tdiff)

    partials={ 'E0': derivTracker.dFdE0, 'E1': derivTracker.dFdE1,
               'a' : derivTracker.dFda, 'b': derivTracker.dFdb }

    # print "partials[E0] = %.7f"%partials['E0']
    # print "partials[E1] = %.7f"%partials['E1']
    # print "partials[a] = %.7f"%partials['a']
    # print "partials[b] = %.7f"%partials['b']
    
    
    error=0.0
    for i in ['E0','E1','a','b']:
        for j in ['E0','E1','a','b']:
            error+=partials[i]*COV[(i, j)]*partials[j]

    return np.sqrt(error)


# DETERMINE THE ERROR OF SIMULTANEOUS FIT TO 2PT & 3PT FNS
def Deriv(t,T,covariance):
    # Make a new instance of Deriv class
    derivTracker=DERIV(options.snkRow,options.srcRow,params)
    derivTracker.populatePartials(T,t)

    partials=np.zeros(numFitParams)
    partials[0]=derivTracker.dRdA
    partials[1]=derivTracker.dRdB
    partials[2]=derivTracker.dRdC
    if options.srcRow!=options.snkRow:
        partials[3]=derivTracker.dRdD
        partials[4]=derivTracker.dRdE0
        partials[5]=derivTracker.dRdE1
        partials[6]=derivTracker.dRda
        partials[7]=derivTracker.dRdb
    else:
        partials[3]=derivTracker.dRdE0
        partials[4]=derivTracker.dRdE1
        partials[5]=derivTracker.dRda
        partials[6]=derivTracker.dRdb
    
    # C2=a+b*np.exp(-(E1-E0)*T)
    # C4=A+B*np.exp(-(E1-E0)*T)+C*np.exp(-(E1-E0)*(T/2.0))*np.cosh((E1-E0)*(t-(T/2.0)))
    # dRdA=C2/(C2**2)
    # dRdB=(C2*np.exp(-(E1-E0)*T))/(C2**2)
    # dRdC=(C2*np.exp(-(E1-E0)*(T/2))*np.cosh((E1-E0)*(t-T/2)))/(C2**2)
    # dRdE0=(C2*(T*B*np.exp(-(E1-E0)*T)+(T/2)*C*np.exp(-(E1-E0)*(T/2))*np.cosh((E1-E0)*(t-T/2))-(t-T/2)*C*np.exp(-(E1-E0)*(T/2))*np.sinh((E1-E0)*(t-T/2)))-C4*b*np.exp(-(E1-E0)*T))/(C2**2)
    # dRdE1=(C2*(-T*B*np.exp(-(E1-E0)*T)-(T/2)*C*np.exp(-(E1-E0)*(T/2))*np.cosh((E1-E0)*(t-T/2))+(t-T/2)*C*np.exp(-(E1-E0)*(T/2))*np.sinh((E1-E0)*(t-T/2)))+C4*b*np.exp(-(E1-E0)*T))/(C2**2)
    # dRda=-C4/(C2**2)
    # dRdb=(-C4*np.exp(-(E1-E0)*T))/(C2**2)
    # partials=np.zeros(7)
    # partials[0]=dRdA
    # partials[1]=dRdB
    # partials[2]=dRdC
    # partials[3]=dRdE0
    # partials[4]=dRdE1
    # partials[5]=dRda
    # partials[6]=dRdb
    
    error=0.0
    for i in range(0,numFitParams):
        for j in range(0,numFitParams):
            # error+=partials[i]*COV[i,j]*partials[j]
            error+=partials[i]*covariance[i][j]*partials[j]
            
    return np.sqrt(error)



# def extChargeError(t,T,COV=[]):
    
#     C2=a
#     C4=A+B*np.exp(-(E1-E0)*T)+C*np.exp(-(E1-E0)*(T/2.0))*np.cosh((E1-E0)*(t-(T/2.0)))
    
#     dRdA=C2/(C2**2)
#     dRdB=(C2*np.exp(-(E1-E0)*T))/(C2**2)
#     dRdC=(C2*np.exp(-(E1-E0)*(T/2))*np.cosh((E1-E0)*(t-T/2)))/(C2**2)
#     dRdE0=(C2*(T*B*np.exp(-(E1-E0)*T)+(T/2)*C*np.exp(-(E1-E0)*(T/2))*np.cosh((E1-E0)*(t-T/2))-(t-T/2)*C*np.exp(-(E1-E0)*(T/2))*np.sinh((E1-E0)*(t-T/2)))-C4*b*np.exp(-(E1-E0)*T))/(C2**2)
#     dRdE1=(C2*(-T*B*np.exp(-(E1-E0)*T)-(T/2)*C*np.exp(-(E1-E0)*(T/2))*np.cosh((E1-E0)*(t-T/2))+(t-T/2)*C*np.exp(-(E1-E0)*(T/2))*np.sinh((E1-E0)*(t-T/2)))+C4*b*np.exp(-(E1-E0)*T))/(C2**2)
#     dRda=-C4/(C2**2)
#     dRdb=(-C4*np.exp(-(E1-E0)*T))/(C2**2)
    
#     partials=np.zeros(7)
#     partials[0]=dRdA
#     partials[1]=dRdB
#     partials[2]=dRdC
#     partials[3]=dRdE0
#     partials[4]=dRdE1
#     partials[5]=dRda
#     partials[6]=dRdb
    
#     error=0.0
#     for i in range(0,7):
#         for j in range(0,7):
#             error+=partials[i]*COV[i,j]*partials[j]
            
#     return np.sqrt(error)



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

def avgJackknifeRatio(data3JksAvg,data2JkAvg,numT):
    ratioAvg=np.zeros(numT)
    for t in range(0,numT):
        dum=0.0
        for n in range(0,gauge_configs):
            # dum+=(data3JksAvg[n,t]/data2JkAvg[n,t])
            dum+=(data3JksAvg[n,t]/data2JkAvg[n,numT])

        ratioAvg[t]=(1.0/(1.0*gauge_configs))*dum

    return ratioAvg

def avgJackknifeRatioErr(data3JkAvg,data2JkAvg,ratioAvg,numT):
    ratioErr=np.zeros(numT)
    for t in range(0,numT):
        dum=0.0
        for n in range(0,gauge_configs):
            # dum+=np.power((data3JkAvg[n,t]/data2JkAvg[n,t])-ratioAvg[t],2)
            dum+=np.power((data3JkAvg[n,t]/data2JkAvg[n,numT])-ratioAvg[t],2)

        ratioErr[t]=np.sqrt(((1.0*(gauge_configs-1))/(1.0*gauge_configs))*dum)

    return ratioErr

def makeTimeArrays(N):
    time=np.zeros(N)
    for i in range(0,N):
        time[i]=i-N/2

    return time


def makeEffCharge(data3Jk,tsep):
    effCh=np.zeros(tsep)
    for t in range(0,tsep):
        for j in range(0,gauge_configs):
            effCh[t]+=data3Jk[j,t]/(1.0*twoState(params['a'],params['b'],params['E0'],params['E1'],tsep))

        effCh[t]*=(1.0/(1.0*gauge_configs))

    return effCh

def makeEffChargeError(data3Jk,effCh,tsep):
    effChErr=np.zeros(tsep)
    for t in range(0,tsep):
        for j in range(0,gauge_configs):
            effChErr[t]+=np.power((data3Jk[j,t]/(1.0*twoState(params['a'],params['b'],params['E0'],params['E1'],tsep)))-effCh[t],2)

        effChErr[t]=np.sqrt(((1.0*(gauge_configs-1))/(1.0*gauge_configs))*effChErr[t])

    return effChErr


# MAKE THE JACKKNIFE SAMPLES OF 2PT & 4PT DATA
print "Making jackknife samples"
data2Jks=jackknifeSamples(data2pt,1,21)
t6Jks=jackknifeSamples(dataT6,complexity,6)
t8Jks=jackknifeSamples(dataT8,complexity,8)
t10Jks=jackknifeSamples(dataT10,complexity,10)
if options.cutNoisy == 0:
    t12Jks=jackknifeSamples(dataT12,complexity,12)
    t14Jks=jackknifeSamples(dataT14,complexity,14)
    t16Jks=jackknifeSamples(dataT16,complexity,16)
elif options.cutNoisy == 2:
    t12Jks=jackknifeSamples(dataT12,complexity,12)
    t14Jks=jackknifeSamples(dataT14,complexity,14)


# MAKE AVERAGES PER JACKKNIFE SAMPLE OF 2PT & 4PT DATA
print "Calculating averages per jackknife sample"
data2JkAvg=avgJackknifeSamples(data2Jks,21)
t6JksAvg=avgJackknifeSamples(t6Jks,6)
t8JksAvg=avgJackknifeSamples(t8Jks,8)
t10JksAvg=avgJackknifeSamples(t10Jks,10)
if options.cutNoisy == 0:
    t12JksAvg=avgJackknifeSamples(t12Jks,12)
    t14JksAvg=avgJackknifeSamples(t14Jks,14)
    t16JksAvg=avgJackknifeSamples(t16Jks,16)
elif options.cutNoisy == 2:
    t12JksAvg=avgJackknifeSamples(t12Jks,12)
    t14JksAvg=avgJackknifeSamples(t14Jks,14)


# GET AN AVERAGE DETERMINATION PER TIME SLICE
print "Calculating average estimates per time slice"
t6RatioAvg=avgJackknifeRatio(t6JksAvg,data2JkAvg,6)
print t6RatioAvg
print t6JksAvg
print data2JkAvg
t8RatioAvg=avgJackknifeRatio(t8JksAvg,data2JkAvg,8)
t10RatioAvg=avgJackknifeRatio(t10JksAvg,data2JkAvg,10)
if options.cutNoisy == 0:
    t12RatioAvg=avgJackknifeRatio(t12JksAvg,data2JkAvg,12)
    t14RatioAvg=avgJackknifeRatio(t14JksAvg,data2JkAvg,14)
    t16RatioAvg=avgJackknifeRatio(t16JksAvg,data2JkAvg,16)
elif options.cutNoisy == 2:
    t12RatioAvg=avgJackknifeRatio(t12JksAvg,data2JkAvg,12)
    t14RatioAvg=avgJackknifeRatio(t14JksAvg,data2JkAvg,14)

# GET AN ERROR DETERMINATION PER TIME SLICE
print "Calculating error estimates per time slice"
t6RatioErr=avgJackknifeRatioErr(t6JksAvg,data2JkAvg,t6RatioAvg,6)
t8RatioErr=avgJackknifeRatioErr(t8JksAvg,data2JkAvg,t8RatioAvg,8)
t10RatioErr=avgJackknifeRatioErr(t10JksAvg,data2JkAvg,t10RatioAvg,10)
if options.cutNoisy == 0:
    t12RatioErr=avgJackknifeRatioErr(t12JksAvg,data2JkAvg,t12RatioAvg,12)
    t14RatioErr=avgJackknifeRatioErr(t14JksAvg,data2JkAvg,t14RatioAvg,14)
    t16RatioErr=avgJackknifeRatioErr(t16JksAvg,data2JkAvg,t16RatioAvg,16)
elif options.cutNoisy == 2:
    t12RatioErr=avgJackknifeRatioErr(t12JksAvg,data2JkAvg,t12RatioAvg,12)
    t14RatioErr=avgJackknifeRatioErr(t14JksAvg,data2JkAvg,t14RatioAvg,14)

# MAKE TIME ARRAYS FOR PLOTTING
D6time=makeTimeArrays(6)
D8time=makeTimeArrays(8)
D10time=makeTimeArrays(10)
if options.cutNoisy == 0:
    D12time=makeTimeArrays(12)
    D14time=makeTimeArrays(14)
    D16time=makeTimeArrays(16)
elif options.cutNoisy == 2:
    D12time=makeTimeArrays(12)
    D14time=makeTimeArrays(14)

bareTime=np.linspace(-9,9,500)
bareCharge=np.zeros(500)
bareChargeError=np.zeros(500)


######################################################################
# DETERMINE PREFACTORS/KINEMATIC FACTORS THAT NEED TO BE INCLUDED
######################################################################
redstarFactor=np.sqrt(2) # to account for redstar isovector/scalar currents
kinefactor=1.0
if srcMom.split('.')[2] < 4:
    gndmass=np.sqrt((E0**2) - ((2*np.pi/32.0)*(momenta["Src"]["val"][0]+momenta["Src"]["val"][1]+momenta["Src"]["val"][2]))**2)
else:
    gndmass=0.4757

momz=momenta["Src"]["val"][2]

print "GND MASS = %.7f"%gndmass
numeratorE=4*E0*(E0+gndmass)# E0 IS FITTED NUCLEON ENERGY, NOT JUST GND MASS
if gamma==0:
    kinefactor=numeratorE/(4.0*gndmass*(gndmass+E0))
if gamma==3:
    kinefactor=numeratorE/(4.0*(gndmass*E0+gndmass**2)) # not so certain on the factor of 4 in denom here (10/12/19)
if gamma==11:
    kinefactor=1.0*numeratorE/(4.0*(gndmass*E0+gndmass**2+(((2*np.pi)/(32.0))*momz)**2))
if gamma==7:
    kinefactor=-1.0*numeratorE/(4.0*(gndmass+E0)*(((2*np.pi)/(32.0))*momz))
    # kinefactor=-1.0/(4*(((2*np.pi)/(1.0*32))*momz)*E0) # produces charges ~0.68
if gamma==8:
    kinefactor=numeratorE/(4.0*E0*(E0+gndmass))
if gamma==4:
    kinefactor=numeratorE/(4.0*(E0+gndmass)*(((2*np.pi)/(32.0))*momz))

print "Kinematic factor applied to ratio = %.7f"%kinefactor
########################################################################


# Set the asymptotic bare/renormalized charge and error
bareCharge[:]=(A/a)*kinefactor*redstarFactor
# if bareCharge[0] < 0:
#     print "WARNING: DOING SOME TOM FOOLERY HERE!!!"
#     bareCharge[:]*=-1.0 #BAD!
#     redstarFactor*=-1.0 #BAD!

# bareChargeError[:]=np.sqrt(np.power(dA/a,2)+np.power((A/(a**2)),2)*np.power(da,2)-((2*A)/(np.power(a,3)))*covAa)*kinefactor*redstarFactor

bareChargeError[:]=np.abs(np.sqrt((1.0/a)*COV[('A', 'A')]*(1.0/a)+(1.0/a)*COV[('A', 'a')]*(-A/(a**2))+(-A/(a**2))*COV[('a', 'A')]*(1.0/a)+(-A/(a**2))*COV[('a', 'a')]*(-A/(a**2)))*kinefactor*redstarFactor)


renormCharge=renormalizations[gamma]*bareCharge
renormChargeError=renormalizations[gamma]*bareChargeError

print "Bare charge = %.7f"%bareCharge[0]
print "Bare charge error = %.7f"%bareChargeError[0]
print "Renorm charge = %.7f"%renormCharge[0]
print "Renorm charge error = %.7f"%renormChargeError[0]

time6=np.linspace(-0.5,5.5,500)
time8=np.linspace(-0.5,7.5,500)
time10=np.linspace(-0.5,9.5,500)
time12=np.linspace(-0.5,11.5,500)
time14=np.linspace(-0.5,13.5,500)
time16=np.linspace(-0.5,15.5,500)
tsep6Fit=np.zeros(500)
tsep6FitError=np.zeros(500)
tsep8Fit=np.zeros(500)
tsep8FitError=np.zeros(500)
tsep10Fit=np.zeros(500)
tsep10FitError=np.zeros(500)
tsep12Fit=np.zeros(500)
tsep12FitError=np.zeros(500)
tsep14Fit=np.zeros(500)
tsep14FitError=np.zeros(500)
tsep16Fit=np.zeros(500)
tsep16FitError=np.zeros(500)
# Define the 3pt/2pt fitted ratios
for i in range(0,500):
    tsep6Fit[i]=fit3pt2ptRatio(params,time6[i],6)*renormalizations[gamma]*kinefactor*redstarFactor
    tsep8Fit[i]=fit3pt2ptRatio(params,time8[i],8)*renormalizations[gamma]*kinefactor*redstarFactor
    tsep10Fit[i]=fit3pt2ptRatio(params,time10[i],10)*renormalizations[gamma]*kinefactor*redstarFactor
    if options.cutNoisy == 0:
        tsep12Fit[i]=fit3pt2ptRatio(params,time12[i],12)*renormalizations[gamma]*kinefactor*redstarFactor
        tsep14Fit[i]=fit3pt2ptRatio(params,time14[i],14)*renormalizations[gamma]*kinefactor*redstarFactor
        tsep16Fit[i]=fit3pt2ptRatio(params,time16[i],16)*renormalizations[gamma]*kinefactor*redstarFactor
    elif options.cutNoisy == 2:
        tsep12Fit[i]=fit3pt2ptRatio(params,time12[i],12)*renormalizations[gamma]*kinefactor*redstarFactor
        tsep14Fit[i]=fit3pt2ptRatio(params,time14[i],14)*renormalizations[gamma]*kinefactor*redstarFactor


# Define the 3pt/2pt fitted ratio errors
for i in range(0,500):
    tsep6FitError[i]=Deriv(time6[i],6,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
    tsep8FitError[i]=Deriv(time8[i],8,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
    tsep10FitError[i]=Deriv(time10[i],10,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
    if options.cutNoisy == 0:
        tsep12FitError[i]=Deriv(time12[i],12,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
        tsep14FitError[i]=Deriv(time14[i],14,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
        tsep16FitError[i]=Deriv(time16[i],16,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
    elif options.cutNoisy == 2:
        tsep12FitError[i]=Deriv(time12[i],12,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
        tsep14FitError[i]=Deriv(time14[i],14,covariance)*renormalizations[gamma]*kinefactor*redstarFactor


# Determine the effective charges (3pt/2pt fit)
effChT6=makeEffCharge(t6JksAvg,6)*kinefactor*renormalizations[gamma]*redstarFactor
effChT8=makeEffCharge(t8JksAvg,8)*kinefactor*renormalizations[gamma]*redstarFactor
effChT10=makeEffCharge(t10JksAvg,10)*kinefactor*renormalizations[gamma]*redstarFactor
if options.cutNoisy == 0:
    effChT12=makeEffCharge(t12JksAvg,12)*kinefactor*renormalizations[gamma]*redstarFactor
    effChT14=makeEffCharge(t14JksAvg,14)*kinefactor*renormalizations[gamma]*redstarFactor
    effChT16=makeEffCharge(t16JksAvg,16)*kinefactor*renormalizations[gamma]*redstarFactor
elif options.cutNoisy == 2:
    effChT12=makeEffCharge(t12JksAvg,12)*kinefactor*renormalizations[gamma]*redstarFactor
    effChT14=makeEffCharge(t14JksAvg,14)*kinefactor*renormalizations[gamma]*redstarFactor

# Determine the error in effective charges (3pt/2pt fit)
# effCh* right above was determined prior to rescaling by all factors
# Thus in determining error on the effective charges, the JkAvgs must be rescaled
effChErrT6=makeEffChargeError(t6JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT6,6)
effChErrT8=makeEffChargeError(t8JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT8,8)
effChErrT10=makeEffChargeError(t10JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT10,10)
if options.cutNoisy == 0:
    effChErrT12=makeEffChargeError(t12JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT12,12)
    effChErrT14=makeEffChargeError(t14JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT14,14)
    effChErrT16=makeEffChargeError(t16JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT16,16)
elif options.cutNoisy == 2:
    effChErrT12=makeEffChargeError(t12JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT12,12)
    effChErrT14=makeEffChargeError(t14JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT14,14)


# textstr=''
# title_str=operator2+" -- "+operator2

#########################################################################
# POTENTIALLY DETERMINE/PLOT THE EFFECTIVE 2PT EFFECTIVE ENERGY AND FIT
#########################################################################
if options.show2ptFitToData == 1:
    time2pt=np.linspace(0,20,3000)
    fitRange2pt=np.linspace(float(options.fitRange2pt.split('.')[0]),\
                            float(options.fitRange2pt.split('.')[2]),3000)

    # Determine the central value of fit to effective energy
    # Also write out the central values
    W=open("effE_m0p2390_fit2pt_p"+str(options.pi.split('.')[0])+str(options.pi.split('.')[1])+str(options.pi.split('.')[2])+".dat",'w')
    fit2ptEffE=np.linspace(0,20,3000)
    for i in range(0,len(time2pt)-1):
        fit2ptEffE[i]=(1.0/(time2pt[1]-time2pt[0]))*np.log(twoState(params['a'],params['b'],params['E0'],params['E1'],time2pt[i])/twoState(params['a'],params['b'],params['E0'],params['E1'],time2pt[i+1]))
        
    fit2ptEffEError=np.linspace(0,20,3000)
    for i in range(0,len(time2pt)-1):
        fit2ptEffEError[i]=Deriv2pt(time2pt[i],time2pt[1]-time2pt[0],covariance)
        W.write(str(time2pt[i])+" "+str(fit2ptEffE[i])+" "+str(fit2ptEffEError[i])+"\n")
    W.close()


    plt.rc('text', usetex=True)
    plt.xlabel(r't_{\rm{sep}}',fontsize=14)
    plt.ylabel(r'E_{\rm{eff}}', fontsize=14)
    plt.plot(time2pt,fit2ptEffE[:],'b',label=r'$2-\text{state Fit}')#,markersize=4,markerfacecolor="Blue",markeredgecolor='b')
    plt.fill_between(time2pt,fit2ptEffE+fit2ptEffEError,fit2ptEffE-fit2ptEffEError,color='#B2E5FF',alpha=0.8)
    plt.ylim([0,1.1])
    plt.show()
    sys.exit()

########################################################
# PLOT THE RENORMALIZED EFFECTIVE CHARGES WITH FITS
########################################################
plt.rc('text', usetex=True)
# title_str = "Nucleon effective "+charge+r" Charge - Isoclover $32^3\times64$    $m_\pi=%d$ MeV    $a=%.3f$ fm"%(356,0.098)+"\n"+title_str
# plt.title(title_str)
plt.xlabel(r'$\tau-T/2$',fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
# plt.figtext(0.5, 0.2, textstr, fontsize=10, verticalalignment='baseline', horizontalalignment='center')

# For plot clarity
offset=0.1
# plt.rc('text', usetex=True)




# Method to manage plotting of either effective charges or 3pt/2pt ratio of raw lattice data
def plotData(pt,effCh,effChE):
    buffMin=2
    params={ 6: {'buffMax': 5, 'Mark': 'o', 'Color': 'c', 'Label': r'$T/a=6$',\
                 'Shift': -3*offset },
             8: {'buffMax': 7, 'Mark': '^', 'Color': 'b', 'Label': r'$T/a=8$',\
                 'Shift': -2*offset},
             10: {'buffMax': 9, 'Mark': '*', 'Color': 'm', 'Label': r'$T/a=10$',\
                  'Shift': -offset},
             12: {'buffMax': 11, 'Mark': 's', 'Color': 'g', 'Label': r'$T/a=12$', 'Shift': 0 },
             14: {'buffMax': 13, 'Mark': 'h', 'Color': 'r', 'Label': r'$T/a=14$',\
                  'Shift': offset },
             16: {'buffMax': 15, 'Mark': 'v', 'Color': 'k', 'Label': r'$T/a=16$',\
                  'Shift': 2*offset } }

    nn=len(pt)
    # Do things based on the times series
    plt.plot(pt[:buffMin]+params[nn]['Shift'],effCh[:buffMin],params[nn]['Mark'],\
             markersize=3,mfc='grey',mec='grey')
    plt.errorbar(pt[:buffMin]+params[nn]['Shift'],effCh[:buffMin],yerr=effChE[:buffMin],\
                 fmt=params[nn]['Mark'],label='',ecolor='grey',mfc='grey',mec='grey')

    plt.plot(pt[buffMin:params[nn]['buffMax']]+params[nn]['Shift'],\
             effCh[buffMin:params[nn]['buffMax']],\
             params[nn]['Mark'],label=params[nn]['Label'],markersize=3,\
             mfc=params[nn]['Color'],mec=params[nn]['Color'])
    plt.errorbar(pt[buffMin:params[nn]['buffMax']]+params[nn]['Shift'],\
                 effCh[buffMin:params[nn]['buffMax']],\
                 yerr=effChE[buffMin:params[nn]['buffMax']],fmt=params[nn]['Mark'],label='',\
                 ecolor=params[nn]['Color'],mfc=params[nn]['Color'],mec=params[nn]['Color'])

    plt.plot(pt[params[nn]['buffMax']:]+params[nn]['Shift'],effCh[params[nn]['buffMax']:],\
             params[nn]['Mark'],markersize=3,mfc='grey',mec='grey')
    plt.errorbar(pt[params[nn]['buffMax']:]+params[nn]['Shift'],effCh[params[nn]['buffMax']:],\
                 yerr=effChE[params[nn]['buffMax']:],fmt=params[nn]['Mark'],label='',\
                 ecolor='grey',mfc='grey',mec='grey')
    

###########################################
# PLOT EITHER RAW 3PT/2PT LATTICE RATIO
# OR 3PT/(2PT FIT) - I.E. EFFECTIVE CHARGE
###########################################
dataTimeList=[D6time,D8time,D10time]
dataList=dataErrList=""
if options.rawRatio:
    dataList=[t6RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
              t8RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
              t10RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma]]
    dataErrList=[t6RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
                 t8RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
                 t10RatioErr*redstarFactor*kinefactor*renormalizations[gamma]]
    if options.cutNoisy == 0:
        dataList.extend([t12RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
                         t14RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
                         t16RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma]])
        dataErrList.extend([t12RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
                            t14RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
                            t16RatioErr*redstarFactor*kinefactor*renormalizations[gamma]])
        dataTimeList.extend([D12time,D14time,D16time])
    elif options.cutNoisy == 2:
        dataList.extend([t12RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
                         t14RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma]])
        dataErrList.extend([t12RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
                            t14RatioErr*redstarFactor*kinefactor*renormalizations[gamma]])
        dataTimeList.extend([D12time,D14time])
else:
    dataList=[effChT6,effChT8,effChT10]
    dataErrList=[effChErrT6,effChErrT8,effChErrT10]
    if options.cutNoisy == 0:
        dataList.extend([effChT12,effChT14,effChT16])
        dataErrList.extend([effChErrT12,effChErrT14,effChErrT16])
        dataTimeList.extend([D12time,D14time,D16time])
    elif options.cutNoisy == 2:
        dataList.extend([effChT12,effChT14])
        dataErrList.extend([effChErrT12,effChErrT14])
        dataTimeList.extend([D12time,D14time]) 


################################
# PLOT THE SELECTED DATA
################################
for n, k in enumerate(dataTimeList):
    plotData(k,dataList[n],dataErrList[n])



YMIN=YMAX=0.0
if gamma==0:
    plt.ylabel(r'$g_S^{u-d}$',fontsize=20)
    YMIN=0.45
    YMAX=1.3
    # # For unphased pz=4 data
    # YMIN=-3
    # YMAX=4
    # # For phased pz=4 data
    YMIN=-1
    YMAX=2.25
if gamma==1:
    plt.ylabel(r'$g_{V_1}^{u-d}$',fontsize=20)
if gamma==2:
    plt.ylabel(r'$g_{V_2}^{u-d}$',fontsize=20)
if gamma==4:
    plt.ylabel(r'$g_{V_3}^{u-d}$',fontsize=20)
    YMIN=0.7
    YMAX=1.10
    # # For the unphased data at pz=4
    # YMIN=-0.2
    # YMAX=2.25
    # # For the phased data at pz=4
    # YMIN=0.94
    # YMAX=1.15
    # # For the phased data at pz=4 (3pt/2pt raw ratio)
    # YMIN=0.25
    # YMAX=1.5
if gamma==8:
    plt.ylabel(r'$g_{V_4}^{u-d}$',fontsize=20)
    YMIN=0.94
    YMAX=1.10
    # # For the unphased data at pz=4
    # YMIN=0.5
    # YMAX=1.3
if gamma==7:
    plt.ylabel(r'$g_{A_4}^{u-d}$',fontsize=20)
    YMIN=0.75
    YMAX=1.35
    # # For the unphased pz=4 data
    # YMIN=0
    # YMAX=2
    # # For the unphased pz=4 data (3pt/2pt raw ratio)
    # YMIN=0.25
    # YMAX=1.75
    # # For the phased pz=4 data
    # YMIN=1.1
    # YMAX=1.5
if gamma==11:
    plt.ylabel(r'$g_{A_3}^{u-d}$',fontsize=20)
    YMIN=1.0
    YMAX=1.35
    # # For unphased pz=4 data
    # YMIN=0
    # YMAX=2
    # # For phased pz=4 data
    # YMIN=1.1
    # YMAX=1.55
if gamma==3:
    plt.ylabel(r'$g_{T_{12}}^{u-d}$',fontsize=20)
    YMIN=0.95
    YMAX=1.3
    # # For unphased pz=4 data
    YMIN=0.25
    YMAX=2.25
if gamma==15:
    plt.ylabel(r'$g_P^{u-d}$',fontsize=20)
    YMIN=-0.4
    YMAX=0.4
if gamma==15 and options.srcRow!=options.snkRow:
    plt.ylabel(r'$g_P^{u-d}$',fontsize=20)
    YMIN=-1.2
    YMAX=1.2

# # Let's modify the vertical range permanently for the poor phased results
# YMIN=0.2
# YMAX=1.5


# Extracted charge was obtained from fitting 3pt data that was already rescaled (sqrt(2) for Redstar generated files)
# plt.plot(bareTime,bareCharge,'black')
# plt.fill_between(bareTime,bareCharge+bareChargeError,bareCharge-bareChargeError,color='#d8dcd6')
plt.plot(bareTime,renormCharge,'black')
plt.fill_between(bareTime,renormCharge+renormChargeError,renormCharge-renormChargeError,color='#d8dcd6')
plt.xlim([-9,9])
plt.xlim([-7.5,7.5])
plt.ylim([YMIN,YMAX])



# Include the fits to the 3pt/2pt data
plt.plot(time6[:]-3-3*offset, tsep6Fit[:],'c')
plt.fill_between(time6[:]-3-3*offset,tsep6Fit[:]+tsep6FitError[:],
                 tsep6Fit[:]-tsep6FitError[:],color='#13EAC9',alpha=0.6,
                 where=np.logical_and(time6[:]>=2,time6[:]<=4))
plt.fill_between(time6[:]-3-3*offset,tsep6Fit[:]+tsep6FitError[:],
                 tsep6Fit[:]-tsep6FitError[:],color='grey',alpha=0.1,
                 where=np.logical_or(time6[:]<2,time6[:]>4))

plt.plot(time8[:]-4-2*offset, tsep8Fit[:],'b')
plt.fill_between(time8[:]-4-2*offset,tsep8Fit[:]+tsep8FitError[:],
                 tsep8Fit[:]-tsep8FitError[:],color='b',alpha=0.2,# '#B2E5FF'
                 where=np.logical_and(time8[:]>=2,time8[:]<=6))
plt.fill_between(time8[:]-4-2*offset,tsep8Fit[:]+tsep8FitError[:],
                 tsep8Fit[:]-tsep8FitError[:],color='grey',alpha=0.1,
                 where=np.logical_or(time8[:]<2,time8[:]>6))

plt.plot(time10[:]-5-offset, tsep10Fit,'m')
plt.fill_between(time10[:]-5-offset,tsep10Fit[:]+tsep10FitError[:],
                 tsep10Fit[:]-tsep10FitError[:],color='#D767AD',alpha=0.4,
                 where=np.logical_and(time10[:]>=2,time10[:]<=8))
plt.fill_between(time10[:]-5-offset,tsep10Fit[:]+tsep10FitError[:],
                 tsep10Fit[:]-tsep10FitError[:],color='grey',alpha=0.1,
                 where=np.logical_or(time10[:]<2,time10[:]>8))

if options.cutNoisy == 0:
    plt.plot(time12[:]-6, tsep12Fit,'g')
    plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
                     tsep12Fit[:]-tsep12FitError[:],color='g',alpha=0.2,#BCFFB8
                     where=np.logical_and(time12[:]>=2,time12[:]<=10))
    plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
                     tsep12Fit[:]-tsep12FitError[:],color='grey',alpha=0.1,
                     where=np.logical_or(time12[:]<2,time12[:]>10))

    plt.plot(time14[:]-7+offset, tsep14Fit,'r')
    plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
                     tsep14Fit[:]-tsep14FitError[:],color='r',alpha=0.2,#FDD2D2
                     where=np.logical_and(time14[:]>=2,time14[:]<=12))
    plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
                     tsep14Fit[:]-tsep14FitError[:],color='grey',alpha=0.1,
                     where=np.logical_or(time14[:]<2,time14[:]>12))

    plt.plot(time16[:]-8+2*offset, tsep16Fit,'k')
    plt.fill_between(time16[:]-8+2*offset,tsep16Fit[:]+tsep16FitError[:],
                     tsep16Fit[:]-tsep16FitError[:],color='#516572',alpha=0.6,
                     where=np.logical_and(time16[:]>=2,time16[:]<=14))
    plt.fill_between(time16[:]-8+2*offset,tsep16Fit[:]+tsep16FitError[:],
                     tsep16Fit[:]-tsep16FitError[:],color='grey',alpha=0.1,
                     where=np.logical_or(time16[:]<2,time16[:]>14))

elif options.cutNoisy == 2:
    plt.plot(time12[:]-6, tsep12Fit,'g')
    plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
                     tsep12Fit[:]-tsep12FitError[:],color='g',alpha=0.2,#BCFFB8
                     where=np.logical_and(time12[:]>=2,time12[:]<=10))
    plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
                     tsep12Fit[:]-tsep12FitError[:],color='grey',alpha=0.1,
                     where=np.logical_or(time12[:]<2,time12[:]>10))

    plt.plot(time14[:]-7+offset, tsep14Fit,'r')
    plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
                     tsep14Fit[:]-tsep14FitError[:],color='r',alpha=0.2,#FDD2D2
                     where=np.logical_and(time14[:]>=2,time14[:]<=12))
    plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
                     tsep14Fit[:]-tsep14FitError[:],color='grey',alpha=0.1,
                     where=np.logical_or(time14[:]<2,time14[:]>12))


# # Lastly, add a text box with the reduced chi^2 of simultaneous fit
# txtbxstr='$\\chi_r^2$ = %.4f'%chi2
# plt.text(0.00,(YMAX-YMIN)*0.05+YMIN,txtbxstr,horizontalalignment='center',fontsize=14)

ax = plt.gca()
legend = ax.legend(fontsize=16,ncol=3,loc=9,markerscale=2,fancybox=True)
# legend.get_frame().set_facecolor('b')
legend.get_frame().set_alpha(0.5)
# ax.legend(framealpha=0.5)
# ax.legend(fontsize=14,ncol=3,loc=9,markerscale=2,fancybox=True,framealpha=0.5)

# plt.legend(fontsize=14,ncol=3,loc=9,markerscale=2)#,fancybox=True,framealpha=0.5)
plt.savefig(output_name,dpi=500)
plt.show()

