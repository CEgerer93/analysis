#!/usr/bin/python

import numpy as np
import sys

########################################################
# Dictionary associating Chroma Gamma with string label
########################################################
gammaDict={0: '$\mathbb{1}$', 1: '$\gamma_x$', 2: '$\gamma_y$', 3: '$\gamma_x\gamma_y$', \
           4: '$\gamma_z$', 5: '$\gamma_x\gamma_z', 6: '$\gamma_y\gamma_z$', 7: '$\gamma_4\gamma_5$', \
           8: '$\gamma_4$', 11: '$\gamma_z\gamma_5$', 13: '$\gamma_y\gamma_5$', 14:'$\gamma_x\gamma_5$'}

##########################################################################
# Dictionary for converting gammas in circular basis into Cartesian basis
##########################################################################
crtsnGamma={1: { 'norm': 1/np.sqrt(2), 'firstRow': 1.0, 'secondRow': -1.0 }, \
            2: { 'norm': 1/np.sqrt(2), 'firstRow': 1.0, 'secondRow': 1.0 }, \
            5: { 'norm': 1/(2*np.sqrt(2)), 'firstRow': 1.0, 'secondRow': 1.0 }, \
            6: { 'norm': 1/(2*np.sqrt(2)), 'firstRow': 1.0, 'secondRow': -1.0}, \
            13: { 'norm': 1/np.sqrt(2), 'firstRow': 1.0, 'secondRow': 1.0}, \
            14: { 'norm': 1/np.sqrt(2), 'firstRow': 1.0, 'secondRow': -1.0} }


##################################
# Master dictionary for plotting
##################################
pltoffset=0.1
pltParams={ 4: { 'mark': 'o', 'color': 'c', 'label': r'$T/a=4', 'shift': -3*pltoffset },
            6: { 'mark': '^', 'color': 'b', 'label': r'$T/a=6', 'shift': -2*pltoffset },
            8: { 'mark': '*', 'color': 'm', 'label': r'$T/a=8', 'shift': -pltoffset },
            10: { 'mark': 's', 'color': 'g', 'label': r'$T/a=10', 'shift': 0 },
            12: { 'mark': 'h', 'color': 'r', 'label': r'$T/a=12', 'shift': pltoffset },
            14: { 'mark': 'v', 'color': 'k', 'label': r'$T/a=14', 'shift': 2*pltoffset } }


class correlator:
    def __init__(self,cfgs,T,corr,complexity,rdstrFact):
        self.cfgs=cfgs
        self.T=T
        self.corr=corr
        self.complexity=complexity
        self.jks=np.zeros((cfgs,cfgs-1,T))
        self.avgjk=np.zeros((cfgs,T))
        self.rdstrFact=rdstrFact


    # Form jackknife samples of corr array
    def makeJks(self):
        for n in range(0,self.cfgs):
            if n > 0:
                for l in range(0,n):
                    for t in range(0,self.T):
                        self.jks[n,l,t]=self.corr[t+l*self.T,self.complexity]*self.rdstrFact
            for l in range(n,self.cfgs-1):
                for t in range(0,self.T):
                    self.jks[n,l,t]=self.corr[t+(l+1)*self.T,self.complexity]*self.rdstrFact

    # Form average per jackknife sample
    def avgJks(self):
        for n in range(0,self.cfgs):
            for t in range(0,self.T):
                dum=0.0
                for g in range(0,self.cfgs-1):
                    dum+=self.jks[n,g,t]
                self.avgjk[n,t]=(1.0/(1.0*(self.cfgs-1)))*dum

    
class ratio:
    def __init__(self,numerator,src2pt,snk2pt):
        self.N=numerator
        self.Src=src2pt
        self.Snk=snk2pt
        self.Avg=np.zeros(numerator.T)
        self.Err=np.zeros(numerator.T)
        self.pltParams=pltParams[numerator.T]
        
        
    # Form ratio per jackknife average
    def avgratio(self):
        for t in range(0,self.N.T):
            dum=0.0
            for n in range(0,self.N.cfgs):
                dum+=(self.N.avgjk[n,t]/self.Snk.avgjk[n,self.N.T])*\
                    np.sqrt((self.Src.avgjk[n,self.N.T-t]*self.Snk.avgjk[n,t]*self.Snk.avgjk[n,self.N.T])\
                    /(self.Snk.avgjk[n,self.N.T-t]*self.Src.avgjk[n,t]*self.Src.avgjk[n,self.N.T]))
            self.Avg[t]=(1.0/(1.0*self.N.cfgs))*dum

    # Form ratio error per jackknife average
    def avgratioErr(self):
        for t in range(0,self.N.T):
            dum=0.0
            for n in range(0,self.N.cfgs):
                dum+=np.power(((self.N.avgjk[n,t]/self.Snk.avgjk[n,self.N.T])*\
                               np.sqrt((self.Src.avgjk[n,self.N.T-t]*self.Snk.avgjk[n,t]*\
                                        self.Snk.avgjk[n,self.N.T])/\
                                       (self.Snk.avgjk[n,self.N.T-t]*\
                                        self.Src.avgjk[n,t]*self.Src.avgjk[n,self.N.T])))-self.Avg[t],2)
            self.Err[t]=np.sqrt(((1.0*(self.N.cfgs-1))/(1.0*self.N.cfgs))*dum)


    # Plot this ratio
    def plotRatio(self,ts,ax):
        ax.errorbar(ts+self.pltParams['shift'],self.Avg,self.Err,fmt=self.pltParams['mark'],
                    label=self.pltParams['label'],
                    markersize=5,mfc=self.pltParams['color'],mec=self.pltParams['color'],
                    ecolor=self.pltParams['color'])





# Utility to merge correlators in a circular basis to extract gammas in cartesian basis
def mergeCorrs(c1, c2, g):

    # Copy one of the correlators for timing info and sizing
    merge=c1

    if len(c1[:,0]) != len(c2[:,0]):
        print "Correlators to merge are not of the same size = cfgs x tsep"
        sys.exit()
        
    for n, d in enumerate(c1):
        
        merge[n,1]=crtsnGamma[g]['norm']*(crtsnGamma[g]['firstRow']*d[1]+crtsnGamma[g]['secondRow']*c2[n,1])
        merge[n,2]=crtsnGamma[g]['norm']*(crtsnGamma[g]['firstRow']*d[2]+crtsnGamma[g]['secondRow']*c2[n,2])

    return merge



    # # Add momenta/displacement to figure
    # def addMomsZs(self,ax,snkMom,srcMom,zsep):
    #     ax.text(-7,1.2*ax.get_ylim()[0],r'$p_f=\left(%d,%d,%d\right)$'\
    #             %(snkMom.split('.')[0],snkMom.split('.')[1],snkMom.split('.')[2]),fontsize=13)
    #     ax.text(-2,1.2*ax.get_ylim()[0],r'$p_i=\left(%d,%d,%d\right)$'\
    #             %(srcMom.split('.')[0],srcMom.split('.')[1],srcMom.split('.')[2]),fontsize=13)
    #     ax.text(3,1.2*ax.get_ylim()[0],r'$\vec{z}=\left(%d,%d,%d\right)$'\
    #             %(zsep.split('.')[0],zsep.split('.')[1],zsep.split('.')[2]),fontsize=13)

    
        

# class tsep:
#     def __init__(self, numT ):
#         self.numT = numT
#         self.times = None
#         self.correlator = correlator
#         self.effMat = effMat
#         self.effMatErr = effMatErr
#         self.avgJks = None
#         self.avgJkRatio = None
#         self.avgJkRatioErr = None
        

#     def makeTimeArrays(self.N):
#         time=np.zeros(N)
#         for i in range(0,N):
#             time[i]=i-N/2
            
#         return time



# class xmbf2ptSimulFit:
#     def __init__(self, modelParams, tau, ensemCorrs):
#         self.K = modelParams[0] # number of model functions to simultaneously fit
#         self.V = modelParams[1] # number of indep variables in data set
#         self.N = modelParams[2] # number of gauge configs in this data set
#         self.comp = modelParams[3] # real/imag to consider
#         self.tau = tau # src/snk separation array
#         self.ensemCorrs = ensemCorrs
#         self.Times = [[],[]]
#         self.Times[0] = [self.tau[i]+1 for i in range(0,len(self.tau))]
#         self.Times[1] = [self.tau[i] for i in range(0,len(self.tau))]
#         self.D = None #[] # Empty array to hold all npt data

#     # Populate data arrays
#     def popDat(self):
#         for f in self.ensemCorrs:
#             # self.D.append(np.loadtxt(f, delimiter=' ', skiprows=1))
#             self.D = np.loadtxt(f, delimiter=' ', skiprows=1)

#     # Write parsed data to file
#     def Write(self,W):
#         W.write(str(self.K)+"\n")
#         W.write(str(self.V)+"\n")
#         W.write(str(len(self.tau))+"\n")
#         W.write(str(self.N)+"\n")
#         # Print time series
#         for i in range(0,len(self.tau)):
#             W.write(str(i+1)+" "+str(int(self.tau[i]))+"\n")
            

#         # Create blocks of data
#         for j in range(0,self.N): # gauge config loop
#             for i in range(0,len(self.tau)):
#                 # Write gauge config, element of time series, real/imaginary correlator to file
#                 # DATA MUST BE WRITTEN IN SCIENTIFIC NOTATION!!!!!!!!!
#                 W.write(str(j+1)+" "+str(i+1)+" "+str("%.8e"%(self.D[j*len(self.tau)+i,self.comp]))+"\n")


# class xmbf3ptSimulFit:
#     def __init__(self, modelParams, tau, ensemCorrs):
#         self.K = modelParams[0] # number of model functions to simultaneously fit
#         self.V = modelParams[1] # number of indep variables in data set
#         self.N = modelParams[2] # number of gauge configs in this data set
#         self.comp = modelParams[3] # real/imag to consider
#         self.tau = tau # current insertion time array
#         self.Max_t = max(self.tau)*len(self.tau)
#         self.ensemCorrs = ensemCorrs
#         self.Times = -np.ones((2,self.Max_t))
#         self.D = [] # Empty array to hold all npt data

#     # Make the running time series for simultaneous fitting
#     def fill_time_series(self):
#         cntr=0
#         # Loop over the src/snk times
#         for ts in self.tau:
#             # Loop over the insertions times
#             for tins in range(0,ts):
#                 self.Times[0,cntr]=tins
#                 self.Times[1,cntr]=ts
#                 cntr+=1 # Increment the line counter
#                 # Skip (max tsep) - t rows in tau array
#                 if tins==ts-1:
#                     cntr+=(max(self.tau)-tins-1)


#     # Populate data arrays
#     def popDat(self):
#         for f in self.ensemCorrs:
#             self.D.append(np.loadtxt(f, delimiter=' ', skiprows=1))

#     # Write parsed data to file
#     def Write(self,W):
#         W.write(str(self.K)+"\n")
#         W.write(str(self.V)+"\n")
#         W.write(str(self.Max_t)+"\n")
#         W.write(str(self.N)+"\n")
#         # Print time series
#         for i in range(1,self.Max_t+1):
#             W.write(str(i)+" "+str(int(self.Times[0,i-1]))+" "+str(int(self.Times[1,i-1]))+"\n")
            

#         # Create blocks of data
#         for j in range(0,self.N): # gauge config loop
#             dum=1
#             for k in range(0,len(self.D)):
#                 while dum < (self.tau[k]+1+k*max(self.tau)):
#                     # Write gauge config, element of time series, real/imaginary correlator to file
#                     # DATA MUST BE WRITTEN IN SCIENTIFIC NOTATION!!!!!!!!!
#                     W.write(str(j+1)+" "+str(dum)+" "+str("%.8e"%(self.D[k][j*self.tau[k]+self.Times[0,dum-1],\
# self.comp]))+"\n")
#                     dum+=1
#                 else:
#                     while dum <= (k+1)*max(self.tau):
#                         W.write(str(j+1)+" "+str(dum)+" "+str("%.8e"%0.00+"\n"))
#                         dum+=1

