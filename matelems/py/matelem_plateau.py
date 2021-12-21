#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import pylab # to save figures to file
import sys,optparse
from util import *

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-c", "--cfgs", type="int", default=1,
                  help='Number of configurations (default = 1)')
parser.add_option("-x", "--complexity", type="int", default=0,
                  help='Component of Npt function to analyze [Real=1 Imag=2] (default = 1)')
parser.add_option("-i", "--pi", type="string", default="0.0.0",
                  help='Source momentum <pix>.<piy>.<piz> (default = 0.0.0)')
parser.add_option("-q", "--Q", type="string", default="0.0.0",
                  help='Momentum transfer <qx>.<qy>.<qz> (default = 0.0.0)')
parser.add_option("-f", "--pf", type="string", default="0.0.0",
                  help='Sink momentum <pfx>.<pfy>.<pfz> (default = 0.0.0)')
parser.add_option("-z", "--vecZ", type="string", default="0.0.0",
                  help='Displacement vector <x>.<y>.<z> (default = 0.0.0)')
parser.add_option("-N", "--nptCorrsFile", type="string", default="",
                  help='List file containing full path/name for Npt correlator(s) (default="")')
parser.add_option("-r", "--twoptSrcCorrsFile", type="string", default="",
                  help='File containing full path/name for Src 2pt correlator (default="")')
parser.add_option("-k", "--twoptSnkCorrsFile", type="string", default="",
                  help='File containing full path/name for Snk 2pt correlator (default="")')
parser.add_option("-g", "--chromaGamma", type="int", default="-1",
                  help='Gamma matrix to consider (default = -1)')
parser.add_option("-o", "--output", type="string", default="",
                  help='Override name of plateau plot\n (default = "pgitd_plateau_pf<snkX><snkY><snkZ>_pi<srcX><srcY><srcZ>_rows<snkRow><srcRow>_g<gamma>")')
parser.add_option("-w", "--srcRow", type="string", default="X",
                  help='Source interpolator row (default = X)')
parser.add_option("-y", "--snkRow", type="string", default="X",
                  help='Sink interpolator row (default = X)')
parser.add_option("-l", "--lightBkgd", type="int", default=1,
                  help='Format figs for light (1) or dark (0) background (default = 1)')
parser.add_option("-s", "--showFig", type="int", default=0,
		  help='Show the figure (default = 0)')


# Parse the input arguments
(options, args) = parser.parse_args()

################
# INITIALIZE GLOBAL PROPERTIES OF FIGURES
################
# Finalize the figures
plt.rc('text', usetex=True)
plt.rcParams["mathtext.fontset"]="stix"
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rc('xtick.major',size=5)
plt.rc('ytick.major',size=5)
plt.rc('xtick',labelsize=18)
plt.rc('ytick',labelsize=18)
plt.rc('axes',labelsize=18)
truthTransparent=False
FrameAlpha=1
# Optionally swap default black labels for white
if options.lightBkgd == 0:
    truthTransparent=True
    plt.rcParams['text.color'] = 'white'
    plt.rcParams['axes.edgecolor'] = 'white'
    plt.rcParams['legend.frameon' ] = False
    plt.rc('axes',edgecolor='white')
    plt.rc('axes',labelcolor='white')
    plt.rc('xtick',color='white')
    plt.rc('ytick',color='white')
    plt.rc('text',color='white')
    FrameAlpha=0


cfgs=options.cfgs
complexity=options.complexity
gamma=options.chromaGamma
srcMom=options.pi
snkMom=options.pf
insMom=options.Q
vecZ=options.vecZ

redstarFactor=np.sqrt(2)

# Grab and read the Src 2pt correlator
with open(options.twoptSrcCorrsFile) as pt:
    twoptSrcCorr=pt.read().rstrip()
corrSrc2pt=np.loadtxt(twoptSrcCorr, delimiter=' ', skiprows=1)
# Grab and read the Snk 2pt correlator
with open(options.twoptSnkCorrsFile) as pt:
    twoptSnkCorr=pt.read().rstrip()
corrSnk2pt=np.loadtxt(twoptSnkCorr, delimiter=' ', skiprows=1)

#############################################
#       Grab and read Npt correlators
#############################################
with open(options.nptCorrsFile) as tpt:
    nptCorrs=tpt.read().splitlines()
data3pt=[]
for F in nptCorrs:

    # If we are dealing with more than one file, parse and combine before proceeding
    if len(F.split(' ')) == 1:
        data3pt.append(np.loadtxt(F, delimiter=' ', skiprows=1))
        continue
    elif len(F.split(' ')) == 2:
        dumCorr1=np.loadtxt(F.split(' ')[0], delimiter=' ', skiprows=1)
        dumCorr2=np.loadtxt(F.split(' ')[1], delimiter=' ', skiprows=1)
        print dumCorr1
        mergedCorr=mergeCorrs(dumCorr1,dumCorr2,gamma)
        data3pt.append(mergedCorr)
    else:
        print "Don't know how to handle %i correlators"%len(F.split(' '))
        sys.exit()




src2pt=correlator(cfgs,len(corrSrc2pt)/cfgs,corrSrc2pt,1,1.0)
snk2pt=correlator(cfgs,len(corrSnk2pt)/cfgs,corrSnk2pt,1,1.0)
src2pt.makeJks()
snk2pt.makeJks()
src2pt.avgJks()
snk2pt.avgJks()



R=[]
for n, c in enumerate(data3pt):
    dum3pt=correlator(cfgs,len(c)/cfgs,c,complexity,redstarFactor)
    dum3pt.makeJks()
    dum3pt.avgJks()
    # Make the ratio for this data3pt entry
    r=ratio(dum3pt,src2pt,snk2pt)
    r.avgratio()
    r.avgratioErr()
    # Append
    R.append(r)




# fig=plt.figure()
ax=plt.gca()
# title_str = "Nucleon effective "+charge+r" Charge - Isoclover $32^3\times64$    $m_\pi=%d$ MeV    $a=%.3f$ fm"%(356,0.098)+"\n"+title_str
# plt.title(title_str)

if complexity == 1:
    ax.set_ylabel(r'Re $M_%s\left(p_f,p_i,z_3\right)$'%gammaDict[gamma],fontsize=20)
if complexity == 2:
    ax.set_ylabel(r'Im $M_%s\left(p_f,p_i,z_3\right)$'%gammaDict[gamma],fontsize=20)
ax.set_xlabel(r'$\left(\tau-T/2\right)a^{-1}$',fontsize=20)

def makeTimeArrays(N):
    time=np.zeros(N)
    for i in range(0,N):
        time[i]=i-N/2

    return time


for n, r in enumerate(R):
    ts=makeTimeArrays(len(r.Avg))
    # ts=[i-len(r.Avg)/2 for i in range(0,len(r.Avg))]
    # print ts
    r.plotRatio(ts,ax)

# Set a suitable vertical range
ymin, ymax = ax.get_ylim()
Ymin = ymin - 0.8*(ymax-ymin)
Ymax = ymax + 0.8*(ymax-ymin)
ax.set_ylim([Ymin,Ymax])

# Add the momenta and displacements to plot
# ax.text(-7,1.05*ax.get_ylim()[0],r'$\vec{p}_i=\left(%s,%s,%s\right)$'\
#         %(srcMom.split('.')[0],srcMom.split('.')[1],srcMom.split('.')[2]),fontsize=18)
# ax.text(-2,1.05*ax.get_ylim()[0],r'$\vec{p}_f=\left(%s,%s,%s\right)$'\
#         %(snkMom.split('.')[0],snkMom.split('.')[1],snkMom.split('.')[2]),fontsize=18)
# ax.text(3,1.05*ax.get_ylim()[0],r'$\vec{z}=\left(%s,%s,%s\right)$'\
#         %(vecZ.split('.')[0],vecZ.split('.')[1],vecZ.split('.')[2]),fontsize=18)

xext=abs(ax.get_xlim()[1]-ax.get_xlim()[0])
yext=abs(ax.get_ylim()[1]-ax.get_ylim()[0])

ax.text(ax.get_xlim()[0]+0.125*xext,ax.get_ylim()[0]+0.1*yext,\
        r'$\vec{p}_i=\left(%s,%s,%s\right)\quad\quad\vec{p}_f=\left(%s,%s,%s\right)\quad\quad\vec{z}=\left(%s,%s,%s\right)$'\
        %(srcMom.split('.')[0],srcMom.split('.')[1],srcMom.split('.')[2],\
          snkMom.split('.')[0],snkMom.split('.')[1],snkMom.split('.')[2],\
          vecZ.split('.')[0],vecZ.split('.')[1],vecZ.split('.')[2]),fontsize=18)
ax.text(ax.get_xlim()[0]+0.125*xext,ax.get_ylim()[0]+0.05*yext,\
        r'$r_i=\mu_{%s}\quad\quad\quad\quad r_f=\mu_{%s}$'%(options.srcRow,options.snkRow),fontsize=18)


ax.legend(fontsize=16,ncol=3,loc=9,markerscale=1,fancybox=True)
# ax.legend(fontsize=16,ncol=3,loc=9,markerscale=1,fancybox=True)
# ax.legend.get_frame().set_alpha(FrameAlpha)

# Name the output
if options.output == "":
    output_name="ppdf_plateau_pf%s%s%s_pi%s%s%s_r%s%s_g%i_z%s%s%s_c%i"\
        %(snkMom.split('.')[0],snkMom.split('.')[1],snkMom.split('.')[2],\
          srcMom.split('.')[0],srcMom.split('.')[1],srcMom.split('.')[2],\
          options.snkRow,options.srcRow,gamma,\
          vecZ.split('.')[0],vecZ.split('.')[1],vecZ.split('.')[2],complexity)
else:
    output_name=options.output
plt.savefig(output_name,transparent=truthTransparent,bbox_inches='tight',pad_inches=0)#,dpi=500)
if options.showFig == 1:
    plt.show()





# # Set the momenta throughout


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


# # Set the generated effective charge plot with fits
# if options.output == "":
#     output_name="p%i%i%i_g%i_rows%s%s"%(momenta["Src"]["val"][0],momenta["Src"]["val"][1],momenta["Src"]["val"][2],gamma,options.snkRow,options.srcRow)
# else:
#     output_name=options.output

# # Just some temp stuff
# dataT6=data3pt[0]
# dataT8=data3pt[1]
# dataT10=data3pt[2]
# if options.cutNoisy == 0:
#     dataT12=data3pt[3]
#     dataT14=data3pt[4]
#     dataT16=data3pt[5]
# elif options.cutNoisy == 2:
#     dataT12=data3pt[3]
#     dataT14=data3pt[4]



# # SET RENORMALIZATION CONSTANTS
# renormalizations = { 0 : 0.800, 1 : 0.834, 2 : 0.834, 4 : 0.834, 8 : 0.834, 7 : 0.887, 11 : 0.887,
#                      13 : 0.887, 14 : 0.887, 3 : 0.925, 12 : 0.925, 15: 1.000 }

# print "Gamma  = %d"%gamma
# print "Renorm = %.7f"%renormalizations[gamma]


# def makeTimeArrays(N):
#     time=np.zeros(N)
#     for i in range(0,N):
#         time[i]=i-N/2

#     return time



# # MAKE THE JACKKNIFE SAMPLES OF 2PT & 4PT DATA
# print "Making jackknife samples"
# data2Jks=jackknifeSamples(data2pt,1,21)
# t6Jks=jackknifeSamples(dataT6,complexity,6)
# t8Jks=jackknifeSamples(dataT8,complexity,8)
# t10Jks=jackknifeSamples(dataT10,complexity,10)
# if options.cutNoisy == 0:
#     t12Jks=jackknifeSamples(dataT12,complexity,12)
#     t14Jks=jackknifeSamples(dataT14,complexity,14)
#     t16Jks=jackknifeSamples(dataT16,complexity,16)
# elif options.cutNoisy == 2:
#     t12Jks=jackknifeSamples(dataT12,complexity,12)
#     t14Jks=jackknifeSamples(dataT14,complexity,14)


# # MAKE AVERAGES PER JACKKNIFE SAMPLE OF 2PT & 4PT DATA
# print "Calculating averages per jackknife sample"
# data2JkAvg=avgJackknifeSamples(data2Jks,21)
# t6JksAvg=avgJackknifeSamples(t6Jks,6)
# t8JksAvg=avgJackknifeSamples(t8Jks,8)
# t10JksAvg=avgJackknifeSamples(t10Jks,10)
# if options.cutNoisy == 0:
#     t12JksAvg=avgJackknifeSamples(t12Jks,12)
#     t14JksAvg=avgJackknifeSamples(t14Jks,14)
#     t16JksAvg=avgJackknifeSamples(t16Jks,16)
# elif options.cutNoisy == 2:
#     t12JksAvg=avgJackknifeSamples(t12Jks,12)
#     t14JksAvg=avgJackknifeSamples(t14Jks,14)


# # GET AN AVERAGE DETERMINATION PER TIME SLICE
# print "Calculating average estimates per time slice"
# t6RatioAvg=avgJackknifeRatio(t6JksAvg,data2JkAvg,6)
# print t6RatioAvg
# print t6JksAvg
# print data2JkAvg
# t8RatioAvg=avgJackknifeRatio(t8JksAvg,data2JkAvg,8)
# t10RatioAvg=avgJackknifeRatio(t10JksAvg,data2JkAvg,10)
# if options.cutNoisy == 0:
#     t12RatioAvg=avgJackknifeRatio(t12JksAvg,data2JkAvg,12)
#     t14RatioAvg=avgJackknifeRatio(t14JksAvg,data2JkAvg,14)
#     t16RatioAvg=avgJackknifeRatio(t16JksAvg,data2JkAvg,16)
# elif options.cutNoisy == 2:
#     t12RatioAvg=avgJackknifeRatio(t12JksAvg,data2JkAvg,12)
#     t14RatioAvg=avgJackknifeRatio(t14JksAvg,data2JkAvg,14)

# # GET AN ERROR DETERMINATION PER TIME SLICE
# print "Calculating error estimates per time slice"
# t6RatioErr=avgJackknifeRatioErr(t6JksAvg,data2JkAvg,t6RatioAvg,6)
# t8RatioErr=avgJackknifeRatioErr(t8JksAvg,data2JkAvg,t8RatioAvg,8)
# t10RatioErr=avgJackknifeRatioErr(t10JksAvg,data2JkAvg,t10RatioAvg,10)
# if options.cutNoisy == 0:
#     t12RatioErr=avgJackknifeRatioErr(t12JksAvg,data2JkAvg,t12RatioAvg,12)
#     t14RatioErr=avgJackknifeRatioErr(t14JksAvg,data2JkAvg,t14RatioAvg,14)
#     t16RatioErr=avgJackknifeRatioErr(t16JksAvg,data2JkAvg,t16RatioAvg,16)
# elif options.cutNoisy == 2:
#     t12RatioErr=avgJackknifeRatioErr(t12JksAvg,data2JkAvg,t12RatioAvg,12)
#     t14RatioErr=avgJackknifeRatioErr(t14JksAvg,data2JkAvg,t14RatioAvg,14)

# # MAKE TIME ARRAYS FOR PLOTTING
# D6time=makeTimeArrays(6)
# D8time=makeTimeArrays(8)
# D10time=makeTimeArrays(10)
# if options.cutNoisy == 0:
#     D12time=makeTimeArrays(12)
#     D14time=makeTimeArrays(14)
#     D16time=makeTimeArrays(16)
# elif options.cutNoisy == 2:
#     D12time=makeTimeArrays(12)
#     D14time=makeTimeArrays(14)

# bareTime=np.linspace(-9,9,500)
# bareCharge=np.zeros(500)
# bareChargeError=np.zeros(500)


# ######################################################################
# # DETERMINE PREFACTORS/KINEMATIC FACTORS THAT NEED TO BE INCLUDED
# ######################################################################
# redstarFactor=np.sqrt(2) # to account for redstar isovector/scalar currents
# kinefactor=1.0
# if srcMom.split('.')[2] < 4:
#     gndmass=np.sqrt((E0**2) - ((2*np.pi/32.0)*(momenta["Src"]["val"][0]+momenta["Src"]["val"][1]+momenta["Src"]["val"][2]))**2)
# else:
#     gndmass=0.4757

# momz=momenta["Src"]["val"][2]

# print "GND MASS = %.7f"%gndmass
# numeratorE=4*E0*(E0+gndmass)# E0 IS FITTED NUCLEON ENERGY, NOT JUST GND MASS
# if gamma==0:
#     kinefactor=numeratorE/(4.0*gndmass*(gndmass+E0))
# if gamma==3:
#     kinefactor=numeratorE/(4.0*(gndmass*E0+gndmass**2)) # not so certain on the factor of 4 in denom here (10/12/19)
# if gamma==11:
#     kinefactor=1.0*numeratorE/(4.0*(gndmass*E0+gndmass**2+(((2*np.pi)/(32.0))*momz)**2))
# if gamma==7:
#     kinefactor=-1.0*numeratorE/(4.0*(gndmass+E0)*(((2*np.pi)/(32.0))*momz))
#     # kinefactor=-1.0/(4*(((2*np.pi)/(1.0*32))*momz)*E0) # produces charges ~0.68
# if gamma==8:
#     kinefactor=numeratorE/(4.0*E0*(E0+gndmass))
# if gamma==4:
#     kinefactor=numeratorE/(4.0*(E0+gndmass)*(((2*np.pi)/(32.0))*momz))

# print "Kinematic factor applied to ratio = %.7f"%kinefactor
# ########################################################################


# # Set the asymptotic bare/renormalized charge and error
# bareCharge[:]=(A/a)*kinefactor*redstarFactor
# # if bareCharge[0] < 0:
# #     print "WARNING: DOING SOME TOM FOOLERY HERE!!!"
# #     bareCharge[:]*=-1.0 #BAD!
# #     redstarFactor*=-1.0 #BAD!

# # bareChargeError[:]=np.sqrt(np.power(dA/a,2)+np.power((A/(a**2)),2)*np.power(da,2)-((2*A)/(np.power(a,3)))*covAa)*kinefactor*redstarFactor

# bareChargeError[:]=np.abs(np.sqrt((1.0/a)*COV[('A', 'A')]*(1.0/a)+(1.0/a)*COV[('A', 'a')]*(-A/(a**2))+(-A/(a**2))*COV[('a', 'A')]*(1.0/a)+(-A/(a**2))*COV[('a', 'a')]*(-A/(a**2)))*kinefactor*redstarFactor)


# renormCharge=renormalizations[gamma]*bareCharge
# renormChargeError=renormalizations[gamma]*bareChargeError

# print "Bare charge = %.7f"%bareCharge[0]
# print "Bare charge error = %.7f"%bareChargeError[0]
# print "Renorm charge = %.7f"%renormCharge[0]
# print "Renorm charge error = %.7f"%renormChargeError[0]

# time6=np.linspace(-0.5,5.5,500)
# time8=np.linspace(-0.5,7.5,500)
# time10=np.linspace(-0.5,9.5,500)
# time12=np.linspace(-0.5,11.5,500)
# time14=np.linspace(-0.5,13.5,500)
# time16=np.linspace(-0.5,15.5,500)
# tsep6Fit=np.zeros(500)
# tsep6FitError=np.zeros(500)
# tsep8Fit=np.zeros(500)
# tsep8FitError=np.zeros(500)
# tsep10Fit=np.zeros(500)
# tsep10FitError=np.zeros(500)
# tsep12Fit=np.zeros(500)
# tsep12FitError=np.zeros(500)
# tsep14Fit=np.zeros(500)
# tsep14FitError=np.zeros(500)
# tsep16Fit=np.zeros(500)
# tsep16FitError=np.zeros(500)
# # Define the 3pt/2pt fitted ratios
# for i in range(0,500):
#     tsep6Fit[i]=fit3pt2ptRatio(params,time6[i],6)*renormalizations[gamma]*kinefactor*redstarFactor
#     tsep8Fit[i]=fit3pt2ptRatio(params,time8[i],8)*renormalizations[gamma]*kinefactor*redstarFactor
#     tsep10Fit[i]=fit3pt2ptRatio(params,time10[i],10)*renormalizations[gamma]*kinefactor*redstarFactor
#     if options.cutNoisy == 0:
#         tsep12Fit[i]=fit3pt2ptRatio(params,time12[i],12)*renormalizations[gamma]*kinefactor*redstarFactor
#         tsep14Fit[i]=fit3pt2ptRatio(params,time14[i],14)*renormalizations[gamma]*kinefactor*redstarFactor
#         tsep16Fit[i]=fit3pt2ptRatio(params,time16[i],16)*renormalizations[gamma]*kinefactor*redstarFactor
#     elif options.cutNoisy == 2:
#         tsep12Fit[i]=fit3pt2ptRatio(params,time12[i],12)*renormalizations[gamma]*kinefactor*redstarFactor
#         tsep14Fit[i]=fit3pt2ptRatio(params,time14[i],14)*renormalizations[gamma]*kinefactor*redstarFactor


# # Define the 3pt/2pt fitted ratio errors
# for i in range(0,500):
#     tsep6FitError[i]=Deriv(time6[i],6,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
#     tsep8FitError[i]=Deriv(time8[i],8,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
#     tsep10FitError[i]=Deriv(time10[i],10,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
#     if options.cutNoisy == 0:
#         tsep12FitError[i]=Deriv(time12[i],12,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
#         tsep14FitError[i]=Deriv(time14[i],14,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
#         tsep16FitError[i]=Deriv(time16[i],16,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
#     elif options.cutNoisy == 2:
#         tsep12FitError[i]=Deriv(time12[i],12,covariance)*renormalizations[gamma]*kinefactor*redstarFactor
#         tsep14FitError[i]=Deriv(time14[i],14,covariance)*renormalizations[gamma]*kinefactor*redstarFactor


# # Determine the effective charges (3pt/2pt fit)
# effChT6=makeEffCharge(t6JksAvg,6)*kinefactor*renormalizations[gamma]*redstarFactor
# effChT8=makeEffCharge(t8JksAvg,8)*kinefactor*renormalizations[gamma]*redstarFactor
# effChT10=makeEffCharge(t10JksAvg,10)*kinefactor*renormalizations[gamma]*redstarFactor
# if options.cutNoisy == 0:
#     effChT12=makeEffCharge(t12JksAvg,12)*kinefactor*renormalizations[gamma]*redstarFactor
#     effChT14=makeEffCharge(t14JksAvg,14)*kinefactor*renormalizations[gamma]*redstarFactor
#     effChT16=makeEffCharge(t16JksAvg,16)*kinefactor*renormalizations[gamma]*redstarFactor
# elif options.cutNoisy == 2:
#     effChT12=makeEffCharge(t12JksAvg,12)*kinefactor*renormalizations[gamma]*redstarFactor
#     effChT14=makeEffCharge(t14JksAvg,14)*kinefactor*renormalizations[gamma]*redstarFactor

# # Determine the error in effective charges (3pt/2pt fit)
# # effCh* right above was determined prior to rescaling by all factors
# # Thus in determining error on the effective charges, the JkAvgs must be rescaled
# effChErrT6=makeEffChargeError(t6JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT6,6)
# effChErrT8=makeEffChargeError(t8JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT8,8)
# effChErrT10=makeEffChargeError(t10JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT10,10)
# if options.cutNoisy == 0:
#     effChErrT12=makeEffChargeError(t12JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT12,12)
#     effChErrT14=makeEffChargeError(t14JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT14,14)
#     effChErrT16=makeEffChargeError(t16JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT16,16)
# elif options.cutNoisy == 2:
#     effChErrT12=makeEffChargeError(t12JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT12,12)
#     effChErrT14=makeEffChargeError(t14JksAvg*kinefactor*renormalizations[gamma]*redstarFactor,effChT14,14)


# # textstr=''
# # title_str=operator2+" -- "+operator2

# #########################################################################
# # POTENTIALLY DETERMINE/PLOT THE EFFECTIVE 2PT EFFECTIVE ENERGY AND FIT
# #########################################################################
# if options.show2ptFitToData == 1:
#     time2pt=np.linspace(0,20,3000)
#     fitRange2pt=np.linspace(float(options.fitRange2pt.split('.')[0]),\
#                             float(options.fitRange2pt.split('.')[2]),3000)

#     # Determine the central value of fit to effective energy
#     # Also write out the central values
#     W=open("effE_m0p2390_fit2pt_p"+str(options.pi.split('.')[0])+str(options.pi.split('.')[1])+str(options.pi.split('.')[2])+".dat",'w')
#     fit2ptEffE=np.linspace(0,20,3000)
#     for i in range(0,len(time2pt)-1):
#         fit2ptEffE[i]=(1.0/(time2pt[1]-time2pt[0]))*np.log(twoState(params['a'],params['b'],params['E0'],params['E1'],time2pt[i])/twoState(params['a'],params['b'],params['E0'],params['E1'],time2pt[i+1]))
        
#     fit2ptEffEError=np.linspace(0,20,3000)
#     for i in range(0,len(time2pt)-1):
#         fit2ptEffEError[i]=Deriv2pt(time2pt[i],time2pt[1]-time2pt[0],covariance)
#         W.write(str(time2pt[i])+" "+str(fit2ptEffE[i])+" "+str(fit2ptEffEError[i])+"\n")
#     W.close()


#     plt.rc('text', usetex=True)
#     plt.xlabel(r't_{\rm{sep}}',fontsize=14)
#     plt.ylabel(r'E_{\rm{eff}}', fontsize=14)
#     plt.plot(time2pt,fit2ptEffE[:],'b',label=r'$2-\text{state Fit}')#,markersize=4,markerfacecolor="Blue",markeredgecolor='b')
#     plt.fill_between(time2pt,fit2ptEffE+fit2ptEffEError,fit2ptEffE-fit2ptEffEError,color='#B2E5FF',alpha=0.8)
#     plt.ylim([0,1.1])
#     plt.show()
#     sys.exit()

# ########################################################
# # PLOT THE RENORMALIZED EFFECTIVE CHARGES WITH FITS
# ########################################################
# plt.rc('text', usetex=True)
# # title_str = "Nucleon effective "+charge+r" Charge - Isoclover $32^3\times64$    $m_\pi=%d$ MeV    $a=%.3f$ fm"%(356,0.098)+"\n"+title_str
# # plt.title(title_str)
# plt.xlabel(r'$\tau-T/2$',fontsize=20)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# # plt.figtext(0.5, 0.2, textstr, fontsize=10, verticalalignment='baseline', horizontalalignment='center')

# # For plot clarity
# offset=0.1
# # plt.rc('text', usetex=True)




# # Method to manage plotting of either effective charges or 3pt/2pt ratio of raw lattice data
# def plotData(pt,effCh,effChE):
#     buffMin=2
#     params={ 6: {'buffMax': 5, 'Mark': 'o', 'Color': 'c', 'Label': r'$T/a=6$',\
#                  'Shift': -3*offset },
#              8: {'buffMax': 7, 'Mark': '^', 'Color': 'b', 'Label': r'$T/a=8$',\
#                  'Shift': -2*offset},
#              10: {'buffMax': 9, 'Mark': '*', 'Color': 'm', 'Label': r'$T/a=10$',\
#                   'Shift': -offset},
#              12: {'buffMax': 11, 'Mark': 's', 'Color': 'g', 'Label': r'$T/a=12$', 'Shift': 0 },
#              14: {'buffMax': 13, 'Mark': 'h', 'Color': 'r', 'Label': r'$T/a=14$',\
#                   'Shift': offset },
#              16: {'buffMax': 15, 'Mark': 'v', 'Color': 'k', 'Label': r'$T/a=16$',\
#                   'Shift': 2*offset } }

#     nn=len(pt)
#     # Do things based on the times series
#     plt.plot(pt[:buffMin]+params[nn]['Shift'],effCh[:buffMin],params[nn]['Mark'],\
#              markersize=3,mfc='grey',mec='grey')
#     plt.errorbar(pt[:buffMin]+params[nn]['Shift'],effCh[:buffMin],yerr=effChE[:buffMin],\
#                  fmt=params[nn]['Mark'],label='',ecolor='grey',mfc='grey',mec='grey')

#     plt.plot(pt[buffMin:params[nn]['buffMax']]+params[nn]['Shift'],\
#              effCh[buffMin:params[nn]['buffMax']],\
#              params[nn]['Mark'],label=params[nn]['Label'],markersize=3,\
#              mfc=params[nn]['Color'],mec=params[nn]['Color'])
#     plt.errorbar(pt[buffMin:params[nn]['buffMax']]+params[nn]['Shift'],\
#                  effCh[buffMin:params[nn]['buffMax']],\
#                  yerr=effChE[buffMin:params[nn]['buffMax']],fmt=params[nn]['Mark'],label='',\
#                  ecolor=params[nn]['Color'],mfc=params[nn]['Color'],mec=params[nn]['Color'])

#     plt.plot(pt[params[nn]['buffMax']:]+params[nn]['Shift'],effCh[params[nn]['buffMax']:],\
#              params[nn]['Mark'],markersize=3,mfc='grey',mec='grey')
#     plt.errorbar(pt[params[nn]['buffMax']:]+params[nn]['Shift'],effCh[params[nn]['buffMax']:],\
#                  yerr=effChE[params[nn]['buffMax']:],fmt=params[nn]['Mark'],label='',\
#                  ecolor='grey',mfc='grey',mec='grey')
    

# ###########################################
# # PLOT EITHER RAW 3PT/2PT LATTICE RATIO
# # OR 3PT/(2PT FIT) - I.E. EFFECTIVE CHARGE
# ###########################################
# dataTimeList=[D6time,D8time,D10time]
# dataList=dataErrList=""
# if options.rawRatio:
#     dataList=[t6RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
#               t8RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
#               t10RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma]]
#     dataErrList=[t6RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
#                  t8RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
#                  t10RatioErr*redstarFactor*kinefactor*renormalizations[gamma]]
#     if options.cutNoisy == 0:
#         dataList.extend([t12RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
#                          t14RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
#                          t16RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma]])
#         dataErrList.extend([t12RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
#                             t14RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
#                             t16RatioErr*redstarFactor*kinefactor*renormalizations[gamma]])
#         dataTimeList.extend([D12time,D14time,D16time])
#     elif options.cutNoisy == 2:
#         dataList.extend([t12RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma],\
#                          t14RatioAvg[:]*redstarFactor*kinefactor*renormalizations[gamma]])
#         dataErrList.extend([t12RatioErr*redstarFactor*kinefactor*renormalizations[gamma],\
#                             t14RatioErr*redstarFactor*kinefactor*renormalizations[gamma]])
#         dataTimeList.extend([D12time,D14time])
# else:
#     dataList=[effChT6,effChT8,effChT10]
#     dataErrList=[effChErrT6,effChErrT8,effChErrT10]
#     if options.cutNoisy == 0:
#         dataList.extend([effChT12,effChT14,effChT16])
#         dataErrList.extend([effChErrT12,effChErrT14,effChErrT16])
#         dataTimeList.extend([D12time,D14time,D16time])
#     elif options.cutNoisy == 2:
#         dataList.extend([effChT12,effChT14])
#         dataErrList.extend([effChErrT12,effChErrT14])
#         dataTimeList.extend([D12time,D14time]) 


# ################################
# # PLOT THE SELECTED DATA
# ################################
# for n, k in enumerate(dataTimeList):
#     plotData(k,dataList[n],dataErrList[n])



# YMIN=YMAX=0.0
# if gamma==0:
#     plt.ylabel(r'$g_S^{u-d}$',fontsize=20)
#     YMIN=0.45
#     YMAX=1.3
#     # # For unphased pz=4 data
#     # YMIN=-3
#     # YMAX=4
#     # # For phased pz=4 data
#     YMIN=-1
#     YMAX=2.25
# if gamma==1:
#     plt.ylabel(r'$g_{V_1}^{u-d}$',fontsize=20)
# if gamma==2:
#     plt.ylabel(r'$g_{V_2}^{u-d}$',fontsize=20)
# if gamma==4:
#     plt.ylabel(r'$g_{V_3}^{u-d}$',fontsize=20)
#     YMIN=0.7
#     YMAX=1.10
#     # # For the unphased data at pz=4
#     # YMIN=-0.2
#     # YMAX=2.25
#     # # For the phased data at pz=4
#     # YMIN=0.94
#     # YMAX=1.15
#     # # For the phased data at pz=4 (3pt/2pt raw ratio)
#     # YMIN=0.25
#     # YMAX=1.5
# if gamma==8:
#     plt.ylabel(r'$g_{V_4}^{u-d}$',fontsize=20)
#     YMIN=0.94
#     YMAX=1.10
#     # # For the unphased data at pz=4
#     # YMIN=0.5
#     # YMAX=1.3
# if gamma==7:
#     plt.ylabel(r'$g_{A_4}^{u-d}$',fontsize=20)
#     YMIN=0.75
#     YMAX=1.35
#     # # For the unphased pz=4 data
#     # YMIN=0
#     # YMAX=2
#     # # For the unphased pz=4 data (3pt/2pt raw ratio)
#     # YMIN=0.25
#     # YMAX=1.75
#     # # For the phased pz=4 data
#     # YMIN=1.1
#     # YMAX=1.5
# if gamma==11:
#     plt.ylabel(r'$g_{A_3}^{u-d}$',fontsize=20)
#     YMIN=1.0
#     YMAX=1.35
#     # # For unphased pz=4 data
#     # YMIN=0
#     # YMAX=2
#     # # For phased pz=4 data
#     # YMIN=1.1
#     # YMAX=1.55
# if gamma==3:
#     plt.ylabel(r'$g_{T_{12}}^{u-d}$',fontsize=20)
#     YMIN=0.95
#     YMAX=1.3
#     # # For unphased pz=4 data
#     YMIN=0.25
#     YMAX=2.25
# if gamma==15:
#     plt.ylabel(r'$g_P^{u-d}$',fontsize=20)
#     YMIN=-0.4
#     YMAX=0.4
# if gamma==15 and options.srcRow!=options.snkRow:
#     plt.ylabel(r'$g_P^{u-d}$',fontsize=20)
#     YMIN=-1.2
#     YMAX=1.2

# # # Let's modify the vertical range permanently for the poor phased results
# # YMIN=0.2
# # YMAX=1.5


# # Extracted charge was obtained from fitting 3pt data that was already rescaled (sqrt(2) for Redstar generated files)
# # plt.plot(bareTime,bareCharge,'black')
# # plt.fill_between(bareTime,bareCharge+bareChargeError,bareCharge-bareChargeError,color='#d8dcd6')
# plt.plot(bareTime,renormCharge,'black')
# plt.fill_between(bareTime,renormCharge+renormChargeError,renormCharge-renormChargeError,color='#d8dcd6')
# plt.xlim([-9,9])
# plt.xlim([-7.5,7.5])
# plt.ylim([YMIN,YMAX])



# # Include the fits to the 3pt/2pt data
# plt.plot(time6[:]-3-3*offset, tsep6Fit[:],'c')
# plt.fill_between(time6[:]-3-3*offset,tsep6Fit[:]+tsep6FitError[:],
#                  tsep6Fit[:]-tsep6FitError[:],color='#13EAC9',alpha=0.6,
#                  where=np.logical_and(time6[:]>=2,time6[:]<=4))
# plt.fill_between(time6[:]-3-3*offset,tsep6Fit[:]+tsep6FitError[:],
#                  tsep6Fit[:]-tsep6FitError[:],color='grey',alpha=0.1,
#                  where=np.logical_or(time6[:]<2,time6[:]>4))

# plt.plot(time8[:]-4-2*offset, tsep8Fit[:],'b')
# plt.fill_between(time8[:]-4-2*offset,tsep8Fit[:]+tsep8FitError[:],
#                  tsep8Fit[:]-tsep8FitError[:],color='b',alpha=0.2,# '#B2E5FF'
#                  where=np.logical_and(time8[:]>=2,time8[:]<=6))
# plt.fill_between(time8[:]-4-2*offset,tsep8Fit[:]+tsep8FitError[:],
#                  tsep8Fit[:]-tsep8FitError[:],color='grey',alpha=0.1,
#                  where=np.logical_or(time8[:]<2,time8[:]>6))

# plt.plot(time10[:]-5-offset, tsep10Fit,'m')
# plt.fill_between(time10[:]-5-offset,tsep10Fit[:]+tsep10FitError[:],
#                  tsep10Fit[:]-tsep10FitError[:],color='#D767AD',alpha=0.4,
#                  where=np.logical_and(time10[:]>=2,time10[:]<=8))
# plt.fill_between(time10[:]-5-offset,tsep10Fit[:]+tsep10FitError[:],
#                  tsep10Fit[:]-tsep10FitError[:],color='grey',alpha=0.1,
#                  where=np.logical_or(time10[:]<2,time10[:]>8))

# if options.cutNoisy == 0:
#     plt.plot(time12[:]-6, tsep12Fit,'g')
#     plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
#                      tsep12Fit[:]-tsep12FitError[:],color='g',alpha=0.2,#BCFFB8
#                      where=np.logical_and(time12[:]>=2,time12[:]<=10))
#     plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
#                      tsep12Fit[:]-tsep12FitError[:],color='grey',alpha=0.1,
#                      where=np.logical_or(time12[:]<2,time12[:]>10))

#     plt.plot(time14[:]-7+offset, tsep14Fit,'r')
#     plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
#                      tsep14Fit[:]-tsep14FitError[:],color='r',alpha=0.2,#FDD2D2
#                      where=np.logical_and(time14[:]>=2,time14[:]<=12))
#     plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
#                      tsep14Fit[:]-tsep14FitError[:],color='grey',alpha=0.1,
#                      where=np.logical_or(time14[:]<2,time14[:]>12))

#     plt.plot(time16[:]-8+2*offset, tsep16Fit,'k')
#     plt.fill_between(time16[:]-8+2*offset,tsep16Fit[:]+tsep16FitError[:],
#                      tsep16Fit[:]-tsep16FitError[:],color='#516572',alpha=0.6,
#                      where=np.logical_and(time16[:]>=2,time16[:]<=14))
#     plt.fill_between(time16[:]-8+2*offset,tsep16Fit[:]+tsep16FitError[:],
#                      tsep16Fit[:]-tsep16FitError[:],color='grey',alpha=0.1,
#                      where=np.logical_or(time16[:]<2,time16[:]>14))

# elif options.cutNoisy == 2:
#     plt.plot(time12[:]-6, tsep12Fit,'g')
#     plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
#                      tsep12Fit[:]-tsep12FitError[:],color='g',alpha=0.2,#BCFFB8
#                      where=np.logical_and(time12[:]>=2,time12[:]<=10))
#     plt.fill_between(time12[:]-6,tsep12Fit[:]+tsep12FitError[:],
#                      tsep12Fit[:]-tsep12FitError[:],color='grey',alpha=0.1,
#                      where=np.logical_or(time12[:]<2,time12[:]>10))

#     plt.plot(time14[:]-7+offset, tsep14Fit,'r')
#     plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
#                      tsep14Fit[:]-tsep14FitError[:],color='r',alpha=0.2,#FDD2D2
#                      where=np.logical_and(time14[:]>=2,time14[:]<=12))
#     plt.fill_between(time14[:]-7+offset,tsep14Fit[:]+tsep14FitError[:],
#                      tsep14Fit[:]-tsep14FitError[:],color='grey',alpha=0.1,
#                      where=np.logical_or(time14[:]<2,time14[:]>12))


# # # Lastly, add a text box with the reduced chi^2 of simultaneous fit
# # txtbxstr='$\\chi_r^2$ = %.4f'%chi2
# # plt.text(0.00,(YMAX-YMIN)*0.05+YMIN,txtbxstr,horizontalalignment='center',fontsize=14)

# ax = plt.gca()
# legend = ax.legend(fontsize=16,ncol=3,loc=9,markerscale=2,fancybox=True)
# # legend.get_frame().set_facecolor('b')
# legend.get_frame().set_alpha(0.5)
# # ax.legend(framealpha=0.5)
# # ax.legend(fontsize=14,ncol=3,loc=9,markerscale=2,fancybox=True,framealpha=0.5)

# # plt.legend(fontsize=14,ncol=3,loc=9,markerscale=2)#,fancybox=True,framealpha=0.5)
# plt.savefig(output_name,dpi=500)
# plt.show()

