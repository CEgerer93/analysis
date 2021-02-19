#!/dist/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
# import scipy as scipy
import pylab # to save figures to file
import sys,optparse
from collections import OrderedDict

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-f", "--pITDdata", type="str", default="",
                  help='pITD data file (default = "")')
parser.add_option("-m", "--pitdArbnormFile", type="str", default="",
                  help='A secondary pITD to plot of a different normalization (default = "")')
parser.add_option("-i", "--showPITDInset", type="int", default=0,
                  help='Make PDF inset visible (default = 0)')
parser.add_option("-z", "--zlim", type="int", default=0,
                  help='Plot pITD up to and including zlim (default = 0)')
parser.add_option("-E", "--pITDEvo", type="string", default="",
                  help='Evolved pITD data (default = "")')
parser.add_option("-o", "--output", type="str", default="",
                  help='Name of output pITD (default = "")')
parser.add_option("-a", "--alpha", type="float", default=0.0,
                  help='Alpha parameter of pheno PDF (default = 0.0)')
parser.add_option("-b", "--beta", type="float", default=0.0,
                  help='Beta parameter of pheno PDF (default = 0.0)')
parser.add_option("-p", "--getITD", type="int", default=1,
                  help='Bool to set whether LC ITD should be obtained (default = 1)')

# Parse the input arguments
(options, args) = parser.parse_args()
pITDFile=options.pITDdata
pitdArbnormFile=options.pitdArbnormFile
output_name=options.output

pITDEvoFile=options.pITDEvo

# Color according the z value
mainColors=['blue','red','green','purple','orange','magenta',(0.1,1,0.1),'black','gray',\
            'teal','indigo','brown','salmon','dodgerblue','darkgreen','crimson','cyan']

# # Initialize the figure
# fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))


# # Initialize a standalone figure
# fig_imag, aximag = plt.figure(2)
# fig_real, axreal = plt.figure(1)


fig_real = plt.figure()
fig_imag = plt.figure()

axreal = fig_real.gca()
aximag = fig_imag.gca()

# Optionally include an inset figure on REAL pITD plot
if options.showPITDInset:
    axreal_inset=fig_real.add_axes([0.6,0.6,0.28,0.28])
    axreal_inset.set_xlim([-0.1,0.9])
    axreal_inset.set_ylim([0.97,1.001])
    

# def evoKernel(val):
#     res = integrate.quad(lambda u: ((1+u**2)/(1.0*(1-u)))*( , 0, 1)
    


# # Apply perturbative evolution to data, evolving to a common scale (lattice spacing) provided by user
# def evolution(pITD,z0,asoverpi):
#     # pITD     = Ioffe-time pseudo-distribution (discrete, unevolved)
#     # z0       = evolution scale
#     # asoverpi = \alpha_s/pi
    
#     # z3       = spatial separation for this entry
    
#     pITD_evo[i] = pITD[i] - (2.0/3.0)*asoverpi*np.log((z0**2)/(z3**2))*evoKernel(pITD[i])

    


#     return pITD_evo




# Allow for different plotting symbols & horizontal shifts if more than 1 pITD is passed
symbs=['o','v']
shifts=[0.0,0.05]
# Open the file and grab its contents
for n, f in enumerate([pITDFile,pitdArbnormFile]):
    with open(f) as ptr:

        for cnt, line in enumerate(ptr):
            
            # Capture the line and remove spaces
            D=line.split(' ')
            # Store elements of D
            ioffeTime=float(D[0])
            comp=int(D[1])
            matelem=float(D[2])
            error=float(D[3])
            zsep=int(D[4].rstrip()) # need to remove newline character off end
            
            # Don't plot if zsep > zlim
            if zsep > options.zlim:
                break
        
            # Now plot
            thisColor=mainColors[abs(zsep)]
            
            # Add to ax1 if REAL color=thisColor,mfc=None,
            if comp == 1:
                axreal.errorbar(ioffeTime+shifts[n],matelem,yerr=error,fmt=symbs[n],\
                                color=thisColor,mec=thisColor,mfc=thisColor,label="z=%s"%abs(zsep))
                if options.showPITDInset:
                   axreal_inset.errorbar(ioffeTime+shifts[n],matelem,yerr=error,fmt=symbs[n],\
                                color=thisColor,mec=thisColor,mfc="None")
                    
            # Add to ax2 if imaginary
            if comp == 2:
                aximag.errorbar(ioffeTime+shifts[n],matelem,yerr=error,fmt=symbs[n],\
                                color=thisColor,mec=thisColor,mfc=thisColor,label="z=%s"%abs(zsep))
                    

# Open the evolved pITD data, if passed
if options.pITDEvo != "":
    with open(pITDEvoFile) as ptr:
        for cnt, line in enumerate(ptr):

            # Capture the line and remove spaces
            E=line.split(' ')
            # Store elements of the evolved data
            ioffeTime=float(E[0])
            comp=int(E[1])
            matelem=float(E[2])
            error=float(E[3])
            zsep=int(E[4].rstrip()) # need to remove newline character off end


            # Now plot
            thisColor=mainColors[abs(zsep)]
            # Add to ax1 if REAL
            if comp == 1:
                axreal.errorbar(ioffeTime,matelem,yerr=error,fmt='^',mfc=thisColor,mec=thisColor,\
                                color=thisColor,label="z=%s"%abs(zsep))

            # Add to ax2 if imaginary
            if comp == 2:
                aximag.errorbar(ioffeTime,matelem,yerr=error,fmt='^',mfc=thisColor,mec=thisColor,\
                                color=thisColor,label="z=%s"%abs(zsep))



plt.rcParams["mathtext.fontset"]="stix"
aximag.set_xlim([-0.1,10])
aximag.set_ylim([-0.1,0.85]) # b_b0
# aximag.set_ylim([-0.1,3.85]) # pion_pion_2
axreal.set_xlim([-0.1,10])
axreal.set_ylim([-0.5,1.05])

axreal.set_xlim([-0.1,6])
axreal.set_ylim([0.35,1.05])

axreal.set_xlabel(r'$\nu$',fontsize=16)
aximag.set_xlabel(r'$\nu$',fontsize=16)
axreal.set_ylabel(r'Re $\mathscr{M}\left(\nu,z^2\right)$',fontsize=16)
aximag.set_ylabel(r'Im $\mathscr{M}\left(\nu,z^2\right)$',fontsize=16)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))

axreal.legend(by_label.values(), by_label.keys(),fontsize=12,loc='lower left')
aximag.legend(by_label.values(), by_label.keys(),fontsize=12,loc='lower center')


# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace=None)
# plt.savefig(output_name+'.pdf')#,dpi=400)
# plt.savefig(output_name+'.png',dpi=200)



def phenoPDF2(x,a,b):
    return ((x**a)*(1-x)**b)/special.beta(a+1,b+1)


def ITD_real(nu,a,b):
    f = lambda x: np.cos(nu*x)*phenoPDF2(x,a,b)
    results = integrate.quad(f, 0, 1)
    return results[0]

def ITD_imag(nu,a,b):
    f = lambda x: np.sin(nu*x)*phenoPDF2(x,a,b)
    results = integrate.quad(f, 0, 1)
    return results[0]


if options.getITD:
    itdNuSamplings=np.linspace(0,10,500)
    realITD=np.zeros(500)
    imagITD=np.zeros(500)
    for n in range(0,len(realITD)):
        realITD[n]=ITD_real(itdNuSamplings[n],options.alpha,options.beta)
        imagITD[n]=ITD_imag(itdNuSamplings[n],options.alpha,options.beta)

        
    # Now that the FT of the determined PDF has been completed, add it to the RE/IM rITD plots
    axreal.plot(itdNuSamplings,realITD,'b',label='PDF fit to avg - no jks')
    aximag.plot(itdNuSamplings,imagITD,'b',label='PDF fit to avg - no jks')
    
axreal.legend(by_label.values(), by_label.keys(),fontsize=12,loc='lower left')
aximag.legend(by_label.values(), by_label.keys(),fontsize=12,loc='lower right')



fig_real.savefig(output_name+'_RE.png',dpi=400)
fig_imag.savefig(output_name+'_IM.png',dpi=400)
plt.show()
