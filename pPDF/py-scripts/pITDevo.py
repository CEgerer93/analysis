#!/dist/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special
import scipy.interpolate as interpol
# import scipy as scipy
import pylab # to save figures to file
import sys,optparse
from collections import OrderedDict

usage = "Usage: %prog [options] "
parser = optparse.OptionParser(usage);

parser.add_option("-f", "--pITDdata", type="str", default="",
                  help='pITD data file (default = "")')
parser.add_option("-o", "--output", type="str", default="",
                  help='Name of output pITD (default = "")')
parser.add_option("-a", "--alpha", type="float", default=0.0,
                  help='Alpha parameter of pheno PDF (default = 0.0)')
parser.add_option("-b", "--beta", type="float", default=0.0,
                  help='Beta parameter of pheno PDF (default = 0.0)')

# Parse the input arguments
(options, args) = parser.parse_args()
pITDFile=options.pITDdata
output_name=options.output


# Color according the z value
mainColors=['blue','red','green','purple','orange','magenta',(0.1,1,0.1),'black','gray','gray']


fig_real = plt.figure()
fig_imag = plt.figure()
axreal = fig_real.gca()
aximag = fig_imag.gca()



realITD=[]
imagITD=[]
# Open the file and grab its contents
with open(pITDFile) as ptr:
    for cnt, line in enumerate(ptr):

        # Capture the line and remove spaces
        D=line.split(' ')
        # Store elements of D
        ioffeTime=float(D[0])
        comp=int(D[1])
        matelem=float(D[2])
        error=float(D[3])
        zsep=int(D[4].rstrip()) # need to remove newline character off end

        
        # Now plot
        thisColor=mainColors[abs(zsep)]
        
        # Add to ax1 if REAL
        if comp == 1:
            # ax1.errorbar(ioffeTime,matelem,yerr=error,fmt='o',mfc="None",mec=thisColor,color=thisColor,label="z=%s"%abs(zsep))
            axreal.errorbar(ioffeTime,matelem,yerr=error,fmt='o',mfc="None",mec=thisColor,color=thisColor,label="z=%s"%abs(zsep))

            realITD.append(ioffeTime,matelem,error,zsep)
        # Add to ax2 if imaginary
        if comp == 2:
            # ax2.errorbar(ioffeTime,matelem,yerr=error,fmt='o',mfc="None",mec=thisColor,color=thisColor,label="z=%s"%abs(zsep))
            aximag.errorbar(ioffeTime,matelem,yerr=error,fmt='o',mfc="None",mec=thisColor,color=thisColor,label="z=%s"%abs(zsep))


plt.rcParams["mathtext.fontset"]="stix"
aximag.set_xlim([-0.1,10])
aximag.set_ylim([-0.1,0.85]) # b_b0
# aximag.set_ylim([-0.1,3.85]) # pion_pion_2
axreal.set_xlim([-0.1,10])
axreal.set_ylim([-0.5,1.05])
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


itdNuSamplings=np.linspace(0,10,500)
realITD=np.zeros(500)
imagITD=np.zeros(500)
for n in range(0,len(realITD)):
    realITD[n]=ITD_real(itdNuSamplings[n],options.alpha,options.beta)
    imagITD[n]=ITD_imag(itdNuSamplings[n],options.alpha,options.beta)



###############################
# FIT THE REDUCED PITD FOR A GIVEN Z^2
# WITH A 6TH DEGREE POLYNOMIAL IN NU^2
###############################
def pITD_real_polyfit(nu,c2,c4,c6):
    return 1+c2*nu**2+c4*nu**4+c6*nu**6

def pITD_imag_polyfit(nu,c1,c3,c5):
    return c1*nu+c3*nu**3+c5*nu**5


print realITD






# Now that the FT of the determined PDF has been completed, add it to the RE/IM rITD plots
axreal.plot(itdNuSamplings,realITD,'b',label='PDF fit to avg - no jks')
aximag.plot(itdNuSamplings,imagITD,'b',label='PDF fit to avg - no jks')

axreal.legend(by_label.values(), by_label.keys(),fontsize=12,loc='lower left')
aximag.legend(by_label.values(), by_label.keys(),fontsize=12,loc='lower center')



fig_real.savefig(output_name+'_RE.png',dpi=400)
fig_imag.savefig(output_name+'_IM.png',dpi=400)
plt.show()
