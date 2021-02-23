#!/dist/anaconda/bin/python

# from /u/home/cegerer/src/sharedfuncs import *
from matplotlib import cm
from sharedfuncs import *

# Color according the z value
# mainColors=['blue','red','green','purple','orange','magenta',(0.1,1,0.1),'black','gray',\
#             'teal','indigo','brown','salmon','dodgerblue','darkgreen','crimson','cyan']
mainColors=['gray']
for dc in [ cm.hsv(d) for d in np.linspace(0,1,17) ]:
    mainColors.append(dc)

# Global function to determine ioffe-time
def ioffeTime(p,z,L):
    return ((2*np.pi)/L)*float(p.split('.')[2])*float(z)


# Global function to make zsep string
def makeZSep(z):
    s=''
    for i in range(0,z):
        s+='3'
    return s

# Global function to cast momentum to string
def mom2Str(x,y,z):
    return str(x)+str(y)+str(z)
    

# Global function to make mom_type
def mom_type(m):
    tmp=[]
    [tmp.append(i) for i in m]
    for n, m in enumerate(tmp):
        if int(m) < 0:
            m = str(int(m)*-1)

    # Now sort and return
    tmp.sort(reverse=True)
    return "%s%s%s"%(tmp[0],tmp[1],tmp[2])            


# Class to easily parse fit param files
class nptOp:
    def __init__(self,t,iso,row,mom,op,disp=''):
        self.t = t
        self.iso = iso
        self.row = 'r'+str(row)
        self.mom = str(mom)
        self.op = op+'__'+mom_type(self.mom)
        self.disp = disp
        self.name = ''

    def setName(self):
        if self.name != '': # reset the name if something non-trivial is read
            self.name = ''
        for s in [self.t, self.iso, self.row, self.mom, self.op, self.disp]:
            self.name += '%s%s'%(s,',')
        self.name = self.name.rstrip(",")


class matelem:
    def __init__(self,data,forceNullM):
        self.data = data  # Passed data
        self.params = {}  # Dict of all fitted parameters
        self.n = None     # Number of unique fitted parameters
        self.cfgs = 0
        self.jkparams = {}
        self.avgjkparams = {}
        self.forceNullM = forceNullM # zero matrix element if comp = 2 and pz = 0

    def parseParams(self):
        uniqParams=[]
        for n, d in enumerate(self.data):
            if d[0] not in uniqParams:
                self.params.update({ d[0]: [] })
                self.jkparams.update({ d[0]: None })
                self.avgjkparams.update({ d[0]: None })

        # Now set number of unique fitted paramers & # cfgs fitted
        self.n = len(self.params)
        self.cfgs = len(self.data)/self.n
        
    def getParams(self):
        for n, d in enumerate(self.data):
            # Fill from fit results, zeroing if comp = 2 and pz = 0
            if self.forceNullM:
                self.params[d[0]].append(0.0)
            else:
                self.params[d[0]].append(d[1])

    # Jackknife the fitted parameters
    def jk(self):
        for k, v in self.params.items():
            self.jkparams[k] = makeJks(self.params[k])

    # Form average of fitted parameters per jackknife sample
    def avgjk(self):
        for k, v in self.params.items():
            self.avgjkparams[k] = makeAvgJks(self.jkparams[k])



class correlator:
    def __init__(self,o1,o2,o3,mType,defSearchDir):
        self.snk=o1
        self.ins=o2
        self.src=o3
        self.mType=mType # res matrix element type - so far, summedRatio
        self.dsd=defSearchDir # default search directory
        self.name=self.snk.name+"."+self.ins.name+"."+self.src.name

        self.amplitude = { 1: None, 2: None }

    def resRead(self):
        for a, m in self.amplitude.items():
            # Force matelem to zero if we are reading imaginary data at pz=0 (i.e. comp = 2, pz = 0)
            zeroOut=False
            # if a == 2 and self.ins.disp == '':
            #     zeroOut=True

            compTag=""
            if a == 1:
                compTag="RE"
            if a == 2:
                compTag="IM"
            
            # Extract matrix element from fit results - potentially zeroing out
            # self.amplitude[a] = matelem(np.genfromtxt("%s%s.%s.comp_%d.RES.dat"%\
            self.amplitude[a] = matelem(np.genfromtxt("%s%s.%s.%s_fit.dat"%\
                                                      (self.dsd,self.name,self.mType,compTag),dtype=None),\
                                        forceNullM=zeroOut)        



class ratio:
    def __init__(self,N,D):
        self.nr=N.amplitude[1].avgjkparams['B']
        self.ni=N.amplitude[2].avgjkparams['B']
        self.dr=D.amplitude[1].avgjkparams['B']
        self.di=D.amplitude[2].avgjkparams['B']
        self.avgjkratio=np.zeros(len(self.nr))
        self.R = { 1: {'ensem': [], 'jks': [], 'avgjks': []},
                   2: {'ensem': [], 'jks': [], 'avgjks': []} }
        # self.R = { 1: {'jks': [], 'avg': 0.0, 'err': 0.0 },
        #            2: {'jks': [], 'avg': 0.0, 'err': 0.0 } }


    # def fillEnsem(self):
    #     for g in range(0,len(self.nr)):
    #         self.R[1]['ensem'].append( self.nr[g]

    # Form the ratio per jackknife ensemble average
    def fillJkAvgs(self):
        for g in range(0,len(self.nr)):
            self.R[1]['jks'].append((self.nr[g]*self.dr[g]+self.ni[g]*self.di[g])/(self.dr[g]**2+self.di[g]**2))
            self.R[2]['jks'].append((self.dr[g]*self.ni[g]-self.nr[g]*self.di[g])/(self.dr[g]**2+self.di[g]**2))

    def avg(self):
        for j in range(0,len(self.R[1]['jks'])):
            self.R[1]['avg'] += self.R[1]['jks'][j]
            self.R[2]['avg'] += self.R[2]['jks'][j]
        self.R[1]['avg']*=(1.0/len(self.R[1]['jks']))
        self.R[2]['avg']*=(1.0/len(self.R[2]['jks']))



class pitd:
    def __init__(self,B,Bnd,N,Nnd):
        self.B_r=B.amplitude[1].params['A']
        self.B_i=B.amplitude[2].params['A']
        self.N_r=N.amplitude[1].params['A']
        self.N_i=N.amplitude[2].params['A']
        self.Bnd_r=Bnd.amplitude[1].params['A']
        self.Bnd_i=Bnd.amplitude[2].params['A']
        self.Nnd_r=Nnd.amplitude[1].params['A']
        self.Nnd_i=Nnd.amplitude[2].params['A']
        self.pITD = { 1: {'jks': [], 'avg': 0.0, 'err': 0.0 },
                      2: {'jks': [], 'avg': 0.0, 'err': 0.0 } }

    def fillJkAvgs(self):
        for g in range(0,len(self.B_r)):
            self.pITD[1]['jks'].append((self.B_r[g]*self.Nnd_r[g])/(self.Bnd_r[g]*self.N_r[g]))
            self.pITD[2]['jks'].append((self.B_i[g]*self.Nnd_r[g])/(self.Bnd_r[g]*self.N_r[g]))

    def avg(self):
        for c, v in self.pITD.items():
            for j in range(0,len(self.B_r)):
                v['avg'] += v['jks'][j]
            v['avg'] *= (1.0/len(self.B_r))

    def err(self):
        for c, v in self.pITD.items():
            dum=0.0
            dumc=len(self.B_r)
            for j in range(0,dumc):
                dum += np.power( v['jks'][j] - v['avg'], 2)
            v['err'] = np.sqrt( ((dumc-1)/(1.0*dumc))*dum)




class pitdPoly:
    def __init__(self,a,b,c,nu):
        # self.a = a
        # self.b = b
        # self.c = c
        self.paramOrder = {0: 'a', 1: 'b', 2: 'c', 3: 'rChi2'}


    # Evaulate the polynomial fit
    def func(self,comp,nu):
        if comp == 0:
            return 1.0+self.avgParams[0]['a']*self.nu**2+self.avgParams[0]['b']*self.nu**4\
                +self.avgParams[0]['c']*self.nu**6
        if comp == 1:
            return self.avgParams[1]['a']*self.nu+self.avgParams[1]['b']*self.nu**3\
                +self.avgParams[1]['c']*self.nu**5


    # Evaluate the derivative of polynomial fit
    def dfunc(self,comp,nu):
        partials={} # local dict for ordering partials
        for k in self.avgParams[0].keys():
            partials.update({k: 0.0})

        if comp == 0:
            partials['a'] = self.nu**2
            partials['b'] = self.nu**4
            partials['c'] = self.nu**6
            partials['nu'] = 2*self.avgParams[0]['a']*self.nu+4*self.avgParams[0]['b']*self.nu**3\
                             +6*self.avgParams[0]['c']*self.nu**5
        if comp == 1:
            partials['a'] = self.nu
            partials['b'] = self.nu**3
            partials['c'] = self.nu**5
            partials['nu'] = self.avgParams[1]['a']+3*self.avgParams[1]['b']*self.nu**2\
                             +5*self.avgParams[1]['c']*self.nu**4

        error=0.0
        for i in self.avgParams[comp]:
            for j in self.avgParams[comp]:
                error+=partials[i]*self.cov[comp][(i,j)]*partials[j]
        return np.sqrt(error)

            
    # Plot the fit
    def plotFit(self,axR,axI):
        axR.plot(self.nu,self.func(0,self.nu),color=mainColors[self.zsep])
        axR.fill_between(self.nu,self.func(0,self.nu)+self.dfunc(0,self.nu),\
                         self.func(0,self.nu)-self.dfunc(0,self.nu),color=mainColors[self.zsep],alpha=0.25)
        axI.plot(self.nu,self.func(1,self.nu),color=mainColors[self.zsep])
        axI.fill_between(self.nu,self.func(1,self.nu)+self.dfunc(1,self.nu),\
                         self.func(1,self.nu)-self.dfunc(1,self.nu),color=mainColors[self.zsep],alpha=0.25)
        




class polyFit(pitdPoly):
    def __init__(self,r,i,zlabel):
        self.realJks=r['jks']
        self.imagJks=i['jks']
        self.dumPoly = pitdPoly(0,0,0,0)
        self.cfgs = len(self.realJks)
        self.avgParams={0: {}, 1: {}}
        self.cov={0: {}, 1: {}}

        self.nu=np.linspace(0,20,800)
        self.zsep=zlabel
        

    def getAvgParams(self):
        for k,v in self.dumPoly.paramOrder.items():
            avgR=0.0
            avgI=0.0
            for g in range(0,self.cfgs):
                avgR += float(self.realJks[g][k])
                avgI += float(self.imagJks[g][k])
            avgR*=(1.0/(1.0*self.cfgs))
            avgI*=(1.0/(1.0*self.cfgs))

            # Add these average values
            self.avgParams[0].update({v: avgR})
            self.avgParams[1].update({v: avgI})

    def getParamCov(self):
        for lk,lv in self.dumPoly.paramOrder.items():
            for rk,rv in self.dumPoly.paramOrder.items():
                key=(lv,rv)
                dumR=0.0
                dumI=0.0
                for g in range(0,self.cfgs):
                    dumR += (float(self.realJks[g][lk])-self.avgParams[0].get(lv))*\
                            (float(self.realJks[g][rk])-self.avgParams[0].get(rv))
                    dumI += (float(self.imagJks[g][lk])-self.avgParams[1].get(lv))*\
                            (float(self.imagJks[g][rk])-self.avgParams[1].get(rv))

                dumR *= ((1.0*(self.cfgs-1))/self.cfgs)
                dumI *= ((1.0*(self.cfgs-1))/self.cfgs)
                
                self.cov[0].update({key: dumR})
                self.cov[1].update({key: dumI})

        
