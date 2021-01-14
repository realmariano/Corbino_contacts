# -*- coding: utf-8 -*-
"""
VTP de un QD
"""

import scipy.integrate 
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt
from time import process_time
from itertools import product
from scipy import constants
from datetime import datetime

np.seterr(divide='ignore', invalid='ignore')


def DFermi(es,mu,T):
    alpha = constants.hbar*constants.e/(2*.067*constants.m_e)
    kT = constants.k*T/alpha
    arg1 = (es-mu)/kT
    # recordar que mu = hbar e B/2meff. 
    # Entonces el mu que entra como variable es B
    if arg1 < 100:
        F0 = 1/(np.exp(arg1)+1)
        DF = F0**2*np.exp(arg1)
    else:
        F0 = np.exp(-arg1)/(np.exp(-arg1)+1)
        DF = F0*(1-F0)
    DF = DF/kT
    if DF != DF : 
        DF = 0
    else:
        DF = DF
    return DF

def Tau(es,g1,g2):
    # return g1/(es**2+g2**2)
    return (g1*g2)/(es**2+.25*(g1+g2)**2)

def Onsager(es,g1,g2,mu,T,idx):
    alpha = constants.hbar*constants.e/(2*.067*constants.m_e)
    return Tau(es, g1, g2)*DFermi(es, mu, T)*((es-mu)*alpha)**idx

def Integrado(g1,g2,mu,T,idx):
    E = np.linspace(-25,25,100)
    In = 0
    for i in range (E.size-1):
      In = scipy.integrate.quad(Onsager,E[i],E[i+1],args=(g1,g2,mu,T,idx))[0]+In  
    # return (scipy.integrate.quad(Onsager,-10,0,args=(g1,g2,mu,T,idx))[0]
    #         +
    #         scipy.integrate.quad(Onsager,0,10,args=(g1,g2,mu,T,idx))[0])
    return In

#%%
muvec = np.linspace(-1,1,100)
# g1,g2,T = .1,.3,1e-3
g1,T = .01,.3
dt = 1e-3
g2 = np.linspace(1e-5,1,5)
L11 = np.zeros((muvec.size,g2.size),dtype=np.double);
L12 = np.zeros((L11.shape),dtype=np.double);
for (m,g) in product(range(muvec.size),range(g2.size)):
    L11[m,g] = Integrado(g1, g2[g], muvec[m], T, 0)
    L12[m,g] = Integrado(g1, g2[g], muvec[m], T, 1)

alpha = constants.hbar*constants.e/(2*.067*constants.m_e)
VTP = L12/L11*dt/T/alpha
#%%
plt.figure()
for g in range(g2.size):
    # plt.plot(muvec,L12[:,g]/L11[:,g],
    plt.plot(muvec,VTP[:,g],
             label='g ='"%8.5f"%g2[g])
plt.xlim(-.5,.5)
plt.legend(loc=1)
plt.show()

#%% Save data in txt file
now = datetime.now()
dd = now.strftime("%Y%m%d %H%M%S")

f = open('QDotdata'+dd+'.txt','w')
f.write('g2;B;VTP')
for m in range(muvec.size):
    for g in range(g2.size):
        f.write( '\n '+ "%8.5f"%g2[g]+";"+ 
                "%8.5f"%muvec[m]+ ";"+
                "%8.16f"%VTP[m,g])
f.close()


#%% plot maxima como funcion del acople
plt.figure()
plt.plot(g2,VTP.max(0),'o')
plt.show()