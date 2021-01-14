# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:40:45 2020

@author: dmgre
@minor modifications: mreal
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
import pandas as pd 
from datetime import datetime

print("import finished")

#%%
# Constants to use ---------------------------------------------------
#============================================================================
mGaAs = 0.067
print("mass GaAs (factor)= {}".format(mGaAs))
# Densidad de portadores rho = 3.06e15m^-2

Ef = (np.pi*3.06e15*constants.hbar**2)/(constants.k*constants.m_e * mGaAs)
print("Fermi energy (K)= {}".format(Ef))
#Energia de Fermi dividido constante de Boltzmann

prefactor = constants.e*Ef*constants.k/(np.pi**2*constants.hbar)
print("prefactor= {}".format(prefactor))
# Combinacion de factores que aparece seguido


# Function definition  ---------------------------------------------------
#============================================================================
def murulo(nu,T,Gamma,l):
    # murulo esta en unidades de Ef
    # T hay que ponerlo en unidades de Kelvin
    # Gamma es la temperatura de Dingle y 
    # tiene la constante de Bolztmann
    # en estas unidades nu = 1/B
    # Ef es la energia de Fermi y es una constante global
    x = 2*np.pi**2*T*nu*l
    g = -2*np.pi*Gamma*l*nu
    s = 2*np.pi*l*(nu+.5)
    x, g = x/Ef,g/Ef
    return -2*np.pi**2*T/np.sinh(x)*np.sin(s)*np.exp(g)

def murulosuma(nu,T,Gamma,lmax):
    a = murulo(nu,T,Gamma,1)
    for l in range(2,lmax):
        a = murulo(nu, T, Gamma, l)+a
    return a

def murulosuma2(nu,T,Gamma):
    a = murulo(nu,T,Gamma,1)
    for l in range(2,200):
        a = murulo(nu, T, Gamma, l)+a
    return a

def Jrulo(nu,T,Gamma,l):
    # Ef = 126.87
#    Ef = 500
    x = 2*np.pi**2*T*nu*l
    g = -2*np.pi*Gamma*l*nu
    s = 2*np.pi*l*(nu+.5)
    x, g = x/Ef,g/Ef
    phi = x*(x*1/np.tanh(x)-1)/np.sinh(x)*1/l
    return phi*np.sin(s)*np.exp(g)

def Jrulosuma(nu,T,Gamma):
    a = Jrulo(nu, T, Gamma, 1)
    for l in range(2,150):
        a = Jrulo(nu, T, Gamma, l)+a
    return a

#%%
nu = np.linspace(1e-1,1,100)
Temp = np.array([.3,.8,1.5])
mu = np.zeros((nu.size,Temp.size))
for t in range(Temp.size):
    for i in range(nu.size):
        mu[i,t] = murulosuma2(nu[i], Temp[t], .1)

plt.figure()
for t in range(3):
    plt.plot(nu,mu[:,t],label = 'T = {}'.format(Temp[t]))
plt.legend(loc=0)
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\mu/E_f$')
plt.show()
#%%  Este es el importante
Temp = np.array([.3,.4,.5,.8,1.5])
# nu = np.linspace(1e-5,19,1000)
nu = np.linspace(1e-5,20,1000)
Corriente = np.zeros((nu.size,Temp.size),dtype=np.double)
Vtp = np.zeros((Corriente.shape),dtype=np.double)
dT = 1e-3
# AA = Ef*constants.k/constants.e
AA = Ef*constants.k/(constants.e*np.pi)
for t in range(Temp.size):
    for i in range(nu.size):
        Corriente[i,t] = Jrulosuma(nu[i], Temp[t], .1)
        # Vtp[i,t] = 0.01086*(1/nu[i])*Corriente[i,t]*dT/Temp[t]
        Vtp[i,t] = AA*(1/nu[i])*Corriente[i,t]*dT/Temp[t]
        
plt.figure()
for t in range(Temp.size):
    plt.plot(nu,Vtp[:,t],label = 'T = {}'.format(Temp[t]))
plt.legend(loc=0)
# plt.xlim((5,10))
plt.xlabel(r'$\nu$')
plt.ylabel(r'$Vtp$')
plt.show()
#%% Save data VTP vs nu
now = datetime.now()
dd = now.strftime("%Y%m%d %H%M%S")
f = open('murulodata'+dd+'.txt','w')
f.write('Temperatura;nu;VTP')
for i in range(nu.size):
    for t in range(Temp.size):
        f.write('\n'+"%8.5f"%Temp[t]+";"+
                "%8.5f"%nu[i]+";"+
                "%8.16f"%Vtp[i,t])
f.close()

#%% plot como funcion de B
plt.figure()
# cte = np.pi.constants.hbar*constants.c
for t in range(Temp.size):
    plt.plot(1/nu,Vtp[:,t],label = 'T = {}'.format(Temp[t]))
# plt.xlim((0,5))
plt.legend(loc=0)
plt.xlabel(r'$B[T]$')
plt.ylabel(r'$VTP$')
plt.show()
#%%
nu = np.linspace(1e-5,50,1000)
CorrienteT = np.zeros((1000,1))
for i in range(Temp.size):
    CorrienteT[i] = Jrulosuma(nu[i],.5,.1)
plt.figure()
plt.plot(nu,CorrienteT)
plt.show()

#%%
columns = [str(str_string) for str_string in Temp]

df = pd.DataFrame(data=Vtp, columns=columns)
df.insert(0, "nu", nu, True) 
filename = 'abc.csv'

df.to_csv(filename, index=False)
