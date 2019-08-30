# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 13:51:35 2018

@author: lawre
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import sympy as sp
import itertools

####################################################################################################
#Functions
def PiaB_SS(a, n, x): #generates stationary state of symmetric Particle in a Box
    if n%2==0:     
        return sp.Piecewise((0, x<=-a/2), (0, x>=a/2), (np.sqrt(2/a)*sp.sin(n*np.pi*x/a), True))
    else:
        return sp.Piecewise((0, x<=-a/2), (0, x>=a/2), (np.sqrt(2/a)*sp.cos(n*np.pi*x/a), True))
    
def HO_Hamiltonian(h_bar, m, k, func, x): #generates Hamiltonian of Harmonic Oscillator
    Kinetic=(-h_bar**2/(2*m))*sp.diff(func, x, 2)
    Potential=((k*x**2)/2)*func
    return Kinetic+Potential
        
def HO_SS(m, w, h_bar, n, x): #generates stationary state of Harmonic Oscillator
    z=np.sqrt(m*w/h_bar)*x
    Hermite=((-1)**n)*sp.exp(x**2)*sp.diff(sp.exp(-x**2), x, n)
    Hermite=Hermite.subs(x, z)
    return 1/np.sqrt(2**n*np.math.factorial(n))*(m*w/(np.pi*h_bar))**(1./4)*sp.exp(-m*w*x**2/(2*h_bar))*Hermite

def System_Hamiltonian(x, n, h_bar, m, k, initial, final, a): #generates system hamiltonian for particle in a box with harmonic oscillator hamiltonian
    Sys_Hamiltonian=np.zeros((len(n), len(n)))
    for i, j in itertools.product(n, n):
        func_i=PiaB_SS(a, i+1, x)
        func_j=PiaB_SS(a, j+1, x)        
        
        integrand=func_j*HO_Hamiltonian(h_bar, m, k, func_i, x) #bra-ket notation is <PiaB_i|HO Hamiltonian|PiaB_j>
        integrand=sp.lambdify(x, integrand, "numpy") #converts symbolic computation to numeric computation
        
        Sys_Hamiltonian_ij, err=integrate.quad(integrand, initial, final)
        Sys_Hamiltonian[i, j]=Sys_Hamiltonian_ij

    return Sys_Hamiltonian

def Approx_Energy(x, n, a, initial, final, HO_state, SysHam): #finds approximate energy of harmonic oscillator state using system hamiltonian
    E_approx=0
    for i, j in itertools.product(n, n):
        il=PiaB_SS(a, i+1, x)*HO_state
        il=sp.lambdify(x, il, "numpy") #converts symbolic computation to numeric computation
        braket_il, err=integrate.quad(il, initial, final)
        
        lj=HO_state*PiaB_SS(a, j+1, x)
        lj=sp.lambdify(x, lj, "numpy") #converts symbolic computation to numeric computation
        braket_lj, err=integrate.quad(lj, initial, final)
        
        E_approx+=braket_il*SysHam[i, j]*braket_lj #bra-ket notation is <HO_0|PiaB_i>SystemHamiltonian(i, j)<PiaB_j|HO_0>
    return E_approx

#######################################################################################################
#Variables
h_bar=1. #reduced planck constant in natural units
m=1.#mass of particle
w=1.#angular frequency
k=m*w**2 #spring constant

a=10. #size of box for Particle in a Box
initial=-2*a#lower bound of integration
final=2*a #upper bound of integration

Q=[2, 3, 4, 5, 6, 7] #number of basis functions

x=sp.Symbol('x') #x as a symbol in symbolic computation

HO_state=HO_SS(m, w, h_bar, 0, x) #ground state for Harmonic Oscillator

##########################################################################################################
#Calculations
Approx_E=[]
for q in range(len(Q)):
    n=np.arange(0, Q[q]+1)
    SysHam=System_Hamiltonian(x, n, h_bar, m, k, initial, final, a)
    E_Approx=Approx_Energy(x, n, a, initial, final, HO_state, SysHam)
    Approx_E+=[E_Approx]

print Approx_E
Actual_E=[5.00e-01]*len(Q)

############################################################################################################
#Plotting
plt.title("Convergence Plot")
plt.plot(Q, Approx_E, label='Approximate Ground State Energy')
plt.plot(Q, Actual_E, label='Actual Ground State Energy')
plt.xlabel('Number of Basis Functions')
plt.ylabel('Energy')
plt.legend()
