# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 23:17:23 2019

@author: lawre
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import sympy as sp
import itertools

##################################################################################################################
#Functions
def HO_Hamiltonian(h_bar, m, k, func, x): #generates Hamiltonian of Harmonic Oscillator
    Kinetic=(-h_bar**2/(2*m))*sp.diff(func, x, 2)
    Potential=((k*x**2)/2)*func
    return Kinetic+Potential

    
def HO_SS(m, w, h_bar, n, x): #generates stationary state of Harmonic Oscillator
    z=np.sqrt(m*w/h_bar)*x
    Hermite=((-1)**n)*sp.exp(x**2)*sp.diff(sp.exp(-x**2), x, n)
    Hermite=Hermite.subs(x, z)
    return 1/np.sqrt(2**n*np.math.factorial(n))*(m*w/(np.pi*h_bar))**(1./4)*sp.exp(-m*w*x**2/(2*h_bar))*Hermite

def Cross_term_x(c, x, func): #one half of cross term 
    return c*x*func

def Cross_term_y(y, func): #second half of cross term
    return y*func

def System_Hamiltonian_direc(r, n, h_bar, m, kr, wr, initial, final): #generates harmonic oscillator system hamiltonian in specified direction
    Sys_Hamiltonian_direc=np.zeros((len(n), len(n)))
    for i, j in itertools.product(n, n):
        func_i=HO_SS(m, wr, h_bar, i, x)
        func_j=HO_SS(m, wr, h_bar, j, x)       
        
        integrand=func_j*HO_Hamiltonian(h_bar, m, kr, func_i, x) #bra-ket notation is <HO_i|HO Hamiltonian|HO_j>
        integrand=sp.lambdify(x, integrand, "numpy") #converts symbolic computation to numeric computation
        
        Sys_Hamiltonian_ij, err=integrate.quad(integrand, initial, final)
        Sys_Hamiltonian_direc[i, j]=Sys_Hamiltonian_ij

    return Sys_Hamiltonian_direc


def Crossterm_Hamiltonian_x(x, n, h_bar, m, initial, final): #generates system hamiltonian for coupling cross term in the x-direction
    Crossterm_Hamiltonian_x=np.zeros((len(n), len(n)))
    for i, j in itertools.product(n, n):
        func_i=HO_SS(m, wx, h_bar, i, x)
        func_j=HO_SS(m, wx, h_bar, j, x)       
        
        integrand=func_j*Cross_term_x(c, x, func_i) #bra-ket notation is <HO_i|Crossterm(x) Hamiltonian|HO_j>
        integrand=sp.lambdify(x, integrand, "numpy") #converts symbolic computation to numeric computation
        
        Crossterm_Hamiltonian_ij, err=integrate.quad(integrand, initial, final)
        Crossterm_Hamiltonian_x[i, j]=Crossterm_Hamiltonian_ij

    return Crossterm_Hamiltonian_x
    
def Crossterm_Hamiltonian_y(y, n, h_bar, m, initial, final): #generates system hamiltonian for coupling cross term in the x-direction
    Crossterm_Hamiltonian_y=np.zeros((len(n), len(n)))
    for i, j in itertools.product(n, n):
        func_i=HO_SS(m, wy, h_bar, i, y)
        func_j=HO_SS(m, wy, h_bar, j, y)       
        
        integrand=func_j*Cross_term_y(y, func_i) #bra-ket notation is <HO_i|Crossterm(y) Hamiltonian|HO_j>
        integrand=sp.lambdify(y, integrand, "numpy") #converts symbolic computation to numeric computation
        
        Crossterm_Hamiltonian_ij, err=integrate.quad(integrand, initial, final)
        Crossterm_Hamiltonian_y[i, j]=Crossterm_Hamiltonian_ij

    return Crossterm_Hamiltonian_y
    
def System_Hamiltonian_total(x, y, n, h_bar, m, kx, ky, initial, final): #generates system hamiltonian for entire 2D operator
    Sys_Hamiltonian_x=System_Hamiltonian_direc(x, n, h_bar, m, kx, wx, initial, final)
    Sys_Hamiltonian_y=System_Hamiltonian_direc(y, n, h_bar, m, ky, wy, initial, final)
    Ct_Hamiltonian_x=Crossterm_Hamiltonian_x(x, n, h_bar, m, initial, final)
    Ct_Hamiltonian_y=Crossterm_Hamiltonian_y(y, n, h_bar, m, initial, final)
    
    #combined separate system hamiltonians
    Sys_Hamiltonian=np.kron(Sys_Hamiltonian_x, np.eye(len(n)))+np.kron(np.eye(len(n)), Sys_Hamiltonian_y)+np.kron(Ct_Hamiltonian_x, Ct_Hamiltonian_y)
    
    return Sys_Hamiltonian

def Approx_Energy(x, n, initial, final): #finds approximate energy of eigenstate using system hamiltonian
    Sys_Hamiltonian=System_Hamiltonian_total(x, y, n, h_bar, m, kx, ky, initial, final)
    eigval, eigvec=np.linalg.eigh(Sys_Hamiltonian)
    
    return eigval.min()

#################################################################################################
#Variables
h_bar=1.#reduced planck constant in natural units
m=1.#mass of particle

wx=1 #angular frequency in x direction
wy=1 #angular frequency in y direction
kx=m*wx**2 #force constant in x direction
ky=m*wy**2 #force constant in y direction
c=1.0 #coupling term

a=10. #size of system
initial=-2*a#lower bound of integration
final=2*a #upper bound of integration

Q=[1, 2, 3, 4, 5, 6, 7, 8, 9] #number of basis functions



x=sp.Symbol('x') #x as a symbol in symbolic computation
y=sp.Symbol('y') #y as a symbol in symbolic computation

###########################################################################################################
#Calculations    
Approx_E=[]
for q in range(len(Q)):
    n=np.arange(0, Q[q]+1)
    E_Approx=Approx_Energy(x, n, initial, final)
    Approx_E+=[E_Approx]

print Approx_E

#################################################################################################################
#Plotting
plt.title("Convergence for Lambda=%f"%(c))
plt.xlabel('Number of Basis Functions')
plt.ylabel("Approximate Ground State Energy")
plt.plot(Q, Approx_E)