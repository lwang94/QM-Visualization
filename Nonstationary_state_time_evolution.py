# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

###############################################################################################################
#Functions
def NonStationaryState(h_bar, m, norm_c, a, n, x, time): #generates real and imaginary parts of nonstationary state as sum of stationary states for particle in a box
    E=(h_bar*np.pi*n)**2/(2*m*a**2) #energy
    NSS=[0]*len(x)
    for i in range(len(norm_c)):
        NSS+=norm_c[i]*np.sqrt(2/a)*np.sin(n[i]*np.pi*x/a)*np.exp(-1j*E[i]*time/h_bar) #returns sum of stationary states multiplied by arbitrary coefficient and time-dependence
    NSS_real=np.real(np.array(NSS)) #real part of wavefunction
    NSS_imag=np.imag(np.array(NSS)) #imaginary part of wavefunction
    return NSS_real, NSS_imag
    
def ProbDensity(h_bar, m, norm_c, a, n, x, time): #generates probability density of wave-function for particle in a box at specific time
    E=(h_bar*np.pi*n)**2/(2*m*a**2) #energy
    NSS=[0]*len(x)
    for i in range(len(norm_c)):
        NSS+=norm_c[i]*np.sqrt(2/a)*np.sin(n[i]*np.pi*x/a)*np.exp(-1j*E[i]*time/h_bar) #returns sum of stationary states multiplied by arbitrary coefficient and time-dependence
    ProbDensity=NSS*np.conj(NSS) #converts wavefunction into probability density function
    return ProbDensity


#######################################################################################################  
#Variables
h_bar=6.626e-34/(2*np.pi) #reduced Planck constant
m=9.11e-31 #mass of electron
a=1. #length of box

x=np.linspace(0, a, 100) #x-axis (position space)

M=2 #number of states
n=np.random.randint(1, 5, size=M) #vector containing random values of n (the order of state)
c=np.random.randint(-10, 10, size=M) #vector containig random coefficients in infinite sum (easier to work with integers)
norm_c=1/np.sqrt(sum(c**2))*c #normalized coefficients

t=np.arange(0, 2500, 500) #time

######################################################################################################
#Calculations
ProbMatrix=np.zeros((len(t), len(x))) #empty matrix 
realMatrix=np.zeros((len(t), len(x))) #empty matrix
imagMatrix=np.zeros((len(t), len(x))) #empty matrix

for i in range(len(t)):
    ProbMatrix[i]=ProbDensity(h_bar, m, norm_c, a, n, x, t[i]) #set values to probability density matrix
    realMatrix[i], imagMatrix[i]=NonStationaryState(h_bar, m, norm_c, a, n, x, t[i]) #set values to real and imaginary matrices

#######################################################################################################
#Plotting
fig1, axes1=plt.subplots(len(t), 1, sharex=True, sharey=True)

for i, ax1 in enumerate(axes1.flatten()):
    ax1.plot(x, ProbMatrix[i], label='Probability Density')
    ax1.plot(x, realMatrix[i], ':', label='Real part')
    ax1.plot(x, imagMatrix[i], ':', label='Imaginary part')
    ax1.set_title('t=%d'%(t[i]))
    ax1.legend()
    
Title='Wavefunction='
for i in range(M):
    Title+=r'%.1E$\Psi_{%d}$+'%(norm_c[i], n[i])

fig1.suptitle(Title[:-1])

plt.show()