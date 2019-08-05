# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.integrate as integrate

###############################################################################################################
#Functions    
def ProbDensity(h_bar, m, norm_c, a, n, x, time): #generates probability density of wave-function for particle in a box at specific time
    E=(h_bar*np.pi*n)**2/(2*m*a**2) #energy
    NSS=[0]*len(x)
    for i in range(len(norm_c)):
        NSS+=norm_c[i]*np.sqrt(2/a)*np.sin(n[i]*np.pi*x/a)*np.exp(-1j*E[i]*time/h_bar) #returns sum of stationary states multiplied by arbitrary coefficient and time-dependence
    ProbDensity=NSS*np.conj(NSS) #converts wavefunction into probability density function
    return ProbDensity

def InitNonStationaryState(x, c, n, a, M): #generates non-stationary state as sum of particle-in-a-box stationary states at initial time, t=0
    NSS=[0]*len(x)
    for i in range(M):
        NSS+=c[i]*np.sqrt(2/a)*np.sin(n[i]*np.pi*x/a)
    return NSS

def wrapper_InitNonStationaryState(x, n, a, M, *args): #wrapper function to generate parameters for InitNonStationaryState when using curve fitting function (opt.curve_fit)
    c= list(args[0][:M])
    return InitNonStationaryState(x, c, n, a, M) 

def func(x, a): #generates the wavefunction
    y=[]
    for i in range(len(x)):
        if x[i]<=a/2:
            y.append(np.sqrt(960/a**5)*x[i]*(a/2-x[i]))
        else:
            y.append(0)
    return y

#######################################################################################################  
#Variables
h_bar=6.626e-34/(2*np.pi) #reduced Planck constant
m=9.11e-31 #mass of electron
a=1. #length of box

x=np.linspace(0, a, 100) #x-axis
t=np.arange(0, 2500, 500) #time

wavefunc=func(x, a) #wavefunction

M=10 #number of states
n=np.arange(1, M+1) #list containing nth stationary state
c=np.random.randint(-50, 50, size=M) #randomized list of coefficients
norm_c=1/np.sqrt(sum(c**2))*c #normalize
norm_c, pcov=opt.curve_fit(lambda x, *norm_c: wrapper_InitNonStationaryState(x, n, a, M, norm_c), x, wavefunc, p0=norm_c) #use curve fitting to obtain correct list of coefficients

####################################################################################################################
#Calculations
Matrix=np.zeros((len(t), len(x))) #empty matrix
for i in range(len(t)):
    Matrix[i]=ProbDensity(h_bar, m, norm_c, a, n, x, t[i]) #set values to probability density matrix

################################################################################################################
#Plotting
fig, axes=plt.subplots(len(t), 1, sharex=True, sharey=True)

for i, ax in enumerate(axes.flatten()):
    ax.plot(x, Matrix[i], label='t=%d'%(t[i]))
    ax.plot(x, np.array(wavefunc)**2, label='original wavefunction')
    ax.legend()
fig.suptitle('Question 3: M=%d'%(M))

plt.show()
