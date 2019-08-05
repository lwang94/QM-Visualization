# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

###############################################################################################################
#Functions
def InitNonStationaryState(x, c, n, a, M): #generates non-stationary state as sum of particle-in-a-box stationary states at initial time, t=0
    NSS=[0]*len(x)
    for i in range(M):
        NSS+=c[i]*np.sqrt(2/a)*np.sin(n[i]*np.pi*x/a)
    return NSS

def wrapper_InitNonStationaryState(x, n, a, M, *args): #wrapper function to generate parameters for InitNonStationaryState when using curve fitting function (opt.curve_fit)
    c= list(args[0][:M])
    return InitNonStationaryState(x, c, n, a, M) 

def func(x, a): #generates the wavefunction to be fitted
    y=[]
    #dummy function
    for i in range(len(x)):
        if x[i]<=a/2:
            y.append(np.sqrt(960/a**5)*x[i]*(a/2-x[i]))
        else:
            y.append(0)
    return y

def RSS(y_fit, y): #generates residual sum of squares for a fit
    return sum(y_fit-y)**2
    
#######################################################################################################  
#Variables
h_bar=6.626e-34/(2*np.pi) #reduced Planck constant
m=9.11e-31 #mass of electron
a=1. #length of box

x=np.linspace(0, a, 100) #x-axis
wavefunc=func(x, a) #wavefunction of question 2

M=10 #number of states
n=list(np.arange(1, M+1)) #list containing values of n from 1 to M

c=np.random.randint(-50, 50, size=M) #coefficients in infinite sum (easier to work with integers)
norm_c=(1/np.sqrt(sum(c**2)))*c #normalization factor

####################################################################################################
#Calculations
Matrix=np.zeros((M, len(x))) #empty matrix
for i in range(M):
    M_iter=i+1 #variable values of M
    n_iter=n[:M_iter] #variably sized list containing values of n
    norm_c_iter=norm_c[:M_iter] #variably sized list containing values of normalized coefficients

    popt, pcov=opt.curve_fit(lambda x, *norm_c_iter: wrapper_InitNonStationaryState(x, n_iter, a, M_iter, norm_c_iter), x, wavefunc, p0=norm_c_iter) #curve fitting to find coefficients in sum of stationary states that make up the wavefunction
    Matrix[i]=InitNonStationaryState(x, popt, n_iter, a, M_iter) #set values of non-stationary state fitting matrix 

###################################################################################################
#Plotting
fig, axes=plt.subplots(M/2, 2, sharex=True, sharey=True)

for i, ax in enumerate(axes.flatten()):
    _RSS=RSS(Matrix[i], np.array(wavefunc))
    
    ax.plot(x, Matrix[i], label='M=%d, RSS=%.2E'%(n[i], _RSS))
    ax.plot(x, wavefunc, label='Wavefunction')
    
    M_plot=i+1
    Title=''
    for j in range(M_plot):
        Title+=r'%.1E$\Psi_{%d}$+'%(popt[j], n[j])
    ax.set_title(Title[:-1], fontsize=9)
    ax.legend()

plt.show()
