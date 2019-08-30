# QM Visualization
This is a collection of Python programs that graphs solutions to quantum mechanics problems
for research or educational purposes.

## Technologies
Project is created with:
Python version: 2.7
### Libraries
NumPy version: 1.15.3
Matplotlib version: 2.2.3
SciPy version: 0.15.1

## Nonstationary_state_time_evolution.py
Program that graphs how the probability density, real part and imaginary part of a 
wavefunction evolves with time. In this case, the Hamiltonian is for the Particle in a Box
and the wavefunction is a quantum superposition of random stationary states with random
coefficients. The variable M defines how many stationary states are used to create the 
wavefunction.

## fit_to_nonstationary_state_with_increasing_num_stationary_states.py
Program that approximates a wavefunction using a quantum superposition of stationary states. 
In this case, the stationary states are for the Particle in a Box. The default function to 
be fitted to is:
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Psi&space;=\sqrt{960/a^{5}}\left&space;[&space;x(\frac{a}{2}-x)&space;\right&space;]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Psi&space;=\sqrt{960/a^{5}}\left&space;[&space;x(\frac{a}{2}-x)&space;\right&space;]" title="\Psi =\sqrt{960/a^{5}}\left [ x(\frac{a}{2}-x) \right ]" /></a>
for x<a/2, otherwise <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Psi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Psi" title="\Psi" /></a>=0. The variable M defines the total number of states used in the
approximation. The program also graphs how the approximation becomes better as the number
of states used increases.

## time_evolution_of_custom_wavefunction.py
Program that graphs how the probability density of a wavefunction evolves with time. In this
case, the Hamiltonian is for the Particle in a Box and crucially, the wavefunction is 
defined specifically by the user. The default function is:
<a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Psi&space;=\sqrt{960/a^{5}}\left&space;[&space;x(\frac{a}{2}-x)&space;\right&space;]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Psi&space;=\sqrt{960/a^{5}}\left&space;[&space;x(\frac{a}{2}-x)&space;\right&space;]" title="\Psi =\sqrt{960/a^{5}}\left [ x(\frac{a}{2}-x) \right ]" /></a>
for x<a/2, otherwise <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Psi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Psi" title="\Psi" /></a>=0. The variable M defines the number of stationary states used
in the approximation

## Status
All three programs can be used and the variables can be adjusted for the users purpose. 
### TO DO:
The variables should be coded as user inputs rather the user having to hard code them in.
