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

## Approximate_ground_state_energy_of_harmonic_oscillator_using_particle_in_box_wavefunctions.py
Program that approximates the energy of the ground state harmonic oscillator using a variational method and particle in a box basis wavefunctions. An explanation of the principle behind the variational method can be found here: https://en.wikipedia.org/wiki/Variational_method_(quantum_mechanics). The output plot is a graph of the approximation to the ground state energy as a function of the number of basis functions. It should converge to the actual ground state energy as the number of basis functions increases.

## Approximation_of_coupled_2D_Harmonic_Oscillator.py
Program that approximates the ground state energy of a 2 dimensional coupled harmonic oscillator. The Hamiltonian in this case is: 
H = H<sub>x</sub> + H<sub>y</sub> + &lambda;xy, where H<sub>x</sub> and H<sub>y</sub> is a 1D harmonic oscillator hamiltonian in the x and y direction respectively. The output plot is a graph of the approximation as a function of the number of basis functions. It should converge as the number of basis functions increases. 

## Status
All programs can be used and the variables can be adjusted for the users purpose. 
### TO DO:
The variables should be coded as user inputs rather the user having to hard code them in.
