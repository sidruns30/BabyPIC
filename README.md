
# BabayPIC

## Overview 
**BabyPIC** is a 2D particle in cell (PIC) code written in Python that solves Maxwell's Equations 
on a Yee grid and moves finite phase fluid elements (FPFEs) over the grid. The 2D Yee grid 
assumes a magnetic field in the z direction and an electric field in the x-y plane. The motion
of FPFEs generates currents in the x and the y directions, which are then used to compute the
electric fields. The current version does not support collisions between the FPFEs. Furthermore, only
conducting and periodic boundary conditions are implemented in the current version.   

### Initial Conditions 
The user needs to enter the **number**, **charges**, **masses**, **positions** and **velocities** of
the FPFEs which are used to compute the electric and magnetic fields at the start of the simulation. 
BabyPIC uses leapfrog integration and so the initial positions must trail the initial velocities by 
half of a timestep (see timestep below).

Furthermore, the user can also enter **divergence-free background electric** and **magnetic fields** as 
an input, which are then added to the fields of the FPFEs at the start of the simulation. Note that the 
electric and magnetic fields are also offset by half of a timestep as the magnetic fields are defined
to be ahead of the electric fields in the intial condition. 

### Equations Solved
The units are set in a manner such that <a href="https://www.codecogs.com/eqnedit.php?latex=c=\epsilon&space;=\mu=1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c=\epsilon&space;=\mu=1" title="c=\epsilon =\mu=1" /></a>

*Initial Integrator*  
<a href="https://www.codecogs.com/eqnedit.php?latex=-\nabla^2&space;V&space;=&space;\rho&space;\rightarrow&space;-\vec{\nabla}&space;V&space;=&space;\vec{E}_{FPFE}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?-\nabla^2&space;V&space;=&space;\rho&space;\rightarrow&space;-\vec{\nabla}&space;V&space;=&space;\vec{E}_{FPFE}" title="-\nabla^2 V = \rho \rightarrow -\vec{\nabla} V = \vec{E}_{FPFE}" /></a>
<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{E}^0&space;=&space;\vec{E}_{FPFE}&space;&plus;&space;\vec{E}_{back}$,&space;$\vec{B}^{1/2}&space;=&space;\vec{B}_{back}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{E}^0&space;=&space;\vec{E}_{FPFE}&space;&plus;&space;\vec{E}_{back}$,&space;$\vec{B}^{1/2}&space;=&space;\vec{B}_{back}" title="\vec{E}^0 = \vec{E}_{FPFE} + \vec{E}_{back}$, $\vec{B}^{1/2} = \vec{B}_{back}" /></a>
<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{v}_{FPFE}^{1/2}&space;=&space;\vec{v}_{FPFE}^0&space;&plus;&space;\frac{qdt}{2m}(\vec{E}^0&space;&plus;&space;\vec{v}_{FPFE}^0&space;\times&space;\vec{B}^{1/2})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{v}_{FPFE}^{1/2}&space;=&space;\vec{v}_{FPFE}^0&space;&plus;&space;\frac{qdt}{2m}(\vec{E}^0&space;&plus;&space;\vec{v}_{FPFE}^0&space;\times&space;\vec{B}^{1/2})" title="\vec{v}_{FPFE}^{1/2} = \vec{v}_{FPFE}^0 + \frac{qdt}{2m}(\vec{E}^0 + \vec{v}_{FPFE}^0 \times \vec{B}^{1/2})" /></a>

*Main Integrator*
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\vec{E}}{\partial&space;t}&space;=&space;\vec{\nabla}&space;\times&space;\vec{B}&space;-&space;4\pi\vec{J}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\vec{E}}{\partial&space;t}&space;=&space;\vec{\nabla}&space;\times&space;\vec{B}&space;-&space;4\pi\vec{J}" title="\frac{\partial \vec{E}}{\partial t} = \vec{\nabla} \times \vec{B} - 4\pi\vec{J}" /></a>
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\vec{B}}{\partial&space;t}&space;=&space;-\vec{\nabla}&space;\times&space;\vec{E}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\vec{B}}{\partial&space;t}&space;=&space;-\vec{\nabla}&space;\times&space;\vec{E}" title="\frac{\partial \vec{B}}{\partial t} = -\vec{\nabla} \times \vec{E}" /></a>
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d(\gamma&space;\vec{v})}{dt}&space;=&space;q(\vec{E}&space;&plus;&space;\vec{v}\times\vec{B})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d(\gamma&space;\vec{v})}{dt}&space;=&space;q(\vec{E}&space;&plus;&space;\vec{v}\times\vec{B})" title="\frac{d(\gamma \vec{v})}{dt} = q(\vec{E} + \vec{v}\times\vec{B})" /></a> 
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d\vec{x}}{dt}&space;=&space;\vec{v}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d\vec{x}}{dt}&space;=&space;\vec{v}" title="\frac{d\vec{x}}{dt} = \vec{v}" /></a>
### Units 


### Grid Structure 

### Interpolation 

## Files 

### boundary_conditions.py

### current_deposition.py

### fields_class.py

### fpfe_class.py

### interpolate_fields.py

### particle_push.py

### poisson.py

### run_simulation.py

## References 
*
*
*
*
