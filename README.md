# BabayPIC

## Overview 
**BabyPIC** is a 2D particle in cell (PIC) code written in Python that solves Maxwell's Equations 
on a Yee grid and moves finite phase fluid elements (FPFEs) over the grid. The 2D Yee grid 
assumes a magnetic field in the $z$ direction and an electric field in the $x-y$ plane. The motion
of FPFEs generates currents in the $x$ and the $y$ directions, which are then used to compute the
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
The units are set in a manner such that $\c=\epsilon=\mu=1$
*Initial Integrator*  
$-\nabla^2 V = \rho \rightarrow -\vec{\nabla} V = \vec{E}_{FPFE}$
$\vec{E}^0 = \vec{E}_{FPFE} + \vec{E}_{back}$, $\vec{B}^{1/2} = \vec{B}_{back}$
$\vec{v}_{FPFE}^{1/2} = \vec{v}_{FPFE}^0 + \frac{qdt}{2m}(\vec{E}^0 + \vec{v}_{FPFE}^0 \times \vec{B}^{1/2})$

*Main Integrator*
$\frac{\partial \vec{E}}{\partial t} = \vec{\nabla} \times \vec{B} - 4\pi\vec{J}$
$\frac{\partial \vec{B}}{\partial t} = -\vec{\nabla} \times \vec{E}$
$\frac{d(\gamma \vec{v})}{dt} = q(\vec{E} + \vec{v}\times\vec{B})$ ($\vec{v}$ is 3-velocity in lab frame)
$\frac{d\vec{x}}{dt} = \vec{v}$
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
