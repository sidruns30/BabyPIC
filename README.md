
# BabyPIC

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
The units are set in a manner such that:
<a href="https://www.codecogs.com/eqnedit.php?latex=c=\epsilon&space;=\mu=1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c=\epsilon&space;=\mu=1" title="c=\epsilon =\mu=1" /></a>

**Initial Integrator**  

<a href="https://www.codecogs.com/eqnedit.php?latex=-\nabla^2&space;V&space;=&space;\rho&space;\rightarrow&space;-\vec{\nabla}&space;V&space;=&space;\vec{E}_{FPFE}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?-\nabla^2&space;V&space;=&space;\rho&space;\rightarrow&space;-\vec{\nabla}&space;V&space;=&space;\vec{E}_{FPFE}" title="-\nabla^2 V = \rho \rightarrow -\vec{\nabla} V = \vec{E}_{FPFE}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{E}^0&space;=&space;\vec{E}_{FPFE}&space;&plus;&space;\vec{E}_{back}$,&space;$\vec{B}^{1/2}&space;=&space;\vec{B}_{back}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{E}^0&space;=&space;\vec{E}_{FPFE}&space;&plus;&space;\vec{E}_{back}$,&space;$\vec{B}^{1/2}&space;=&space;\vec{B}_{back}" title="\vec{E}^0 = \vec{E}_{FPFE} + \vec{E}_{back}$, $\vec{B}^{1/2} = \vec{B}_{back}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\vec{v}_{FPFE}^{1/2}&space;=&space;\vec{v}_{FPFE}^0&space;&plus;&space;\frac{qdt}{2m}(\vec{E}^0&space;&plus;&space;\vec{v}_{FPFE}^0&space;\times&space;\vec{B}^{1/2})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{v}_{FPFE}^{1/2}&space;=&space;\vec{v}_{FPFE}^0&space;&plus;&space;\frac{qdt}{2m}(\vec{E}^0&space;&plus;&space;\vec{v}_{FPFE}^0&space;\times&space;\vec{B}^{1/2})" title="\vec{v}_{FPFE}^{1/2} = \vec{v}_{FPFE}^0 + \frac{qdt}{2m}(\vec{E}^0 + \vec{v}_{FPFE}^0 \times \vec{B}^{1/2})" /></a>

**Main Integrator**

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\vec{E}}{\partial&space;t}&space;=&space;\vec{\nabla}&space;\times&space;\vec{B}&space;-&space;4\pi\vec{J}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\vec{E}}{\partial&space;t}&space;=&space;\vec{\nabla}&space;\times&space;\vec{B}&space;-&space;4\pi\vec{J}" title="\frac{\partial \vec{E}}{\partial t} = \vec{\nabla} \times \vec{B} - 4\pi\vec{J}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\vec{B}}{\partial&space;t}&space;=&space;-\vec{\nabla}&space;\times&space;\vec{E}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\vec{B}}{\partial&space;t}&space;=&space;-\vec{\nabla}&space;\times&space;\vec{E}" title="\frac{\partial \vec{B}}{\partial t} = -\vec{\nabla} \times \vec{E}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d(\gamma&space;\vec{v})}{dt}&space;=&space;q(\vec{E}&space;&plus;&space;\vec{v}\times\vec{B})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d(\gamma&space;\vec{v})}{dt}&space;=&space;q(\vec{E}&space;&plus;&space;\vec{v}\times\vec{B})" title="\frac{d(\gamma \vec{v})}{dt} = q(\vec{E} + \vec{v}\times\vec{B})" /></a> 

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{d\vec{x}}{dt}&space;=&space;\vec{v}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{d\vec{x}}{dt}&space;=&space;\vec{v}" title="\frac{d\vec{x}}{dt} = \vec{v}" /></a>
### Units 


### Grid Structure 

### Interpolation 

## Files 

### boundary_conditions.py
Enforces the boundary conditions at each timestep. The file contains two boundary condition functions, for class FIELDS and FPFE respectively. The current options for the fields boundary conditions include *periodic* and *conducting*. The periodic boundary conditions copy the values from the cell of the last active zones into the ghost cells. Conducting boundary conditions set the electric fields parallel to the axes to be 0. The fpfe boundary conditions (not yet functional) work similarly, with *periodic* and *conducting* options. In the conducting option, the FPFE is lost from the simulation domain and its charge is deposited over the simulation boundary.   
### current_deposition.py
Deposit current on the Yee Grid from the motion of the FPFEs. The current from each FPFE is estimated as its charge times change in position per timestep. The current is then split into 4 parts weighted based on the distance of the FPFE center from each of the 4 sides of the Yee Grid cell. If the FPFE crosses a grid cell within a timestep then the displacement of the FPFE is brokendown into 2 steps such that the FPFE lies within a cell during each of the partial displacements. This is the *zigzag* scheme used in some of the PIC codes.    
### fields_class.py

### fpfe_class.py

### interpolate_fields.py

### particle_push.py

### poisson.py

### run_simulation.py

## References 
* https://static.ias.edu/pitp/2016/sites/pitp/files/
* https://cds.cern.ch/record/2057438/files/1418884\_181-206.pdf
* https://static.ias.edu/pitp/2016/sites/pitp/files/spitkovsky\_slides.pdf
* https://reader.elsevier.com/reader/sd/pii/S0010465503004375?token=60B44E98CFBA28852E0E446C948DB3F10E7775D710652F5AEBB7932AC8312F97A230B07EF7FEB37AE1EECB778DC06776
* https://iopscience.iop.org/article/10.3847/1538-4365/aab114/pdf
