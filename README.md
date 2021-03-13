# Boundary Layer Flow Thermal Solver
This repository solves heat transfer problems in 2D external flows over a flat plate. The user may specify a variety of flow conditions and thermal boundary conditions. The output includes the wall temperature and heat flux distributions on the top and bottom surfaces of the plate. A description of the main settings is given below.

Setting | Options | Description
------- | ------- | -----------
fluid | 'air', 'water', 'user defined' | Defines the fluid object for property look-up in Cantera.<br>Fluid properties may be specified manually if 'user defined' is chosen.
BC | 'prescribed temp', 'conjugate' | Defines the type of boundary condition used.
wall_conduction | 'thin', 'thick' | Defines the type of discretization used in the solid domain.<br>This setting only applies when BC = 'conjugate'.
transient | true, false | Toggles steady/unsteady analysis.<br>This setting only applies when BC = 'conjugate'.
radiation | true, false | Toggles radiation effects on the outer plate surface.<br>This setting only applies when BC = 'conjugate'.
write | true, false | Toggles option to write data to text file.

## To run this code

You will need MATLAB to run this code as well as a compatible version of [Cantera](https://cantera.org/install/index.html) in order to support fluid property look-ups. Alternatively, you may skip the Cantera installation, but be sure to set the "fluid" variable to 'user defined' or else the code will not run. In this case you will be prompted to enter the fluid properties manually. Once all of the code in the repository have been downloaded to a single directory, the "main.m" script should be executed.

## General

Based on graduate course at Stanford University, ME352C: Convective Heat Transfer.

The conjugate methods implemented in this code are based on the following paper:

Hoffman, D. W., and Eaton, J. K. "Conjugate Heat Transfer Analysis Using the Discrete Green's Function." ASME. J. Heat Transfer. 143(3). (2021)

Author: Davis W. Hoffman (email: dwhoff@stanford.edu)\
Last modified: 03/12/2021\
Developed and tested in MATLAB R2020b
