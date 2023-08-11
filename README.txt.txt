UPDATE RULES:

Each time an update takes place they need to be noted in this document. First the date
the update took place is noted. Then one should write the files and corresponding lines 
in the last version of file/code, and a summary of the changes. Also each time an update
has been made to the code, a comment should be placed above the lines added/changed/deleted
stating the date and name of the person who made the change. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UPDATES:
Version 1 Last updated May 1, 2007 by Abigail Wacher.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EXPLANATION OF WHAT IS INSIDE EACH FILE:

->FPG.m
  This is the file containing variables from the Float and Piston/Pipe Geomerty: The values
  for the following are assigned in this file:
  Density of water, gravity, diameter of float, water plain area, float draft, hydrostatic
  stiffness, length of tube, diameters (inner and outer) of tube, mass of float, mass of tube,
  friction or linear damping and spring stiffness. 

->HD.m
This file contains variables for the Hydrodynamic Data. The data in this file is from the software 
Wamit, provided by Dominique Roddier. Currently the hydrodynamic data (in "HD.m") and geometry 
(in "FPG.m") corresponding to the 1:50 HMRC model scaled up to Prototype scale.

->EF.m
This file is called for regular waves, it returns the exciting force Fw, for a given: t,A,T,F
for regular waves. 

->retardation.m
This file calculates the retardation function integral called for irregular waves

->SeatState.m
Called for irregular waves returns the Seat State S(f), for a given: Tz (average wave period),
Hs (significant wave height). This is equivalent to the file EF.m for regular waves, it returns 
the exciting force for irregular waves as well as the times the
exciting force is calculated.

->f.m
This is the heart of the code, where the equations are implemented. This returns the right hand 
side of the system of Ordinary Differential Equations to be integrated by an ode solver in the
main program AQBA.m. For irregular waves this function also produces a time.mat and velocity.mat
files for tracking the solutions at intermediate times for caculating the retardation function. 
These files are only for the simulation, it reads and writes to these files during the simulations. 

->AQBA.m
This is the main program where the system of odes is integrated and results are output in 
graphic form and a data file for the given inputs:
 A: is the wave height for regular waves
 T: is the wave period for regular waves
 Hs: is the significant wave height for irregular waves
 Tz: is the average wave period for irregular waves
 reg is parameter to choose if modeling regular or irregular waves, reg = 1 results in regular waves, 
      otherwise the model is for irregular waves.
 linear is set to one if the linear damper equations are modeled, nonlinear if the nonlinear damping is modeled
 cd is the damping ratio if linear damping is modeled
 b2_in = Ffric friction between the float and the piston for non-linear damping.
    or b2_in = b2 linear damping.
 c_in is the spring stiffness
 Cd is the fluid friction coefficient

 Also input into AQBA.m one inputs a data_filename and figure_filename to 
 produces the files of data and graph for the corresponding simulation.

->driver.m
This file runs the main program, in this sample file you can define your inputs (as above under AQBA.m)
and call the program AQBA.m. The input file names for the outputs of the currently set up simulation 
are "Data.txt" and "Plots.txt" found as example outputs in this directory already.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RUNNING SIMULATIONS

To run the driver all one has to do is type the word driver at the Matlab prompt.
Changing the driver file does not constitute a change to the code. The driver file
drives the program, and the one contained in this directory is a sample of a driver
file. The files 'Data.txt' and 'Plots.txt' found in this directory are sample outputs
of what is produced when the sample driver in this directory is used.

One can change the inputs in the driver file to specification for different wave 
periods and frictions, stiffness, name of output files etc. to get their corresponding 
Data and Plot files. If the geometry/hydrodynamic data changes one needs to update HD.m and 
FPG.m accordingly. Also the extra mass added to the added mass in the file "FPG.m" needs to 
be updated accordingly to correspond to the geometry and hydrodynamic data. If there is no 
extra added mass, simply assign a value 0 to the variable am_extra in the file FPG.m
