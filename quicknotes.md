*Variables*
Note to self: should rename variables to differentiate user-inputted vs. constants


**f.m**
time_temp: existing time vector with current velocity(time) [as y(3)] data appended to it 
velocity_temp: existing velocity vector with current time [as t]data appended to it  --> are those datapoints a single one or vectors?
new_time: vector that contains unique values of time_temp vector
I: vector with indices of unique values in time_temp vector
new_velocity: elements from velocity_temp vector in indices from I

**driver.m**
A: is the wave height for regular waves {float}
T: is the wave period for regular waves {float}
Hs: is the significant wave height for irregular waves {float}
Tz: is the average wave period for irregular waves {float}
reg: option for modelling regular and irregular waves {bool}
linear: option for modelling either linear or nonlinear damping {bool}

cd: damping ratio if using linear damping {bool}{dependent on linear == True}
b2_in: friction between the float and the piston for non-linear damping. {float}, b2_in = b2 whenever linear damping {float}{constant}
c_in: spring stiffness {float}
Cd: fluid friction coefficient {float}{change name of this one}

**Integrals**
use scipy.quad for functions
use scipy.trapz for data

**WAMIT**
--Gotta look into this--