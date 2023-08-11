% File FPG.m is the file containing variables from the Float and
% Piston/Pipe Geomerty
function [M1,M2,c,S,Aw,b2,D,rho,g,Dt_in,Dt_out,Lt,am_extra] = FPG(linear,cd,b2_in,c_in);

rho = 1000;            % Water density (1015 in sea, 1000 in HMRC model)

g = 9.81;              % Gravitational constant

D = 7;                 % Float diameter (Full scale D = 7m)

Aw = (pi/4)*D^2;       % Float water plane Area 

hf = 4.5;              % Float draft

S  = Aw*rho*g;         % Hydrodynamic stiffness parameter
 
Lt = 30;               % Tube length (Full scale Lt = 30m)

Dt_in = 4.7;           % Tube inner diameter

Dt_out = 5;            % Tube outer diameter

% Mass of float
%M1_dry = Aw*hf*rho;   % Dry weight of buoy    
%am_extra = (Dt_out^2-Dt_in^2)*(pi/4)*Lt*rho; % Extra added mass 
%M1 = M1_dry + am_extra;

am_extra = 149000;     % Extra added mass, this value is for HMRC tests
M1 = 237000;           % M1 corresponding to HMRC tests

%M2 = (pi/4)*Lt*rho*Dt_in^2;   % Mass of pipe
M2 = 363000;

if(linear == 1)
                         % cd is the linear damping ratio (optimizing parameter)  
  b2 = cd*sqrt(M1*S);    % Linear damping of relative motion between float and piston
  c = 0;                 % PTO spring stiffness between float and piston
else
  b2 = b2_in;            % Friction between float and piston:  b2 = Ffric
  c =  c_in;             % PTO spring stiffness between float and piston
end


