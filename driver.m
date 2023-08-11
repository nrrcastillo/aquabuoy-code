% A: is the wave height for regular waves
% T: is the wave period for regular waves
% Hs: is the significant wave height for irregular waves
% Tz: is the average wave period for irregular waves
% reg is parameter to choose if modeling regular or irregular waves, 
%     reg = 1 results in regular waves, otherwise the model is for 
%     irregular waves. 
% linear is set to one if the linear damper equations are modeled, 
%     nonlinear if the nonlinear damping is modeled
% cd is the damping ratio if linear damping is modeled
% b2_in = Ffric friction between the float and the piston for
%     non-linear damping.
% OR b2_in = b2 linear damping.
% c_in is the spring stiffness
% Cd is the fluid friction coefficient
 
clear all
close all

reg = 1;
H = 2.5;
A = H/2;            
T = 8.1;
Hs = H;   % doesn't matter for regular waves
Tz = T;   % doesn't matter for regular waves
linear = 0;
cd = 1;
% This call to FPG.m is only called to get the hydrostatic stiffness, hence
% the inputs to this call do not matter.
[M1,M2,c,S,Aw,b2,D,rho,g,Dt_in,Dt_out,Lt,am_extra] = FPG(1,1,1,1);
b2_in = 250000;
c_in = .5*S;
Cd = 3;
AQBA(A,T,Hs,Tz,reg,linear,cd,b2_in,c_in,Cd,'Data.txt','Plots.pdf')

clear all

