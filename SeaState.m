% File SeaState.m containing returns the Seat State S(f), for a given:
% Tz (average wave period),Hs (significant wave height).
function [t_Fa Fa] = SeaState(Tz,Hs,S,T_vec,F_vec,end_time,dtau)

Tp = 1.4*Tz;           % Peak wave period.   
fp = 1/Tp;              % Corresponding peak wave frequency
As = ((5*Hs^2)/16)*fp^4;
Bs = (5/4)*fp^4;

% If one wants to set state for random number generation to be different each time, using
% clock, uncomment the following line.
%rand('state',sum(100*clock));

N = 40; 

% need to go from 1 to N+1 instead of 0 to N because Matlab stores first element of array at 1 position 
i = (0:1:N)';
  
% each component of r is random number between 0 and 1.
r = rand(N+1,1);

% frequency step
df = 0.01;

% frequency components
f = 0.5*df*r + i*df;

% period components
T = 1./f;
  
% PM spectrum
Sf = (As*T.^5).*exp(-Bs*T.^4);

% wave amplitude of each component
a = realsqrt(2*Sf*df);

% random wave phase for each component
fi = 2*pi*rand(N+1,1);

% interpolate F on the periods
fF = interp1(T_vec,F_vec,T,'cubic');
  
t_Fa = (-250:dtau:end_time);
Fa = zeros(1,length(t_Fa));
% Wave exciting force as a function of time
for j = 1:length(t_Fa),
  Fa(j) = sum(a(1:end).*sin(2*pi*f(1:end)*(j-1) + fi(1:end)).*fF(1:end));
end
