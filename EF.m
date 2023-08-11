% File EF.m containing returns the exciting force Fw, for a given: t,A,T,F
function Fw = EF(t,A,T,F)

omega = 2*pi/T;              % For wave period, corresponding wave frequency
gamma = 0;                   % Arbitrary phase, choose 0
Fw = A*F*cos(omega*t+gamma); % exciting force function
