% File retardation.m to calculate the retardation function integral
function hw = retardation(b,omega,t)

temp(1:length(b)-1) = 0.5*(b(2:end)+b(1:end-1)).*cos(.5*(omega(2:end)+omega(1:end-1))*t).*(omega(1:end-1)-omega(2:end));

hw = (2/pi)*sum(temp);