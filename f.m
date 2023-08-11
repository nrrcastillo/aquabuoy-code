% This returns the right hand side of the system of Ordinary Differential
% Equations
function func = f(t,y,b,b_i,omega_i,c,b2,S,A,T,F,M1,M2,a,t_Fa,Fa,reg,linear,Cd,Aw,rho,am_extra)

if(reg == 1)
  % Assign the exciting force for regular waves
  Fw = EF(t,A,T,F);
else
  % Assign the exciting force for irregular waves
  Fw = interp1(t_Fa,Fa,t,'cubic');
end

if(reg == 0)
  % Read velocity data from file
  load('velocity.mat', '-regexp') 

  % add current velocity(time) to velocity vector
  velocity_temp = [velocity; y(3)];
  clear velocity

  % Read time data from file
  load('time.mat', '-regexp') 

  % add current time to time vector
  time_temp = [time; t];
  clear time

  % Remove any of the repetition values in time from time and velocity
  % vectors, by creating new time and new velocity vectors with unique values
  % only
  [new_time, I] = unique(time_temp);
  new_velocity = velocity_temp(I);

  velocity = new_velocity;
  time = new_time;

  if(length(time) > 200)
      time = new_time(end-200:end);
      velocity = new_velocity(end-200:end);
  else
      time = new_time;
      velocity = new_velocity;
  end

  % Write ammended time data to file for next integration point
  save time.mat time
  save velocity.mat velocity

  % For regular waves, along with the retardation function we use
  % the average added mass. If using the Wamit prototype data then
  am = 46000 + am_extra;

  % The retardation function times velocity integrated over time, using the
  % trapezium rule. The last trailing terms are irrelevant as the retardation
  % function goes to zero.
  dt = 0.5;
  int_1toN = 0;
  for n = 1:13,
     if(length(new_time) == 1) % if there is only one point in new_time cannot interpolate
        vel = new_velocity(1);
     else
       t_n = t-n*dt;
       vel = interp1(new_time,new_velocity,t_n,'cubic');
     end
     int_1toN = int_1toN + retardation(b_i,omega_i,n*dt)*vel*dt;
  end

  % The 0th contribution from the retardation function times velocity
  if(length(new_time) == 1) % if there is only one point in new_time cannot interpolate
     vel = new_velocity(1);
  else 
    t_n = t;
    vel = interp1(new_time,new_velocity,t_n,'cubic');
  end
  int_0 = 0.5*retardation(b_i,omega_i,0)*vel*dt;

end

% The fluid friction term:
Fff = 0.5*Cd*Aw*rho*y(3)*abs(y(3));

if(reg == 1)
% If one wants to test the ODEs without the retardation function 
  if(linear == 1)
    func = [y(3);
            y(4);
           (Fw-b*y(3)-S*y(1)-b2*(y(3)-y(4))-c*(y(1)-y(2)) - Fff)/(M1 + a + am_extra);
           (b2*(y(3)-y(4))+c*(y(1)-y(2)))/M2];
  else
    term = b2*tanh(100*(y(3)-y(4)));
    func = [y(3);
            y(4);
           (Fw-b*y(3)-S*y(1)-term-c*(y(1)-y(2)) - Fff)/(M1 + a + am_extra);
           (term+c*(y(1)-y(2)))/M2];
  end
else
  % Right hand side of system of ODEs with non-linear/or linear damping, with retardation function, and inf added mass 
  if(linear == 1)
    func = [y(3);
          y(4);
          (Fw-int_1toN-int_0*y(3)-S*y(1)-b2*(y(3)-y(4))-c*(y(1)-y(2)) - Fff)/(M1 + am);
          (b2*(y(3)-y(4))+c*(y(1)-y(2)))/M2];
  else
    term = b2*tanh(100*(y(3)-y(4)));
    func = [y(3);
            y(4);
           (Fw-int_1toN-int_0*y(3)-S*y(1)-term-c*(y(1)-y(2)) - Fff)/(M1 + am);
           (term+c*(y(1)-y(2)))/M2];
  end
end
