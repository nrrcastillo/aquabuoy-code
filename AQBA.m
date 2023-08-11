% Main program for solving two mass system which models AquaBuOY
% configuration A when modeling with regular or irregular waves 

%Inputs:
% Choose regular waves or irregular waves, reg = 1 for regular waves, or
% reg = 0 for irregular waves. 

% For regular waves need to set period T and amplitude A. Example:
% T = 4;
% A = 1;            

% For irregular waves need average wave period Tz and significant wave height
% Hs. Example:
% Tz = 8.5;
% Hs = 1; 

function AQBA(A,T,Hs,Tz,reg,linear,cd,b2_in,c_in,Cd,data_filename,figure_filename)

% Notation for unknown variables:
% t0 = 0       Initial time of simulation
% y(1) = z1    Displacement of float mass center of gravity
% y(2) = z2    Displacement of piston center of gravity
% y(3) = v1    velocity of float mass center of gravity
% y(4) = v2    velocity of piston center of gravity

%Get the following variables from the Float and Tube Geometry
[M1,M2,c,S,Aw,b2,D,rho,g,Dt_in,Dt_out,Lt,am_extra] = FPG(linear,cd,b2_in,c_in);   

% Get the following variables from the Hyrdodynamic data file HD.m
[T_vec,F_vec,a_vec,b_vec] = HD;

% Interpolate the given vectors T, F, a, b to make smooth functions of period T
% also calculate omegas as functions of the periods
T_i = (T_vec(1):.45:T_vec(end));
omega_i = 2*pi./T_i;
% F_i = interp1(T_vec,F_vec,T_i,'cubic');
% a_i = interp1(T_vec,a_vec,T_i,'cubic');
b_i = interp1(T_vec,b_vec,T_i,'cubic');
        
    % Initial positions and initial velocities: y0 = [z1(t0) z2(t0) v1(t0) v2(t0)]
    y0 = [0 0 0 0];

    % The time span where the equations are solved: tspan = [t0 t_final] in
    % time steps specified, start few periods prior to t = 0 so graph
    % initiating at t = 0 is already the developed solution
    if(reg == 1)
      dtau = .1;  
      T = T; 
      end_time = 100;
      tspan = (-250:dtau:end_time);
      t_Fa = 0;
      Fa = 0;
    else
      dtau = 1; %should be a number so that tspan(1)/dtau is an integer see few lines below.
      T = Tz;
      end_time = 1200;
      tspan = (-250:dtau:end_time);
      [t_Fa Fa] = SeaState(Tz,Hs,S,T_vec,F_vec,end_time,dtau); 
    end
    t0 = round(1-tspan(1)/dtau);
    
    % Get hydrodynamic data for particular T
    F = interp1(T_vec,F_vec,T,'cubic');
    a = interp1(T_vec,a_vec,T,'cubic');
    b = interp1(T_vec,b_vec,T,'cubic');

     % Save initial time and velocity in data file for accessing in function f.m  
     delete time.mat
     delete velocity.mat
     time = tspan(1);
     save time.mat time
     velocity = y0(3);
     save velocity.mat velocity
 
     % Describe the sparsity pattern of the jacobian matrix
     jpat = [0 0 1 0; 0 0 0 1; 1 1 1 1; 1 1 1 1];
     options = odeset('JPattern',jpat);
     % Can adjust the accuracy of the computations by setting the following values   
        options = odeset('RelTol',1e-5,...
                     'AbsTol',[1e-10 1e-10 1e-10 1e-10]);
    [time,y] = ode23tb(@f,tspan,y0,[],b,b_i,omega_i,c,b2,S,A,T,F,M1,M2,a,t_Fa,Fa,reg,linear,Cd,Aw,rho,am_extra);
    
    % Only want solutions after time t=0 once solutions are settled.
    tp = tspan(t0:end);
    X1 = y(t0:end,1);
    X2 = y(t0:end,2);
    V1 = y(t0:end,3);
    V2 = y(t0:end,4);

    % limit the results to only desplay solutions within the excursion
    % limit, for the full size prototype this would be 2m. Denis
    % Leterneau's idea, not used so far
%    limit = 1;
%    if(limit == 1)
%       for i = 1:length(X1),
%         if(X1(i)-X2(i) > 2) 
%           X2(i) = X1(i) - 2;
%           V2(i) = 0;
%         elseif(X1(i)-X2(i) < -2)
%           X2(i) = X1(i) + 2;
%           V2(i) = 0;
%         end
%       end 
%    end

    xamp = X1-X2;
    max_xamp = max(abs(X1-X2));
    rel_V = V1-V2;
 
    % Graph position of float and piston and amplitude difference of them.
    figure(1),axes('position',[.3  .3  .4  .4])
    subplot(3,2,1),plot(tp,X1,'r-','LineWidth',.1); hold on;
    subplot(3,2,1),plot(tp,X2,'b-','LineWidth',.1); hold on;
    subplot(3,2,1),plot(tp,xamp,'g-','LineWidth',.1); hold on;
    xlabel('time t [s]');
    ylabel('Displacement [m]'); 
    if(reg == 1) 
        if(linear == 1)
           title(['                           Prototype scale. T = ', num2str(T),', H = ',num2str(2*A),', b2 = ',num2str(b2),', c = ',num2str(c),', max(xamp) = ', num2str(max_xamp)]);
        else
           title(['                           Prototype scale. T = ', num2str(T),', H = ',num2str(2*A),', F_{fric} = ',num2str(b2),', c = ',num2str(c),', max(xamp) = ', num2str(max_xamp)]);
        end
    else
        if(linear == 1)
           title(['                           Prototype scale. Tz = ', num2str(T),', Hs = ',num2str(Hs),', b2 = ',num2str(b2),', c = ',num2str(c),', max(xamp) = ', num2str(max_xamp)]);
        else
           title(['                           Prototype scale. Tz = ', num2str(T),', Hs = ',num2str(Hs),', F_{fric} = ',num2str(b2),', c = ',num2str(c),', max(xamp) = ', num2str(max_xamp)]);
        end
    end
    legend('Float displacement: z_1','Piston displacement: z_2','Relative displacement: z_1-z_2','Location','BestOutside');
  
    if(linear == 1)
      force = b2*(rel_V) + c*(xamp);
    else
      force = b2*tanh(100*(rel_V)) + c*(xamp);        
    end
     
    ave_force = zeros(1,length(tp));
    for j = 1:length(tp),
      temp = force(1:j);
      ave_force(j) = sum(temp)/length(temp);
    end
    % Look at the graph force
    figure(1), axes('position',[.3  .3  .4  .4]);
    subplot(3,2,3), 
    plot(tp,force/1000,'k-'); axis tight;
    xlabel('time t [s]');

    figure(1), axes('position',[.3  .3  .4  .4]);
    subplot(3,2,5), 
    plot(xamp,force/1000,'g*'); grid on; axis tight;
    xlabel('Extension [m]');
    if(linear == 1)
      ylabel('                                            F_{pto}=b2(v_1-v_2) [kN]'); 
    else
      ylabel('                                            F_{pto} = F_{fric}tanh(100(v_1-v_2)) + c(z_1-z_2) [kN]'); 
    end
            
    power = force.*(rel_V);
    ave_power = zeros(1,length(tp));
    for j = 1:length(tp),
      temp = power(1:j);
      ave_power(j) = sum(temp)/length(temp);
    end
    % Look at the graph power 
    figure(1), axes('position',[.3  .3  .4  .4])
    subplot(3,2,4), 
    plot(tp,power/1000,'k-'); axis tight;
    xlabel('time t [s]');
    ylabel('P_{pto} = F_{pto}(v_1-v_2) [kW]'); 
    
    subplot(3,2,6), plot(tp,ave_power/1000,'k-'); axis tight;
    xlabel('time t [s]');
    ylabel('Average P_{pto} [kW]');
          
    saveas(figure(1),figure_filename) 
    
    fid = fopen(data_filename,'w+');
    fprintf(fid,'INPUTS for Prototype scale model:\n\n');
    if(reg == 1)
      fprintf(fid,'Regular Waves.\n');        
      fprintf(fid,'Wave Height H [m]: %4.2f\n',2*A);
      fprintf(fid,'Wave Period T [s]: %4.2f\n',T);
    else
      fprintf(fid,'Irregular Waves.\n');        
      fprintf(fid,'Significant Wave Height Hs [m]: %4.2f\n',Hs);
      fprintf(fid,'Average Wave Period Tz [s]: %4.2f\n',Tz);
    end
    fprintf(fid,'Float diameter [m]: %4.2f\n',D);
    fprintf(fid,'Tube inner diameter [m]: %4.2f\n',Dt_in);
    fprintf(fid,'Tube outer diameter [m]: %4.2f\n',Dt_out);
    fprintf(fid,'Tube length [m]: %4.2f\n',Lt);
    if(linear == 1)
      fprintf(fid,'cd: %4.2f\n',cd);
    else
      fprintf(fid,'Ffric [N]: %4.2f\n',b2);   
    end
    fprintf(fid,'PTO Spring Stiffness c [N/m]: %4.2f\n\n\n',c);
      
    fprintf(fid,'OUTPUTS:\n\n');
    fprintf(fid,'Peak Float Amplitude [m]: %6.3f\n',max(X1));
    fprintf(fid,'Average Float Amplitude [m]: %6.3f\n',sum(abs(X1))/length(X1));
    fprintf(fid,'Root Mean Square Float Amplitude [m]: %6.3f\n\n',realsqrt(sum(X1.^2)/length(X1)));

    fprintf(fid,'Peak Relative Amplitude |Float - Piston| [m]: %6.3f\n',max_xamp);
    fprintf(fid,'Average Relative Amplitude |Float - Piston| [m]: %6.3f\n',sum(abs(xamp))/length(xamp));
    fprintf(fid,'Root Mean Square Relative Amplitude |Float - Piston| [m]: %6.3f\n\n',realsqrt(sum(xamp.^2)/length(xamp)));
      
    fprintf(fid,'Peak PTO Force [kN]: %4.2f\n',max(force)/1000);
    fprintf(fid,'Average PTO Force [kN]: %4.2f\n',sum(abs(force))/length(force)/1000);
    fprintf(fid,'Root Mean Square PTO Force [kN]: %4.2f\n\n',realsqrt(sum(force.^2)/length(force))/1000);
      
    fprintf(fid,'Peak Power [kW]: %4.2f\n',max(power)/1000);
    fprintf(fid,'Average Power [kW]: %4.2f\n',sum(power)/length(power)/1000);
    fprintf(fid,'Root Mean Sqare Power [kW]: %4.2f\n\n',realsqrt(sum(power.^2)/length(power))/1000);
            
    fprintf(fid,'Peak Velocity of Float [m/s]: %6.3f\n',max(V1));
    fprintf(fid,'Average Velocity of Float [m/s]: %6.3f\n',sum(abs(V1))/length(V1));
    fprintf(fid,'Root Mean Square Velocity of Float [m/s]: %6.3f\n\n',realsqrt(sum(V1.^2)/length(V1)));
      
    fprintf(fid,'Peak Velocity of Piston [m/s]: %6.4f\n',max(V2));
    fprintf(fid,'Average Velocity of Piston [m/s]: %6.4f\n',sum(abs(V2))/length(V2));
    fprintf(fid,'Root Mean Square Velocity of Piston [m/s]: %6.4f\n',realsqrt(sum(V2.^2)/length(V2)));
    fclose(fid);
    

    