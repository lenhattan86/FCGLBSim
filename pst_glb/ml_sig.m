function f = ml_sig(t,k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
% function desription:
% modulates the active load at a bus

global lmod_sig n_lmod
global bus_v  OLC_gain  OLC_mod OLCtime disturbance_size disturbance_mod controlled_load OLC_bus DELAY
global OLC_capacity
global GLB_bus
global METHOD RUNNING_MODE INCIDENT_START
global load_freq bus dc_buses non_dc_buses big_L dc_beta load_bus_index Q_glb L d d_set
global gamma a c Phi mu fcp_lambda fcp_D sumOfD_j fcp_gamma
global g

control_d = zeros(length(OLC_bus),1);
changedLoad = 0*(disturbance_size);

if t>INCIDENT_START
    if RUNNING_MODE == RunningMode.LoadChange    
        % TODO: change the impact time.    
        changedLoad = disturbance_size;            
    end
end

if METHOD==Method.proposed
   f=0; % dummy variable
   if n_lmod>0
      timestep=0.01;
      OLCstep = OLCtime/timestep;
      %OLCstep=1;     %step used for discrete time control 
      lmod_sig(:,k)=0;
      
      if t>INCIDENT_START
        lmod_sig(disturbance_mod,k)=changedLoad;  %step increase of load at some of the buses
      else
        mu(k) = 0;
      end
      
      if k>1
        temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-1));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        freq_deviation = freq_deviation/timestep/(120*pi);
        load_freq(:, k) = freq_deviation + 1; % convert from rad to normalized Hz (1 = 60hz)
      else
        load_freq(:, k) = 1;        
      end
      
      delay_step = DELAY/timestep;
      if k>1 
        if mod(k, OLCstep) == 0
%           if mod(k, 1000) == 0
%             disp('debug');
%           end
          % algorithm parts   
          
          delay_mu = 0;
          if k>delay_step
            delay_mu = mu(k-delay_step);            
          end
          
          control_d = (ones(size(c))./c).*(freq_deviation - a*delay_mu);
          control_d = max(min(control_d,OLC_capacity(:,2)),OLC_capacity(:,1)); 
          
%           if k> (delay_step-OLCstep)
          mu(k+OLCstep) = mu(k) + fcp_lambda * ((sum(a.*control_d)*(fcp_gamma) - mu(k))./(fcp_gamma)); %To fix problems with small gamma          
%           end
          
          %mu(k+OLCstep) = mu(k) + fcp_lambda * (sum(a.*control_d) - mu(k)/(fcp_gamma));           
          controlled_load(:,k) = control_d;
%           lmod_sig(OLC_mod,k) = lmod_sig(OLC_mod,k) + controlled_load(:,k);
        else
          control_d = controlled_load(:,k-1);
          controlled_load(:,k) = control_d;
%           lmod_sig(OLC_mod,k) = lmod_sig(OLC_mod,k) + control_d;     
        end
%         if k>delay_step
%           lmod_sig(OLC_mod,k) = lmod_sig(OLC_mod,k) + controlled_load(:,k-delay_step);
%         end 
          lmod_sig(OLC_mod,k) = lmod_sig(OLC_mod,k) + controlled_load(:,k);
      else
        mu(k+1) = -9999;
      end    
   end
elseif METHOD == Method.optimal
  f=0; % dummy variable
  if n_lmod>0
    timestep=0.01;
    OLCstep = OLCtime/timestep;
    %OLCstep=1;     %step used for discrete time control 
    lmod_sig(:,k)=0;

    if t>INCIDENT_START
      lmod_sig(disturbance_mod,k)=changedLoad;  %step increase of load at some of the buses
    else         
      mu(k) = 0;
    end

    if k>=OLCstep+1
      temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
      freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
          (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
          (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
      freq_deviation = freq_deviation/OLCtime/(120*pi);
      load_freq(:, k) = freq_deviation + 1; % convert from rad to normalized Hz (1 = 60hz)
    else
      load_freq(:, k) = 1;        
    end

    if k>=OLCstep+1 
      N =  length(OLC_bus);      
      % use CVX to compute control_d.
      warning off;
      cvx_begin quiet
        variables d_j(N) w        
        minimize( (fcp_gamma)/2* ((sum(a.*d_j))^2) + ...
              sum((c'/2).*(d_j.^2)) + sumOfD_j*w^2 ...
              )        
        subject to
            d_j >= OLC_capacity(:,1);
            d_j <= OLC_capacity(:,2);
            sum(d_j) + sumOfD_j*w == sum(disturbance_size);
      cvx_end         
      control_d  = d_j;   
      if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
          error('This optimization is unsolved');
      end  
      warning on;
      % bound
      lmod_sig(OLC_mod,k)= lmod_sig(OLC_mod,k) + control_d;
      controlled_load(:,k) = control_d;
    end      

end
elseif METHOD == Method.offline
   f=0; % dummy variable
   if n_lmod>0
      timestep=0.01;
      OLCstep = OLCtime/timestep;
      %OLCstlmod_sigep=1;     %step used for discrete time control 
      lmod_sig(:,k)=0;
      
      if t>INCIDENT_START
        lmod_sig(disturbance_mod,k)=changedLoad;  %step increase of load at some of the buses
      else         
        mu(k) = 0;
      end
      
      if k>=OLCstep+1
        temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-1));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        freq_deviation = freq_deviation/OLCtime/(120*pi);
        load_freq(:, k) = freq_deviation + 1; % convert from rad to normalized Hz (1 = 60hz)
      else
        load_freq(:, k) = 1;        
      end
      
      if k>=OLCstep+1 
          control_d = zeros(length(OLC_bus),1);
        % controlled_load is preloaded
%         if t>INCIDENT_START+1.45 % 1.45 is the best
        if t>INCIDENT_START+timestep*2
            control_d= controlled_load(:,length(controlled_load(1,:)));
        end
        lmod_sig(OLC_mod,k)= lmod_sig(OLC_mod,k) + control_d;
        controlled_load(:,k) = control_d;
      end      
lmod_sig
   end
elseif METHOD == Method.OLC
    f=0; %dummy variable   
   if n_lmod>0
      timestep=0.01;
      OLCstep = OLCtime/timestep;
      delay_step = DELAY/timestep;
      %OLCstep=1;     %step used for discrete time control 
      lmod_sig(:,k)=0;
      if t>INCIDENT_START+timestep
      lmod_sig(disturbance_mod,k)=changedLoad;  %step increase of load at some of the buses
      end

      %if k>=OLCstep+1 && mod(k-1,OLCstep)==0
      if k >1
        temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-1));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);          
        freq_deviation = freq_deviation/timestep/(120*pi);        
        load_freq(:, k) = freq_deviation+1; % convert from rad to normalized Hz (1 = 60hz) 
        
        if mod(k, OLCstep) == 0          
          control_d = (ones(size(c))./c).*(freq_deviation);                   
          control_d=max(min(control_d,OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu                
          controlled_load(:,k) = control_d;
        else
          control_d=controlled_load(:,k-1);    %frequency-pu                
          controlled_load(:,k) = control_d;
        end
        if k>delay_step
          lmod_sig(OLC_mod,k) = lmod_sig(OLC_mod,k) + controlled_load(:,k-delay_step);
        end 
      end      
       
   end
elseif METHOD == Method.NONE
    f=0; %dummy variable
    lmod_sig(:,k)=zeros(n_lmod,1);
    
    % added the following code to measure the bus frequency of the default
    % system
    timestep=0.01;
    OLCstep = OLCtime/timestep;
    
    lmod_sig(:,k)=0;
    
    if t>INCIDENT_START
      lmod_sig(disturbance_mod,k)=changedLoad;  %step increase of load at some of the buses
    end
      
    if k>=OLCstep+1 && mod(k-1,OLCstep)==0
        temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
         (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
         (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        freq_deviation = freq_deviation/OLCtime/(120*pi);
        load_freq(:, k) = freq_deviation +1; % convert from rad to normalized Hz (1 = 60hz)  
    end
else % default
    f=0; %dummy variable
    lmod_sig(:,k)=zeros(n_lmod,1);
%     if n_lmod~=0
%        lmod_sig(:,k)=zeros(n_lmod,1);
%         lmod_sig(:,k) = 0.1*randn(n_lmod,1);
%        if t>=0.0
%           lmod_sig(1,k)=.1;
%        elseif t>0.2
%            lmod_sig(:,k)=zeros(n_lmod,1);
%        end
%     end
    % added the following code to measure the bus frequency of the default
    % system
    timestep=0.01;
    OLCstep = OLCtime/timestep;
    if k>=OLCstep+1 && mod(k-1,OLCstep)==0
        temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
         (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
         (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        freq_deviation = freq_deviation/OLCtime/(120*pi);
        load_freq(:, k) = freq_deviation +1; % convert from rad to normalized Hz (1 = 60hz)  
    end
end
%****************************************************************