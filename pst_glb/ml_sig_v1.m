function f = ml_sig(t,k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
% function desription:
% modulates the active load at a bus

global lmod_sig n_lmod
global bus_v  OLC_gain  OLC_mod OLCtime disturbance_size disturbance_mod controlled_load OLC_bus  
global OLC_capacity
global GLB_bus
global METHOD
global load_freq bus dc_buses non_dc_buses big_L dc_alpha dc_beta load_bus_index Q_glb L d d_set
if strcmp(METHOD, 'DGLB')
   f=0; % dummy variable
   if n_lmod>0
      timestep=0.01;
      OLCstep = OLCtime/timestep;
      %OLCstep=1;     %step used for discrete time control 
      lmod_sig(:,k)=0;
      
      if t>0.5
        lmod_sig(disturbance_mod,k)=disturbance_size;  %step increase of load at some of the buses
      end
      
      d(:,k) = d_set;
      
      if k>=OLCstep+1
      temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        freq_deviation = freq_deviation/OLCtime/(120*pi);
        load_freq(:, k) = freq_deviation+1; % convert from rad to normalized Hz (1 = 60hz)
      else
        load_freq(:, k) = 1;
      end
      
      if k>=OLCstep+1 
        if k==OLCstep+1
         save('glb_func/debug.mat');
        end
        %% algorithm parts        
        N_bus = length(lmod_sig(:,1)); % TODO: should be the number of all buses. legnth(bus(:,1)).
        % Step 1: Initialize y        
        
        %simp01: alpha = [0.8147 0.9058 0.1270  0.9134 0.6324  0.0975  0.2785  0.5469 0.9575  0.9649 0.1576 0.9706  0.9572 0.4854 0.8003  0.1419 0.4218  0.9157  0.7922];
        %simp01: alpha = 1/2*alpha';
        % basic rules: OLC_gain cannot be too large
        %% Working versions
        % case 1: small w(Q) quadratic + no summing queue + no random alpha
%        L(:,k) = OLC_gain.*(freq_deviation + 0.00065*Q_glb(k-1)); %w(Q) quadratic % + 0.0005*Q_glb(k-1)
%         Q_glb(k) = max(0, 0 - sum(L(:,k))) ; % if nothing queue, sum(L) == 0;
%         controlled_load = 1*L(:,k)+0;
%controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu

        % case 2: w(Q) linear + no random alpha
%          L(:,k) = OLC_gain.*(freq_deviation + 0.005); %w(Q) linear (converges with summing queue)
%          Q_glb(k) = 0 - sum(L(:,k)) + Q_glb(k-1);
%          controlled_load = 1*L(:,k)+0;
%controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu

        % case 3: dropping some workload~ + no random alpha
%         L(:,k) = OLC_gain.*(freq_deviation + 0.001*Q_glb(k-1)); %w(Q) quadratic (does not converge with summng queue)
%          Q_glb(k) = -0.121*big_L - sum(L(:,k)) + Q_glb(k-1); % summing queue & dropping the workload.
%          Q_glb(k) = max(0,Q_glb(k));freq_deviation
%          controlled_load = 1*L(:,k)+0;
%controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu

        % case 4: small w(Q) quadratic + no summing queue + random alpha
%         L(:,k) = OLC_gain.*(alpha.*freq_deviation + 0.0005*Q_glb(k-1)); %w(Q) quadratic (does not converge with summng queue)
%         Q_glb(k) = max(0, 0 - sum(L(:,k))) ; % if nothing queue, sum(L(:,k)) == 0;
%         controlled_load = alpha.*L(:,k)+0;
%controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
        
        % case 5 increase the range of controlled load
%         B = 0.0;
%         A = 0.001; % A=0.001 is not small as Q_glb is the sum.
%         L(:,k) = OLC_gain.*(freq_deviation + A*Q_glb(k-1) + B); 
%         Q_glb(k) = max(0, 0 - sum(L(:,k))) ; % if nothing queue, sum(L) == 0;
%         controlled_load = 1*L(:,k)+0;
%         controlled_load=max(min(controlled_load, 4*OLC_capacity(:,2)),4*OLC_capacity(:,1));    %frequency-pu

        
        %% non-working versions
        % case 1: quite large w(Q) quadratic + no summing queue + no random alpha
%         L(:,k) = OLC_gain.*(freq_deviation + 0.001*Q_glb(k-1)); %w(Q) quadratic (does not converge with summng queue)
%         Q_glb(k) = max(0,0 - sum(L(:,k))) ; % if nothing queue, sum(L(:,k)) == 0;
%         controlled_load = 1*L(:,k)+0;
%controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu

         % case 2: w(Q) quadratic + summing
%          L(:,k) = OLC_gain.*(freq_deviation + 0.0001*Q_glb(k-1)); %w(Q) quadratic (does not converge with summng queue)
%          Q_glb(k) = 0 - sum(L(:,k)) + Q_glb(k-1); % summing queue
%          Q_glb(k) = max(0, Q_glb(k));
%          controlled_load = 1*L(:,k)+0;
%controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
        
        % case 2: w(Q) quadratic + summing + random alpha
%          L(:,k) = OLC_gain.*((ones(size(alpha))./alpha).*freq_deviation + 0.0001*Q_glb(k-1)); %w(Q) quadratic (does not converge with summng queue)
%          Q_glb(k) = 0 - sum(L(:,k)) + Q_glb(k-1); % summing queue
%          Q_glb(k) = max(0,Q_glb(k));
%          controlled_load = alpha.*L(:,k)+0;
%controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu

        
        % case 3
        %L(:,k) = OLC_gain.*(freq_deviation + 0.01*Q_glb(k-1)+0.1); % Josh playing with various functions
        %The above setting converges in Fig 1 but not 2.  It looks like the
        %queue keeps oscillating between 0 and 60 but is a bounded
        %oscillation.
        %Q_glb(k) = 0 - sum(L(:,k)) + Q_glb(k-1); % if we need to consider the previous Q_glb queue, sum(L(:,k)) == 0;
        %Q_glb(k) = max(0,0 - sum(L(:,k)) + Q_glb(k-1)); % max(0,Q_glb) because a negative queue can not be used;
        
%         % small alpha
%          L(:,k) = OLC_gain.*((ones(size(alpha))./alpha).*freq_deviation + 0.0001*Q_glb(k-1)); %w(Q) quadratic (does not converge with summng queue)
%          Q_glb(k) = 0 - sum(L(:,k)) + Q_glb(k-1); % summing queue
%          Q_glb(k) = max(0,Q_glb(k));
%          controlled_load = alpha.*L(:,k)+0;

        %% Trying
        
         % case 1: quite large w(Q) quadratic + no summing queue + no random alpha
%          avg_range = 100;
%          Q_temp = 0;
%          iCount = 0;
%          for i=1:avg_range
%              if (k-i>0)             
%                  Q_temp = Q_temp + Q_glb(k-i);
%                  iCount = iCount+1;
%              end
%          end
%          Q_temp = Q_temp/iCount;
%          
%         L(:,k) = OLC_gain.*(freq_deviation + 0.01*(Q_temp)); %w(Q) quadratic % + 0.0005*Q_glb(k-1)
%         Q_glb(k) = max(0, 0 - sum(L(:,k))) ; % if nothing queue, sum(L) == 0;
%         controlled_load = 1*L(:,k)+0;

      % case 2: cut down the constant in Q
%        idle_power = 3; % sum(beta)
%        W_workload = sum(bus(load_bus_index,6)) - idle_power;       
%        B = 0.0;
%        A = 0.001;       
%        over_load = 2.8;
%        L(:,k) = OLC_gain.*(freq_deviation + A*Q_glb(k-1) + B); 
%        Q_glb(k) = max(0, 0 - sum(L(:,k)) + Q_glb(k-1)) ; % if nothing queue, sum(L) == 0;
%        controlled_load = 1*L(:,k)+0;
%        controlled_load=max(min(controlled_load, 4*OLC_capacity(:,2)),4*OLC_capacity(:,1));    %frequency-pu

      % case 3: smaller OLC_gain
       %simp05: temp = OLC_gain;
       %temp = 83.4*ones(6,1); %simp05: make all a's the same 83 works, 83.4 does not (kappa = 0.002)
       temp = 160*ones(6,1);
       alpha = ones(6,1); %simp01: Make all energy linear functions = 1
       kappa = 0.001; %simp02: Add a control parameter in front of w'(Q)
       A = 1/2; %0.0005/2; %simp04
        %L_opt = temp.*(freq_deviation + 0.0005*Q_glb(k-1)); %w(Q) quadratic % + 0.0005*Q_glb(k-1)
        L_opt = temp.*(alpha.*(load_freq(:,k-1)-1) + kappa*2*A*Q_glb(k-1)); %w(Q) quadratic % + 0.0005*Q_glb(k-1)
        L(:,k) = max(min(L_opt, OLC_capacity(:,2)),OLC_capacity(:,1));
        %Q_glb(k) = max(0, 0 - sum(L(:,k))) ; % if nothing queue, sum(L) == 0;
        Q_glb(k) = 0 - sum(L(:,k)) ; %To match the signal "Q" from the paper
        controlled_load = 1*L(:,k)+0;
        
        %controlled_load=max(min(controlled_load, OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
        
        %%        
        d(:,k) = d_set + controlled_load;
        lmod_sig(OLC_mod,k)= lmod_sig(OLC_mod,k) + controlled_load;
      end
      
   end
elseif strcmp(METHOD, 'DGLB_admm')
   f=0; % dummy variable
   if n_lmod>0
      timestep=0.01;
      OLCstep = OLCtime/timestep;
      %OLCstep=1;     %step used for discrete time control 
      lmod_sig(:,k)=0;
      if t>0.5
        lmod_sig(disturbance_mod,k)=disturbance_size;  %step increase of load at some of the buses
      end
      
      if k>=OLCstep+1
      temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        freq_deviation = freq_deviation/OLCtime/(120*pi);
      end
      
      if k>=OLCstep+1 && mod(k-1,OLCstep)==0
        if k==OLCstep+1
         save('glb_func/debug.mat');
        end
        %% algorithm parts
        omega  = freq_deviation;        
        N_bus = length(lmod_sig(:,1)); % TODO: should be the number of all buses. legnth(bus(:,1)).        
        D=0.001;
        numIter = 10;
        cost_gain = 1;
        % Step 1: Initialize y
        
        y_left      = zeros(N_bus,numIter); % copy for {(j-1,j) @ j}
        y_right      = zeros(N_bus,numIter); % copy for{(j,j+1) @ j}

        
        mu_left = zeros(N_bus,numIter);
        mu_right = zeros(N_bus,numIter);
        
        L          = zeros(N_bus,1);
        d          = zeros(N_bus,1);
        v          = zeros(N_bus,1);
        
        alpha = 1;
        beta = 0;
        p_r_j = 0;
        rho = 0.1;
        
        for j=1:N_bus            
            y_left(j, 1)        = (j-1)/N_bus*big_L;                        
            y_right(j, 1)       = (j)/N_bus*big_L;
        end        
        
        for iCount=2:numIter            
            for j=1:N_bus % TODO: do it in parrallel 
                % Step 2: solve (39)
                if(j>1)
                  [y_j isSuccessful] = solve_cost_eq( omega(j), mu_left(j, iCount-1), rho, ...
                            y_right(j, iCount-1), ...
                            y_right(j-1, iCount-1), ...
                            cost_gain, alpha, beta, 1);                   
                  y_left(j, iCount)    = y_j;
                  y_right(j-1, iCount) = y_j;
                else
                  y_left(N_bus, iCount)  = 0;
                end
                % y_left(j, iCount) = y_j; %??
                
                % Step 3: solve (40)
                if(j<N_bus)
                  [y_j isSuccessful] = solve_cost_eq( omega(j), mu_right(j, iCount-1), rho, ...
                           y_left(j, iCount), y_left(j+1, iCount), ...
                            cost_gain, alpha, beta, -1);                  
                  y_right(j, iCount)  = y_j;                
                  y_left(j+1, iCount) = y_j;
                else
                  y_right(N_bus, iCount)  = big_L;  
                end
                
                % Step 4: Dual update
                if(j>1)
                  mu_left(j, iCount)  = mu_left(j, iCount-1) + rho*(y_right(j-1, iCount) - y_left(j, iCount));
                end
                if(j<N_bus)
                  mu_right(j, iCount) = mu_right(j, iCount-1) + rho*(y_right(j, iCount)  - y_left(j+1, iCount));                
                end
                % Stop condition.
            end
        end    
        % step 5: obtain the final solution:
        for j=1:N_bus
            L(j) = y_right(j,iCount)-y_left(j,iCount);
            d(j) = alpha*L(j)+beta;
        end             
        %controlled_load = (d(j) - bus(load_bus_index,6));
        %controlled_load = 0.00001*randn(size(OLC_mod));
        controlled_load = OLC_gain.*(freq_deviation) - randn(size(OLC_mod));
        %controlled_load = 0*(freq_deviation);
        load_freq(:, k) = freq_deviation +1; % convert from rad to normalized Hz (1 = 60hz)          
        controlled_load=max(min(controlled_load,OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
      end
      lmod_sig(OLC_mod,k)= lmod_sig(OLC_mod,k)+ controlled_load;
   end
elseif strcmp(METHOD, 'GLB')
    %% TODO: only control the data center loads
   f=0; %dummy variable
   if n_lmod>0
   %     if t <30
   %         timestep = 0.01;
   %     else
   %         timestep =0.025;
   %     end;
      timestep=0.01;
      OLCstep = OLCtime/timestep;
      %OLCstep=1;     %step used for discrete time control 
      lmod_sig(:,k)=0;
      if t>0.5
      lmod_sig(disturbance_mod,k)=disturbance_size;  %step increase of load at some of the buses
      end

      if k>=OLCstep+1 && mod(k-1,OLCstep)==0
        temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        controlled_load = OLC_gain.*(freq_deviation);
        load_freq(:, k) = freq_deviation /(2*pi)+1; % convert from rad to normalized Hz (1 = 60hz)  
        controlled_load=max(min(controlled_load,OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
        
         J = 19;
         p = bus(1:J,6); % where can we get these values?
         p(dc_buses) = zeros(1,length(dc_buses)); % use all active load for data centers.
         r = bus(1:J,4);                  

         Alpha = dc_alpha * ones(size(dc_buses));
         Beta = dc_beta * ones(size(dc_buses));

         cost_gain = 1 * ones(J,1);
         D = -0.001*ones(J,1);
         
         omega = load_freq(:, k);
         save('glb_func/debug.mat');
         [isSuccessful, d, L] = centralized_glb(p, r, big_L, D, J, Alpha, Beta, omega, cost_gain, dc_buses, non_dc_buses);
         
         lmod_sig(GLB_bus,k)= lmod_sig(GLB_bus,k) + controlled_load(GLB_bus);
      end
      %lmod_sig(OLC_mod,k)= lmod_sig(OLC_mod,k)+ controlled_load;
      
   end
elseif strcmp(METHOD, 'OLC')
   f=0; %dummy variable
   if n_lmod>0
      timestep=0.01;
      OLCstep = OLCtime/timestep;
      %OLCstep=1;     %step used for discrete time control 
      lmod_sig(:,k)=0;
      if t>0.5
      lmod_sig(disturbance_mod,k)=disturbance_size;  %step increase of load at some of the buses
      end

      if k>=OLCstep+1 && mod(k-1,OLCstep)==0
        temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
        freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
        freq_deviation = freq_deviation/OLCtime/(120*pi);
        controlled_load = OLC_gain.*(freq_deviation);
        load_freq(:, k) = freq_deviation+1; % convert from rad to normalized Hz (1 = 60hz)  
        controlled_load=max(min(controlled_load,OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
      end      
      lmod_sig(OLC_mod,k)= lmod_sig(OLC_mod,k)+ controlled_load;
   end
elseif strcmp(METHOD, 'OLC_few_loads') % do not work well. ==> change size of lmod_sig
    f=0; %dummy variable 
    num_loads = 18;
    if n_lmod>0
       timestep=0.01;
       OLCstep = OLCtime/timestep;
       lmod_sig(:,k)=0;
       if t>0.5
         lmod_sig(disturbance_mod,k)=disturbance_size;  %step increase of load at some of the buses
       end

       if k>=OLCstep+1 && mod(k-1,OLCstep)==0
           temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
           freq_deviation = temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
            (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
            (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi);
            freq_deviation = freq_deviation/OLCtime/(120*pi);
           controlled_load = OLC_gain.*(freq_deviation);
           load_freq(:, k) = freq_deviation +1; % convert from rad to normalized Hz (1 = 60hz)  
           controlled_load=max(min(controlled_load,OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
           controlled_load=max(min(controlled_load,OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
       end
      lmod_sig(1:num_loads,k)= lmod_sig(1:num_loads,k)+ controlled_load(1:num_loads);
    end
else
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