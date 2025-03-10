function f = mtg_sig(t,k)
% Syntax: f = mtg_sig(t,k)
% 12:37 PM 7/0/98
% defines modulation signal for turbine power reference
global tg_sig n_tg n_tgh
global mac_spd  governor_gain  remove_governor 
global generator_capacity
global mac_con 
global METHOD RUNNING_MODE INCIDENT_START METHOD
global bus_v  OLC_gain  OLC_mod OLCtime disturbance_size disturbance_gen_mod mac_scale control_gen_mod

changedLoad = 0*(disturbance_size);
if t>= INCIDENT_START
    if RUNNING_MODE == RunningMode.GenLoss
        changedLoad = disturbance_size;     
    elseif RUNNING_MODE == RunningMode.RenGen
            % TODO: add the renewable generation in 
    end
end

%% the below code is the original.
f=0; %dummy variable
if n_tg~=0|n_tgh~=0
    tg_sig(:,k) = zeros(n_tg+n_tgh,1);
    if t<=0.0
      tg_sig(:,k) = zeros(n_tg+n_tgh,1);
    else
      tg_sig(:,k) = zeros(n_tg+n_tgh,1);
      %tg_sig(1,k) = -1.0*t;
      %tg_sig(1,k) = -0.01;
      if t>= INCIDENT_START
          % power lost on disturbance_gen_mod
          tg_sig(disturbance_gen_mod,k) = -changedLoad'./mac_scale(disturbance_gen_mod);  
      end
      
      % Zhao used the following code in his paper: control on control_gen_mod
      % Droop control (refer tg_con matrix) is already enabled we don't need the following code.
%       if  METHOD==Method.proposed_dc || METHOD==Method.dc || METHOD==Method.OLC
%           tg_sig(control_gen_mod,k) = ...
%               max(min(-governor_gain(control_gen_mod).*(mac_spd(control_gen_mod,k)-1), ...
%                     generator_capacity(control_gen_mod,2)),generator_capacity(control_gen_mod,1));  
%       end    
    end
end
return
