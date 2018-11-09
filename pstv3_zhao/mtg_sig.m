function f = mtg_sig(t,k)
% Syntax: f = mtg_sig(t,k)
% 12:37 PM 7/0/98
% defines modulation signal for turbine power reference
global tg_sig n_tg n_tgh
global mac_spd  governor_gain  remove_governor 
global generator_capacity
f=0; %dummy variable
if n_tg~=0|n_tgh~=0
  tg_sig(:,k) = zeros(n_tg+n_tgh,1);
  if t<=0.0
     tg_sig(:,k) = zeros(n_tg+n_tgh,1);
  else
     tg_sig(:,k) = zeros(n_tg+n_tgh,1); 
     if true
        tg_sig(:,k) = max(min(-governor_gain.*(mac_spd(:,k)-1),generator_capacity(:,2)),generator_capacity(:,1));  
%      tg_sig(remove_governor,k) = 0;
     end
     
     if t>1 % power losss
        tg_sig(1,k)=-[1];  %step increase of load at some of the buses
     end 
  end
end
return
