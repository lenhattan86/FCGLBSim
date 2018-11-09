function f = ml_sig(t,k)
% Syntax: f = ml_sig(t,k)
%4:40 PM 15/08/97
% defines modulation signal for lmod control
% function desription:
% modulates the active load at a bus

global lmod_sig n_lmod
global bus_v  OLC_gain  OLC_mod OLCtime disturbance_size disturbance_mod controlled_load OLC_bus  
global OLC_capacity
f=0; %dummy variable
%lmod_sig(:,k)=zeros(n_lmod,1);   %to return to original, uncomment this line
%if n_lmod~=0
%    lmod_sig(:,k)=zeros(n_lmod,1);
    %lmod_sig(:,k) = 0.1*randn(n_lmod,1);
%    if t>=0.0
%       lmod_sig(1,k)=.1;
%    elseif t>0.2
%        lmod_sig(:,k)=zeros(n_lmod,1);
%    end

%***********Load modulation by Changhong, just a test
% 
% if n_lmod~=0
%    if t<=0
%       lmod_sig(:,k)=zeros(n_lmod,1);
%    else
%       lmod_sig(:,k)=zeros(n_lmod,1);
%       lmod_sig(1,k)=0.5;
%    end
% end;

%***********

%**********Load modulation by Changhong, used for OLC************
%**********Feb 4 Monday 2013
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
   if false
       if t>0.5
         lmod_sig(disturbance_mod,k)=disturbance_size;  %step increase of load at some of the buses
       end
   % load control:
   end
   if false
       if k>=OLCstep+1 && mod(k-1,OLCstep)==0
           temp_delta_theta=angle(bus_v(OLC_bus,k))-angle(bus_v(OLC_bus,k-OLCstep));
           controlled_load=max(min(OLC_gain.*((temp_delta_theta.*(abs(temp_delta_theta)<=1.9*pi)+...
               (temp_delta_theta-2*pi).*(temp_delta_theta>1.9*pi)+...
               (temp_delta_theta+2*pi).*(temp_delta_theta<-1.9*pi))/OLCtime/(120*pi)),OLC_capacity(:,2)),OLC_capacity(:,1));    %frequency-pu
       end
      lmod_sig(OLC_mod,k)= lmod_sig(OLC_mod,k)+ controlled_load;
   end
end;
%****************************************************************