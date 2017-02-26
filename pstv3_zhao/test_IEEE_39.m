% load data for OFC

generator_capacity_factor = 0.1;    % ratio with setppoint load

generator_bus = bus(abs(bus(:,4))>0.001,1); 
generator_setpoint = bus(generator_bus, 4);
generator_base = mac_con(:,3); 
generator_setpoint_macpu = basmva*generator_setpoint./generator_base;

governor_gain = 25*generator_setpoint_macpu;
generator_capacity = [-generator_setpoint_macpu*generator_capacity_factor, generator_capacity_factor*generator_setpoint_macpu];

 remove_governor = [2;4;6;8;10];
% remaining_governor = [1; 3; 5; 7;];
% remove_governor = [];
%remaining_governor = transpose(1:1:10);

governor_gain(remove_governor)=0;


if (~isempty(remove_governor))
remove_ratio = sum(generator_setpoint(remove_governor))/sum(bus(generator_bus,4));   % with respect to loads
OLC_capacity_factor = remove_ratio*generator_capacity_factor;
else
OLC_capacity_factor = 0;
end;


load_bus_index = bus(abs(bus(:,6))>0.0001,1); 
number_load_bus = length(load_bus_index);

disturbance_size = [1; 1 ; 1];         % in pu
disturbance_bus = [ 4; 15; 16];     % buses
disturbance_mod = [2; 6; 7];     % modules


load_con = zeros(number_load_bus,5);
load_con(:,1) = load_bus_index;
load_con(:,2:5) = repmat([0 0 0.5 0],number_load_bus,1);
% sample: load_con = [...
   % 1  0 0 .5 0;
  % ...]


lmod_con = zeros(number_load_bus,7);
lmod_con(:,1) = transpose(1:1:number_load_bus);
lmod_con(:,2) = load_bus_index;
lmod_con(:,3:7) = repmat([100 3 -3 1 0.5],number_load_bus,1);


OLC_mod=lmod_con(:,1);         % all loads can do OLC
%OLC_mod(disturbance_mod)=[];
OLC_bus = lmod_con(OLC_mod,2);



controlled_load = zeros(length(OLC_mod),1);
OLCtime=0.01;

load_setpoint = bus(load_bus_index,6);
OLC_setpoint = load_setpoint(OLC_mod);
OLC_capacity = [-OLC_setpoint*OLC_capacity_factor, OLC_setpoint*OLC_capacity_factor];  % OLC_load_number * 2 vector 


OLC_gain = 25* OLC_setpoint ; % uniform cost function 1/2*k*p^2, but for load need to rescale by base power




% % Add governor model for OFC simulation (Zhao, Low, CDC2014)
tg_con = [... 
1  1  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
1  2  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
1  3  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
1  4  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
1  5  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
1  6  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
1  7  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
1  8  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
1  9  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0; 
1  10  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0;
]; 

tg_con(:,4) = 0;    %redesign governor control at mtg_sig.m
tg_con(:,5) = 2;   
%tg_con(:,5) = generator_setpoint_macpu*(1+generator_capacity_factor);
% if (~isempty(remove_governor))
%     tg_con(remove_governor,4) = 0;
% end;

