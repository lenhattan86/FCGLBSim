global RUNNING_MODE METHOD INCIDENT_START
global gamma a c Phi mu fcp_lambda fcp_alpha sumOfD_j disturbance_gen_mod mac_scale

SYS_FREQ=60;
DEBUG = false;

data = [true, false, false];
if data(1)
    dfile = 'datane_glb.m';% New England
    pmech_max = 0.2500;
%     disturbance_size = 1*[1, 1 , 1];    sumOfD_j =  0.0;     % in pu
%     disturbance_bus = [ 4, 8, 20];     % buses % for datane_glb
    disturbance_size = 0.15*[1]; sumOfD_j = 0.0;     % in pu
%     disturbance_bus = [4];     % buses % for datane_glb    
    disturbance_bus = [39];     % buses % for datane_glb
    num_generators = 10;
    disturbance_gen_bus = [39];
    remove_governor = [1:9];
elseif data(2)
    dfile= 'data16m_glb.m'; 
    pmech_max = 0.8333;
    disturbance_size = 60*[1, 1 , 1];         % in pu
    disturbance_bus = [37, 42, 52];     % buses for data16m_glb 
    num_generators  = 16;
elseif data(3)
    dfile= 'datanp48_glb.m'; % % NY & 
    pmech_max = 0.8333;
    disturbance_size = 20*[1, 1, 1];         % in pu
    disturbance_bus = [15, 30, 59];     % buses for data16m_glb 
    num_generators = 48;       
else    
   
end

lfile =length(dfile);
dfile = lower(dfile(1:lfile-2));
eval(dfile);

BASE_POWER = 100; % 100 MVA
DC_CAPACITY = 30; % 50 MW % datane_glb
DC_DEMAND = 25; % 50 MW % datane_glb
DC_CAPACITY = DC_CAPACITY/BASE_POWER;
DC_DEMAND = DC_DEMAND/BASE_POWER;
DC_num = 6;

%% 

% ibus_con = zeros(length(mac_con(:,1)),1);
%ibus_con(10) = 1;

generator_capacity_factor = 0.1;    % ratio with setppoint load

%generator_bus = bus(abs(bus(:,4))>0.0001,1); 
generator_bus = mac_con(:,2);
generator_setpoint = bus(generator_bus(1:num_generators), 4);
generator_base = mac_con(1:num_generators,3); 
generator_setpoint_macpu = BASE_POWER*generator_setpoint./generator_base;
mac_scale = generator_base/BASE_POWER;
% generator_setpoint_macpu = generator_setpoint;

% governor_gain = 50*generator_setpoint_macpu;
% generator_capacity = [-generator_setpoint_macpu*generator_capacity_factor, generator_capacity_factor*generator_setpoint_macpu];


% remaining_governor = [1; 3; 5; 7;];
% remove_governor = [];
%remaining_governor = transpose(1:1:10);

% governor_gain(remove_governor)=0;


% if (~isempty(remove_governor))
%    remove_ratio = sum(generator_setpoint(remove_governor))/sum(bus(generator_bus,4));   % with respect to loads
%    OLC_capacity_factor = remove_ratio*generator_capacity_factor;
% else
%    OLC_capacity_factor = 0;
% end;

OLC_capacity_factor = FLEX;

load_bus_index = bus(abs(bus(:,6))>0,1);    

gen_bus_index = bus(abs(bus(:,4))>0,1);    

% OLC_bus = bus(abs(bus(:,6))>0,1);
temp = bus(abs(bus(:,6))>DC_CAPACITY,1);
OLC_bus = temp(1:DC_num);

number_load_bus = length(load_bus_index);


disturbance_mod = zeros(size(disturbance_bus));


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


%OLC_mod=lmod_con(:,1);         % all loads can do OLC
OLC_mod = zeros(size(OLC_bus));
iCount = 1;
for i=1:length(load_bus_index)
   for j=1:length(OLC_bus)
      if OLC_bus(j) == load_bus_index(i)
         OLC_mod(j) = i;
      end
   end
   for j=1:length(disturbance_bus)
      if disturbance_bus(j) == load_bus_index(i)
         disturbance_mod(j) = i;
      end
   end
end

for i=1:length(gen_bus_index)
  for j=1:length(disturbance_gen_bus)
    if disturbance_gen_bus(j) == gen_bus_index(i)
       disturbance_gen_mod(j) = i;
    end
  end
end
%OLC_mod(disturbance_mod)=[];
%OLC_bus = lmod_con(OLC_mod,2);

controlled_load = zeros(length(OLC_bus),1);
if METHOD == Method.offline
    load(['output/' char(RUNNING_MODE) '_' char(Method.proposed) '_' num2str(FLEX) '_' num2str(DELAY) '_' num2str(WEIGHT) '.mat'], 'controlled_load');
    if(length(controlled_load(:,1))~=length(OLC_bus))
        error('load the wrong file');
    end
end

OLCtime=DELAY;

OLC_setpoint = bus(OLC_bus,6);

DC_RANGE = [DC_CAPACITY*(1 - OLC_capacity_factor), DC_CAPACITY];

OLC_capacity = DC_RANGE - DC_DEMAND;


OLC_gain = 25* OLC_setpoint ; % uniform cost function 1/2*k*p^2, but for load need to rescale by base power
% OLC_gain = 10*ones(size(OLC_setpoint));

 % Add governor model 
% governor model
% tg_con matrix format
%column data unit 
% 1 turbine model number (=1)
% 2 machine number
% 3 speed set point wf pu
% 4 steady state gain 1/R pu
% 5 maximum power order Tmax pu on generator base
% 6 servo time constant Ts sec
% 7 governor time constant Tc sec
% 8 transient gain time constant T3 sec
% 9 HP section time constant T4 sec
% 10 reheater time constant T5 sec


tg_con = ones(num_generators,1)* [ 1  1  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0]; 
tg_con(:,4) = 0;    %redesign governor control at mtg_sig.m
tg_con(:,5) = generator_setpoint_macpu*(1+generator_capacity_factor);
% tg_con(:,5) = generator_setpoint_macpu;
if (~isempty(remove_governor))
   tg_con(remove_governor,4) = 0;
end;

Q_glb = zeros(1,1);
N = length(OLC_bus);
non_dc_buses = [];
dc_alpha = 1;
dc_beta = 0.2;
big_L = sum(OLC_setpoint); % assuming workload are equivalent to power
load_freq = ones(length(OLC_bus),1);
L = zeros(length(OLC_bus),1);
d = zeros(length(OLC_bus),1);


%% control parameters

fcp_lambda  = 0.00001;
% fcp_alpha = 0.001;
fcp_alpha = 0.005;
mu = zeros(1,1);


c = ones(size(OLC_gain))./OLC_gain; c = c';

%     a = ones(size(c)); a = 0.35*2.*rand(size(c)) + 0;
% 19
% a = [1.8492    0.4475    0.7471    0.1750    1.2802    0.3612    0.0901    1.4463    0.6949    1.3212    0.7677    1.2547    0.0433    1.8211    1.6011 1.4917    1.6262    0.7666    1.2346];
% 6
a = [0.15 0.2 0.25 0.4 0.45 0.7];
% Phi
temp1 =  gamma .* (sum(a./c));
temp2 =  gamma .* (sum(a.^2./c));
Phi = temp1./(1+temp2);    

% scale up/down the costs
fcp_alpha = fcp_alpha/WEIGHT;
c = c/WEIGHT;

%% linear model
global g;
randTemp = [0.1270  0.8147    0.9058  0.9134    0.6324    0.0975];
randTemp = randTemp*DC_num/(sum(randTemp));
% g = - randTemp;
% g = randTemp;
g = -0.5*randTemp;


%% 

fprintf('===========================================\n');
fprintf('===========================================\n');
fprintf('Number of buses in the power network: %d\n',length(bus(:,1)));
fprintf('Number of data center buses: %d\n',length(OLC_bus));
fprintf('Frequency Base: %d (Hz) \n',SYS_FREQ);
fprintf('Power Base: %d (MVA) \n',BASE_POWER);
fprintf('===========================================\n');
fprintf('===========================================\n');


%%
OUTPUT = '';