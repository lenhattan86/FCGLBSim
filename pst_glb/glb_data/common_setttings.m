global RUNNING_MODE METHOD INCIDENT_START
global gamma a c Phi mu fcp_lambda fcp_alpha sumOfD_j disturbance_gen_mod mac_scale fcp_beta fcp_gamma NEW_ENG_BASE control_gen_mod
global control_d TIME_STEP DELAY POWER_LOSS 


SYS_FREQ=60;
DEBUG = false;

data = [true, false, false];
if data(1)
    dfile = 'datane_glb.m';% New England
%     dfile = 'datane.m';% New England
    pmech_max = 0.2500;
%     disturbance_size = 1*[1, 1 , 1];    sumOfD_j =  0.0;     % in pu
%     disturbance_bus = [ 4, 8, 20];     % buses % for datane_glb
     if(~exist('SINGLE_RUN','var') || SINGLE_RUN)
         POWER_LOSS = 400;        
     end
     disturbance_size = POWER_LOSS*[1]; sumOfD_j = 0.0;     % 1 GW
     
%    disturbance_size = 50*[1]; sumOfD_j = 0.0;     % 50 MW in CDC and TCNS Submission #1
%     disturbance_bus = [4];     % buses % for datane_glb    
    disturbance_bus = [];     % buses % for datane_glb
    num_generators  = 10; %10
    disturbance_gen_bus = [39]; %39
    control_gen_bus = [30:38];
elseif data(2)
    dfile= 'data16m_glb.m'; 
    pmech_max = 0.8333;
    disturbance_size = 60*[1, 1 , 1];         % in pu
    disturbance_bus  = [37, 42, 52];     % buses for data16m_glb 
    num_generators   = 16;
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

% BASE_POWER = 14000/sum(bus(:,6)); % MW
NEW_ENG_BASE = 14000/sum(bus(:,6)); %14 GW scaled down as Ratio of total power of simulator
BASE_POWER = 100; % in MW
DC_SCALEUP = 1;
DC_CAPACITY = 30*DC_SCALEUP; % 50 MW % datane_glb
DC_CAPACITY = DC_CAPACITY/NEW_ENG_BASE*BASE_POWER; %scale down DC to simulator total power
DC_DEMAND = 25*DC_SCALEUP; % 50 MW % datane_glb 
DC_DEMAND = DC_DEMAND/NEW_ENG_BASE*BASE_POWER; %scale down DC to simulator total power
DC_CAPACITY = DC_CAPACITY/BASE_POWER; %scale down to pu
DC_DEMAND = DC_DEMAND/BASE_POWER; %scale down to pu
DC_num = 10;
disturbance_size = disturbance_size/NEW_ENG_BASE*BASE_POWER;

%%


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

governor_gain = 25*generator_setpoint_macpu;
generator_capacity = [-generator_setpoint_macpu*generator_capacity_factor, generator_capacity_factor*generator_setpoint_macpu];


% remaining_governor = [1; 3; 5; 7;];
remove_governor = [];
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
  for j=1:length(control_gen_bus)
    if control_gen_bus(j) == gen_bus_index(i)
       control_gen_mod(j) = i;
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

% OLCtime = DELAY;
OLCtime = TIME_STEP;

OLC_setpoint = bus(OLC_bus,6);

DC_RANGE = [DC_DEMAND*(1- OLC_capacity_factor), DC_CAPACITY];

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
% 4 steady state gain 1/R pu %% TODO: reduce impact from Power generator.
% 5 maximum power order Tmax pu on generator base
% 6 servo time constant Ts sec
% 7 governor time constant Tc sec
% 8 transient gain time constant T3 sec
% 9 HP section time constant T4 sec
% 10 reheater time constant T5 sec


% tg_con = ones(num_generators,1)* [ 1  1  1  25.0  1.0  0.1  0.5 0.0 1.25 5.0]; % old
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

tg_con(:,4) = 20*ones(size(tg_con(:,4)));

% tg_con(:,4) = 0;    %redesign governor control at mtg_sig.m
tg_con(:,5) = generator_setpoint_macpu*(1+generator_capacity_factor);
% if (~isempty(remove_governor))
%    tg_con(remove_governor,4) = 0;
% end;

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

GAMMA = GAMMA*(NEW_ENG_BASE/BASE_POWER)^2; %scale from 14 GW to total power of simulator
GAMMA = GAMMA*BASE_POWER^2; %scale down to pu
WEIGHT = WEIGHT*NEW_ENG_BASE/BASE_POWER; %scale from 14 GW to total power of simulator
WEIGHT = WEIGHT*BASE_POWER*SYS_FREQ; %Scaling down to pu for power and frequency
GAMMA = GAMMA/WEIGHT; %cost in terms of freq dev.

fcp_gamma  = GAMMA;
fcp_beta  = WEIGHT;
mu = zeros(1,END_TIME*100);
% mu = -2*(10^-3)*ones(1,END_TIME*100);

% c is eta, 0.11*2.*rand(1,DC_num) + 0
% mean(c) =  0.11/POWER_BASE^2
c = [ 0.0877    0.0644    0.1763    0.0762    0.0583    0.1124    0.0807    0.1627    0.1154    0.1770]';
%c = [ 0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11    0.11]';
% c = c*(BASE_POWER^2);
% c  = c/ 10;
% 1/mean(c)
c = c*(NEW_ENG_BASE/BASE_POWER)^2;
c = c*BASE_POWER^2;

% aRange = [1/2 1/1.1];
% aMean = 1/1.8;
% a = (aMean-mean(aRange))/3*randn(1,DC_num) + aMean
% a = [0.5724    0.5293    0.5031    0.5548    0.7054    0.5261    0.6548    0.5327    0.5655    0.5141]'; 
%a = [ 1/1.2  1/1.65  1/1.75  1/1.8  1/1.91  1/1.92  1/1.93  1/1.97  1/1.98  1/2]';
a = [ 1/1.1  1/1.5  1/1.7  1/1.8  1/1.9  1/1.95  1/1.95  1/2 1/2  1/2.1]';
% a = a*(BASE_POWER^2);
% a = a/10;

% scale up/down the costs
%a = a/fcp_beta;
c = c/fcp_beta; %cost in terms of freq dev

%control_d = zeros(length(OLC_bus),1);

%%
if RUNNING_MODE == RunningMode.GenLoss && METHOD ~= Method.optimal
%     disturbance_size = disturbance_size/mac_con(1,3);    
    disturbance_size = disturbance_size/BASE_POWER;  
else
    disturbance_size = disturbance_size/BASE_POWER;
end


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