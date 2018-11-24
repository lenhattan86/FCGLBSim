%% global OLC_bus sumOfD_j OLC_capacity fcp_alpha a c
clear; close all; clc;

addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');
default_settings;
METHOD = Method.optimal;
strScenario = [char(RUNNING_MODE) '_' char(METHOD) '_' num2str(FLEX)];
extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
        num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];

common_setttings;

% vars ={'load_freq','controlled_load','NEW_ENG_BASE','BASE_POWER','bus'};
load(['output/GenLoss_proposed_dc' extra '.mat'] );
sumOfD_j = 0;


w = 1- mean(load_freq(:,length(load_freq(1,:))));
N =  length(OLC_bus);      
% use CVX to compute control_d.
total_load_dev = sum(controlled_load(:,size(controlled_load,2)) );
warning off;
cvx_begin %quiet
  variables delta_j(N)        
  minimize( (fcp_gamma)/2* ((sum(a.*delta_j))^2) + ...
        sum((c/2).*(delta_j.^2)) + sumOfD_j/2*w^2 ...
        )        
  subject to
      delta_j >= ones(size(delta_j))*OLC_capacity(:,1);
      delta_j <= ones(size(delta_j))*OLC_capacity(:,2);
%       sum(delta_j) + sumOfD_j*w + w*sum(tg_con(:,4).*bus(30:39,4) ) == sum(disturbance_size);
%     sum(delta_j) + sumOfD_j*w - sum(pelect(:, length(pelect))-bus(30:39,4)) == 0;
     sum(delta_j)  == total_load_dev;
cvx_end
control_d  = delta_j;      
if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
    error('This optimization is unsolved');
end  
warning on;
% controlled_load(:,1) = control_d;
optVal = cvx_optval;
optCost = cvx_optval - sumOfD_j/2*w^2;

%%
strScenario
matFile = [ strScenario '.mat'];
save(['output/' matFile]);