clear global; clear; clc; close all; 
global RUNNING_MODE METHOD INCIDENT_START DELAY FLEX END_TIME WEIGHT TIME_STEP fcp_lambda isSave
global POWER_LOSS
isSave = true;
% TIME_STEP = 0.01;
%%
IS_MULTIPLE_RUN = true;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');
%%
% Figure 4 & 5
% runs = [true false false false false false false false false false];
% Figure 6
% runs = [false false false false true false false false false false];
% runs = [false false false false false true false false false false];
% runs = [false false false false false false true false false false];
% runs = [false false false false false false false true false false];
%% trade off
runs = [false false false false false false false false true false ];
%% common
% RUNNING_MODE = RunningMode.LoadChange;
iRun = 0;
%% proposed
iRun = iRun+1;
if runs(iRun)
  default_settings;
  METHOD = Method.proposed_dc;
  s_simu_glb;
end
%% OLC
iRun = iRun+1;
if runs(iRun)
  default_settings;
  METHOD = Method.OLC;
  s_simu_glb;
end
%% Droop control
iRun = iRun+1;
if runs(iRun)
  default_settings;
  METHOD = Method.dc;
  s_simu_glb;
end
%% NONE
iRun = iRun+1;
if runs(iRun)
  default_settings;  
  METHOD = Method.NONE;
  s_simu_glb;
end
%% Importance of global costs
iRun = iRun+1;
if runs(iRun)
  
  GAMMA_ARRAY = [0.005 0.015:0.015:0.3];
  for i_GAMMA_ARRAY=1:length(GAMMA_ARRAY)      
    default_settings;
    END_TIME = 100;
    METHOD = Method.proposed_dc;
    GAMMA = GAMMA_ARRAY(i_GAMMA_ARRAY);
    s_simu_glb;
  end
end

if runs(iRun)  
  GAMMA_ARRAY = [0.005 0.015:0.015:0.3];
  for i_GAMMA_ARRAY=1:length(GAMMA_ARRAY) 
    default_settings;
    END_TIME = 100;
    METHOD = Method.OLC;
    GAMMA = GAMMA_ARRAY(i_GAMMA_ARRAY);
    s_simu_glb;
  end
end

%% Demand flexibility
iRun = iRun+1;

FLEX_ARRAY = 0:0.05:1.0;
if runs(iRun)
  for i_FLEX_ARRAY=1:length(FLEX_ARRAY)  
    default_settings;
    METHOD = Method.proposed_dc;  
    FLEX = FLEX_ARRAY(i_FLEX_ARRAY);
    s_simu_glb;
  end
end

if runs(iRun)    
  for i_FLEX_ARRAY=1:length(FLEX_ARRAY)    
    default_settings;
    METHOD = Method.OLC;
    FLEX = FLEX_ARRAY(i_FLEX_ARRAY);
    s_simu_glb;
  end
end

%% Delay 
iRun = iRun+1;
if runs(iRun)
  DELAY_ARRAY =  [0 0.1 0.2   0.5  1.0]; 
  scale_lambda = [1 1/2 1/2.5 1/5 1/10];
  for i_DELAY_ARRAY=1:length(DELAY_ARRAY)
    default_settings;
    END_TIME = 500;
    METHOD = Method.proposed_dc;
    FLEX = 0.4;
    TIME_STEP = 0.1;
    DELAY = DELAY_ARRAY(i_DELAY_ARRAY);
    fcp_lambda = fcp_lambda*scale_lambda(i_DELAY_ARRAY);
    s_simu_glb;
  end
end
%% impact of powerloss
iRun = iRun+1;
% POWER_LOSSES = [0 50 100 200 300 400 500 600 700 800 900 1000];
POWER_LOSSES = [0];
if runs(iRun)
    for i=1:length(POWER_LOSSES)
        default_settings;
        POWER_LOSS = POWER_LOSSES(i);
        METHOD = Method.proposed_dc;
        s_simu_glb;
    end
end

if runs(iRun)    
    for i=1:length(POWER_LOSSES)
        default_settings;
        POWER_LOSS = POWER_LOSSES(i);
        METHOD = Method.OLC;
        s_simu_glb;
    end
end

%% Importance of frequency deviation. trade-off frequency dev. & costs.
iRun = iRun+1;
if runs(iRun)  
%   WEIGHT_ARRAY = [0:4 5:10:155 160:40:300];
  WEIGHT_ARRAY = [320:40:400];
  for i_WEIGHT_ARRAY=1:length(WEIGHT_ARRAY) 
    default_settings;
    METHOD = Method.proposed_dc;
    WEIGHT = WEIGHT_ARRAY(i_WEIGHT_ARRAY);
    s_simu_glb;
  end
  
  for i_WEIGHT_ARRAY=1:length(WEIGHT_ARRAY) 
    default_settings;
    METHOD = Method.OLC;
    WEIGHT = WEIGHT_ARRAY(i_WEIGHT_ARRAY);
    s_simu_glb;
  end
end

%% TIME STEP 
iRun = iRun+1;
if runs(iRun)
  TIME_STEP_ARRAY = 0.01*2.^(5:10);
  for i_TIME_STEP_ARRAY=1:length(TIME_STEP_ARRAY)
    default_settings;
    END_TIME = 1000;
    METHOD = Method.proposed_dc;
    FLEX = 0.4;
    DELAY = 0.0;
    TIME_STEP = TIME_STEP_ARRAY(i_TIME_STEP_ARRAY);
    s_simu_glb;
  end
end