clear global; clear; clc; close all; 
global RUNNING_MODE METHOD INCIDENT_START DELAY FLEX END_TIME WEIGHT
%%
IS_MULTIPLE_RUN = true;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');
%%
% runs = [true true true true true true true true];
% runs = [false false false false true true true true];
runs = [false false false false false false false false]; runs(6) = true;
%% common
% RUNNING_MODE = RunningMode.LoadChange;
iRun = 0;
%% proposed
iRun = iRun+1;
if runs(iRun)
  %
  default_settings;    
  METHOD = Method.proposed;    
  s_simu_glb;  
end
%% OLC
iRun = iRun+1;
if runs(iRun)
  default_settings;  
  METHOD = Method.OLC;
  s_simu_glb;
end
%% NONE
iRun = iRun+1;
if runs(iRun)
  default_settings;  
  METHOD = Method.NONE;
  s_simu_glb;
end
%% optimal
iRun = iRun+1;
if runs(iRun)
  default_settings;
  RUNNING_MODE = RunningMode.GenLoss;
  METHOD = Method.optimal;
  computeOptimalCosts;
end

%% Demand flexibility
iRun = iRun+1;
if runs(iRun)
  default_settings;
  END_TIME = 150;
  METHOD = Method.proposed;
  FLEX_ARRAY = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
  %   FLEX_ARRAY = [0.1 0.5 1.0];
  for i_FLEX_ARRAY=1:length(FLEX_ARRAY)      
    FLEX = FLEX_ARRAY(i_FLEX_ARRAY);
    s_simu_glb;
  end
end

%% Delay 
iRun = iRun+1;
if runs(iRun)
  default_settings;
  END_TIME = 1000;
  METHOD = Method.proposed;
  FLEX = 0.5;
%   DELAY_ARRAY = 0.01*2.^(0:9);
  DELAY_ARRAY = 0.01*2.^(0:4);
  for i_DELAY_ARRAY=1:length(DELAY_ARRAY)      
    DELAY = DELAY_ARRAY(i_DELAY_ARRAY);
    s_simu_glb;
  end
end
%% Importance of frequency deviation.
iRun = iRun+1;
if runs(iRun)
  default_settings;
  END_TIME = 150;
  METHOD = Method.proposed;
  FLEX = 0.5;
  WEIGHT_ARRAY = [0.1:0.2:1 2.^(0:3)];
  for i_WEIGHT_ARRAY=1:length(WEIGHT_ARRAY)      
    WEIGHT = WEIGHT_ARRAY(i_WEIGHT_ARRAY);
    s_simu_glb;
  end
end

%% Importance of global costs
iRun = iRun+1;
if runs(iRun)
  default_settings;
  END_TIME = 150;
  METHOD = Method.proposed;
  FLEX = 0.5;
  GAMMA_ARRAY = [0:0.1:5];
  for i_GAMMA_ARRAY=1:length(GAMMA_ARRAY)      
    GAMMA = GAMMA_ARRAY(i_GAMMA_ARRAY);
    s_simu_glb;
  end
end