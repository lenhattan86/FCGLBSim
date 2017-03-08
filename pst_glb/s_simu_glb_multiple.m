clear global; clear; clc; close all; 
global RUNNING_MODE METHOD INCIDENT_START DELAY FLEX END_TIME WEIGHT
%%
IS_MULTIPLE_RUN = true;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');
%%
% runs = [true true true true true true true true];
% runs = [true true true true false false false false];
% runs = [false false false false true true true true];
runs = [false false false false false false false false]; runs(5) = true;
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
  METHOD = Method.optimal;
  computeOptimalCosts;
end

%% Demand flexibility
iRun = iRun+1;

%FLEX_ARRAY = 0.1:0.05:1.0;
FLEX_ARRAY = [0.21 0.22];
if runs(iRun)
  for i_FLEX_ARRAY=1:length(FLEX_ARRAY)  
    default_settings;
    END_TIME = 200;
    METHOD = Method.proposed;  
    FLEX = FLEX_ARRAY(i_FLEX_ARRAY);
    s_simu_glb;
  end
end

if runs(iRun)    
  for i_FLEX_ARRAY=1:length(FLEX_ARRAY)    
    default_settings;
    END_TIME = 200;
    METHOD = Method.OLC;
    FLEX = FLEX_ARRAY(i_FLEX_ARRAY);
    s_simu_glb;
  end
end

%% Delay 
iRun = iRun+1;
if false %runs(iRun)
  DELAY_ARRAY = 0.01*2.^(0:4);
  for i_DELAY_ARRAY=1:length(DELAY_ARRAY)
    default_settings;
    END_TIME = 1000;
    METHOD = Method.proposed;
    FLEX = 0.5;
    DELAY = DELAY_ARRAY(i_DELAY_ARRAY);
    s_simu_glb;
  end
end
%% Importance of frequency deviation.
iRun = iRun+1;
if runs(iRun)
  
  WEIGHT_ARRAY = [5:10:100];
  for i_WEIGHT_ARRAY=1:length(WEIGHT_ARRAY) 
    default_settings;
    END_TIME = 1000;
    METHOD = Method.proposed;
    WEIGHT = WEIGHT_ARRAY(i_WEIGHT_ARRAY);
    s_simu_glb;
  end
  
  WEIGHT_ARRAY = [5:10:100];
  for i_WEIGHT_ARRAY=1:length(WEIGHT_ARRAY) 
    default_settings;
    END_TIME = 1000;
    METHOD = Method.OLC;
    WEIGHT = WEIGHT_ARRAY(i_WEIGHT_ARRAY);
    s_simu_glb;
  end
end

%% Importance of global costs
iRun = iRun+1;
if runs(iRun)
  
  GAMMA_ARRAY = [0.005 0.015:0.015:0.3];
  for i_GAMMA_ARRAY=1:length(GAMMA_ARRAY)      
    default_settings;
    END_TIME = 100;
    METHOD = Method.proposed;
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