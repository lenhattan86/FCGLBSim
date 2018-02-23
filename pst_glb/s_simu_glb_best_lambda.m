clear global; clear; clc; close all; 
global RUNNING_MODE METHOD INCIDENT_START DELAY FLEX END_TIME WEIGHT TIME_STEP fcp_lambda isSave
%%
IS_MULTIPLE_RUN = true;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');

%%

if true %runs(iRun) 
%   DELAY_ARRAY = 0.01*2.^(5:10);
  TIME_STEP = 0.1;
%   DELAY_ARRAY = [0 0.1 0.2 0.4 0.8 1.0 2.0 3.0 4.0 5.0]; 
  DELAY_ARRAY = [0 0.1 1.0 5.0]; 
  lambdas = [0.001 0.0001];
  bestLambdas = zeros(size(DELAY_ARRAY));
  for i_DELAY_ARRAY=1:length(DELAY_ARRAY)
    shortestConvTime = inf;
    for iLambda = 1:length(lambdas)
      isSave = false;
      default_settings;
      END_TIME = 200;
      METHOD = Method.proposed;
      FLEX = 0.4;
      
      DELAY = DELAY_ARRAY(i_DELAY_ARRAY);
      fcp_lambda = lambdas(iLambda);
      s_simu_glb;
      convergentTime;
      if shortestConvTime > conv_time
        shortestConvTime = conv_time;
        bestLambdas(i_DELAY_ARRAY) = fcp_lambda;
      end
    end
  end
  bestLambdaFile = sprintf('output/bestLambdas_%s.mat', num2str(TIME_STEP));
  save('output/bestLambdas.mat','bestLambdas');
end

