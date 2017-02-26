clear; close all; clc;
addpath('glb_func');
addpath('glb_data');

additional_setttings
%%
% 1 - a_j*Ph

%%

IS_MULTIPLE_RUN = true;
IS_PLOT = false;

% phi = 2^(-9); 2(-8); ... 2^0.
for powerOfPhi = -9:1:0;
    Phi = 1-2^powerOfPhi;
    a = ones(1,DC_num);  % assume that a_j = 1. => compute gamma according to Phi
    gamma = Phi/((1-Phi)/sum(c))
    
    if Phi <= min(ones(size(a))./a) % Phi <= all 1/a_j
        if ASSUMPTION~=1
            error('[INFO] Phi must be > one of 1/a_j');
        end
    else
        if ASSUMPTION==1
            error('[INFO] Phi  must be <= min(ones(size(a))./a)');
        end
    end
    1-mean(Phi*a)
    EXTRA=['_powerOfPhi_' int2str(powerOfPhi)]; 
    s_simu_glb;    
end