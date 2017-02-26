% global OLC_bus sumOfD_j OLC_capacity fcp_alpha a c
common_setttings;
load(['output/' char(RUNNING_MODE) '_' char(Method.proposed) '_' num2str(FLEX) '_' num2str(DELAY) '_' num2str(WEIGHT) '.mat'], 'load_freq');
load(['output/' char(RUNNING_MODE) '_' char(Method.proposed) '_' num2str(FLEX) '_' num2str(DELAY) '_' num2str(WEIGHT) '.mat'], 'sumOfD_j');
sumOfD_j = 0;
METHOD = Method.optimal;

w = 1- mean(load_freq(:,length(load_freq(1,:))));
N =  length(OLC_bus);      
% use CVX to compute control_d.
warning off;
cvx_begin quiet
  variables delta_j(N)        
  minimize( fcp_alpha/2* ((sum(a'.*delta_j))^2) + ...
        sum((c'/2).*(delta_j.^2)) + sumOfD_j/2*w^2 ...
        )        
  subject to
      delta_j >= OLC_capacity(:,1);
      delta_j <= OLC_capacity(:,2);
      sum(delta_j) + sumOfD_j*w == sum(disturbance_size);
cvx_end      
control_d  = delta_j;      
if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
    error('This optimization is unsolved');
end  
warning on;
controlled_load(:,1) = control_d;
optVal = cvx_optval

optCost = cvx_optval - sumOfD_j/2*w^2;

%%
strScenario = [char(RUNNING_MODE) '_' char(METHOD) '_' num2str(FLEX) '_' num2str(DELAY)]
matFile = [ strScenario '.mat'];
save(['output/' matFile]);