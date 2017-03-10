%%
clear all; clear global; close; clc;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');

global RUNNING_MODE METHOD INCIDENT_START DELAY FLEX END_TIME WEIGHT
FLEX = 0.5;
WEIGHT = 1;
END_TIME = 35;
%METHOD = Method.optimal;
 METHOD = Method.OLC;
RUNNING_MODE  = RunningMode.Motivation
% common_setttings

%%
strScenario = [char(RUNNING_MODE) '_' char(METHOD)];% 

GAMMA_ARRAY = 0:0.1:10;
N = 2;
gamma = 0;

DC_DEMAND = 20;
DC_CAPACITY = 30;
fcp_alpha = 1;
a = [0.9 0.5];
c = [1 1];
b = -28; % b = sum(a) * DC_DEMAND.
OLC_capacity = [-10 10];
disturbance = 10;

for i = 1:length(GAMMA_ARRAY)  
  if METHOD == Method.optimal
    gamma = GAMMA_ARRAY(i);
  end
  warning off;
  cvx_begin quiet
    variables delta_j(N)            
    minimize( gamma * fcp_alpha/2* ((sum(a'.*delta_j))^2) + ...
          sum((c'/2).*(delta_j.^2)) ...
          )        
    subject to
        delta_j >= OLC_capacity(:,1);
        delta_j <= OLC_capacity(:,2);
        sum(delta_j) == disturbance;
  cvx_end      
  control_d  = delta_j;      
  if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
      error('This optimization is unsolved');
  end  
  warning on;
  controlled_load(:,1) = control_d;
  localCosts(i) = sum((c'/2).*(delta_j.^2));
  globalCosts(i) = GAMMA_ARRAY(i) * fcp_alpha/2* ((sum(a'.*delta_j))^2) ;
  totalCosts(i)= localCosts(i) + globalCosts(i);
end

%%
matFile = [ strScenario '.mat'];
save(['output/' matFile]);

%%
figure_settings;
figIdx = 0;
figure; 

figSize = figOneCol;
yVals =  totalCosts;
plot(GAMMA_ARRAY, yVals, 'b-', 'linewidth',lineWidth);

%     legend(hBar, strLegends,'Location','southeast','FontSize', fontLegend,'Orientation','horizontal');
ylim([0, max(yVals)]);    
ylabel('normalized cost','fontname', 'Arial','fontsize', fontAxis);    
xlabel('weight','fontname', 'Arial','fontsize', fontAxis);
set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

