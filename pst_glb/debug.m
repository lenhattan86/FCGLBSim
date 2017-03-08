clear; close all; clc;
addpath('glb_data');
addpath('glb_func');
figure_settings;
warning off;
folder = 'output/';
dataFiles = {%'GenLoss_proposed_0.4_0.01_75_0.005'...
            %,'GenLoss_proposed_0.4_0.01_75_0.015'...
            'GenLoss_proposed_0.4_0.01_75_0.03'...
            ,'GenLoss_proposed_0.4_0.01_75_0.045'...
            ,'GenLoss_proposed_0.4_0.01_75_0.06'...
            ,'GenLoss_proposed_0.4_0.01_75_0.075'...
            ,'GenLoss_proposed_0.4_0.01_75_0.09'...
            };

figure; 
for i=1:length(dataFiles)      
  load([folder dataFiles{i} '.mat']);
  hPlot(i) = plot(t,load_freq(1,:)*60, 'linewidth', lineWidth);  hold on;
  hold on;
end    

figSize = figOneCol;

xlabel('Time (sec)','fontname', 'Arial','fontsize', fontAxis);
ylabel('Frequency (Hz)','fontname', 'Arial','fontsize', fontAxis);
set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
