clear; close all; clc; clear global;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');
figure_settings;
warning off;
%%

PLOTS = [true true true false false false false];

X_LIM = 30;

% folder = 'output01/';
folder = 'output/';
default_settings;
extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
        num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];
% dataFiles = {
%     ['GenLoss_OLC_0.4_0.01_75_0.16_0_0.001'];
%     ['GenLoss_proposed_0.4_0.01_75_0.16_0_0.001'];
%     ['GenLoss_NONE_0.4_0.01_75_0.16_0_0.001'];
%     };
% strLegends = {strOLC, strProposed, strNone, strOptimal};

lines = { lineOLC, lineProposed, lineNone, lineOptimal};
colors = {colorOLC, colorProposed, colorNone, colorOptimal};

dataFiles = {
    ['GenLoss_OLC' extra];
    ['GenLoss_proposed_dc' extra];
    ['GenLoss_dc' extra];
    };
strLegends = {strOLC, strProposed, strDroopControl, strOptimal};

optimalFile = ['GenLoss_optimal_0.4'];
  
figIdx = 0;

%%
if PLOTS(1)
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'load_freq';
    
    figSize = figTwoThirdCol;
    for i=1:length(dataFiles)%-1
      load([folder dataFiles{i} '.mat']);
      hPlot(i) = plot(t,load_freq(1,:)*60, lines{i}, 'linewidth', lineWidth,'Color',colors{i});  hold on;
      plot(t,load_freq*60, lines{i}, 'linewidth',1,'Color',colors{i});
      hold on;
    end
    legend(hPlot, strLegends, 'Location','northeast','FontSize', fontLegend,'Orientation','vertical');
    xlabel(strTime, 'fontname', 'Arial','fontsize', fontAxis);
    ylabel(strFrequency,'fontname', 'Arial','fontsize', fontAxis);
%     ylim([59.9, 60.0]); 
    xlim([0, X_LIM]);
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
end

if false
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'load_freq_2';
    
    figSize = figTwoThirdCol;
    for i=3:length(dataFiles)
      fileToLoad = [folder dataFiles{i} '.mat']
      if(exist(fileToLoad,'file'))
          load(fileToLoad);
          hPlot(i) = plot(t,load_freq(1,:)*60, lines{i}, 'linewidth', lineWidth,'Color',colors{i});  hold on;
          plot(t,load_freq*60, lines{i}, 'linewidth',lineWidth,'Color',colors{i});
          hold on;
      end
    end
%     legend(hPlot, strLegends, 'Location','best','FontSize', fontLegend,'Orientation','horizontal');
    
    xlabel(strTime, 'fontname', 'Arial','fontsize', fontAxis);
    ylabel(strFrequency,'fontname', 'Arial','fontsize', fontAxis);
    % ylim([59.90, 60.0]); 
    xlim([0, X_LIM]);
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
end
%%
if PLOTS(2) 
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'load';
    
    figSize = figTwoThirdCol;
    for i=2:length(dataFiles)-1
      load([folder dataFiles{i} '.mat']);
%       hPlot(i) = plot(t, controlled_load(1,:), lines{i}, 'linewidth', lineWidth,'Color',colors{i});
%       hold on;
%       plot(t, controlled_load*NEW_ENG_BASE/BASE_POWER, lines{i}, 'linewidth',lineWidth,'Color',colors{i});
      plot(t, controlled_load*NEW_ENG_BASE, '-.','linewidth',lineWidth);
      hold on;
    end
    xlim([0, X_LIM]);
    a = 1:10 ;
%     strLegends = strread(num2str(a),'%s');
%     legend(hPlot, strLegends,'Location','best','FontSize', fontLegend,'Orientation','horizontal');
%     legend(strLegends,'Location','best','FontSize', fontLegend-2,'Orientation','horizontal','Right');
%     xlim([0, X_LIM]);
    xlabel(strTime,'fontname', 'Arial','fontsize', fontAxis);
    ylabel(strPower,'fontname', 'Arial','fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

end

%%
if PLOTS(3)
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'cost_time';
    
    figSize = figTwoThirdCol;
        
    for i=1:length(dataFiles)-1
      load([folder dataFiles{i} '.mat']);
      
      d_j = controlled_load;
      independent_cost = WEIGHT*(c'/2) * (d_j.^2)/2;  
      interdependent_cost = WEIGHT*(fcp_gamma)/2*(a'*d_j).^2;  
      
      hPlot(i) = plot(t,independent_cost  + interdependent_cost, lines{i}, 'linewidth', lineWidth,'Color',colors{i});      
      hold on;
    end       

    
    legend(strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
    xlim([0, X_LIM]);
    xlabel(strTime,'fontname', 'Arial','fontsize', fontAxis);
    ylabel(strTotalCost,'fontname', 'Arial','fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
end

%%
if PLOTS(4)
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'cost';
    
    figSize = figHalfCol;
    
    for i=1:length(dataFiles)
      load([folder dataFiles{i} '.mat']);
      costs = fcp_alpha/2*(a*controlled_load).^2 + (ones(size(OLC_gain))./OLC_gain)' * (controlled_load.^2)/2;
      hPlot(i) = plot(t,costs, lines{i}, 'linewidth', lineWidth,'Color',colors{i});      
      hold on;
    end
    
    load([folder optimalFile '.mat'],'optCost');
    for t_idx = 1:length(t)
      if t(t_idx) < INCIDENT_START
        costs(t_idx) = 0;
      else
        costs(t_idx) = optCost;
      end
    end
    
    hPlot(i+1) = plot(t, costs, lines{i+1}, 'linewidth', lineWidth,'Color',colors{i+1});
    hold on;
    
    legend(hPlot, strLegends,'Location','southeast','FontSize', fontLegend,'Orientation','horizontal');
    xlim([0, X_LIM]);
    xlabel(strTime,'fontname', 'Arial','fontsize', fontAxis);
    ylabel('global + local costs','fontname', 'Arial','fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
end
%% demand flexiblity
if PLOTS(5)
    dataFiles = {
       'GenLoss_proposed_0.1_0.01_1'...
      ,'GenLoss_proposed_0.2_0.01_1'...
      ,'GenLoss_proposed_0.3_0.01_1'...
      ,'GenLoss_proposed_0.4_0.01_1'...
      ,'GenLoss_proposed_0.5_0.01_1'...        
      ,'GenLoss_proposed_0.6_0.01_1'...
      ,'GenLoss_proposed_0.7_0.01_1'...
      ,'GenLoss_proposed_0.8_0.01_1'...
      ,'GenLoss_proposed_0.9_0.01_1'...
      ,'GenLoss_proposed_1_0.01_1'...   
      };
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'load_freq_flex';
    
    for i=1:length(dataFiles)      
      load([folder dataFiles{i} '.mat']);
      hPlot(i) = plot(t,load_freq(1,:)*60, 'linewidth', lineWidth);  hold on;
      hold on;
    end    
    
    figSize = figHalfCol;
    
    xlabel(strTime,'fontname', 'Arial','fontsize', fontAxis);
    ylabel(strFrequency,'fontname', 'Arial','fontsize', fontAxis);
    strLegends={'10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'};
    legend(hPlot, strLegends,'Location','eastoutside','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
end

%% weight of frequency deviation (beta)
if PLOTS(6)
    dataFiles = {
       'GenLoss_proposed_0.5_0.01_0.1'...
      ,'GenLoss_proposed_0.5_0.01_0.3'...
      ,'GenLoss_proposed_0.5_0.01_0.5'...
      ,'GenLoss_proposed_0.5_0.01_0.7'...
      ,'GenLoss_proposed_0.5_0.01_0.9'...        
      ,'GenLoss_proposed_0.5_0.01_1'...
      ,'GenLoss_proposed_0.5_0.01_2'...
      ,'GenLoss_proposed_0.5_0.01_4'...
      ,'GenLoss_proposed_0.5_0.01_8'...
      };
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'load_freq_flex';
    
    for i=1:length(dataFiles)      
      load([folder dataFiles{i} '.mat']);
      hPlot(i) = plot(t,load_freq(1,:)*60, 'linewidth', lineWidth);  hold on;
      hold on;
    end    
    
    figSize = figHalfCol;
    
    xlabel(strTime,'fontname', 'Arial','fontsize', fontAxis);
    ylabel(strFrequency,'fontname', 'Arial','fontsize', fontAxis);
%     strLegends={'10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'};
%     legend(hPlot, strLegends,'Location','eastoutside','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

end
%% weight of frequency deviation (beta)
if PLOTS(7)
    dataFiles = {'GenLoss_proposed_0.5_0.01_1_1'...
            ,'GenLoss_proposed_0.5_0.02_1_1'...
            ,'GenLoss_proposed_0.5_0.04_1_1'...
            ,'GenLoss_proposed_0.5_0.08_1_1'...
            ,'GenLoss_proposed_0.5_0.16_1_1'...
            }; 
    figIdx=figIdx +1;
    figures{figIdx} = figure; 
    fileNames{figIdx} = 'load_freq_flex';
    
    for i=1:length(dataFiles)      
      load([folder dataFiles{i} '.mat']);
      hPlot(i) = plot(t,load_freq(1,:)*60, 'linewidth', lineWidth);  hold on;
      hold on;
    end    
    
    figSize = figHalfCol;
    
    xlabel(strTime,'fontname', 'Arial','fontsize', fontAxis);
    ylabel(strFrequency,'fontname', 'Arial','fontsize', fontAxis);
%     strLegends={'10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'};
%     legend(hPlot, strLegends,'Location','eastoutside','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
end
%%
warning on;
return;
%%
for i=1:length(fileNames)
    fileName = fileNames{i};
    epsFile = [ LOCAL_FIG fileName '.eps'];     
    print (figures{i}, '-depsc', epsFile);
    
    pdfFile = [ fig_path fileName '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end
