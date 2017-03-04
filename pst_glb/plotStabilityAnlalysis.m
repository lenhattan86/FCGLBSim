clear; close all; clc; clear global;
addpath('glb_data');
addpath('glb_func');
figure_settings;
warning off;
%%

PLOTS = [true true false false false false false];

X_LIM = 35;

% folder = 'output01/';
folder = 'output/';
GAMMA = 0.08;
DELAY = 0.01;
FLEX = 0.4;
WEIGHT = 75;
extra = ['_' num2str(FLEX) '_' num2str(DELAY) '_' num2str(WEIGHT) '_' num2str(GAMMA)];
dataFiles = {
    ['GenLoss_OLC' extra];
    ['GenLoss_proposed' extra];   
    ['GenLoss_NONE' extra]; 
    };  
optimalFile = ['GenLoss_optimal' extra];

  
figIdx = 0;
  strLegends = {strOLC, strProposed, strNone, strOptimal};
  lines = { lineOLC, lineProposed, lineNone, lineOptimal};
  colors = {colorOLC, colorProposed, colorNone, colorOptimal};
%%
if PLOTS(1)
    figure; 
    figSize = figOneCol;
    for i=1:length(dataFiles)
      load([folder dataFiles{i} '.mat']);
      hPlot(i) = plot(t,load_freq(1,:)*60, lines{i}, 'linewidth', lineWidth,'Color',colors{i});  hold on;
      plot(t,load_freq*60, lines{i}, 'linewidth',lineWidth,'Color',colors{i});
      hold on;
    end
    legend(hPlot, strLegends, 'Location','best','FontSize', fontLegend,'Orientation','horizontal');
    
    xlabel('Time (sec)', 'fontname', 'Arial','fontsize', fontAxis);
    ylabel('Frequency (Hz)','fontname', 'Arial','fontsize', fontAxis);
    % ylim([59.90, 60.0]); 
    xlim([0, X_LIM]);
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'load_freq';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end
%%
if PLOTS(2) 
    figure; 
    figSize = figOneCol;
    for i=2:length(dataFiles)-1
      load([folder dataFiles{i} '.mat']);
%       hPlot(i) = plot(t, controlled_load(1,:), lines{i}, 'linewidth', lineWidth,'Color',colors{i});
%       hold on;
%       plot(t, controlled_load*NEW_ENG_BASE/BASE_POWER, lines{i}, 'linewidth',lineWidth,'Color',colors{i});
      plot(t, controlled_load*NEW_ENG_BASE/BASE_POWER, '-.','linewidth',lineWidth);
      hold on;
    end
    a = 1:10 ;
%     strLegends = strread(num2str(a),'%s');
%     legend(hPlot, strLegends,'Location','best','FontSize', fontLegend,'Orientation','horizontal');
%     legend(strLegends,'Location','best','FontSize', fontLegend-2,'Orientation','horizontal','Right');
    xlim([0, X_LIM]);
    xlabel('Time (sec)','fontname', 'Arial','fontsize', fontAxis);
    ylabel('MW','fontname', 'Arial','fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'load';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end
%%
if PLOTS(3)
    figure; 
    figSize = figOneCol;
    
    for i=1:length(dataFiles)
      load([folder dataFiles{i} '.mat']);
      costs = fcp_alpha/2*(a*controlled_load).^2;
      hPlot(i) = plot(t,costs, lines{i}, 'linewidth', lineWidth,'Color',colors{i});      
      hold on;
    end       

    
    legend(hPlot, strLegends,'Location','best','FontSize', fontLegend,'Orientation','horizontal');
    xlim([0, X_LIM]);
    xlabel('Time (sec)','fontname', 'Arial','fontsize', fontAxis);
    ylabel('global','fontname', 'Arial','fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'global_cost';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end

%%
if PLOTS(4)
    figure; 
    figSize = figOneCol;
    
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
    xlabel('Time (sec)','fontname', 'Arial','fontsize', fontAxis);
    ylabel('global + local costs','fontname', 'Arial','fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'cost';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
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
    figure; 
    for i=1:length(dataFiles)      
      load([folder dataFiles{i} '.mat']);
      hPlot(i) = plot(t,load_freq(1,:)*60, 'linewidth', lineWidth);  hold on;
      hold on;
    end    
    
    figSize = figOneCol;
    
    xlabel('Time (sec)','fontname', 'Arial','fontsize', fontAxis);
    ylabel('Frequency (Hz)','fontname', 'Arial','fontsize', fontAxis);
    strLegends={'10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'};
    legend(hPlot, strLegends,'Location','eastoutside','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'load_freq_flex';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%         print ('-depsc', epsFile);
%     end
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
    figure; 
    for i=1:length(dataFiles)      
      load([folder dataFiles{i} '.mat']);
      hPlot(i) = plot(t,load_freq(1,:)*60, 'linewidth', lineWidth);  hold on;
      hold on;
    end    
    
    figSize = figOneCol;
    
    xlabel('Time (sec)','fontname', 'Arial','fontsize', fontAxis);
    ylabel('Frequency (Hz)','fontname', 'Arial','fontsize', fontAxis);
%     strLegends={'10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'};
%     legend(hPlot, strLegends,'Location','eastoutside','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'load_freq_flex';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%         print ('-depsc', epsFile);
%     end
end
%% weight of frequency deviation (beta)
if PLOTS(7)
    dataFiles = {'GenLoss_proposed_0.5_0.01_1_1'...
            ,'GenLoss_proposed_0.5_0.02_1_1'...
            ,'GenLoss_proposed_0.5_0.04_1_1'...
            ,'GenLoss_proposed_0.5_0.08_1_1'...
            ,'GenLoss_proposed_0.5_0.16_1_1'...
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
%     strLegends={'10%','20%','30%','40%','50%','60%','70%','80%','90%','100%'};
%     legend(hPlot, strLegends,'Location','eastoutside','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'load_freq_flex';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%         print ('-depsc', epsFile);
%     end
end
%%
warning on;
return;
%%
for i=1:length(fileNames)
    fileName = fileNames{i};
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end
