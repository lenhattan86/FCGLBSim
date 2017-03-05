clear; close all; clc; clear global;
addpath('glb_data');
addpath('glb_func');
figure_settings;
warning off;
%%

PLOTS = [true true true true];
X_LIM = 35;

% folder = 'output01/';
folder = 'output/';

GAMMA = 0.08;
DELAY = 0.01;
FLEX = 0.4;
WEIGHT = 75;
extra = ['_' num2str(FLEX) '_' num2str(DELAY) '_' num2str(WEIGHT) '_' num2str(GAMMA)];
  
optimalFile = ['GenLoss_optimal' extra];
  
figIdx = 0;
strLegends = {strOLC, strProposed, strOptimal};
lines = { lineOLC, lineProposed, lineOptimal};
colors = {colorOLC, colorProposed, colorOptimal};

dataFiles = {
  ['GenLoss_OLC' extra];
  ['GenLoss_proposed' extra];   
};

%%
if PLOTS(1)
    figure; 
    figSize = figTwoThirdCol;
    
    for iFile=1:length(dataFiles)
      load([folder dataFiles{iFile} '.mat']);
      d_j = controlled_load(:,length(controlled_load(1,:)));
      costs(iFile,1) = sum((c/2).* (d_j.^2)/2);  
      costs(iFile,2) = (fcp_gamma)/2*(sum(a.*d_j)).^2;  
    end
    
    load([folder optimalFile '.mat']);    
    costs(iFile+1,1) = sum((c/2).* (d_j.^2)/2);   
    costs(iFile+1,2) = (fcp_gamma)/2*(sum(a.*d_j)).^2;       
    
    costs = costs*(NEW_ENG_BASE^2);
    
    hBar = bar(costs, 0.2, 'stacked');
    set(gca,'XTickLabel', strLegends,'FontSize',fontAxis);
    hold on;
    colorBreakDown = {colorLocal, colorGlobal};
    for i=1:length(hBar)
       hBar(i).FaceColor = colorBreakDown{i};
    end
    
    legend(hBar, {strLocalCost,strGlobalCost},'Location','best','FontSize', fontLegend,'Orientation','horizontal');
    xlim([0.5, 3.5]);   
    ylim([0 max(sum(costs,2)*1.3)]);
    ylabel(strTotalCost,'fontname', fontName,'fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'cost';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end

%%
if PLOTS(2)
    figure; 
    figSize = figTwoThirdCol;
    
    for iFile=1:length(dataFiles)
      load([folder dataFiles{iFile} '.mat']);
      d_j = controlled_load(:,length(controlled_load(1,:)));
      COSTS(iFile,:) = (c/2).* (d_j.^2)/2;  
%       COSTS(iFile,2) = (fcp_gamma)/2*(sum(a.*d_j)).^2;  
    end
    
    COSTS = COSTS*(NEW_ENG_BASE^2);
    
    hBar = bar(COSTS', 0.8, 'grouped');
%     set(gca,'XTickLabel', strLegends,'FontSize',fontAxis);
    hold on;
    
    colorBreakDown = {colorOLC, colorProposed};
    for i=1:length(hBar)
       hBar(i).FaceColor = colorBreakDown{i};
    end
    
    legend(hBar, {strOLC, strProposed},'Location','best','FontSize', fontLegend,'Orientation','horizontal');
    xlim([0.5, length(COSTS(1,:))+0.5]); 
    ylim([0 max(max(COSTS))*1.3]);
    ylabel(strTotalCost,'fontname', fontName,'fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'breakdown_cost';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end

%%
if PLOTS(3)

    figure; 
    figSize = figTwoThirdCol;
    
    for iFile=1:length(dataFiles)
      load([folder dataFiles{iFile} '.mat']);
      d_j = controlled_load(:,length(controlled_load(1,:)));
      LOADS(iFile,:) = d_j;  
    end
    
    LOADS = LOADS*(NEW_ENG_BASE);
    
    hBar = bar(LOADS', 0.8, 'grouped');
%     set(gca,'XTickLabel', strLegends,'FontSize',fontAxis);
    hold on;
    
        colorBreakDown = {colorOLC, colorProposed};
    for i=1:length(hBar)
       hBar(i).FaceColor = colorBreakDown{i};
    end
    
    legend(hBar, {strOLC, strProposed},'Location','best','FontSize', fontLegend,'Orientation','horizontal');
    xlim([0.5, length(LOADS(1,:))+0.5]); 
    ylim([min(min(LOADS))*1.3 max(max(LOADS))*1.3]);
    ylabel(strPower,'fontname', fontName,'fontsize', fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'load_deviation';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end

%%

%%
return;
%%
for iFile=1:length(fileNames)
    fileName = fileNames{iFile};
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end