clear; close all; clc; clear global;
addpath('glb_data');
addpath('glb_func');
figure_settings;
warning off;
%%

PLOTS = [true true];
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
  
%%
if PLOTS(1)
    dataFiles = {
      ['GenLoss_OLC' extra];
      ['GenLoss_proposed' extra];   
    };
    figure; 
    figSize = figOneCol;
    
    for iFile=1:length(dataFiles)
      load([folder dataFiles{iFile} '.mat']);
      d_j = controlled_load(:,length(controlled_load(1,:)));
      costs(iFile,1) = sum((c/2).* (d_j.^2)/2);  
      costs(iFile,2) = (fcp_gamma)/2*(sum(a.*d_j)).^2;  
    end
    
    load([folder optimalFile '.mat']);    
    costs(iFile+1,1) = sum((c/2).* (d_j.^2)/2);   
    costs(iFile+1,2) = (fcp_gamma)/2*(sum(a.*d_j)).^2;       
    
    costs = costs*(BASE_POWER^2);
    
    hBar = bar(costs, 0.2, 'stacked');
    set(gca,'XTickLabel', strLegends,'FontSize',fontAxis);
    hold on;
    
    legend(hBar, {strLocalCost,strGlobalCost},'Location','northeast','FontSize', fontLegend,'Orientation','horizontal');
    xlim([0.5, 3.5]);    
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