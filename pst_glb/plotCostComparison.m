clear; close all; clc; clear global;
addpath('glb_data');
addpath('glb_func');
figure_settings;
%%

PLOTS = [true true];
X_LIM = 35;

% folder = 'output01/';
folder = 'output/';


  
optimalFile = 'LoadChange_centralized_0.5_0.01';
  
figIdx = 0;
strLegends = {strOLC, strProposed, strOptimal};
lines = { lineOLC, lineProposed, lineOptimal};
colors = {colorOLC, colorProposed, colorOptimal};
  
%%
if PLOTS(1)
    dataFiles = {
      'LoadChange_OLC_0.5_0.01';
      'LoadChange_proposed_0.5_0.01';   
    };
    figure; 
    figSize = figOneCol;
    
    for i=1:length(dataFiles)
      load([folder dataFiles{i} '.mat']);
      d_j = controlled_load(:,length(controlled_load(1,:)));
      costs(i) = fcp_alpha/2*(sum(a.*d_j')).^2 + (ones(size(OLC_gain))./OLC_gain)' * (d_j.^2)/2;      
    end
    
    load([folder optimalFile '.mat'],'optCost');    
    costs(i+1) = optCost;    
    
    hBar = bar(costs, 0.2);
    set(gca,'XTickLabel',strLegends,'FontSize',fontAxis);
    hold on;
    
%     legend(hBar, strLegends,'Location','southeast','FontSize', fontLegend,'Orientation','horizontal');
    xlim([0.5, 3.5]);    
    ylabel('global + local costs','fontname', fontName,'fontsize', fontAxis);    
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
for i=1:length(fileNames)
    fileName = fileNames{i};
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end