
clear; close all; clc; clear global;
addpath('glb_data');
addpath('glb_func');
figure_settings;

PLOTS = [false true false];
%%
FILES = {'Motivation_OLC'...
            ,'Motivation_optimal'...
            };
          
%           FILES = {'Motivation_OLC_0.5_1'
%             };
%%
figure_settings;
figIdx = 0;
folder = 'output/';
figSize = figOneCol;

strLegends = {strOLC, strProposed};
colors = {colorOLC, colorProposed};
lines = { lineOLC, lineProposed};

%% 
if PLOTS(1)
  figure; 
  
  for iFile = 1:length(FILES)
    fileToLoad = [folder FILES{iFile} '.mat'];
    load(fileToLoad, 'totalCosts'); 
    load(fileToLoad, 'GAMMA_ARRAY'); 
    yVals(iFile,:) = totalCosts;
  %   plot(GAMMA_ARRAY, yVals, lines{iFile},'linewidth',lineWidth,'Color', colors{iFile});
  %   hold on;
  end

  plot(GAMMA_ARRAY, yVals(1,:)./yVals(2,:),'linewidth', lineWidth);
  xlim([0 5]);
  ylim([0 1.5]);

  % legend(strLegends,'Location','northwest','FontSize', fontLegend,'Orientation','horizontal');
  ylabel('normalized cost','fontname', 'Arial','fontsize', fontAxis);    
  xlabel('\gamma','fontname', 'Arial','fontsize', fontAxis);
  set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

  if is_printed
    figIdx=figIdx +1;
    fileNames{figIdx} = 'motivation';
    epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
  end
end
%%

if PLOTS(2)
    figSize = figHalfCol;
  strLegends = {strFCLocal,strFCOptimal};
  
  
  iFile=1;
  fileToLoad = [folder FILES{iFile} '.mat'];
  load(fileToLoad, 'globalCosts'); 
  load(fileToLoad, 'localCosts'); 
  load(fileToLoad, 'totalCosts');
  load(fileToLoad, 'GAMMA_ARRAY'); 
  
  Y1(1,:) = localCosts;
  Y1(2,:) = globalCosts;
  
  iFile=2;
  fileToLoad = [folder FILES{iFile} '.mat'];
  load(fileToLoad, 'globalCosts'); 
  load(fileToLoad, 'localCosts'); 
  load(fileToLoad, 'totalCosts');
  load(fileToLoad, 'GAMMA_ARRAY'); 
  Y2(1,:) = localCosts;
  Y2(2,:) = globalCosts;
  
  Y_MAX = max(max(max(Y1)), max(max(Y2)));
  figure;  
  plot(GAMMA_ARRAY, Y1(1,:), lines{1},'linewidth',lineWidth,'Color', colors{1});
  hold on;
  plot(GAMMA_ARRAY, Y2(1,:), lines{2},'linewidth',lineWidth,'Color', colors{2});
  ylim([0 Y_MAX]);
  xlim([0 7]);
  legend(strLegends,'Location','northwest','FontSize', fontLegend,'Orientation','vertical');
  ylabel(strCost,'fontname', fontName,'fontsize', fontAxis);    
  xlabel('\gamma','fontname', fontName,'fontsize', fontAxis);
  set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

  if is_printed
    figIdx=figIdx +1;
    fileNames{figIdx} = 'mov_local';
    epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
  end

  
  figure; 
  plot(GAMMA_ARRAY, Y1(2,:), lines{1},'linewidth',lineWidth,'Color', colors{1});
  hold on;
  plot(GAMMA_ARRAY, Y2(2,:), lines{2},'linewidth',lineWidth,'Color', colors{2});
  ylim([0 Y_MAX]);
  xlim([0 7]);
  legend(strLegends,'Location','northwest','FontSize', fontLegend,'Orientation','vertical');
  ylabel(strCost,'fontname', fontName,'fontsize', fontAxis);    
  xlabel('\gamma','fontname', fontName,'fontsize', fontAxis);
  set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

  if is_printed
    figIdx=figIdx +1;
    fileNames{figIdx} = 'mov_global';
    epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
  end
end

%%
if PLOTS(3)    
  strLegends = {strLocalCost, strGlobalCost, strFCOptimal, strFCLocal};
  
  
  iFile=2;
  fileToLoad = [folder FILES{iFile} '.mat'];
  load(fileToLoad, 'globalCosts'); 
  load(fileToLoad, 'localCosts'); 
  load(fileToLoad, 'totalCosts');
  load(fileToLoad, 'GAMMA_ARRAY'); 
  
  yVals(1,:) = localCosts;
  yVals(2,:) = globalCosts;
  
 
  bar(GAMMA_ARRAY, yVals',1,'stacked','EdgeColor','none');
  hold on;
  plot(GAMMA_ARRAY, totalCosts, lines{iFile},'linewidth',lineWidth,'Color', colors{iFile});
   hold on;
 
  
  iFile=1;
  fileToLoad = [folder FILES{iFile} '.mat'];
  load(fileToLoad, 'totalCosts');
  load(fileToLoad, 'GAMMA_ARRAY'); 
  
  plot(GAMMA_ARRAY, totalCosts, lines{iFile},'linewidth',lineWidth,'Color', colors{iFile});
 
  xlim([0 10]);
  legend(strLegends,'Location','northwest','FontSize', fontLegend,'Orientation','vertical');
  ylabel(strCost,'fontname', fontName,'fontsize', fontAxis);    
  xlabel('\gamma','fontname', fontName,'fontsize', fontAxis);
  set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

  if is_printed
    figIdx=figIdx +1;
    fileNames{figIdx} = 'motivation';
    epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
  end
end

return;
%%
for i=1:length(fileNames)
    fileName = fileNames{i};
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName '.pdf'];
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end