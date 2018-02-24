clear; close all; clc;
addpath('glb_data');
addpath('glb_func');
figure_settings;
warning off;
%%
is_printed = true;
figIdx = 0;
delta_frequency = 10^-6;
figSize = figOneCol;
 PLOTS = [false true false false false false];
folder = 'output/';
strLegends = {strOLC, strProposed, strNone, strOptimal};
  lines = { lineOLC, lineProposed, lineNone, lineOptimal};
  colors = {colorOLC, colorProposed, colorNone, colorOptimal};
DEBUG_PLOTS = [false false];  
PLOTS = [true false];  
%%

if DEBUG_PLOTS(1)      
    TIME_STEP = 0.1;
%     DELAYS = [0 0.1 0.2 0.4 0.8 1.0 2.0 3.0 4.0 5.0];
    DELAYS = [0 0.1 0.2 1.0 3.0 5.0];
%     lambdas = {zeros(1,10), 0.0001*ones(1,10)};
    lambdas =   {zeros(1,10), 0.002*ones(1,10), 0.001*ones(1,10), 5e-04*ones(1,10), 0.0002*ones(1,10), 0.0001*ones(1,10), 5e-05*ones(1,10), 1e-06*ones(1,10)};
    legendStr = {'\beta=0', num2str(0.002),num2str(0.001), num2str(5e-04), num2str(0.0002), num2str(0.0001),num2str(5e-05),num2str(1e-06)};
    INCIDENT_START = 5;
    conv_times = zeros(length(lambdas),length(DELAYS));
    figure;
    for iLambda = 1:length(lambdas)
      for iDelay = 1:length(DELAYS)
        fileToLoad = sprintf('%sGenLoss_proposed_0.4_%s_75_0.16_%s_%s.mat', folder, num2str(TIME_STEP), num2str(DELAYS(iDelay)), num2str(lambdas{iLambda}(iDelay)));
%         fileToLoad
        [conv_times(iLambda,iDelay), iConverge, converged_freq] = computeConvergentTime( fileToLoad, delta_frequency );
      end
      plot(DELAYS,conv_times(iLambda,:), 'linewidth',lineWidth);
      hold on;
    end
    legend(legendStr);
    xlim([0 max(DELAYS)]);
%     ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('delay (secs)','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);        
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end    
end
%%
if PLOTS(1)
  TIME_STEP = 0.1;
  figSize = figTwoThirdCol;
%   DELAYS =       [0     0.1   0.2   1.0     3.0   5.0];
%   best_lambdas = [0.002 0.002 0.001 0.0002  0.0001  0.0001];   
%   DELAYS =       [0     0.1   0.2   1.0];
%   best_lambdas = [0.002 0.002 0.001 0.0002];     
   DELAYS =       [0     0.5   1.0];
   best_lambdas = [0.002 0.0005 0.0002];     
  
  legendStr = cell(1,length(DELAYS));
  for i=1:length(DELAYS)
    if DELAYS(i)==0 
      legendStr{i} = sprintf('No delay',num2str(DELAYS(i)));
    else
      legendStr{i} = sprintf('\\Delta=%s second', num2str(DELAYS(i)));
    end
  end
  for iDelay = 1:length(DELAYS)    
    fileToLoad = sprintf('%sGenLoss_proposed_0.4_%s_75_0.16_%s_%s.mat', folder, num2str(TIME_STEP), ...
      num2str(DELAYS(iDelay)), num2str(best_lambdas(iDelay)));
%     fileToLoad
    [conv_times(iDelay), iConverge, converged_freq] = computeConvergentTime(fileToLoad, delta_frequency );   
    load(fileToLoad);
    plot(t(501:length(t))-5,load_freq(1,501:length(t))*60, 'LineWidth',lineWidth);
    hold on;
  end  
  conv_times = conv_times - 5;
  legend(legendStr,'Location','best','FontSize', fontLegend,'Orientation','vertial');
  xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
  ylabel('Frequency (Hz)','fontname', 'Arial','fontsize',fontAxis);
  ylim([59.96 60]);
  xlim([0 100]);
  set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
  
  if is_printed
    figIdx=figIdx +1;
    fileNames{figIdx} = 'delay_freq';
    epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
  end      
end
%%
if DEBUG_PLOTS(2)
      fileToLoad = [folder FILES{1} '.mat'];
    load(fileToLoad);
    plot(t,load_freq(1,:)*60, 'LineWidth',lineWidth);
    xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('Frequency (Hz)','fontname', 'Arial','fontsize',fontAxis)    
    figSize = figOneCol;
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'freq_conv';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end 
end
%% 

for iFile=1:length(fileNames)
    fileName = fileNames{iFile};
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end