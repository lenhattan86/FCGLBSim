clear; close all; clc;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');
figure_settings;
warning off;
%%
is_printed = true;
figIdx = 0;
delta_frequency = 10^-5;
figSize = figTwoThirdCol;
 PLOTS = [true false false false true false true true];
folder = 'output/';
strLegends = {strOLC, strProposed, strNone, strOptimal};
lines = { lineOLC, lineProposed, lineNone, lineOptimal};
colors = {colorOLC, colorProposed, colorNone, colorOptimal};
%% demand flexibility
if PLOTS(1)    
    default_settings;
    FLEXES = [0.1:0.05:1.0];
    for iFlex = 1:length(FLEXES)
        FLEX = FLEXES(iFlex);
        extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
            num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];
        FILES{iFlex}   = ['GenLoss_OLC' extra];
        FILES_2{iFlex} = ['GenLoss_proposed_dc' extra];
    end
    FLEXES = 100*FLEXES;
    range = 1:length(FLEXES);
    converge_time = zeros(2,length(FILES));
    converged_freq = zeros(2,length(FILES));
    costs = zeros(2,length(FILES));   
    [converge_time(1,:), converged_freq(1,:), costs(1,:)] = obtainResult(folder, FILES, delta_frequency);
    [converge_time(2,:), converged_freq(2,:), costs(2,:)] = obtainResult(folder, FILES_2, delta_frequency);
    
%     figure
%     for iFig=1:2
%       plot(FLEXES(range),converge_time(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
%       hold on;
%     end
%     
%     xlim([FLEXES(1) max(FLEXES(range))]);
%     ylim([0 max(max(ceil(converge_time(range)/5)*5),1)]);    
%     xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
%     ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);    
%     legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');    
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
%     
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'flex_converge_times';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%         print ('-depsc', epsFile);
%     end
    
    % frequency
%     figure
%     for iFig=1:2
%       plot(FLEXES(range),60 - converged_freq(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
%       hold on;
%     end
% 
%     xlim([FLEXES(1) max(FLEXES(range))]);
%     xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
%     ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
%      legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
%    
%     
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'flex_converged_freq_dev';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%         print ('-depsc', epsFile);
%     end   
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'flex_costs';
    figures{figIdx} = figure; 
    
    plot(FLEXES(range),(costs(1,range)-costs(2,range))./costs(1,range)*100,lines{2}, 'linewidth', lineWidth,'Color',colors{2});
    %for iFig=1:2
    %  plot(FLEXES(range),costs(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
    %  hold on;
    %end
    xlim([FLEXES(1) max(FLEXES(range))]);
%     ylim([0 40]);
%     bar(FLEXES(range), costs(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(FLEXES(range))+5]); 
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
    ylabel('cost savings (%)','fontname', fontName,'fontsize',fontAxis);    
    %legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
end
%% control delay
if PLOTS(2)    
    DELAYS       = [0 0.1 0.2 1.0 3.0 5.0]; 
    scale_lambda = [1 1/2 1/2 1/2.5 1/5 1/10];
    for iFile = 1:length(DELAYS)
        default_settings;
        TIME_STEP=0.1;
        DELAY = DELAYS(iFile);        
        fcp_lambda = fcp_lambda*scale_lambda(iFile);
        extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
        num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];        
        FILES{iFlex} = ['GenLoss_proposed_dc' extra];
    end
    
    range = 1:length(DELAYS);
    INCIDENT_START = 5;
    converge_time = zeros(size(FILES));
    converged_freq = zeros(size(FILES));
    
    for iFile=1:length(FILES)
       fileToLoad = [folder FILES{iFile} '.mat'];
       [converge_time(iFile), iConverge, converged_freq(iFile) ] = computeConvergentTime( fileToLoad, delta_frequency );
    end
%     costs = costs.*DELAYS;    
%     costs = costs*(BASE_POWER^2);
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figIdx=figIdx +1;
    fileNames{figIdx} = 'delay_converge_times';
    figures{figIdx} = figure; 
    
    plot(DELAYS(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(DELAYS(range))]);
%     ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('delay (secs)','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    
    % frequency
    figIdx=figIdx +1;
    fileNames{figIdx} = 'delay_converged_freq_dev';
    figures{figIdx} = figure; 
    
    bar(DELAYS(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]);    
    xlabel('delay (secs)','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    
    % cost
%     figure      
%     plot(DELAYS(range), costs(range));
%     xlim([min(WEIGHT_ARRAY(range))-5, max(WEIGHT_ARRAY(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]); 
%     xlabel(strDelay,'fontname', fontName,'fontsize',fontAxis);
%     ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'delay_costs';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%         print ('-depsc', epsFile);
%     end     
%     

%     fileToLoad = [folder FILES{1} '.mat'];
%     load(fileToLoad);
%     plot(t,load_freq*60, 'LineWidth',lineWidth);
%     xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
%     ylabel('Frequency (Hz)','fontname', 'Arial','fontsize',fontAxis)    
%     figSize = figOneCol;
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
%     
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'freq_conv';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%         print ('-depsc', epsFile);
%     end 
    
end
%% step size
if PLOTS(3)    
    FILES = {'GenLoss_proposed_dc_0.4_0.01_75_0.16_0'...
            ,'GenLoss_proposed_dc_0.4_0.02_75_0.16_0'...
            ,'GenLoss_proposed_dc_0.4_0.04_75_0.16_0'...
            ,'GenLoss_proposed_dc_0.4_0.08_75_0.16_0'...
            ,'GenLoss_proposed_dc_0.4_0.16_75_0.16_0'...
            };    
    TIME_STEPS = 0.01*2.^(0:4);
    range = 1:length(TIME_STEPS);
    converge_time = zeros(size(FILES));
    converged_freq = zeros(size(FILES));
    
    for iFile=1:length(FILES)
       fileToLoad = [folder FILES{iFile} '.mat'];
        if exist(fileToLoad,'file')                     
          load(fileToLoad);
          freq_dev = 1 - load_freq(:,:);
          delta = inf;
          start = INCIDENT_START/0.01 + 100;
          iConverge = length(load_freq);

          convergedVal = min(min(load_freq));
          convergedLag = 1 - convergedVal;
          for i=start:length(load_freq)      
            freq_lag = abs(load_freq(:,i)-convergedVal);
            freq_dev_lag = abs(freq_dev(:,i) - freq_dev(:,i-1));
            delta = max(freq_lag, freq_dev_lag);
            if (max(delta/convergedLag) < delta_frequency)
                iConverge = i;
                break;
            end
          end
          converge_time(iFile) = t(iConverge);
          converged_freq(iFile) = convergedVal*60;         
              
          d_j = controlled_load(:,length(controlled_load(1,:)));
%           costs(iFile) = (fcp_alpha*fcp_gamma)/2*(sum(a.*d_j)).^2 + sum((c/2).* (d_j.^2)/2);    
        end
    end
%     costs = costs.*DELAYS;    
%     costs = costs*(BASE_POWER^2);
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
        figIdx=figIdx +1;
    fileNames{figIdx} = 'delay_converge_times';
    figures{figIdx} = figure; 
    
    plot(TIME_STEPS(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(TIME_STEPS(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('time step (secs)','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
    
    % frequency
     figIdx=figIdx +1;
    fileNames{figIdx} = 'delay_converged_freq_dev';
    figures{figIdx} = figure; 
    bar(TIME_STEPS(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]);    
    xlabel('time step (secs)','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_converged_freq_dev';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end   
    
    % cost
%     figIdx=figIdx +1;
%     fileNames{figIdx} = 'delay_costs';
%     figures{figIdx} = figure;       
%     plot(DELAYS(range), costs(range));
%     xlim([min(WEIGHT_ARRAY(range))-5, max(WEIGHT_ARRAY(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]); 
%     xlabel(strDelay,'fontname', fontName,'fontsize',fontAxis);
%     ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    

end
%% importance of frequency deviation.
if PLOTS(4)    
    FILES = {'GenLoss_OLC_0.5_0.01_0'...
            ,'GenLoss_OLC_0.5_0.01_0.2'...
            ,'GenLoss_OLC_0.5_0.01_0.4'...
            ,'GenLoss_OLC_0.5_0.01_0.6'...
            ,'GenLoss_OLC_0.5_0.01_0.8'...
            ,'GenLoss_OLC_0.5_0.01_1'...            
            ,'GenLoss_OLC_0.5_0.01_2'...
            ,'GenLoss_OLC_0.5_0.01_4'...
            ,'GenLoss_OLC_0.5_0.01_8'...            
            };
    
    WEIGHT_ARRAY = [0:0.2:1 2.^(1:3)];
    range = 1:length(WEIGHT_ARRAY);
    converge_time = zeros(size(FILES));
    converged_freq = zeros(size(FILES));
    
    for iFile=1:length(FILES)
       fileToLoad = [folder FILES{iFile} '.mat'];
        if exist(fileToLoad,'file')                     
          load(fileToLoad);
          freq_dev = 1 - load_freq(:,:);
          delta = inf;
          start = INCIDENT_START/0.01 + 100;
          iConverge = start;

          convergedVal = min(min(load_freq));
          convergedLag = 1 - convergedVal;
          for i=length(load_freq):-1:start      
            freq_lag = abs(load_freq(:,i)-convergedVal);
            freq_dev_lag = abs(freq_dev(:,i) - freq_dev(:,i-1));
            delta = max(freq_lag, freq_dev_lag);
            if (max(delta/convergedLag) >= delta_frequency)
                iConverge = i;
                break;
            end
          end
          converge_time(iFile) = t(iConverge);
          converged_freq(iFile) = convergedVal*60;         
              
          d_j = controlled_load(:,length(controlled_load(1,:)));
          costs(iFile) = (fcp_alpha*fcp_gamma)/2*(sum(a.*d_j')).^2 + sum((c/2).* (d_j.^2)/2);    
        end
    end
    costs = costs.*WEIGHT_ARRAY;    
    costs = costs*(BASE_POWER^2);
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figIdx=figIdx +1;
    fileNames{figIdx} = 'flex_converge_times';
    figures{figIdx} = figure; 
    
    plot(WEIGHT_ARRAY(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(WEIGHT_ARRAY(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('weight','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    

    
    % frequency
    figIdx=figIdx +1;
    fileNames{figIdx} = 'flex_converged_freq_dev';
    figures{figIdx} = figure; 
    
    bar(WEIGHT_ARRAY(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]);    
    xlabel('weight','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
 
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'flex_costs';
    figures{figIdx} = figure; 
    
    plot(WEIGHT_ARRAY(range), costs(range), barWidth);
%     xlim([min(WEIGHT_ARRAY(range))-5, max(WEIGHT_ARRAY(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]); 
    xlabel('weight','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    

end
%% importance of interdepence cost
if PLOTS(5)       
    
    default_settings;
    GAMMA_ARRAY = [0.03:0.015:0.3];
    for iFlex = 1:length(GAMMA_ARRAY)
        GAMMA = GAMMA_ARRAY(iFlex);
        extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
            num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];
        FILES{iFlex}   = ['GenLoss_OLC' extra];
        FILES_2{iFlex} = ['GenLoss_proposed_dc' extra];
    end
    
    range = 1:length(GAMMA_ARRAY);
    converge_time = zeros(2,length(FILES));
    converged_freq = zeros(2,length(FILES));
    costs = zeros(2,length(FILES));   
    [converge_time(1,:), converged_freq(1,:), costs(1,:)] = obtainResult(folder, FILES, delta_frequency);
    [converge_time(2,:), converged_freq(2,:), costs(2,:)] = obtainResult(folder, FILES_2, delta_frequency);
    
%     figIdx=figIdx +1;
%     fileNames{figIdx} = 'gamma_converge_times';
%     figures{figIdx} = figure;
%     
%     for iFig=1:2
%       plot(GAMMA_ARRAY(range),converge_time(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
%       hold on;
%     end
%     
%     xlim([GAMMA_ARRAY(1) max(GAMMA_ARRAY(range))]);
%     %ylim([0 max(ceil(converge_time(range)/5)*5)]); 
%     ylim([0 max(converge_time(2,range))]); 
%     xlabel('\gamma','fontname', fontName,'fontsize',fontAxis);
%     ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);    
%     legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');    
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);  
    
    
     % frequency
%     figIdx=figIdx +1;
%     fileNames{figIdx} = 'gamma_converged_freq_dev';
%     figures{figIdx} = figure;
%     
%     for iFig=1:2
%       plot(GAMMA_ARRAY(range),60 - converged_freq(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
%       hold on;
%     end
% 
%     xlim([GAMMA_ARRAY(1) max(GAMMA_ARRAY(range))]);
%     xlabel('\gamma','fontname', fontName,'fontsize',fontAxis);
%     ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
%      legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'gamma_costs';
    figures{figIdx} = figure;
    
    plot(GAMMA_ARRAY(range),(costs(1,range)-costs(2,range))./costs(1,range)*100,lines{2}, 'linewidth', lineWidth,'Color',colors{2});
    %for iFig=1:2
    %  plot(GAMMA_ARRAY(range),costs(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
    %  hold on;
    %end
    xlim([GAMMA_ARRAY(1) max(GAMMA_ARRAY(range))]);
%     ylim([0 40]);
    xlabel(strGamma,'fontname', fontName,'fontsize',fontAxis);
    ylabel('cost savings (%)','fontname', fontName,'fontsize',fontAxis);    
    %legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
end
%% importance of costs
if PLOTS(6)    
    FILES = {'GenLoss_OLC_0.4_0.01_75_0.01'...
            ,'GenLoss_OLC_0.4_0.01_75_0.02'...
            ,'GenLoss_OLC_0.4_0.01_75_0.03'...
            ,'GenLoss_OLC_0.4_0.01_75_0.04'...
            ,'GenLoss_OLC_0.4_0.01_75_0.05'...
            ,'GenLoss_OLC_0.4_0.01_75_0.06'...
            ,'GenLoss_OLC_0.4_0.01_75_0.07'...
            ,'GenLoss_OLC_0.4_0.01_75_0.08'...
            ,'GenLoss_OLC_0.4_0.01_75_0.09'...
            ,'GenLoss_OLC_0.4_0.01_75_0.1'...
            ,'GenLoss_OLC_0.4_0.01_75_0.11'...
            ,'GenLoss_OLC_0.4_0.01_75_0.12'...
            ,'GenLoss_OLC_0.4_0.01_75_0.13'...
            ,'GenLoss_OLC_0.4_0.01_75_0.14'...
            ,'GenLoss_OLC_0.4_0.01_75_0.15'...
            ,'GenLoss_OLC_0.4_0.01_75_0.16'...
            ,'GenLoss_OLC_0.4_0.01_75_0.17'...
            ,'GenLoss_OLC_0.4_0.01_75_0.18'...
            ,'GenLoss_OLC_0.4_0.01_75_0.19'...
            ,'GenLoss_OLC_0.4_0.01_75_0.2'...
            };
          
     FILES_2 = {'GenLoss_proposed_dc_0.4_0.01_75_0.01'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.02'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.03'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.04'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.05'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.06'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.07'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.08'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.09'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.1'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.11'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.12'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.13'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.14'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.15'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.16'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.17'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.18'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.19'...
            ,'GenLoss_proposed_dc_0.4_0.01_75_0.2'...
            };
    
    WEIGHT_ARRAY = [0.1:0.2:1 2.^(0:3)];
    range = 2:length(WEIGHT_ARRAY);
    converge_time = zeros(size(FILES));
    converged_freq = zeros(size(FILES));
    costs = zeros(size(FILES));
    for iFile=1:length(FILES)
       fileToLoad = [folder FILES{iFile} '.mat'];
        if exist(fileToLoad,'file')                     
          load(fileToLoad);
          freq_dev = 1 - load_freq(:,:);
          delta = inf;
          start = INCIDENT_START/0.01 + 100;
          iConverge = start;

          convergedVal = min(min(load_freq));
          convergedLag = 1 - convergedVal;
          for i=start:length(load_freq)      
            freq_lag = abs(load_freq(:,i)-convergedVal);
            freq_dev_lag = abs(freq_dev(:,i) - freq_dev(:,i-1));
            delta = max(freq_lag, freq_dev_lag);
            if (max(delta/convergedLag) < delta_frequency)
                iConverge = i;
                break;
            end
          end
          converge_time(iFile) = t(iConverge);
          converged_freq(iFile) = convergedVal*60;         
              
          d_j = controlled_load(:,length(controlled_load(1,:)));
          costs(iFile) = (fcp_alpha*fcp_gamma)/2*(sum(a.*d_j')).^2 + sum((c/2).* (d_j.^2)/2);    
        end
    end
    WEIGHT_ARRAY =  [0.1:0.2:1 2.^(0:3)];
    costs = costs.*WEIGHT_ARRAY;    
    costs = costs*(BASE_POWER^2);
    deta  = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figIdx=figIdx +1;
    fileNames{figIdx} = 'beta_converge_times';
    figures{figIdx} = figure;
    
    plot(WEIGHT_ARRAY(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(WEIGHT_ARRAY(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('\beta','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    
    % frequency
    figIdx=figIdx +1;
    fileNames{figIdx} = 'beta_converged_freq_dev';
    figures{figIdx} = figure;
    
    bar(WEIGHT_ARRAY(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]);    
    xlabel('\beta','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'beta_costs';
    figures{figIdx} = figure;
    
    plot(WEIGHT_ARRAY(range), costs(range));
    xlabel('\beta','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost, 'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    

    
end

%% impact of power losses
if PLOTS(7) 
    default_settings;
    POWER_LOSSES = [0 50 100 200 300 400 500 600 700 800 900 1000];
    FILES={};
    for iFile = 1:length(POWER_LOSSES)
        POWER_LOSS = POWER_LOSSES(iFile);
        extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
            num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];
        FILES{iFile}   = ['GenLoss_proposed_dc' extra];
    end
    POWER_LOSSES = POWER_LOSSES/1000;
    range = 1:length(POWER_LOSSES);    
    costs = zeros(1,length(FILES)); 
    last_freqs = zeros(1,length(FILES));
    [temp, temp, costs(:), last_freqs] = obtainResult(folder, FILES, delta_frequency); 
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'power_freq_dev';
    figures{figIdx} = figure; 
    
    plot(POWER_LOSSES(range), 60-last_freqs(1,range),lines{2}, 'linewidth', lineWidth, 'Color',colors{2});
    xlim([POWER_LOSSES(1) max(POWER_LOSSES(range))]);
    xlabel('power loss (GW)','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'power_costs';
    figures{figIdx} = figure; 
    
    plot(POWER_LOSSES(range), costs(1,range),lines{2}, 'linewidth', lineWidth, 'Color',colors{2});
    xlim([POWER_LOSSES(1) max(POWER_LOSSES(range))]);
    xlabel('power loss (GW)','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
end

%% trade-offs of frequency & costs
if PLOTS(8) 
    default_settings;
    WEIGHTS = [1:4 5:10:155 160:40:320];
    FILES={};
    for iFile = 1:length(WEIGHTS)
        WEIGHT = WEIGHTS(iFile);
        extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
            num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];
        FILES{iFile}   = ['GenLoss_proposed_dc' extra];
    end
    
    range = 1:length(WEIGHTS);    
    costs = zeros(1,length(FILES)); 
    last_freqs = zeros(1,length(FILES));
    [temp, temp, costs(:), last_freqs] = obtainResult(folder, FILES, delta_frequency); 
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'weight_freq_dev';
    figures{figIdx} = figure;     
    strLegends={strProposed, strOLC};
    
    for iFile = 1:length(WEIGHTS)
        WEIGHT = WEIGHTS(iFile);
        extra = ['_' num2str(FLEX) '_' num2str(TIME_STEP) '_' num2str(WEIGHT) '_' ...
            num2str(GAMMA) '_' num2str(DELAY) '_' num2str(fcp_lambda) '_' num2str(POWER_LOSS)];
        FILES{iFile}   = ['GenLoss_OLC' extra];
    end
    
%     [temp, temp, temp, OLC_freq] = obtainResult(folder, {'GenLoss_OLC_0.4_0.01_75_0.16_0_0.0001_400'}, delta_frequency); 
    [temp, temp, OLC_costs, OLC_freq] = obtainResult(folder, FILES, delta_frequency); 
    
    plot(WEIGHTS(range), 60-last_freqs(1,range),lines{2}, 'linewidth', lineWidth, 'Color',colors{2});
    hold on;
    plot(WEIGHTS(range), (60-OLC_freq).*ones(size(range)), lines{1}, 'linewidth', lineWidth, 'Color',colors{1});
    legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');    
    xlim([WEIGHTS(1) max(WEIGHTS(range))]);
    ylim([0 max(60-last_freqs(1,range)) * 1.1]);
    xlabel('\alpha ($/MW-Hz)','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
    
    % cost
    figIdx=figIdx +1;
    fileNames{figIdx} = 'weight_costs';
    figures{figIdx} = figure; 
    
    plot(WEIGHTS(range), costs(1,range),lines{2}, 'linewidth', lineWidth, 'Color',colors{2});
    xlim([WEIGHTS(1) max(WEIGHTS(range))]);
    ylim([0 max(costs(1,range)) * 1.1]);
    xlabel('\alpha ($/MW-Hz)','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
    

    
    %% trade-offs
    figIdx=figIdx +1;
    fileNames{figIdx} = 'weight_tradoff';
    figures{figIdx} = figure;
    
%     strLegends={strOLC, strProposed, strOLC, strProposed};
    strLegends={strOLC, strProposed};
    yyaxis left
    plot(WEIGHTS(range), (60-OLC_freq).*ones(size(range)), lines{1}, 'linewidth', lineWidth,'Color','black');
    hold on;
    plot(WEIGHTS(range), 60-last_freqs(1,range), lines{2}, 'linewidth', lineWidth,'Color','black');
    
    hold on;
    plot(WEIGHTS(range), (60-OLC_freq).*ones(size(range)), lines{1}, 'linewidth', lineWidth);
    hold on;
    plot(WEIGHTS(range), 60-last_freqs(1,range),lines{2}, 'linewidth', lineWidth);
      
    
    xlim([WEIGHTS(1) max(WEIGHTS(range))]);
    ylim([0 max(60-last_freqs(1,range)) * 1.1]);
    xlabel('\alpha ($/MW-Hz)','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
    
    yyaxis right  
    plot(WEIGHTS(range), OLC_costs.*ones(size(range)), lines{1}, 'linewidth', lineWidth);
    hold on;
    plot(WEIGHTS(range), costs(1,range),lines{2}, 'linewidth', lineWidth);
    xlim([WEIGHTS(1) max(WEIGHTS(range))]);
%     ylim([0 max(costs(1,range)) * 1.1]);
    xlabel('\alpha ($/MW-Hz)','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
    
    legend( strLegends, 'Location','southeast','FontSize', fontLegend,'Orientation','vertical'); 
    
% trade-offs 1
%     figIdx=figIdx +1;
%     fileNames{figIdx} = 'weight_tradoff';
%     figures{figIdx} = figure;
%     
%     scatter(60-last_freqs(1,range), costs(1,range));
%     xlim([0 max(60-last_freqs(1,range)) * 1.1]);
%     ylim([0 max(costs(1,range)) * 1.1]);
%     xlabel(strFrequency,'fontname', fontName,'fontsize',fontAxis);
%     ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);   
    
end
%%
return
%% convert to pdf FILES
for i=1:length(fileNames)    
    fileName = fileNames{i};
    
    epsFile = [ LOCAL_FIG fileNames{i} '.eps'];
    print (figures{i},'-depsc', epsFile);     
    
    pdfFile = [ fig_path fileName '.pdf']; 
    pdfFile
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end
