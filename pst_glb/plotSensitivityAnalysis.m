clear; close all; clc;
addpath('glb_data');
addpath('glb_func');
figure_settings;
warning off;
%%
is_printed = true;
figIdx = 0;
delta_frequency = 1/100;
figSize = figTwoThirdCol;
% PLOTS = [true false false true false];
 PLOTS = [false false false false false]; PLOTS(1) = true;
folder = 'output/';
strLegends = {strOLC, strProposed, strNone, strOptimal};
  lines = { lineOLC, lineProposed, lineNone, lineOptimal};
  colors = {colorOLC, colorProposed, colorNone, colorOptimal};
%% demand flexibility
if PLOTS(1)    
    FILES = {%'GenLoss_OLC_0.1_0.01_75_0.16'...
            %,'GenLoss_OLC_0.15_0.01_75_0.16'...
            %,'GenLoss_OLC_0.2_0.01_75_0.16'...
            'GenLoss_OLC_0.21_0.01_75_0.16'...
            ,'GenLoss_OLC_0.22_0.01_75_0.16'...
            ,'GenLoss_OLC_0.23_0.01_75_0.16'...
            ,'GenLoss_OLC_0.24_0.01_75_0.16'...
            ,'GenLoss_OLC_0.25_0.01_75_0.16'...
            ,'GenLoss_OLC_0.26_0.01_75_0.16'...
            ,'GenLoss_OLC_0.27_0.01_75_0.16'...
            ,'GenLoss_OLC_0.3_0.01_75_0.16'...
            ,'GenLoss_OLC_0.35_0.01_75_0.16'...
            ,'GenLoss_OLC_0.4_0.01_75_0.16'...
            ,'GenLoss_OLC_0.45_0.01_75_0.16'...
            ,'GenLoss_OLC_0.5_0.01_75_0.16'...
            ,'GenLoss_OLC_0.55_0.01_75_0.16'...
            ,'GenLoss_OLC_0.6_0.01_75_0.16'...
            ,'GenLoss_OLC_0.65_0.01_75_0.16'...
            ,'GenLoss_OLC_0.7_0.01_75_0.16'...
            ,'GenLoss_OLC_0.75_0.01_75_0.16'...
            ,'GenLoss_OLC_0.8_0.01_75_0.16'...
            ,'GenLoss_OLC_0.85_0.01_75_0.16'...
            ,'GenLoss_OLC_0.9_0.01_75_0.16'...
            ,'GenLoss_OLC_0.95_0.01_75_0.16'...
            ,'GenLoss_OLC_1_0.01_75_0.16'...
            };
          
    FILES_2 = {%'GenLoss_proposed_0.1_0.01_75_0.16'...
            %,'GenLoss_proposed_0.15_0.01_75_0.16'...
            %,'GenLoss_proposed_0.2_0.01_75_0.16'...
            'GenLoss_proposed_0.21_0.01_75_0.16'...
            ,'GenLoss_proposed_0.22_0.01_75_0.16'...
            ,'GenLoss_proposed_0.23_0.01_75_0.16'...
            ,'GenLoss_proposed_0.24_0.01_75_0.16'...
            ,'GenLoss_proposed_0.25_0.01_75_0.16'...
            ,'GenLoss_proposed_0.26_0.01_75_0.16'...
            ,'GenLoss_proposed_0.27_0.01_75_0.16'...
            ,'GenLoss_proposed_0.3_0.01_75_0.16'...
            ,'GenLoss_proposed_0.35_0.01_75_0.16'...
            ,'GenLoss_proposed_0.4_0.01_75_0.16'...
            ,'GenLoss_proposed_0.45_0.01_75_0.16'...
            ,'GenLoss_proposed_0.5_0.01_75_0.16'...
            ,'GenLoss_proposed_0.55_0.01_75_0.16'...
            ,'GenLoss_proposed_0.6_0.01_75_0.16'...
            ,'GenLoss_proposed_0.65_0.01_75_0.16'...
            ,'GenLoss_proposed_0.7_0.01_75_0.16'...
            ,'GenLoss_proposed_0.75_0.01_75_0.16'...
            ,'GenLoss_proposed_0.8_0.01_75_0.16'...
            ,'GenLoss_proposed_0.85_0.01_75_0.16'...
            ,'GenLoss_proposed_0.9_0.01_75_0.16'...
            ,'GenLoss_proposed_0.95_0.01_75_0.16'...
            ,'GenLoss_proposed_1_0.01_75_0.16'...
            };
    
    FLEXES = 100*[(0.21:0.01:0.27)  (0.3:0.05:1.0)];
    range = 1:length(FLEXES);
    converge_time = zeros(2,length(FILES));
    converged_freq = zeros(2,length(FILES));
    costs = zeros(2,length(FILES));   
    [converge_time(1,:), converged_freq(1,:), costs(1,:)] = obtainResult(folder, FILES, delta_frequency);
    [converge_time(2,:), converged_freq(2,:), costs(2,:)] = obtainResult(folder, FILES_2, delta_frequency);
    
    figure
    for iFig=1:2
      plot(FLEXES(range),converge_time(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
      hold on;
    end
    
    xlim([FLEXES(1) max(FLEXES(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);    
    legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
    
    % frequency
    figure
    for iFig=1:2
      plot(FLEXES(range),60 - converged_freq(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
      hold on;
    end

    xlim([FLEXES(1) max(FLEXES(range))]);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(FLEXES(range))+5]);    
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
     legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
   
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_converged_freq_dev';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end   
    
    % cost
    figure      
    for iFig=1:2
      plot(FLEXES(range),costs(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
      hold on;
    end
    xlim([FLEXES(1) max(FLEXES(range))]);
    ylim([0 150]);
%     bar(FLEXES(range), costs(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(FLEXES(range))+5]); 
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_costs';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end  
end
%% control delay
if PLOTS(2)    
    FILES = {'GenLoss_proposed_0.5_0.01_1_1'...
            ,'GenLoss_proposed_0.5_0.02_1_1'...
            ,'GenLoss_proposed_0.5_0.04_1_1'...
            ,'GenLoss_proposed_0.5_0.08_1_1'...
            ,'GenLoss_proposed_0.5_0.16_1_1'...
            };    
    DELAYS = 0.01*2.^(0:4);
    range = 1:length(DELAYS);
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
          costs(iFile) = (fcp_alpha*fcp_gamma)/2*(sum(a.*d_j')).^2 + sum((c/2).* (d_j.^2)/2);    
        end
    end
    costs = costs.*DELAYS;    
    costs = costs*(BASE_POWER^2);
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figure
    plot(DELAYS(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(DELAYS(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('delay','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
    
    % frequency
    figure
    bar(DELAYS(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]);    
    xlabel('delay','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_converged_freq_dev';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end   
    
    % cost
    figure      
    plot(DELAYS(range), costs(range));
%     xlim([min(WEIGHT_ARRAY(range))-5, max(WEIGHT_ARRAY(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]); 
    xlabel(strDelay,'fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_costs';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end     
    
end
%% importance of frequency deviation.
if PLOTS(3)    
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
    costs = costs.*WEIGHT_ARRAY;    
    costs = costs*(BASE_POWER^2);
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figure
    plot(WEIGHT_ARRAY(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(WEIGHT_ARRAY(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('weight','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
    
    % frequency
    figure
    bar(WEIGHT_ARRAY(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]);    
    xlabel('weight','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_converged_freq_dev';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end   
    
    % cost
    figure      
    plot(WEIGHT_ARRAY(range), costs(range), barWidth);
%     xlim([min(WEIGHT_ARRAY(range))-5, max(WEIGHT_ARRAY(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]); 
    xlabel('weight','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_costs';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end  
end
%% importance of interdepence cost
if PLOTS(4)    
    FILES = {%'GenLoss_OLC_0.4_0.01_75_0.005'...
            %,'GenLoss_OLC_0.4_0.01_75_0.015'...
            ,'GenLoss_OLC_0.4_0.01_75_0.03'...
            ,'GenLoss_OLC_0.4_0.01_75_0.045'...
            ,'GenLoss_OLC_0.4_0.01_75_0.06'...
            ,'GenLoss_OLC_0.4_0.01_75_0.075'...
            ,'GenLoss_OLC_0.4_0.01_75_0.09'...
            ,'GenLoss_OLC_0.4_0.01_75_0.105'...
            ,'GenLoss_OLC_0.4_0.01_75_0.12'...
            ,'GenLoss_OLC_0.4_0.01_75_0.135'...
            ,'GenLoss_OLC_0.4_0.01_75_0.15'...
            ,'GenLoss_OLC_0.4_0.01_75_0.165'...
            ,'GenLoss_OLC_0.4_0.01_75_0.18'...
            ,'GenLoss_OLC_0.4_0.01_75_0.195'...
            ,'GenLoss_OLC_0.4_0.01_75_0.21'...
            ,'GenLoss_OLC_0.4_0.01_75_0.225'...
            ,'GenLoss_OLC_0.4_0.01_75_0.24'...
            ,'GenLoss_OLC_0.4_0.01_75_0.255'...
            ,'GenLoss_OLC_0.4_0.01_75_0.27'...
            ,'GenLoss_OLC_0.4_0.01_75_0.285'...
            ,'GenLoss_OLC_0.4_0.01_75_0.3'...
            };
          
     FILES_2 = {%'GenLoss_proposed_0.4_0.01_75_0.005'...
            %,'GenLoss_proposed_0.4_0.01_75_0.015'...
            ,'GenLoss_proposed_0.4_0.01_75_0.03'...
            ,'GenLoss_proposed_0.4_0.01_75_0.045'...
            ,'GenLoss_proposed_0.4_0.01_75_0.06'...
            ,'GenLoss_proposed_0.4_0.01_75_0.075'...
            ,'GenLoss_proposed_0.4_0.01_75_0.09'...
            ,'GenLoss_proposed_0.4_0.01_75_0.105'...
            ,'GenLoss_proposed_0.4_0.01_75_0.12'...
            ,'GenLoss_proposed_0.4_0.01_75_0.135'...
            ,'GenLoss_proposed_0.4_0.01_75_0.15'...
            ,'GenLoss_proposed_0.4_0.01_75_0.165'...
            ,'GenLoss_proposed_0.4_0.01_75_0.18'...
            ,'GenLoss_proposed_0.4_0.01_75_0.195'...
            ,'GenLoss_proposed_0.4_0.01_75_0.21'...
            ,'GenLoss_proposed_0.4_0.01_75_0.225'...
            ,'GenLoss_proposed_0.4_0.01_75_0.24'...
            ,'GenLoss_proposed_0.4_0.01_75_0.255'...
            ,'GenLoss_proposed_0.4_0.01_75_0.27'...
            ,'GenLoss_proposed_0.4_0.01_75_0.285'...
            ,'GenLoss_proposed_0.4_0.01_75_0.3'...
            };
    
    GAMMA_ARRAY = [0.03:0.015:0.3];
    range = 1:length(GAMMA_ARRAY);
    converge_time = zeros(2,length(FILES));
    converged_freq = zeros(2,length(FILES));
    costs = zeros(2,length(FILES));   
    [converge_time(1,:), converged_freq(1,:), costs(1,:)] = obtainResult(folder, FILES, delta_frequency);
    [converge_time(2,:), converged_freq(2,:), costs(2,:)] = obtainResult(folder, FILES_2, delta_frequency);
    
    figure
    for iFig=1:2
      plot(GAMMA_ARRAY(range),converge_time(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
      hold on;
    end
    
    xlim([GAMMA_ARRAY(1) max(GAMMA_ARRAY(range))]);
    %ylim([0 max(ceil(converge_time(range)/5)*5)]); 
    ylim([0 max(converge_time(2,range))]); 
    xlabel('\gamma','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);    
    legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);  
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'gamma_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
    
     % frequency
    figure
    for iFig=1:2
      plot(GAMMA_ARRAY(range),60 - converged_freq(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
      hold on;
    end

    xlim([GAMMA_ARRAY(1) max(GAMMA_ARRAY(range))]);
    xlabel('\gamma','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
     legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'gamma_converged_freq_dev';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end   
    
    % cost
    figure      
    for iFig=1:2
      plot(GAMMA_ARRAY(range),costs(iFig,range),lines{iFig}, 'linewidth', lineWidth,'Color',colors{iFig});
      hold on;
    end
    xlim([GAMMA_ARRAY(1) max(GAMMA_ARRAY(range))]);
    xlabel(strGamma,'fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    legend( strLegends, 'Location','best','FontSize', fontLegend,'Orientation','vertical');
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);      
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'gamma_costs';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end  
end
%% importance of costs
if PLOTS(5)    
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
          
     FILES_2 = {'GenLoss_proposed_0.4_0.01_75_0.01'...
            ,'GenLoss_proposed_0.4_0.01_75_0.02'...
            ,'GenLoss_proposed_0.4_0.01_75_0.03'...
            ,'GenLoss_proposed_0.4_0.01_75_0.04'...
            ,'GenLoss_proposed_0.4_0.01_75_0.05'...
            ,'GenLoss_proposed_0.4_0.01_75_0.06'...
            ,'GenLoss_proposed_0.4_0.01_75_0.07'...
            ,'GenLoss_proposed_0.4_0.01_75_0.08'...
            ,'GenLoss_proposed_0.4_0.01_75_0.09'...
            ,'GenLoss_proposed_0.4_0.01_75_0.1'...
            ,'GenLoss_proposed_0.4_0.01_75_0.11'...
            ,'GenLoss_proposed_0.4_0.01_75_0.12'...
            ,'GenLoss_proposed_0.4_0.01_75_0.13'...
            ,'GenLoss_proposed_0.4_0.01_75_0.14'...
            ,'GenLoss_proposed_0.4_0.01_75_0.15'...
            ,'GenLoss_proposed_0.4_0.01_75_0.16'...
            ,'GenLoss_proposed_0.4_0.01_75_0.17'...
            ,'GenLoss_proposed_0.4_0.01_75_0.18'...
            ,'GenLoss_proposed_0.4_0.01_75_0.19'...
            ,'GenLoss_proposed_0.4_0.01_75_0.2'...
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
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figure
    plot(WEIGHT_ARRAY(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(WEIGHT_ARRAY(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('\beta','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'beta_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
    
    % frequency
    figure
    bar(WEIGHT_ARRAY(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]);    
    xlabel('\beta','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize); 
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'beta_converged_freq_dev';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end   
    
    % cost
    figure      
    plot(WEIGHT_ARRAY(range), costs(range));
%     xlim([min(WEIGHT_ARRAY(range))-5, max(WEIGHT_ARRAY(range))+5]);
%     xlim([-5, max(WEIGHT_ARRAY(range))+5]); 
    xlabel('\beta','fontname', fontName,'fontsize',fontAxis);
    ylabel(strTotalCost,'fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'beta_costs';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end  
end
%%
return
%% convert to pdf FILES
for i=1:length(fileNames)
    fileName = fileNames{i};
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName '.pdf'];    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end
