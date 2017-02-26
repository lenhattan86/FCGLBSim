clear; close all; clc;
addpath('glb_data');
addpath('glb_func');
figure_settings;
%%
is_printed = true;
figIdx = 0;
delta_frequency = 0.0005;
figSize = [0.0 0 5.0 3.0];
PLOTS = [false false true];
folder = 'output/';
%% convergence time.


if PLOTS(1)
    
    FILES = {'LoadChange_proposed_0.1_0.01'...
            ,'LoadChange_proposed_0.2_0.01'...
            ,'LoadChange_proposed_0.3_0.01'...
            ,'LoadChange_proposed_0.4_0.01'...
            ,'LoadChange_proposed_0.5_0.01'...        
            ,'LoadChange_proposed_0.6_0.01'...
            ,'LoadChange_proposed_0.7_0.01'...
            ,'LoadChange_proposed_0.8_0.01'...
            ,'LoadChange_proposed_0.9_0.01'...
            ,'LoadChange_proposed_1_0.01'...
            };
    
    FLEXES = [10 20 30 40 50 60 70 80 90 100];
    range = 1:length(FLEXES);
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
          for i=start:length(load_freq)      
            freq_lag = abs(load_freq(:,i)-convergedVal);
            freq_dev_lag = abs(freq_dev(:,i) - freq_dev(:,i-1));
            delta = max(freq_lag, freq_dev_lag);
            if (max(delta) < delta_frequency)
                iConverge = i;   
                break;
            end
          end
          converge_time(iFile) = t(iConverge);
          converged_freq(iFile) = convergedVal*60;
          
              
          d_j = controlled_load(:,length(controlled_load(1,:)));
          costs(iFile) = fcp_alpha/2*(sum(a.*d_j')).^2 + (ones(size(OLC_gain))./OLC_gain)' * (d_j.^2)/2;    
        end
    end
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figure
    plot(FLEXES(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(FLEXES(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
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
    bar(FLEXES(range), 60 - converged_freq(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
    xlim([-5, max(FLEXES(range))+5]);    
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
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
    bar(FLEXES(range), costs(range), barWidth);
%     xlim([min(FLEXES(range))-5, max(FLEXES(range))+5]);
    xlim([-5, max(FLEXES(range))+5]); 
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
    ylabel('cost','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_costs';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end  
end
%%
if PLOTS(2)    
    FILES = {'LoadChange_proposed_0.1_0.01'...
            ,'LoadChange_proposed_0.2_0.02'...
            ,'LoadChange_proposed_0.3_0.04'...
            ,'LoadChange_proposed_0.4_0.08'...
            ,'LoadChange_proposed_0.5_0.16'...
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
          iConverge = start;

          convergedVal = min(min(load_freq));
          for i=start:length(load_freq)      
            freq_lag = abs(load_freq(:,i)-convergedVal);
            freq_dev_lag = abs(freq_dev(:,i) - freq_dev(:,i-1));
            delta = max(freq_lag, freq_dev_lag);
            if (max(delta) < delta_frequency)
                iConverge = i;   
                break;
            end
          end
          converge_time(iFile) = t(iConverge);
          converged_freq(iFile) = convergedVal*60;
        end
    end
    deta = delta_frequency*60;
    converge_time = converge_time - INCIDENT_START;
    
    figure
    plot(DELAYS(range),converge_time(range),'b-x', 'linewidth',lineWidth);
    xlim([0 max(DELAYS(range))]);
    ylim([0 max(ceil(converge_time(range)/5)*5)]);    
    xlabel('flexibility (%)','fontname', fontName,'fontsize',fontAxis);
    ylabel('Convergence time (secs)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_converge_times';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
    %
    figure
    bar(DELAYS(range), 60 - converged_freq(range), barWidth);
%     xlim([min(DELAYS(range))-5, max(DELAYS(range))+5]);
    xlim([-5, max(DELAYS(range))+5]);
    
    xlabel('delay (secs)','fontname', fontName,'fontsize',fontAxis);
    ylabel('frequency deviation (Hz)','fontname', fontName,'fontsize',fontAxis);
    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'delay_converged_freq_dev';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end   
    
end
%% importance of frequency deviation.
if PLOTS(3)
    
    FILES = {'LoadChange_OLC_0.5_0.01_0'...
            ,'LoadChange_OLC_0.5_0.01_0.2'...
            ,'LoadChange_OLC_0.5_0.01_0.4'...
            ,'LoadChange_OLC_0.5_0.01_0.6'...
            ,'LoadChange_OLC_0.5_0.01_0.8'...
            ,'LoadChange_OLC_0.5_0.01_1'...            
            ,'LoadChange_OLC_0.5_0.01_2'...
            ,'LoadChange_OLC_0.5_0.01_4'...
            ,'LoadChange_OLC_0.5_0.01_8'...            
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
          for i=start:length(load_freq)      
            freq_lag = abs(load_freq(:,i)-convergedVal);
            freq_dev_lag = abs(freq_dev(:,i) - freq_dev(:,i-1));
            delta = max(freq_lag, freq_dev_lag);
            if (max(delta) < delta_frequency)
                iConverge = i;   
                break;
            end
          end
          converge_time(iFile) = t(iConverge);
          converged_freq(iFile) = convergedVal*60;
          
              
          d_j = controlled_load(:,length(controlled_load(1,:)));
          costs(iFile) = fcp_alpha/2*(sum(a.*d_j')).^2 + (ones(size(OLC_gain))./OLC_gain)' * (d_j.^2)/2;    
        end
    end
    costs = costs.*WEIGHT_ARRAY;
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
    ylabel('cost','fontname', fontName,'fontsize',fontAxis);    
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);    
    
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'flex_costs';
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
