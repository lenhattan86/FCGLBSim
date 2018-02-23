close all;
addpath('glb_data');
addpath('glb_func');
addpath('glb_classes');

common_setttings
figure_settings
fig_path = '../../FCGLB/full_paper_3/images/'
figSize = [0.0 0 5.0 3.0];
extra='_negative2';
figIdx = 0;
folder = 'output/';
filename = 'GenLoss_proposed_0.4_0.01_75_0.16_0_0.001.mat';
filePath = [folder filename];

%% plot the frequencies of generators
plots = [false true true true];

if plots(1)
    figure; 
    load(['output/DGLB' int2str(100*dc_flex) '.mat']);
    h3 = plot(t,mac_spd(1,:)*60,'b--');  hold on;
    plot(t,mac_spd*60,'b--');
    hold on;

%     load('output/OLC.mat');
%     h4= plot(t,mac_spd(1,:)*60,'k-');  hold on;
%     plot(t,mac_spd*60,'k-');  

    %legend([h1 h2 h3], {'PST', 'GLB', 'OLC'});
    % legend([h1 h2 h3 h4], {'PST', 'GLB','DGLB', 'OLC'});
    % legend([h1 h3 h4], {'PST', 'DGLB', 'OLC'});
%     legend([h3 h4], {'DGLB', 'OLC'});
    legend([h3], {'DGLB'});

    title('Turbine Machine Speed (Hz)','fontname', 'Arial','fontsize',fontAxis);
    xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('Frequency (Hz)','fontname', 'Arial','fontsize',fontAxis);
    % ylim([59.90, 60]); xlim([0, 50]);
    %ylim([59, 60.2]); xlim([0, 30]);

    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'gen_freq';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end
%% plot the frequencies of loads
if plots(2)
    figure; 

    load(filePath);
    h3 = plot(t,load_freq(1,:)*60,'b--');  hold on;
    plot(t,load_freq*60,'b--');
    hold on;

%     load('output/OLC.mat');
%     h4= plot(t,load_freq(1,:)*60,'k-');  hold on;
%     plot(t,load_freq*60,'k-');

    %legend([h1 h2 h3], {'PST', 'GLB', 'OLC'});
    % legend([h1 h2 h3 h4], {'PST', 'GLB','DGLB', 'OLC'});
    % legend([h1 h3 h4], {'PST', 'DGLB', 'OLC'});
%     legend([h3 h4], {'DGLB', 'OLC'});
    legend([h3], {'DGLB'});

    % title('Frequency at load buses (Hz)','fontname', 'Arial','fontsize',fontAxis);
    xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('Frequency (Hz)','fontname', 'Arial','fontsize',fontAxis);
    % ylim([59.90, 60.0]); xlim([0, 50]);


    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'load_freq';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end
%%
% figure; 
% load(['output/DGLB' int2str(100*dc_flex) '.mat']);
% h3 = plot(t,Q_glb*BASE_POWER,'b-');
% title('Queue Size (in terms of power)','fontname', 'Arial','fontsize',fontAxis);
% xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
% ylabel('MW','fontname', 'Arial','fontsize',fontAxis);

%%
if plots(3)
    
    % figure; 
    % load(['output/DGLB' int2str(100*dc_flex) '.mat']);
    % h3 = plot(t, L*BASE_POWER,'b-');
    % % title('changed L (in terms of power)','fontname', 'Arial','fontsize',fontAxis);
    % xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    % ylabel('MW','fontname', 'Arial','fontsize',fontAxis);

    figure; 
    load(['output/DGLB' int2str(100*dc_flex) '.mat']);
    h3 = plot(t, d*BASE_POWER);
    legend({'DC 1','DC 2','DC 3','DC 4','DC 5','DC 6'});
%     ylim([0 DC_CAPACITY*BASE_POWER]);
    % title('d (in terms of power)','fontname', 'Arial','fontsize',fontAxis);
    xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('MW','fontname', 'Arial','fontsize',fontAxis);

    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'dc_power';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
        print ('-depsc', epsFile);
    end
end
%% Plot each bus comparing frequency, Q, and power consumption
if plots(4)
    load(['output/DGLB' int2str(100*dc_flex) '.mat']);
    k0 = 90; kF = 95;
    %k0 = 135; kF = 145;
    % for dc = 1:6
    %     figure;
    %     subplot(3,1,1);
    %     plot(t(k0:kF), d(dc,k0:kF)*BASE_POWER);
    %     ylabel('Consump (MW)','fontname', 'Arial','fontsize',fontAxis);
    %     title(['Data Center ', num2str(dc)]);
    %     subplot(3,1,2);
    %     plot(t(k0:kF),load_freq(dc,k0:kF)*60);
    %     ylabel('Frq (Hz)','fontname', 'Arial','fontsize',fontAxis);
    %     subplot(3,1,3);
    %     plot(t(k0:kF),Q_glb(1,k0:kF)*BASE_POWER)
    %     ylabel('Q (MW)','fontname', 'Arial','fontsize',fontAxis);
    %     xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    %     print ('-depsc', [ LOCAL_FIG 'DC_num', num2str(dc), '.eps']);
    % end

    % figure;
    % plot(t,load_freq*60);
    % ylabel('Frq (Hz)','fontname', 'Arial','fontsize',fontAxis);
    % xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    % legend({'DC 1','DC 2','DC 3','DC 4','DC 5','DC 6'});
    % print ('-depsc', [LOCAL_FIG 'AllDC_freq.eps']);

%     figure;
%     plot(t,d*BASE_POWER);
%     ylabel('Consump (MW)','fontname', 'Arial','fontsize',fontAxis);
%     xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
%     legend({'DC 1','DC 2','DC 3','DC 4','DC 5','DC 6'});
%     set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
%     if is_printed
%       figIdx=figIdx +1;
%       fileNames{figIdx} = 'AllDC_Power';
%       epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
%       print ('-depsc', epsFile);
%     end

    figure;
    plot(t,Q_glb*BASE_POWER);
    ylabel('q (MW)','fontname', 'Arial','fontsize',fontAxis);
    xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    if is_printed
      figIdx=figIdx +1;
      fileNames{figIdx} = 'AllDC_Q';
      epsFile = [ LOCAL_FIG fileNames{figIdx} '.eps'];
      print ('-depsc', epsFile);
    end
end
%%
return;
%%
for i=1:length(fileNames)
    fileName = fileNames{i};
    epsFile = [ LOCAL_FIG fileName '.eps'];
    pdfFile = [ fig_path fileName extra '.pdf']    
    cmd = sprintf(PS_CMD_FORMAT, epsFile, pdfFile);
    status = system(cmd);
end