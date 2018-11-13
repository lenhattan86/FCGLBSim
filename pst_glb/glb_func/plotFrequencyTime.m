function [ success ] = plotFrequencyTime( matFile, isPrintted)
  success = false;
  figure_settings;
  filePath = [matFile];
  if exist(filePath, 'file')
    load(filePath);
    figure;
    
    %h3 = plot(t,load_freq(1,:)*60, 'b-','LineWidth',lineWidth);  hold on;
    %plot(t,load_freq*60,'b-', 'LineWidth',lineWidth);
    plot(t,load_freq*60, 'LineWidth',lineWidth);
    %legend([h3], {'frequencies of load buses'});
    xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('Frequency (Hz)','fontname', 'Arial','fontsize',fontAxis)
    
    ylim([59.8 60]);
    
    figSize = figOneCol;
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

    if isPrintted
      epsFile = [ LOCAL_FIG 'test_freq' '.eps'];
        print ('-depsc', epsFile);
    end
    
    figure
    %h3 = plot(t,load_freq(1,:)*60, 'b-','LineWidth',lineWidth);  hold on;
    %plot(t,load_freq*60,'b-', 'LineWidth',lineWidth);
    plot(t,mac_spd*60, 'LineWidth',lineWidth);
     ylim([59.8 60]);
    %legend([h3], {'frequencies of load buses'});
    xlabel('Time (sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('mac speed (Hz)','fontname', 'Arial','fontsize',fontAxis)
    
%     ylim([59.4 60]);
    
    figSize = figOneCol;
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

    if isPrintted
      epsFile = [ LOCAL_FIG 'test_mac_spd' '.eps'];
        print ('-depsc', epsFile);
    end

    return
    
    
    figure 
    temp = mu(1:length(t));
%     muIds = temp~=0;
%     plot(t(muIds), temp(muIds), 'LineWidth',lineWidth);
    plot(t, temp, 'LineWidth',lineWidth);
    figSize = figOneCol;
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

    xlabel('time(sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('\mu','fontname', 'Arial','fontsize',fontAxis);
    
        
    figure 
    d_j = controlled_load(:,length(controlled_load(1,:)));
    costs(1) = WEIGHT*sum((c/2).* (d_j.^2)/2);  
    costs(2) = WEIGHT*(fcp_gamma)/2*(sum(a.*d_j)).^2;  
    temp = mu(1:length(t));
%     muIds = temp~=0;
%     plot(t(muIds), temp(muIds), 'LineWidth',lineWidth);
    hBar = bar(costs, 0.2, 'stacked');
    figSize = figOneCol;
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
%     xlabel('time(sec)','fontname', 'Arial','fontsize',fontAxis);
%     ylabel('\mu','fontname', 'Arial','fontsize',fontAxis);
    
  else
    error('file doest not exists');
  end
  success = true;
end

