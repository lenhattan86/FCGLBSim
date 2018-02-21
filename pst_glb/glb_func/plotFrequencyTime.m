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
    
%     ylim([59.4 60]);
    
    figSize = figOneCol;
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);

%     if isPrintted
%       epsFile = [ LOCAL_FIG matFile '_freq' '.eps'];
%         print ('-depsc', epsFile);
%     end

    figure 
    temp = mu(1:length(t));
%     muIds = temp~=0;
%     plot(t(muIds), temp(muIds), 'LineWidth',lineWidth);
    plot(t, temp, 'LineWidth',lineWidth);
    figSize = figOneCol;
    set (gcf, 'Units', 'Inches', 'Position', figSize, 'PaperUnits', 'inches', 'PaperPosition', figSize);
    
    xlabel('time(sec)','fontname', 'Arial','fontsize',fontAxis);
    ylabel('\mu','fontname', 'Arial','fontsize',fontAxis)
    
  else
    error('file doest not exists');
  end
  success = true;
end

