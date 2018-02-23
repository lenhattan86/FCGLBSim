function [conv_time, conv_idx, conv_freq] = computeConvergentTime( fileToLoad, delta_frequency )
  if exist(fileToLoad,'file')                     
    load(fileToLoad);
    convergentTime
  else
    fprintf('%s does not exist \n',fileToLoad);
    conv_time = 0;
    conv_idx = 0;
    conv_freq = 0;
  end
end

