function [ converge_time, converged_freq, costs ] = obtainResult(folder, FILES, delta_frequency)
%OBTAINRESULT Summary of this function goes here
%   Detailed explanation goes here
converge_time = zeros(1,length(FILES));
converged_freq = zeros(1,length(FILES));
costs = zeros(1,length(FILES));
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
        converge_time(iFile) = t(iConverge) - INCIDENT_START;
        converged_freq(iFile) = convergedVal*60;  
        
        d_j = controlled_load(:,length(controlled_load(1,:)));
        costs(iFile) = (fcp_gamma)/2*(sum(a.*d_j)).^2 + sum((c/2).* (d_j.^2)/2);
        costs(iFile) = costs(iFile)*(NEW_ENG_BASE^2);
      end
end
 end

