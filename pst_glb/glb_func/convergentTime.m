freq_dev = 1 - load_freq(:,:);
delta = inf;
start = INCIDENT_START/0.01 + 100;
iConverge = length(load_freq);

    convergedVal = min(min(load_freq));
% convergedVal = load_freq(length(load_freq));
convergedLag = 1 - convergedVal;
for i=length(load_freq):-1:start
  freq_lag = abs(load_freq(:,i)-convergedVal);
  freq_dev_lag = abs(freq_dev(:,i) - freq_dev(:,i-1));
  delta = max(freq_lag, freq_dev_lag);
  if (max(freq_lag) >= delta_frequency)
      iConverge = i;
      break;
  end
end
conv_time = t(iConverge);
conv_idx =   iConverge;
conv_freq = convergedVal;