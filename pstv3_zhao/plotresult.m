close all;
%clear all;
%load('load_normal.mat');
% 
% figure; hold on;
% plot(t,mac_spd_1(1,:)*60,'r:','linewidth',2.0);
% plot(t,mac_spd_2(1,:)*60,'b-','linewidth',2.0);
% xlabel('Time (sec)'); ylabel('Frequency (Hz)');


figure; hold on;
plot(t,mac_spd*60, 'k-');
xlabel('Time (sec)','fontname', 'Arial','fontsize',16);
ylabel('Frequency (Hz)','fontname', 'Arial','fontsize',16);
% ylim([59.76, 60]); 
xlim([0, 30]);
set(gca, 'fontname','Arial','fontsize', 14)

datane; test_IEEE_39;

%cost_1 = load_cost = sum((repmat(OLC_gain,1,t).*)/2,1)
cost_1 = zeros(1,length(t));

return;
load('generation_normal.mat');
plot(t,mac_spd*60,'r--');   


cost_2 = CalculateCost(mac_spd, OLC_gain, governor_gain, OLC_capacity, generator_capacity, 10);


figure; hold on;
plot(t,cost_2, 'r--', 'linewidth',2.0);
plot(t,cost_1, 'k-', 'linewidth', 1.5);
xlabel('Time (sec)','fontname', 'Arial','fontsize',16);
ylabel('Frequency (Hz)','fontmagnitudename', 'Arial','fontsize',16);
ylim([59.76, 60]); xlim([0, 30]);
set(gca, 'fontname','Arial','fontsize', 14)