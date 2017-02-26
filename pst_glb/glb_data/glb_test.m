load debug.mat;
J = 3;
p = bus(1:J,6); % where can we get these values?
r = bus(1:J,4);
dc_buses = [1:2];
dc_num = length(dc_buses);
non_dc_buses = [3:J];

Alpha = 1 * ones(J,1);
Beta = 0 * ones(J,1);

cost_gain = 1 * ones(J,1);
D = -0.001*ones(J,1);

big_L = 3;
omega = load_freq(:, k);
omega = omega(1:J);

[isSuccessful, d, L] = centralized_glb(p, r, big_L, D, J, Alpha, Beta, omega, cost_gain, dc_buses, non_dc_buses);