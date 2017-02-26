final_delta = controlled_load(:,length(controlled_load(1,:)));
%% sum(a * delta) - s = 0
mu_final = mu(length(mu));
s = mu_final*(1/fcp_alpha);
demand = sum(a'.*final_delta);
s_equality_constraint = demand - s
%% g'(s) - mu = 0
s = sum(a'.*final_delta);
s_equality_constraint = fcp_alpha*s - mu_final

%% sum(delta_j)
delta_sum = sum(final_delta);
delta_sum
%%
sumOfD_j = (sum(disturbance_size)+delta_sum)/(1-mean(load_freq(:,length(load_freq(1,:)))));
