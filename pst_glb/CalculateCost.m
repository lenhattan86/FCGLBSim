function cost = CalculateCost(mac_spd, OLC_gain, governor_gain, OLC_capacity, generator_capacity,, g_to_basva_ratio );

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[L,t] = size(lmod_sig);

setpoint = zeros(L,1); setpoint(disturb_mod) = disturb_size;



gen_cost = sum((((g_to_basva_ratio*mtg_sig).^2).*repmat(governor_gain,1,t)*g_to_basva_ratio)/2,1);
cost = load_cost+ gen_cost;
end

