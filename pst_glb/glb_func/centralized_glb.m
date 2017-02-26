function [isSuccessful, d, L] = centralized_glb(p, r, big_L, D, J, Alpha, Beta, omega, cost_gain, dc_buses, non_dc_buses)
    % p : bus power
    % r : generator power
    % d : dc load
    N = length(dc_buses);
    cvx_begin
        variables d(J) L(N)
        %dual variable q_lMultipliers
        minimize(sum(cost_gain.*d + 0.5 * D.* omega.^2)+sum(Beta))        
        subject to
            L >= 0;
            sum(L)>=big_L;
            d(non_dc_buses)==0;  
            d >= 0;
            d(dc_buses)>=(Alpha'.*(L) + Beta');
%             sum(-p + r - d - D .* omega) == 0;
            %q_lMultipliers: -q_r+m-w_r-q_l <= 0;
    cvx_end
    isSuccessful = 1;
    if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
        error('This optimization is unsolved');
        isSuccessful = -1;
        cvx_status
    end   
end

