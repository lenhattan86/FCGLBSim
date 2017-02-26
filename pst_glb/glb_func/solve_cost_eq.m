function [y_j isSuccessful] = solve_cost_eq( omega, mu, rho, y_1, y_2 , ...
                            cost_gain, alpha, beta, sign)
    %syms y
    %eqn = sign*cost_gain*(y_j_right-y) - alpha*omega + mu + sign*(y_j_left-y) == 0;
    %y_j = solve(eqn,y);
    %sign*cost_gain*(y_1-y) - alpha*omega + mu + sign*rho*(y_2-y)
    
    if cost_gain == 0
        isSuccessful=0;
    else
        term = (alpha*omega - mu) - (sign*(cost_gain*y_1+rho*y_2));
        y_j = sign*term/(-rho-cost_gain);
        isSuccessful = 1;
    end
    %y_j = max(0,y_j);
    %y_j = min(y_j,100);
end

