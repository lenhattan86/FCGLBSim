function [x_j , y_j, L_j, d_j, v_j, isSuccessful] = solve_local_opt( x_prev, y_prev, lamda, ...
        mu, D_j, p_r_j, cost_gain, rho, alpha, beta, sign)
    
    cvx_begin
        variables d_j v_j L_j x_j y_j ;
        minimize(sum(cost_gain.*d_j + 0.5 * D_j.* v_j.^2) - sign * lamda*x_j ...
            - sign * mu*y_j + rho/2 *((x_prev - x_j)^2+(y_prev-y_j)^2) )
        subject to            
            p_r_j-d_j-v_j == sign * (x_prev-x_j);
            L_j == sign *(y_prev-y_j);
            d_j == alpha*L_j+beta;
            y_j >= 0;
            L_j >= 0;
            d_j >= 0;
    cvx_end
    
    isSuccessful = 1;
    if (~(or(strcmp(cvx_status,'Solved'),strcmp(cvx_status,'Inaccurate/Solved'))))        
        error('This optimization is unsolved');
        isSuccessful = -1;
        cvx_status
    end
end

