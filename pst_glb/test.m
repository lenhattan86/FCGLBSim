load('glb_func/debug.mat');
%% algorithm parts
        omega  = measured_freq;        
        N_bus = length(lmod_sig(:,1)); % TODO: should be the number of all buses. legnth(bus(:,1)).        
        D=0.001;
        numIter = 10;
        cost_gain = 1;
        % Step 1: Initialize y
        
        y_left      = zeros(N_bus,numIter); % copy for {(j-1,j) @ j}
        y_right      = zeros(N_bus,numIter); % copy for{(j,j+1) @ j}

        
        mu_left = zeros(N_bus,numIter);
        mu_right = zeros(N_bus,numIter);
        
        L          = zeros(N_bus,1);
        d          = zeros(N_bus,1);
        v          = zeros(N_bus,1);
        
        alpha = 1;
        beta = 0;
        p_r_j = 0;
        rho = 0.1;
        
        for j=1:N_bus            
            y_left(j, 1)        = (j-1)/N_bus*big_L;                        
            y_right(j, 1)       = (j)/N_bus*big_L;
        end        
        
        for iCount=2:numIter            
            for j=1:N_bus % TODO: do it in parrallel 
                % Step 2: solve (39)
                if(j>1)
                  [y_j isSuccessful] = solve_cost_eq( omega(j), mu_left(j, iCount-1), rho, ...
                            y_right(j, iCount-1), ...
                            y_right(j-1, iCount-1), ...
                            cost_gain, alpha, beta, 1);                   
                  y_left(j, iCount)    = y_j;
                  y_right(j-1, iCount) = y_j;
                else
                  y_left(N_bus, iCount)  = 0;
                end
                % y_left(j, iCount) = y_j; %??
                
                % Step 3: solve (40)
                if(j<N_bus)
                  [y_j isSuccessful] = solve_cost_eq( omega(j), mu_right(j, iCount-1), rho, ...
                           y_left(j, iCount), y_left(j+1, iCount), ...
                            cost_gain, alpha, beta, -1);                  
                  y_right(j, iCount)  = y_j;                
                  y_left(j+1, iCount) = y_j;
                else
                  y_right(N_bus, iCount)  = big_L;  
                end
                
                % Step 4: Dual update
                if(j>1)
                  mu_left(j, iCount)  = mu_left(j, iCount-1) + rho*(y_right(j-1, iCount) - y_left(j, iCount));
                end
                if(j<N_bus)
                  mu_right(j, iCount) = mu_right(j, iCount-1) + rho*(y_right(j, iCount)  - y_left(j+1, iCount));                
                end
                % Stop condition.
            end
        end    
        % step 5: obtain the final solution:
        for j=1:N_bus
            L(j) = y_right(j,iCount)-y_left(j,iCount);
            d(j) = alpha*L(j)+beta;
        end           
        
        y_left