W = 5000;  % Workload demand - comming 
I = length(W);    % Number of workload sources
J = 2;            % Number of data centers
L = zeros(I,J);   % Distributed load from I sources to J data centers.
L_big = sum(W);

alpha = 1;
Alpha = alpha * ones(1,J);
beta = 0.2
Beta = beta * ones(1,J);

gamma = 10; % data center operational cost = gamma * L_j

Gamma = gamma * ones(1,J);

D_j = ones(1,J);

global GLB_bus
%GLB_bus = [3; 13];