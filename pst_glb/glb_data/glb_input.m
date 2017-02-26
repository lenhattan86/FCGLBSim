if ~exist('N')
    N = 3; % Number of Data centers
    J = 40; % Number of sources
end

datacenter_location = [37.390073 122.081852;
                      47.237054	119.852629;
                      45.600287	121.143165;
                      39.791127	89.647865;
                      31.515337	82.850018;
                      37.554376	77.433815;
                      29.784641	95.364304;
                      27.976211	82.455368;
                      35.920474	81.540928;
                      33.165145	80.002899];

state_location = [32 87; 65 150; 34 111; 35 92; 37 120; ...
              39 105; 41.5 72.5; 39 75.5; 39 77; 28 81; ...
              33 83; 21 158; 44 114; 40 89; 40 86; ...
              34 81; 44 100; 36 86; 31 99; 40 112; ...
              44 73; 38 79; 47 121; 39 81; 44 89; 43 107; ...
              42 93; 38 98; 38 85; 31 92; 45 69; ...
              39 76.5; 42.5 72; 43 84; 46 94.5; 33 90; ...
              38 92; 47 110; 41 99; 39 116; 43.5 71.5; ...
              40 74.5; 34 106; 43 76; 35 79; 47 100; ...
              40 83; 35 97; 44 120; 41 78; 41.5 71.5; ...
              ];
%     N = length(datacenter_location);
%     J = length(state_location);
 closestDC = zeros(J,1);
 pi_ij = ones(N,J); % network delay
 for j = 1:1:J
     minDist = inf;
     for i = 1:1:N
         pi_ij(i,j)=sqrt((datacenter_location(i,1)-state_location(j,1))^2+(datacenter_location(i,2)-state_location(j,2))^2);
         if minDist > pi_ij(i,j)
             minDist = pi_ij(i,j);
             closestDC(j) = i;            
         end
     end
 end
 
LOAD_SCALE = 1;
dcLoads = hist(closestDC,N);
L_mean = 5000;    
% setup capacity
MU = 1;
M = zeros (N,1); % Capacity: Number of servers in a data center.  
for i=1:N
  M(i) = dcLoads(i)*L_mean/LOAD_SCALE/MU;               
end