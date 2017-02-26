
% if strcmp(METHOD, 'DGLB')   
%    OLC_capacity = DC_flex*[-DC_DEMAND, DC_DEMAND];  
% elseif strcmp(METHOD, 'OLC')
%    OLC_capacity = [-OLC_setpoint*OLC_capacity_factor, OLC_setpoint*OLC_capacity_factor];
% else
%    error(['Running ' METHOD]);
% end
% 
% OLC_gain = 25* OLC_setpoint ; % uniform cost function 1/2*k*p^2, but for load need to rescale by base power
% 
% if strcmp(METHOD, 'GLB')
%    OLC_bus = bus(abs(bus(:,6))>4,1);    
% elseif strcmp(METHOD, 'DGLB')
%    %OLC_bus = bus(abs(bus(:,6))>5,1);
%    [temp OLC_bus ] = multipe_max( bus(:,6), DC_num);
%    d_set = ones(size(OLC_bus'))*DC_DEMAND;
% elseif strcmp(METHOD, 'OLC')
%    OLC_bus = load_bus_index;
% elseif strcmp(METHOD, 'OLC_few_loads')
%    OLC_bus = bus(abs(bus(:,6))>2,1);
% else
%    disp(['Running ' METHOD]);
% end