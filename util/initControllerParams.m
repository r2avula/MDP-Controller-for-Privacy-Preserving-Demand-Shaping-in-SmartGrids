function [params] = initControllerParams(config)
h_num = 2; %H=1: away, H=2: present
if(h_num ~= 2) 
    error('Not implemented!');
end
C_HgHh = eye(h_num); % Bayesian cost assignment corresponding to detection probability

p_pu = config.powerQuantPU; % in W

priorMaxPowerConstraint = config.priorMaxPowerConstraint;
if(ischar(priorMaxPowerConstraint))
    priorMaxPowerConstraint = str2double(priorMaxPowerConstraint);
end

max_power_in_hmm = min(priorMaxPowerConstraint,config.trainingSMdata_max); % in W

x_max_pu = floor(max_power_in_hmm/p_pu);
x_min_pu = 0;
x_num = x_max_pu-x_min_pu+1;
x_offset = x_min_pu - 1;

slotIntervalInSeconds = config.slotIntervalInSeconds;
slotIntervalInHours = slotIntervalInSeconds/3600; %in h
evalStartHourIndex = 1;
evalEndHourIndex = 24;
k_num_in_day = (evalEndHourIndex-evalStartHourIndex+1)/slotIntervalInHours; %Measurement slots
if(k_num_in_day~=floor(k_num_in_day))
    error('Wrong slotIntervalInSeconds setting!');
end

timeHorizonsPerDay= config.timeHorizonsPerDay;

k_num_in_horizon = k_num_in_day/timeHorizonsPerDay;
if(k_num_in_horizon~=floor(k_num_in_horizon))
    error('Wrong timeHorizonsPerDay setting!');
end

timeHorizonSlots = zeros(timeHorizonsPerDay,2);
temp = 1;
for horizonIdx = 1:timeHorizonsPerDay
    timeHorizonSlots(horizonIdx,1) = temp;
    timeHorizonSlots(horizonIdx,2) = temp+k_num_in_horizon-1;
    temp = temp + k_num_in_horizon;
end

batteryNominalVoltage = config.batteryNominalVoltage;
cellNominalVoltage = 4; %in V
cellsInSeries = ceil(batteryNominalVoltage/cellNominalVoltage);
if(cellsInSeries~=floor(cellsInSeries))
    warning('Battery voltage is modified!');
end
batteryNominalVoltage = cellsInSeries*cellNominalVoltage;% in V
converterEfficiency = (config.converterEfficiency)/100;     
batteryRatedCapacityInAh = config.batteryRatedCapacityInAh; %in Ah
cell_SOC_high = config.cell_SOC_high;
cell_SOC_low = config.cell_SOC_low;
energyCostPer_Wh = (config.energyCostPerkWh)/1000;
capacityCostPer_Ah = (config.capacityCostPerkWh)*batteryNominalVoltage/1000;
z_cap = (batteryRatedCapacityInAh*batteryNominalVoltage); % in Wh

cell_1C_power = config.cell_1C_power; %in W
cell_1C_capacityInAh = config.cell_1C_capacityInAh; %in Ah
legsInParallel = round(batteryRatedCapacityInAh/cell_1C_capacityInAh);
d_max_ch_ess = cell_1C_power*legsInParallel*cellsInSeries/converterEfficiency;
d_max_disch_ess = -cell_1C_power*legsInParallel*cellsInSeries*converterEfficiency;

d_rated = config.converterRatedPower;
d_max_ch = min(d_rated,d_max_ch_ess);
d_max_disch = max(-d_rated,d_max_disch_ess);

d_max_ch_pu = floor(d_max_ch/p_pu);
d_max_disch_pu = ceil(d_max_disch/p_pu);
out_pow_set = (d_max_disch_pu:d_max_ch_pu)*p_pu;
d_num = length(out_pow_set);
d_offset = d_max_disch_pu-1;

y_max_pu = x_max_pu + d_max_ch_pu;
y_min_pu = x_min_pu;
y_num = y_max_pu-y_min_pu+1;
y_offset = y_min_pu - 1;

e_pu = p_pu*slotIntervalInHours; % in Wh
z_min_pu = floor(cell_SOC_low*z_cap/e_pu);
z_max_pu = floor(cell_SOC_high*z_cap/e_pu);
z_num = z_max_pu-z_min_pu+1;
z_offset = z_min_pu - 1;
z_grid = (z_min_pu:z_max_pu)*e_pu;
soc_grid_boundaries = z_grid/z_cap;
if(soc_grid_boundaries(1)<cell_SOC_low && soc_grid_boundaries(2) >cell_SOC_low)
    soc_grid_boundaries(1) = cell_SOC_low;
end
if(soc_grid_boundaries(end)>=cell_SOC_high)
    soc_grid_boundaries(end)= cell_SOC_high;
else
    soc_grid_boundaries = [soc_grid_boundaries cell_SOC_high];
end

soc_grid_bin_mean = zeros(z_num,1);
for bin_idx = 1:z_num
    soc_grid_bin_mean(bin_idx) = (soc_grid_boundaries(bin_idx) + soc_grid_boundaries(bin_idx+1))/2;
end

% Prepare params
params = struct;
params.k_num_in_day = k_num_in_day;
params.h_num = h_num;
params.x_num = x_num;
params.x_offset = x_offset;
params.C_HgHh = C_HgHh;
params.p_pu = p_pu;
params.slotIntervalInHours = slotIntervalInHours;
params.evalStartHourIndex = evalStartHourIndex;
params.paramsPrecisionDigits = config.paramsPrecisionDigits;
params.beliefSpacePrecisionDigits = config.beliefSpacePrecisionDigits;
params.timeHorizonsPerDay = timeHorizonsPerDay;
params.timeHorizonSlots = timeHorizonSlots;

params.d_num = d_num;
params.y_num = y_num;
params.d_offset = d_offset;
params.y_offset = y_offset;

params.out_pow_set = out_pow_set;
params.cellsInSeries = cellsInSeries;
params.legsInParallel = legsInParallel;
params.batteryNominalVoltage = batteryNominalVoltage;
params.batteryRatedCapacityInAh = batteryRatedCapacityInAh;
params.deglifePartitions_num = config.deglifePartitions;
params.energyCostPer_Wh = energyCostPer_Wh;
params.capacityCostPer_Ah = capacityCostPer_Ah;
params.z_num = z_num;
params.z_offset = z_offset;
params.e_pu = e_pu;
params.soc_grid_boundaries = soc_grid_boundaries;
params.soc_grid_bin_mean = soc_grid_bin_mean;
end