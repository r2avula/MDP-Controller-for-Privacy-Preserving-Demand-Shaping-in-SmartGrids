function [params] = initAdversaryParams(config)
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
end