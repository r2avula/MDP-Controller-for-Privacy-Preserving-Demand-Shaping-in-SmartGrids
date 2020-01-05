function [controller_essParams] = getESSParams(controller_Params,controller_config)
cell_SOC_high = controller_config.cell_SOC_high;
cell_SOC_low = controller_config.cell_SOC_low;
cell_1C_power = controller_config.cell_1C_power; %in W

cellCapacityInAh = controller_config.cell_1C_capacityInAh;

soc_grid_boundaries_sim = linspace(cell_SOC_low,cell_SOC_high,9);
cell_pow_set_sim = linspace(-cell_1C_power,cell_1C_power,15);
slotIntervalInSeconds = controller_config.slotIntervalInSeconds;

deglifePartitions_num =  controller_config.deglifePartitions;
deglifePartitions = linspace(1,0.8,deglifePartitions_num+1);
allowedRelativeCapacityChange = (deglifePartitions(1)-deglifePartitions(2))*100;

cellSimData_all_partitions = cell(deglifePartitions_num,1);

for partition_idx = 1:deglifePartitions_num
    cellSimParams = struct;
    cellSimParams.soc_grid_boundaries = soc_grid_boundaries_sim;
    cellSimParams.cell_pow_set = cell_pow_set_sim;
    cellSimParams.initialRelCap = deglifePartitions(partition_idx)*100;
    cellSimParams.allowedRelativeCapacityChange = allowedRelativeCapacityChange; % 20
    cellSimParams.sample_num = controller_config.degSampleNum; % 20
    cellSimParams.slotIntervalInSeconds = slotIntervalInSeconds;
    cellSimParams.SOC_low = controller_config.cell_SOC_low;
    cellSimParams.SOC_high = controller_config.cell_SOC_high;
    cellSimParams.SOC_init = controller_config.cell_SOC_init;
    cellSimParams.cell_voltage_high = controller_config.cell_voltage_high;
    cellSimParams.cell_voltage_low = controller_config.cell_voltage_low;
    cellSimParams.deglifePartitions = controller_config.deglifePartitions;
    cellSimParams.driveToSOH_timeAccelerationFactor = controller_config.driveToSOH_timeAccelerationFactor;
    cellSimParams.driveToSOC_timeAccelerationFactor = controller_config.driveToSOC_timeAccelerationFactor;
    driveToSOC_attempts_max = controller_config.driveToSOC_attempts_max;
    if(ischar(driveToSOC_attempts_max))
        driveToSOC_attempts_max = str2double(driveToSOC_attempts_max);
    end
    cellSimParams.driveToSOC_attempts_max = driveToSOC_attempts_max;
    fileNamePrefix = [controller_config.path_to_degradation_data filesep 'cellSimData'];
    
    [filename,fileExists] = findFileName(cellSimParams,fileNamePrefix,'cellSimParams');
    if(fileExists)
        load(filename,'cellSimData');
        disp(strcat({'cellSimData loaded from '},filename,' .'));
    else
        error('Requires COMSOL simulation!');
    end
    cellSimData_all_partitions{partition_idx} = cellSimData;
end

%% Process cell degradation data
deglifePartitions_num = length(cellSimData_all_partitions);

[soc_num_sim,pow_num_sim,~] = size(cellSimData_all_partitions{1}.simTimeRatio_samples);

soc_grid_bin_mean_sim = zeros(soc_num_sim,1);
for soc_bin_idx_sim = 1:soc_num_sim
    soc_grid_bin_mean_sim(soc_bin_idx_sim) = (soc_grid_boundaries_sim(soc_bin_idx_sim) + soc_grid_boundaries_sim(soc_bin_idx_sim+1))/2;
end

energy_loss_rate_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in W
capacity_loss_rate_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in A
simTimeRatio_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
soc_kp1_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);

slotIntervalInHours = controller_config.slotIntervalInSeconds/3600;
paramsPrecisionDigits = controller_config.paramsPrecisionDigits;

for partition_idx = 1:deglifePartitions_num
    cellSimData = cellSimData_all_partitions{partition_idx};
    
    cell_energy_loss_samples = cellSimData.cell_energy_loss_samples; % in Wh
    cell_capacity_loss_rate_samples = cellSimData.capacity_loss_factor_samples*cellCapacityInAh/slotIntervalInHours; % in A
    simTimeRatio_samples = cellSimData.simTimeRatio_samples;    
    z_kp1_idx_samples = cellSimData.z_kp1_idx_samples;
    
    simTimeThreshold_ene = 0;
    simTimeThreshold_cap = 0;
    
    for soc_bin_idx_sim = 1:soc_num_sim
        for pow_idx_sim = 1:pow_num_sim
            dataSamples = reshape(cell_energy_loss_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            simTimeRatioSamples = reshape(simTimeRatio_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            dataSamples = dataSamples./(simTimeRatioSamples*slotIntervalInHours); %in W
            dataSamples = dataSamples(~isinf(dataSamples)&~isnan(dataSamples)&simTimeRatioSamples>=simTimeThreshold_ene);
            if(~isempty(dataSamples))
                energy_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = abs(mean(dataSamples));
            end
            
            dataSamples = reshape(cell_capacity_loss_rate_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            dataSamples = dataSamples(~isnan(dataSamples)&simTimeRatioSamples>=simTimeThreshold_cap);
            if(~isempty(dataSamples))
                capacity_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = abs(mean(dataSamples));
            end 
                        
            dataSamples = reshape(z_kp1_idx_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            dataSamples = dataSamples(~isnan(dataSamples)&simTimeRatioSamples>=simTimeThreshold_cap);
            if(~isempty(dataSamples))
                soc_kp1_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = mean(soc_grid_bin_mean_sim(dataSamples));
            end 
            
            simTimeRatioSamples(isinf(simTimeRatioSamples)|isnan(simTimeRatioSamples)) = 0;
            simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = abs(mean(simTimeRatioSamples));
        end
    end
    
    energy_loss_rate_mean_cell_partitions_sim(:,:,partition_idx) = round(max(energy_loss_rate_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
    capacity_loss_rate_mean_cell_partitions_sim(:,:,partition_idx) = round(max(capacity_loss_rate_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
    simTimeRatio_mean_cell_partitions_sim(:,:,partition_idx) = round(max(simTimeRatio_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
    soc_kp1_mean_cell_partitions_sim (:,:,partition_idx) = round(max(soc_kp1_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
end

res_energy_loss_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Wh
capacity_loss_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Ah

simTimeRatioThreshold = 0.5;

for partition_idx = 1:deglifePartitions_num
    for soc_bin_idx_sim = 1:soc_num_sim
        for pow_idx_sim = 1:pow_num_sim
            soc_next = soc_kp1_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx);
            if(soc_next > soc_grid_boundaries_sim(end) || soc_next< soc_grid_boundaries_sim(1) || ...
                    simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) <=simTimeRatioThreshold)
                possible_state_transition = 0;
            else
                possible_state_transition = 1;
            end
            if(possible_state_transition)   
                res_energy_loss_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    energy_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;                
                capacity_loss_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    capacity_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;                
            end
        end
    end
end

%% compute ESS losses  in time slot duration

cellsInSeries = controller_Params.cellsInSeries;
legsInParallel = controller_Params.legsInParallel;
slotIntervalInHours = controller_Params.slotIntervalInHours;
converterEfficiency = (controller_config.converterEfficiency)/100;     
out_pow_set = controller_Params.out_pow_set;  
paramsPrecisionDigits = controller_Params.paramsPrecisionDigits;
deglifePartitions_num =  controller_config.deglifePartitions;
d_num = controller_Params.d_num;
z_num = controller_Params.z_num;
d_offset = controller_Params.d_offset;
z_offset = controller_Params.z_offset;
p_pu = controller_Params.p_pu;
e_pu = controller_Params.e_pu;
batteryRatedCapacityInAh = controller_Params.batteryRatedCapacityInAh;
energyCostPer_Wh = controller_Params.energyCostPer_Wh;
capacityCostPer_Ah = controller_Params.capacityCostPer_Ah;
soc_grid_bin_mean = controller_Params.soc_grid_bin_mean;

bat_pow_set_sim = cell_pow_set_sim*cellsInSeries*legsInParallel;
out_pow_set_sim = zeros(1,pow_num_sim);
for pow_idx_sim = 1:pow_num_sim
    if(cell_pow_set_sim(pow_idx_sim)<0)
        out_pow_set_sim(pow_idx_sim) = bat_pow_set_sim(pow_idx_sim)*converterEfficiency;
    else
        out_pow_set_sim(pow_idx_sim) = bat_pow_set_sim(pow_idx_sim)/converterEfficiency;
    end
end

energy_loss_rate_mean_bat_partitions_sim = energy_loss_rate_mean_cell_partitions_sim*cellsInSeries*legsInParallel;
capacity_loss_rate_mean_ess_partitions_sim = capacity_loss_rate_mean_cell_partitions_sim*legsInParallel;

energy_loss_rate_mean_ess_partitions_sim = zeros(size(energy_loss_rate_mean_bat_partitions_sim));
energy_loss_ess_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Wh
capacity_loss_ess_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Ah

for partition_idx = 1:deglifePartitions_num
    for soc_bin_idx_sim = 1:soc_num_sim
        for pow_idx_sim = 1:pow_num_sim
            if(simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) <=simTimeRatioThreshold)
                possible_state_transition = 0;
            else
                possible_state_transition = 1;
            end
            if(possible_state_transition)                
                out_pow = out_pow_set_sim(pow_idx_sim);
                bat_pow = bat_pow_set_sim(pow_idx_sim);                
                
                energy_loss_rate_mean_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    energy_loss_rate_mean_bat_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) + ...
                    abs(out_pow-bat_pow);
                
                energy_loss_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    energy_loss_rate_mean_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;
                
                capacity_loss_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    capacity_loss_rate_mean_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;                
            end
        end
    end
end


%% interpolation
energy_loss_ess_partitions_interpolated_1 = zeros(soc_num_sim,d_num,deglifePartitions_num);
energy_loss_ess_partitions_interpolated = zeros(z_num,d_num,deglifePartitions_num);
capacity_loss_ess_partitions_interpolated_1 = zeros(soc_num_sim,d_num,deglifePartitions_num);
capacity_loss_ess_partitions_interpolated = zeros(z_num,d_num,deglifePartitions_num);

for partition_idx = 1:deglifePartitions_num
    for soc_bin_idx_sim = 1:soc_num_sim
        Yt = out_pow_set_sim';
        Zt = reshape(energy_loss_ess_partitions_sim(soc_bin_idx_sim,:,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3 && ~isempty(find(Zt>0, 1)))
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,pow_idx_sim_temp]=findpeaks(K);
            K(pow_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),out_pow_set','linear','extrap');
            energy_loss_ess_partitions_interpolated_1(soc_bin_idx_sim,:,partition_idx) = Zt2;
        end
        Yt = out_pow_set_sim';
        Zt = reshape(capacity_loss_ess_partitions_sim(soc_bin_idx_sim,:,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3)
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,pow_idx_sim_temp]=findpeaks(K);
            K(pow_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),out_pow_set','linear','extrap');
            capacity_loss_ess_partitions_interpolated_1(soc_bin_idx_sim,:,partition_idx) = Zt2;
        end
    end
    
    for pow_idx = 1:d_num
        Yt = soc_grid_bin_mean_sim;
        Zt = reshape(energy_loss_ess_partitions_interpolated_1(:,pow_idx,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3 && ~isempty(find(Zt>0, 1)))
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,soc_bin_idx_sim_temp]=findpeaks(K);
            K(soc_bin_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),soc_grid_bin_mean,'linear','extrap');
            energy_loss_ess_partitions_interpolated(:,pow_idx,partition_idx) = max(Zt2,0);
        end
        Yt = soc_grid_bin_mean_sim;
        Zt = reshape(capacity_loss_ess_partitions_interpolated_1(:,pow_idx,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3)
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,soc_bin_idx_sim_temp]=findpeaks(K);
            K(soc_bin_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),soc_grid_bin_mean,'linear','extrap');
            capacity_loss_ess_partitions_interpolated(:,pow_idx,partition_idx) = max(Zt2,0);
        end
    end
end

essUsageCost_ess_partitions_interpolated = energyCostPer_Wh*energy_loss_ess_partitions_interpolated + 5*capacityCostPer_Ah*capacity_loss_ess_partitions_interpolated;
essUsageCost_ess_partitions_interpolated = round(max(essUsageCost_ess_partitions_interpolated,0),paramsPrecisionDigits);
energy_loss_ess_partitions_interpolated = round(max(energy_loss_ess_partitions_interpolated,0),paramsPrecisionDigits);
capacity_loss_ess_partitions_interpolated = round(max(capacity_loss_ess_partitions_interpolated,0),paramsPrecisionDigits);

z_kp1_idxs_map = zeros(z_num,d_num,deglifePartitions_num);
for partition_idx = 1:deglifePartitions_num      
    for z_idx=1:z_num
        z_cur = (z_offset+z_idx)*e_pu;
        for d_idx=1:d_num
            out_pow = (d_idx+d_offset)*p_pu;
            en_loss_ess = energy_loss_ess_partitions_interpolated(z_idx,d_idx,partition_idx);
            z_next = z_cur + out_pow*slotIntervalInHours - en_loss_ess;
            z_next_idx = round(z_next/e_pu)-z_offset;
            if(z_next_idx>z_num)
                z_dif = z_next-(z_num+z_offset)*e_pu;
                if(z_dif <e_pu/2)
                    possible_state_transition = 1;
                    z_next_idx = z_num;
                else
                    possible_state_transition = 0;
                end
            elseif(z_next_idx<1)
                z_dif = (z_offset+1)*e_pu - z_next;
                if(z_dif <e_pu/2)
                    possible_state_transition = 1;
                    z_next_idx = 1;
                else
                    possible_state_transition = 0;
                end
            else
                possible_state_transition = 1;
            end
            if(possible_state_transition)
                z_kp1_idxs_map(z_idx,d_idx,partition_idx) = z_next_idx;
            else
                z_kp1_idxs_map(z_idx,d_idx,partition_idx) = nan;
                capacity_loss_ess_partitions_interpolated(z_idx,d_idx,partition_idx) = nan;
                energy_loss_ess_partitions_interpolated(z_idx,d_idx,partition_idx) = nan;
                essUsageCost_ess_partitions_interpolated(z_idx,d_idx,partition_idx) = nan;
            end            
        end
    end    
end

% Prepare params
controller_essParams{deglifePartitions_num} = struct;
for partition_idx = 1:deglifePartitions_num
    controller_essParams{partition_idx}.batteryRatedCapacityInAh = batteryRatedCapacityInAh;
    controller_essParams{partition_idx}.cellsInSeries = cellsInSeries;
    controller_essParams{partition_idx}.legsInParallel = legsInParallel;
    controller_essParams{partition_idx}.capacityLossInAh_map = capacity_loss_ess_partitions_interpolated(:,:,partition_idx);
    controller_essParams{partition_idx}.energyLossInWh_map = energy_loss_ess_partitions_interpolated(:,:,partition_idx);
    controller_essParams{partition_idx}.z_kp1_idxs_map = z_kp1_idxs_map(:,:,partition_idx);
    controller_essParams{partition_idx}.essUsageCost_map = essUsageCost_ess_partitions_interpolated(:,:,partition_idx);
end
end
