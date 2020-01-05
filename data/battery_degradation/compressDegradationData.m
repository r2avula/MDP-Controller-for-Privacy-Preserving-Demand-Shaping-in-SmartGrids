clear;
path_to_original_degradation_data = '$PATH_TO_Li-Ion-Battery-Degradation-ToolKit$';

[path_to_parent,~,~]=fileparts(pwd);
deg_config_filename = [path_to_parent filesep 'config' filesep 'degradation_20percent_1hour.yaml'];
deg_config = ReadYaml(deg_config_filename);
deg_config.path_to_original_degradation_data = path_to_original_degradation_data;

compression_factor = 3;

cell_SOC_high = deg_config.cell_SOC_high;
cell_SOC_low = deg_config.cell_SOC_low;
cell_1C_power = deg_config.cell_1C_power; %in W

old_deglifePartitions_num =  deg_config.deglifePartitions;
old_deglifePartitions = linspace(1,0.8,old_deglifePartitions_num+1);
old_allowedRelativeCapacityChange = (old_deglifePartitions(1)-old_deglifePartitions(2))*100;

soc_grid_boundaries_sim = linspace(cell_SOC_low,cell_SOC_high,9);
cell_pow_set_sim = linspace(-cell_1C_power,cell_1C_power,15);

new_deglifePartitions_num = deg_config.deglifePartitions/compression_factor;

new_deglifePartitions = linspace(1,0.8,new_deglifePartitions_num+1);
new_allowedRelativeCapacityChange = (new_deglifePartitions(1)-new_deglifePartitions(2))*100;
soc_num = length(soc_grid_boundaries_sim)-1;
pow_num = length(cell_pow_set_sim);
new_sample_num = deg_config.degSampleNum*compression_factor;

old_partition_idx_offset = 0;

for partition_idx = 1:new_deglifePartitions_num
    cellSimParams = struct;
    cellSimParams.soc_grid_boundaries = soc_grid_boundaries_sim;
    cellSimParams.cell_pow_set = cell_pow_set_sim;
    cellSimParams.initialRelCap = new_deglifePartitions(partition_idx)*100;
    cellSimParams.allowedRelativeCapacityChange = new_allowedRelativeCapacityChange; 
    cellSimParams.sample_num = deg_config.degSampleNum*compression_factor; 
    cellSimParams.slotIntervalInSeconds = deg_config.slotIntervalInSeconds;
    cellSimParams.SOC_low = deg_config.cell_SOC_low;
    cellSimParams.SOC_high = deg_config.cell_SOC_high;
    cellSimParams.SOC_init = deg_config.cell_SOC_init;
    cellSimParams.cell_voltage_high = deg_config.cell_voltage_high;
    cellSimParams.cell_voltage_low = deg_config.cell_voltage_low;
    cellSimParams.deglifePartitions = deg_config.deglifePartitions/compression_factor;
    cellSimParams.driveToSOH_timeAccelerationFactor = deg_config.driveToSOH_timeAccelerationFactor;
    cellSimParams.driveToSOC_timeAccelerationFactor = deg_config.driveToSOC_timeAccelerationFactor;
    driveToSOC_attempts_max = deg_config.driveToSOC_attempts_max;
    if(ischar(driveToSOC_attempts_max))
        driveToSOC_attempts_max = str2double(driveToSOC_attempts_max);
    end
    cellSimParams.driveToSOC_attempts_max = driveToSOC_attempts_max;
    fileNamePrefix = 'cellSimData';
    [new_fileName,fileExists] = findFileName(cellSimParams,fileNamePrefix,'cellSimParams');
    if(fileExists)
        break;
    end  
    new_cellSimParams = cellSimParams;
    cell_energy_loss_samples = Inf([soc_num,pow_num,new_sample_num]);
    capacity_loss_factor_samples = Inf([soc_num,pow_num,new_sample_num]);
    z_kp1_idx_samples = Inf([soc_num,pow_num,new_sample_num]);
    simTimeRatio_samples = Inf([soc_num,pow_num,new_sample_num]);
    cellsIterated = 0;
    
    range_old_parition_idx = old_partition_idx_offset+1:old_partition_idx_offset+compression_factor;
    old_partition_idx_offset = old_partition_idx_offset+compression_factor;
    sample_num_offset = 0;
    for old_parition_idx = range_old_parition_idx
        cellSimParams = struct;
        cellSimParams.soc_grid_boundaries = soc_grid_boundaries_sim;
        cellSimParams.cell_pow_set = cell_pow_set_sim;
        cellSimParams.initialRelCap = old_deglifePartitions(old_parition_idx)*100;
        cellSimParams.allowedRelativeCapacityChange = old_allowedRelativeCapacityChange; 
        cellSimParams.sample_num = deg_config.degSampleNum;
        cellSimParams.slotIntervalInSeconds = deg_config.slotIntervalInSeconds;
        cellSimParams.SOC_low = deg_config.cell_SOC_low;
        cellSimParams.SOC_high = deg_config.cell_SOC_high;
        cellSimParams.SOC_init = deg_config.cell_SOC_init;
        cellSimParams.cell_voltage_high = deg_config.cell_voltage_high;
        cellSimParams.cell_voltage_low = deg_config.cell_voltage_low;
        cellSimParams.deglifePartitions = deg_config.deglifePartitions;
        cellSimParams.driveToSOH_timeAccelerationFactor = deg_config.driveToSOH_timeAccelerationFactor;
        cellSimParams.driveToSOC_timeAccelerationFactor = deg_config.driveToSOC_timeAccelerationFactor;
        driveToSOC_attempts_max = deg_config.driveToSOC_attempts_max;
        if(ischar(driveToSOC_attempts_max))
            driveToSOC_attempts_max = str2double(driveToSOC_attempts_max);
        end
        cellSimParams.driveToSOC_attempts_max = driveToSOC_attempts_max;
        
        fileNamePrefix = [deg_config.path_to_original_degradation_data filesep 'cellSimData'];
        [filename,fileExists] = findFileName(cellSimParams,fileNamePrefix,'cellSimParams');
        if(fileExists)
            load(filename,'cellSimData');
        else
            error('Something is wrong!');
        end
        simTimeRatio_samples(:,:,sample_num_offset+1:sample_num_offset+deg_config.degSampleNum) = cellSimData.simTimeRatio_samples;
        cell_energy_loss_samples(:,:,sample_num_offset+1:sample_num_offset+deg_config.degSampleNum) = cellSimData.cell_energy_loss_samples;
        capacity_loss_factor_samples(:,:,sample_num_offset+1:sample_num_offset+deg_config.degSampleNum) = cellSimData.capacity_loss_factor_samples;
        z_kp1_idx_samples(:,:,sample_num_offset+1:sample_num_offset+deg_config.degSampleNum) = cellSimData.z_kp1_idx_samples;
        cellsIterated = cellSimData.cellsIterated;
        sample_num_offset = sample_num_offset + deg_config.degSampleNum;
    end    

    cellSimParams = new_cellSimParams;
    cellSimData = struct;
    cellSimData.simTimeRatio_samples = simTimeRatio_samples;
    cellSimData.cell_energy_loss_samples = cell_energy_loss_samples;
    cellSimData.capacity_loss_factor_samples = capacity_loss_factor_samples;
    cellSimData.z_kp1_idx_samples = z_kp1_idx_samples;
    cellSimData.cellsIterated = cellsIterated;
    save(new_fileName,'cellSimData','cellSimParams')    
end