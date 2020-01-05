function [performanceData] = simulateMDP_privacy_cost_tradeoff(evalCacheParams,evaluationSMdata,evaluationGTdata)
controller_Params = evalCacheParams.controller_Params;
adversary_Params = evalCacheParams.adversary_Params;
controller_essParams = evalCacheParams.controller_essParams;
adversary_HMM_params = evalCacheParams.adversary_HMM_params;
deglifePartitions_num = controller_Params.deglifePartitions_num; 
batteryRatedCapacityInAh = controller_Params.batteryRatedCapacityInAh;
timeHorizonsPerDay = controller_Params.timeHorizonsPerDay;

%% Hypothesis-aware EMU simulation with estimated ESS model

z_num = controller_Params.z_num;
slotIntervalInHours = controller_Params.slotIntervalInHours; 

% if(controller_Params.x_num~= adversary_Params.x_num)
%     error('Not implemented!');
% end
controller_x_num = controller_Params.x_num;
controller_y_num = controller_Params.y_num;

adversary_M_b = adversary_HMM_params.M_b;
adversary_P_HgHn1 = adversary_HMM_params.P_HgHn1;
adversary_P_Hk = adversary_HMM_params.P_Hk;
adversary_x_num = adversary_Params.x_num;

p_pu = controller_Params.p_pu;
y_offset = controller_Params.y_offset;
x_offset = controller_Params.x_offset;
d_offset = controller_Params.d_offset;
energyCostPer_Wh = controller_Params.energyCostPer_Wh;
capacityCostPer_Ah = controller_Params.capacityCostPer_Ah;
k_num_in_day = controller_Params.k_num_in_day;
k_num_in_horizon = k_num_in_day/timeHorizonsPerDay;

numDays = size(evaluationSMdata,2);


deglifePartitions = linspace(1,0.8,deglifePartitions_num+1);
soh_partitions = round(5*deglifePartitions-4,2);

paramsPrecisionDigits = controller_Params.paramsPrecisionDigits;
beliefSpacePrecisionDigits = controller_Params.beliefSpacePrecisionDigits;
paramsPrecision = 10^(-paramsPrecisionDigits);

y_k_idxs = zeros(k_num_in_day,numDays);
modifiedSMdata = zeros(k_num_in_day,numDays);
ESSpower = zeros(k_num_in_day,numDays);
z_k_idxs = zeros(k_num_in_day,numDays);
dailyRelativeCapacity = zeros(k_num_in_day,numDays);
dailySOH = zeros(k_num_in_day,numDays);
ESSusage_cost = zeros(k_num_in_day,numDays);

SOH_partition_changed_day_idxs = ones(deglifePartitions_num,1);
EOL_reached = 0;
cumulative_scaled_squared_deviation = 0;
totalCapacityLoss  = 0;
totalEnergyLoss  = 0;
totalEnergyWastage  = 0;
totalESSUsageInkWh = 0;

x_k_idxs = min(round(evaluationSMdata/p_pu)-x_offset,controller_x_num);
h_k_idxs = round(evaluationGTdata);

z_k_idxs(1,1) = z_num;

controllerPolicy_all_partitions = evalCacheParams.controllerPolicy;

current_degradation_partition_idx = 1;

priorBeliefSpacePartitionFirstElements = controllerPolicy_all_partitions{current_degradation_partition_idx}.priorBeliefSpacePartitionFirstElements;
controllerDecisions = controllerPolicy_all_partitions{current_degradation_partition_idx}.controllerDecisions;
z_kp1_idx_est_map_3c = controller_essParams{current_degradation_partition_idx}.z_kp1_idxs_map;
capacityLossInAh_3c = controller_essParams{current_degradation_partition_idx}.capacityLossInAh_map;
energyLossInWh_3c = controller_essParams{current_degradation_partition_idx}.energyLossInWh_map;
essUsageCost_3c = controller_essParams{current_degradation_partition_idx}.essUsageCost_map;

y_kn1_idx = 1;

for day_idx = 1:numDays
    disp(strcat({' Simulating EMU: '}, num2str(day_idx),{'/'},num2str(numDays),{'...'}));
    
    for horizon_idx = 1:timeHorizonsPerDay   
        if(horizon_idx>1)
            belief_kn1 = adversary_HMM_params.P_Hk(:,k_in_day-1);           
        else
            belief_kn1 = adversary_HMM_params.P_Hk(:,end);            
        end                

        for k_in_horizon=1:k_num_in_horizon   
            k_in_day = (horizon_idx-1)*k_num_in_horizon + k_in_horizon;                                         
            
            z_k_idx_sys = z_k_idxs(k_in_day,day_idx);
            
            p_idx_k = find(priorBeliefSpacePartitionFirstElements{k_in_day}>=belief_kn1(1),1)-1;
            if(p_idx_k==0)
                p_idx_k = 1;
            end
            
            d_k_idx_star = controllerDecisions{k_in_day}(z_k_idx_sys,x_k_idxs(k_in_day,day_idx),h_k_idxs(k_in_day,day_idx),p_idx_k);
            
            z_kp1_idx = z_kp1_idx_est_map_3c(z_k_idx_sys,d_k_idx_star);
            y_k_idx = (x_k_idxs(k_in_day,day_idx)+x_offset) + (d_k_idx_star+d_offset) - y_offset;
            y_k_idxs(k_in_day,day_idx) = y_k_idx;
            modifiedSMdata(k_in_day,day_idx) = (y_k_idxs(k_in_day,day_idx) + y_offset)*p_pu;
            
            if(k_in_day<k_num_in_day)
                z_k_idxs(k_in_day+1,day_idx) = z_kp1_idx;
            elseif(day_idx<numDays)
                z_k_idxs(1,day_idx+1) = z_kp1_idx;
            end
            
            cumulative_scaled_squared_deviation = cumulative_scaled_squared_deviation + (y_k_idx - y_kn1_idx)^2;    
            
            totalCapacityLoss = totalCapacityLoss + capacityLossInAh_3c(z_k_idx_sys,d_k_idx_star);
            totalEnergyLoss = totalEnergyLoss + energyLossInWh_3c(z_k_idx_sys,d_k_idx_star);
            dailyRelativeCapacity(k_in_day,day_idx) = (batteryRatedCapacityInAh-totalCapacityLoss)/batteryRatedCapacityInAh;
            dailySOH(k_in_day,day_idx) = 5*dailyRelativeCapacity(k_in_day,day_idx)-4;
            
            ESSpower(k_in_day,day_idx) = p_pu*(d_k_idx_star+d_offset);
            ESSusage_cost(k_in_day,day_idx) = essUsageCost_3c(z_k_idx_sys,d_k_idx_star);
            totalESSUsageInkWh = totalESSUsageInkWh + abs(ESSpower(k_in_day,day_idx))*slotIntervalInHours/1000;
            
            y_kn1_idx = y_k_idx;
            
            if(dailySOH(k_in_day,day_idx) <= soh_partitions(current_degradation_partition_idx + 1))
                disp(strcat({'ESS 3C estimate of SOH partition changed from '},num2str(current_degradation_partition_idx),{' to '},num2str(current_degradation_partition_idx+1),{'...'}));
                current_degradation_partition_idx = current_degradation_partition_idx + 1;
                if(current_degradation_partition_idx == deglifePartitions_num+1)
                    numDays = day_idx-1;
                    EOL_reached = 1;
                    disp('EOL reached as per 3C!');
                    break;
                else
                    SOH_partition_changed_day_idxs(current_degradation_partition_idx) = day_idx+1;        
                    priorBeliefSpacePartitionFirstElements = controllerPolicy_all_partitions{current_degradation_partition_idx}.priorBeliefSpacePartitionFirstElements;
                    controllerDecisions = controllerPolicy_all_partitions{current_degradation_partition_idx}.controllerDecisions;
                    z_kp1_idx_est_map_3c = controller_essParams{current_degradation_partition_idx}.z_kp1_idxs_map;
                    capacityLossInAh_3c = controller_essParams{current_degradation_partition_idx}.capacityLossInAh_map;
                    energyLossInWh_3c = controller_essParams{current_degradation_partition_idx}.energyLossInWh_map;
                    essUsageCost_3c = controller_essParams{current_degradation_partition_idx}.essUsageCost_map;                    
                end
            end           
            
            belief_k = adversary_M_b(:,:,min(y_k_idxs(k_in_day,day_idx),adversary_x_num),k_in_day)*belief_kn1;
            if(sum(belief_k)<paramsPrecision)
                belief_k = adversary_P_HgHn1(:,:,k_in_day)*belief_kn1;
                if(sum(belief_k)<paramsPrecision)
                    belief_k = adversary_P_Hk(:,k_in_day);
                end
            end
            belief_k = belief_k/sum(belief_k);            
            belief_kn1 = round(belief_k,beliefSpacePrecisionDigits);
        end        
        if(EOL_reached)
            break;
        end
    end
    
    if(EOL_reached)
        break;
    end        
end

totalESSusageCost = energyCostPer_Wh*(totalEnergyLoss+totalEnergyWastage) + 5*capacityCostPer_Ah*totalCapacityLoss;
rms_deviation = sqrt(cumulative_scaled_squared_deviation/(k_num_in_day*numDays))*p_pu;

if(EOL_reached)
    y_k_idxs(:,numDays+1:end) = [];
    modifiedSMdata(:,numDays+1:end) = [];
    ESSpower(:,numDays+1:end) = [];
    z_k_idxs(:,numDays+1:end) = [];
    dailyRelativeCapacity(:,numDays+1:end) = [];
    dailySOH(:,numDays+1:end) = [];
    ESSusage_cost(:,numDays+1:end) = [];
end

batteryLifeIndaysExtrapolated = batteryRatedCapacityInAh*0.2*numDays/totalCapacityLoss;

performanceData = struct;
performanceData.ESSusage_cost = ESSusage_cost;
performanceData.y_k_idxs = y_k_idxs;
performanceData.EOL_reached = EOL_reached;
performanceData.totalCapacityLoss = totalCapacityLoss;
performanceData.totalEnergyLoss = totalEnergyLoss;
performanceData.totalEnergyWastage = totalEnergyWastage;
performanceData.totalESSusageCost = totalESSusageCost;
performanceData.modifiedSMdata = modifiedSMdata;
performanceData.ESSpower = ESSpower;
performanceData.z_k_idxs = z_k_idxs;
performanceData.dailyRelativeCapacity = dailyRelativeCapacity;
performanceData.dailySOH = dailySOH;
performanceData.SOH_partition_changed_day_idxs = SOH_partition_changed_day_idxs;
performanceData.batteryLifeIndaysExtrapolated = batteryLifeIndaysExtrapolated;
performanceData.totalESSUsageInkWh = totalESSUsageInkWh;
performanceData.rms_deviation = rms_deviation;
end

