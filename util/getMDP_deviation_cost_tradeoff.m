function [policy_all_partitions] = getMDP_deviation_cost_tradeoff(evalCacheParams)
%% Hypothesis-aware EMU optimal strategy design
controller_Params = evalCacheParams.controller_Params;
controller_essParams_all_partitions = evalCacheParams.controller_essParams;
controller_HMM_params = evalCacheParams.controller_HMM_params;
deglifePartitions_num = controller_Params.deglifePartitions_num;
batteryRatedCapacityInAh = controller_Params.batteryRatedCapacityInAh;
deviationWeight = evalCacheParams.deviationWeight;
costWeight = evalCacheParams.costWeight;
if(evalCacheParams.privacyWeight~=0)
    error('Executing wrong function!');
end

tradeOff_omega_withinUtility_idx = evalCacheParams.tradeOff_omega_withinUtility_idx;
tradeOff_sigma_forcost_idx = evalCacheParams.tradeOff_sigma_forcost_idx;

policyParams = struct;
policyParams.controller_Params = controller_Params;
policyParams.controller_HMM_params = controller_HMM_params;
policyParams.deviationWeight = evalCacheParams.deviationWeight;
policyParams.costWeight = evalCacheParams.costWeight;

policy_all_partitions = cell(deglifePartitions_num,1);
for partition_idx = 1:deglifePartitions_num
    controller_essParams = controller_essParams_all_partitions{partition_idx};
    policyParams.controller_essParams = controller_essParams;
    
    policyFileNamePrefix = strcat('cache/policy_',...
        num2str(batteryRatedCapacityInAh,'%02d'),'_',...
        num2str(partition_idx),'_',...
        num2str(tradeOff_omega_withinUtility_idx),'_',...
        num2str(tradeOff_sigma_forcost_idx),'_');
    
    [policyFileName,fileExists] = findFileName(policyParams,policyFileNamePrefix,'policyParams');
    if(fileExists)
        disp(strcat({'MDP control policy found in: '},policyFileName));
        load(policyFileName,'policy');
    else                        
        k_num_in_day = controller_Params.k_num_in_day;
        h_num = controller_Params.h_num;
        x_num = controller_Params.x_num;
        z_num = controller_Params.z_num;
        d_num = controller_Params.d_num;
        y_num = controller_Params.y_num;
        
        timeHorizonsPerDay = controller_Params.timeHorizonsPerDay;
        
        k_num_in_horizon = k_num_in_day/timeHorizonsPerDay;
        
        d_offset = controller_Params.d_offset;
        x_offset = controller_Params.x_offset;
        y_offset = controller_Params.y_offset;
                
        paramsPrecisionDigits = controller_Params.paramsPrecisionDigits;
        paramsPrecision = 10^(-paramsPrecisionDigits);
        essUsageCost = controller_essParams.essUsageCost_map;
        
        controller_M_b = controller_HMM_params.M_b;
        p_pu = controller_Params.p_pu;
                
        z_kp1_idx_map = controller_essParams.z_kp1_idxs_map;
                
        controllerDecisions = cell(k_num_in_day,1);
        
        valueFuncation_kp1 = zeros(z_num,x_num,h_num,y_num);
                
        disp('Computing policy...');        
        for horizon_idx = timeHorizonsPerDay:-1:1            
            for k_in_horizon = k_num_in_horizon:-1:1
                k_in_day = (horizon_idx-1)*k_num_in_horizon + k_in_horizon;
                               
                controllerDecisions_k = zeros(z_num,x_num,h_num,y_num);
                valueFunction_k = zeros(z_num,x_num,h_num,y_num);
                opt_start = tic;                
                
                for y_kn1_idx = 1:y_num
                    for h_k_idx = 1:h_num
                        for x_k_idx = 1:x_num
                            for z_k_idx = 1:z_num
                                alphaVector_k = zeros(d_num,1);
                                for d_k_idx = 1:d_num
                                    L_bar = essUsageCost(z_k_idx,d_k_idx);
                                    z_kp1_idx_t = z_kp1_idx_map(z_k_idx,d_k_idx);
                                    if(isnan(L_bar)||isnan(z_kp1_idx_t)||z_kp1_idx_t<=0)
                                        alphaVector_k(d_k_idx) = inf;
                                    else
                                        y_k_idx = min((x_k_idx+x_offset) + (d_k_idx+d_offset) - y_offset,y_num);
                                        if(y_k_idx<1)
                                            alphaVector_k(d_k_idx) = inf;
                                        else
                                            temp = deviationWeight*abs((y_kn1_idx - y_k_idx)*p_pu) + costWeight*L_bar;
                                            if(k_in_day<k_num_in_day)
                                                for x_kp1_idx = 1:x_num
                                                    for h_kp1_idx = 1:h_num
                                                        temp = temp + valueFuncation_kp1(z_kp1_idx_t,x_kp1_idx,h_kp1_idx,y_k_idx)*...
                                                            controller_M_b(h_kp1_idx,h_k_idx,x_kp1_idx,k_in_day+1);
                                                    end
                                                end
                                            end
                                            alphaVector_k(d_k_idx) = temp;
                                        end
                                    end
                                end
                                
                                [opt_val,d_k_idx_star] = min(alphaVector_k);
                                controllerDecisions_k(z_k_idx,x_k_idx,h_k_idx,y_kn1_idx) = d_k_idx_star;
                                valueFunction_k(z_k_idx,x_k_idx,h_k_idx,y_kn1_idx) = opt_val;
                            end
                        end
                    end
                end                
                
                opt_time = toc(opt_start);
                
                controllerDecisions{k_in_day} = controllerDecisions_k;
                valueFuncation_kp1 = valueFunction_k;
                
                disp(strcat('Time index: ',num2str(k_in_day),...
                    ', Opt:',num2str(opt_time)));
            end                
        end
        
        %% save designed policy
        
        policy = struct;
        policy.controllerDecisions = controllerDecisions;    
        
        save(policyFileName,'policy','policyParams')
        disp(strcat({'MDP control policy saved in: '},policyFileName));
    end
    policy_all_partitions{partition_idx} = policy;
end
end


