function [estimatedOccupancyData] = runSequentialBayesDetection(detectionParams,HMM_params,evaluationSMdata)
numDays = size(evaluationSMdata,2);
h_num = detectionParams.h_num;
x_num = detectionParams.x_num;
C_HgHh = eye(h_num);
p_pu = detectionParams.p_pu;
x_offset = detectionParams.x_offset;
x_k_obs_idxs = min(round(evaluationSMdata/p_pu)-x_offset,x_num);

timeHorizonsPerDay = detectionParams.timeHorizonsPerDay;
k_num_in_day = detectionParams.k_num_in_day;

k_num_in_horizon = k_num_in_day/timeHorizonsPerDay;

estimatedOccupancyData = zeros(k_num_in_day,numDays);

paramsPrecisionDigits = detectionParams.paramsPrecisionDigits;
paramsPrecision = 10^(-paramsPrecisionDigits);


P_HgHn1 = HMM_params.P_HgHn1;
M_b = HMM_params.M_b;
P_Hk = HMM_params.P_Hk;

for day_idx = 1:numDays    
    for horizon_idx = 1:timeHorizonsPerDay
        for k_in_horizon = 1:k_num_in_horizon
            k_in_day = (horizon_idx-1)*k_num_in_horizon + k_in_horizon;
            
            if(k_in_horizon==1)
                % Reset belief state
                if(k_in_day==1)
                    belief_kn1 = P_Hk(:,k_num_in_day);
                else
                    belief_kn1 = P_Hk(:,k_in_day-1);
                end
            end            
            
            x_k_idx = x_k_obs_idxs(k_in_day,day_idx);
            
            belief_k = M_b(:,:,min(x_k_idx,x_num),k_in_day)*belief_kn1;
            if(sum(belief_k)<paramsPrecision)
                belief_k = P_HgHn1(:,:,k_in_day)*belief_kn1;
                if(sum(belief_k)<paramsPrecision)
                    belief_k = P_Hk(:,k_in_day);
                end
            end
            belief_k = belief_k/sum(belief_k);
            
            [~,detectedState] = max(C_HgHh*belief_k);
            estimatedOccupancyData(k_in_day,day_idx) = detectedState(1);
            belief_kn1 = belief_k;
        end
    end
end
end

