function [bayesRiskAveragedInHorizon,overallBayesRisk] = computeBayesRisk(detectionParams,evaluationGTdata,estimatedOccupancyData)
numDays = size(estimatedOccupancyData,2);
GTData = evaluationGTdata(:,1:numDays);
h_num = detectionParams.h_num;
C_HgHh = eye(h_num);

timeHorizonsPerDay = detectionParams.timeHorizonsPerDay;
k_num_in_day = detectionParams.k_num_in_day;

k_num_in_horizon = k_num_in_day/timeHorizonsPerDay;

bayesRisk  = zeros(k_num_in_day,numDays);

for k_in_day = 1:k_num_in_day
    for h_idx = 1:h_num
        for hh_idx = 1:h_num
            bayesRisk(k_in_day,:) = bayesRisk(k_in_day,:) + C_HgHh(h_idx,hh_idx)*(GTData(k_in_day,1:numDays)==h_idx & estimatedOccupancyData(k_in_day,:)==hh_idx);
        end
    end
end

bayesRiskAveragedInHorizon = zeros(size(bayesRisk));

for day_idx = 1:numDays
    temp = mean(reshape(bayesRisk(:,day_idx),k_num_in_horizon,[]),1);
    bayesRiskAveragedInHorizon(:,day_idx) = repelem(temp,k_num_in_horizon);
end

overallBayesRisk = mean(bayesRisk(:));
end

