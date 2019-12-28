clear;
simStartup();

reDoSimulationEvenIfCacheExists = 0;
numMCevalDays = 6000;

[path_to_parent,cur_folder,~]=fileparts(pwd);
if(ispc)
    path_to_sm_data = strcat(path_to_parent,'\',cur_folder,'\data');
    path_to_degradation_data = strcat(path_to_parent,'\',cur_folder,'\data\battery_degradation');
else
    path_to_sm_data = strcat(path_to_parent,'/',cur_folder,'/data');
    path_to_degradation_data = strcat(path_to_parent,'/',cur_folder,'/data/battery_degradation');
end

controller_config_filename = 'mdp_control_occupancy_data_house_2.yaml';
controller_config = ReadYaml(controller_config_filename);
controller_config.path_to_sm_data = path_to_sm_data;
controller_processedData = fetchData(controller_config);
controller_trainingSMdata = controller_processedData.trainingSMdata;
controller_trainingGTdata = controller_processedData.trainingGTdata;
controller_config.trainingSMdata_max = max(controller_trainingSMdata(:));
controller_Params = initControllerParams(controller_config);
controller_HMM_params = getHMMParams_tnh(controller_Params,controller_trainingSMdata,controller_trainingGTdata);

testSMdata = controller_processedData.testSMdata;
testGTdata = controller_processedData.testGTdata;
availableEvalDays = size(testSMdata,2);

eval_day_idxs = zeros(numMCevalDays,1);

availableEvalDay_idx = 1;
for idx = 1:numMCevalDays
    eval_day_idxs(idx) = availableEvalDay_idx;
    availableEvalDay_idx = availableEvalDay_idx + 1;
    if(availableEvalDay_idx == availableEvalDays+1)
        availableEvalDay_idx = 1;
    end
end

evaluationSMdata = testSMdata(:,eval_day_idxs);
evaluationGTdata = testGTdata(:,eval_day_idxs);

%% optimal demand-shaping with MDP controller

controller_config.path_to_degradation_data = path_to_degradation_data;
controller_essParams = getESSParams(controller_Params,controller_config);


num_tradeoffs = 6;
tradeOffWithinUtility = round(linspace(0,1,num_tradeoffs),2);
tradeOffForcost = round(linspace(0,1,num_tradeoffs),2);
tradeOff_omega_withinUtility_idx = num_tradeoffs;
tradeOff_sigma_forcost_idx = num_tradeoffs;

evalCacheParams = struct;
evalCacheParams.numMCevalDays = numMCevalDays;
evalCacheParams.controller_Params = controller_Params;
evalCacheParams.controller_HMM_params = controller_HMM_params;
evalCacheParams.controller_essParams = controller_essParams;
evalCacheParams.privacyWeight = 0;
evalCacheParams.deviationWeight = 1;
evalCacheParams.costWeight = 0;
evalCacheParams.tradeOff_omega_withinUtility_idx = tradeOff_omega_withinUtility_idx;
evalCacheParams.tradeOff_sigma_forcost_idx = tradeOff_sigma_forcost_idx;
batteryRatedCapacityInAh = controller_Params.batteryRatedCapacityInAh;

evalCacheParams.controllerPolicy = getMDP_deviation_cost_tradeoff(evalCacheParams);

fileNamePrefix = strcat('cache/evaluationData_',num2str(batteryRatedCapacityInAh,'%02d'),'_',...
    num2str(tradeOff_omega_withinUtility_idx,'%02d'),'_',...
    num2str(tradeOff_sigma_forcost_idx,'%02d'),'_');

[filename,fileExists] = findFileName(evalCacheParams,fileNamePrefix,'evalCacheParams');
if fileExists && ~reDoSimulationEvenIfCacheExists
    load(filename,'simulatedControllerData','estimatedOccupancyData','bayesRiskAveragedInHorizon','overallBayesRisk');
    disp(strcat({'evaluationData loaded from '},filename,' .'));
else
    simulatedControllerData = simulateMDP_deviation_cost_tradeoff(evalCacheParams,evaluationSMdata,evaluationGTdata);
    estimatedOccupancyData = runSequentialBayesDetection(adversary_Params,adversary_HMM_params,simulatedControllerData.modifiedSMdata);
    [bayesRiskAveragedInHorizon,overallBayesRisk] = computeBayesRisk(adversary_Params,evaluationGTdata,estimatedOccupancyData);
    save(filename,'simulatedControllerData','estimatedOccupancyData','bayesRiskAveragedInHorizon','overallBayesRisk','evalCacheParams');
    disp(strcat({'evaluationData saved to '},filename,' .'));
end






