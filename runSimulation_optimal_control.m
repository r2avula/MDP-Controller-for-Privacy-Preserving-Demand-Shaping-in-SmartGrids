clear;
[path_to_sm_data,path_to_degradation_data] = simStartup();

reDoSimulationEvenIfCacheExists = 0;
numMCevalDays = 4000;
num_tradeoffs = 6;

adversary_config_filename = 'occupancy_test_house_2_1hour.yaml';
adversary_config = ReadYaml(adversary_config_filename);
adversary_config.path_to_sm_data = path_to_sm_data;
adversary_processedData = fetchData(adversary_config);
adversary_trainingSMdata = adversary_processedData.trainingSMdata;
adversary_trainingGTdata = adversary_processedData.trainingGTdata;
adversary_config.trainingSMdata_max = max(adversary_trainingSMdata(:));
adversary_Params = initAdversaryParams(adversary_config);
adversary_HMM_params = getHMMParams_th(adversary_Params,adversary_trainingSMdata,adversary_trainingGTdata);

controller_config_filename = 'mdp_control_occupancy_data_house_2_1hour.yaml';
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
for idx = 1:numMCevalDays
    eval_day_idxs(idx) = randi(availableEvalDays);
end

evaluationSMdata = testSMdata(:,eval_day_idxs);
evaluationGTdata = testGTdata(:,eval_day_idxs);

plotData = struct;
%% optimal demand-shaping with MDP controller

tradeOff_omega_withinUtility_idx = num_tradeoffs;
tradeOff_sigma_forcost_idx = num_tradeoffs;

controller_config.path_to_degradation_data = path_to_degradation_data;
controller_essParams = getESSParams(controller_Params,controller_config);

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

evalCacheParams.controllerPolicy = getMDP_deviation_cost_tradeoff_policy(evalCacheParams);

fileNamePrefix = strcat('cache/evaluationData_',num2str(batteryRatedCapacityInAh,'%02d'),'_',...
    num2str(tradeOff_omega_withinUtility_idx),'_',...
    num2str(tradeOff_sigma_forcost_idx),'_');

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
plotData.simulatedControllerData_demand_shaping = simulatedControllerData;
plotData.bayesRiskAveragedInHorizon_demand_shaping = bayesRiskAveragedInHorizon;

%% optimal privacy control

tradeOff_omega_withinUtility_idx = 1;
tradeOff_sigma_forcost_idx = num_tradeoffs;

evalCacheParams = struct;
evalCacheParams.numMCevalDays = numMCevalDays;
evalCacheParams.adversary_Params = adversary_Params;
evalCacheParams.adversary_HMM_params = adversary_HMM_params;
evalCacheParams.controller_Params = controller_Params;
evalCacheParams.controller_HMM_params = controller_HMM_params;
evalCacheParams.controller_essParams = controller_essParams;
evalCacheParams.privacyWeight = 1;
evalCacheParams.deviationWeight = 0;
evalCacheParams.costWeight = 0;
evalCacheParams.tradeOff_omega_withinUtility_idx = tradeOff_omega_withinUtility_idx;
evalCacheParams.tradeOff_sigma_forcost_idx = tradeOff_sigma_forcost_idx;

evalCacheParams.controllerPolicy = getMDP_privacy_cost_tradeoff_policy(evalCacheParams);

fileNamePrefix = strcat('cache/evaluationData_',num2str(batteryRatedCapacityInAh,'%02d'),'_',...
    num2str(tradeOff_omega_withinUtility_idx),'_',...
    num2str(tradeOff_sigma_forcost_idx),'_');

[filename,fileExists] = findFileName(evalCacheParams,fileNamePrefix,'evalCacheParams');
if fileExists && ~reDoSimulationEvenIfCacheExists
    load(filename,'simulatedControllerData','estimatedOccupancyData','bayesRiskAveragedInHorizon','overallBayesRisk');
    disp(strcat({'evaluationData loaded from '},filename,' .'));
else
    simulatedControllerData = simulateMDP_privacy_cost_tradeoff(evalCacheParams,evaluationSMdata,evaluationGTdata);
    estimatedOccupancyData = runSequentialBayesDetection(adversary_Params,adversary_HMM_params,simulatedControllerData.modifiedSMdata);
    [bayesRiskAveragedInHorizon,overallBayesRisk] = computeBayesRisk(adversary_Params,evaluationGTdata,estimatedOccupancyData);
    save(filename,'simulatedControllerData','estimatedOccupancyData','bayesRiskAveragedInHorizon','overallBayesRisk','evalCacheParams');
    disp(strcat({'evaluationData saved to '},filename,' .'));
end
plotData.simulatedControllerData_privacy_control = simulatedControllerData;
plotData.bayesRiskAveragedInHorizon_privacy_control = bayesRiskAveragedInHorizon;
plotData.controller_Params = controller_Params;
plotData.evaluationSMdata = evaluationSMdata;

plotFigures__optimal_control(plotData);