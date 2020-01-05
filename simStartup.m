function [path_to_sm_data,path_to_degradation_data] = simStartup()
pathCell = regexp(path, pathsep, 'split');
test_dir = [pwd filesep 'util'];
onPath = any(strcmpi(test_dir, pathCell));

if (~onPath)        
    path(pathdef);
    addpath(genpath('util'));
    addpath(genpath('config'));    
end
path_to_sm_data = [pwd filesep 'data'];
path_to_degradation_data = [pwd filesep 'data' filesep 'battery_degradation'];
rng(1,'twister');
rng(1,'twister');
end
