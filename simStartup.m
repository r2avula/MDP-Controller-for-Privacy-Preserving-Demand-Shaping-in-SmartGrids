function simStartup()
pathCell = regexp(path, pathsep, 'split');
if ispc      
    test_dir = strcat(pwd,'\util');
    onPath = any(strcmpi(test_dir, pathCell));
elseif(isunix)   
    test_dir = strcat(pwd,'/util');
    onPath = any(strcmp(test_dir, pathCell));
else
    error('Unsupported OS!');
end

if (~onPath)        
    path(pathdef);
    addpath(genpath('util'));
    addpath(genpath('config'));    
end
rng(1,'twister');
end
