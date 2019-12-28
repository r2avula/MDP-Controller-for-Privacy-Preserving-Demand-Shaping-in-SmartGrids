function [processedData] = fetchData(config)
trainingHouseIndices = config.trainingHouseIndices;
dataset = config.dataset;
testHouseIndices = config.testHouseIndices;
if(length(testHouseIndices) > 1)
    error('Not implemented!');
end
if(~isequal(dataset,'eco'))
    error('Not implemented!');
end

if(iscell(trainingHouseIndices))
    trainingHouseIndices =  cell2mat(trainingHouseIndices);
    trainingHouseIndices = num2str(reshape(trainingHouseIndices,[],1), '%02d');
end

if(iscell(testHouseIndices))
    testHouseIndices =  cell2mat(testHouseIndices);
    testHouseIndices = num2str(reshape(testHouseIndices,[],1), '%02d');
end



path_to_sm_data = config.path_to_sm_data;

evalStartHourIndex = 1;
evalEndHourIndex = 24;

slotIntervalInSeconds = config.slotIntervalInSeconds;
slotIntervalInHours = slotIntervalInSeconds/3600; %in hours
slot_num_in_day = (evalEndHourIndex-evalStartHourIndex+1)/slotIntervalInHours;

fileNamePrefix = 'cache/smartMeterData_';
smDataParams = struct;
smDataParams.dataset = dataset;
smDataParams.slotIntervalInSeconds = slotIntervalInSeconds;
smDataParams.testHouseIndices = testHouseIndices(1,:);
smDataParams.trainingHouseIndices = trainingHouseIndices;

[filename,fileExists] = findFileName(smDataParams,fileNamePrefix,'smDataParams');
if(fileExists)
    load(filename,'processedData');
else
    houseIndex = testHouseIndices(1,:);
    filename_occupancy_mat = strcat(path_to_sm_data, '/', dataset, '/occupancy/', houseIndex,'_summer.mat');
    if exist(filename_occupancy_mat, 'file') == 0
        filename_occupancy_csv = strcat(path_to_sm_data, '/', dataset, '/occupancy/', houseIndex,'_summer.csv');
        occupancy_raw_data = csvimport(filename_occupancy_csv);
        dateStrings = transpose(string(occupancy_raw_data(2:end,1)));
        occupancyData = transpose(cell2mat(occupancy_raw_data(2:end,2:end)));
        save(filename_occupancy_mat,'dateStrings','occupancyData');
    else
        load(filename_occupancy_mat,'dateStrings','occupancyData');
    end
    
    dateStrings = string(datetime(dateStrings, 'InputFormat', 'dd-MMM-yyyy', 'Format', 'yyyy-MM-dd'));
    
    testDateStrings = dateStrings';
    num_days = length(testDateStrings);
    
    testSMdata = zeros(slot_num_in_day,num_days);
    testGTdata = zeros(slot_num_in_day,num_days);
    invalid_days = zeros(num_days,1);
    
    for day_idx = 1:num_days
        day_str = testDateStrings(day_idx);
        filename_sm = strcat(path_to_sm_data, '/', dataset, '/smartmeter/', houseIndex, '/',day_str,'.mat');
        if exist(filename_sm, 'file')
            sm_consumption = readSmData_ECO(path_to_sm_data,dataset,houseIndex,day_str,slotIntervalInSeconds, 'powerallphases');
            sm_consumption = max(sm_consumption((evalStartHourIndex-1)/slotIntervalInHours+1:evalEndHourIndex/slotIntervalInHours),0);
            
            testSMdata(:,day_idx) = sm_consumption;
            
            temp_occupancyData = mean(reshape(occupancyData(:,day_idx),[],24/slotIntervalInHours),1);
            temp_occupancyData = round(max(temp_occupancyData((evalStartHourIndex-1)/slotIntervalInHours+1:evalEndHourIndex/slotIntervalInHours),0));
            
            away_hourFlag = temp_occupancyData'==0;
            testGTdata(away_hourFlag,day_idx)=1;
            testGTdata(~away_hourFlag,day_idx)=2;
            
        else
            invalid_days(day_idx) = 1;
            continue;
        end   
    end
    
    testSMdata(:,invalid_days==1) = [];
    testGTdata(:,invalid_days==1) = [];
    testDateStrings(invalid_days==1) = [];
    
    
    trainingHouses_num = size(trainingHouseIndices,1);
    temp_SMdata = cell(trainingHouses_num,1);
    temp_GTdata = cell(trainingHouses_num,1);
    temp_num_days = zeros(trainingHouses_num,1);
    for idx = 1:trainingHouses_num
        houseIndex = trainingHouseIndices(idx,:);
        filename_occupancy_mat = strcat(path_to_sm_data, '/', dataset, '/occupancy/', houseIndex,'_summer.mat');
        if exist(filename_occupancy_mat, 'file') == 0
            filename_occupancy_csv = strcat(path_to_sm_data, '/', dataset, '/occupancy/', houseIndex,'_summer.csv');
            occupancy_raw_data = csvimport(filename_occupancy_csv);
            dateStrings = transpose(string(occupancy_raw_data(2:end,1)));
            occupancyData = transpose(cell2mat(occupancy_raw_data(2:end,2:end)));
            save(filename_occupancy_mat,'dateStrings','occupancyData');
        else
            load(filename_occupancy_mat,'dateStrings','occupancyData');
        end
        
        dateStrings = string(datetime(dateStrings, 'InputFormat', 'dd-MMM-yyyy', 'Format', 'yyyy-MM-dd'));
        dateStrings = dateStrings';
        num_days = length(dateStrings);
        
        result_SMdata = zeros(slot_num_in_day,num_days);
        result_GTdata = zeros(slot_num_in_day,num_days);
        invalid_days = zeros(num_days,1);
        
        for day_idx = 1:num_days
            day_str = dateStrings(day_idx);
            filename_sm = strcat(path_to_sm_data, '/', dataset, '/smartmeter/', houseIndex, '/',day_str,'.mat');
            if exist(filename_sm, 'file')
                sm_consumption = readSmData_ECO(path_to_sm_data,dataset,houseIndex,day_str,slotIntervalInSeconds, 'powerallphases');
                sm_consumption = max(sm_consumption((evalStartHourIndex-1)/slotIntervalInHours+1:evalEndHourIndex/slotIntervalInHours),0);
                
                result_SMdata(:,day_idx) = sm_consumption;
                
                temp_occupancyData = mean(reshape(occupancyData(:,day_idx),[],24/slotIntervalInHours),1);
                temp_occupancyData = round(max(temp_occupancyData((evalStartHourIndex-1)/slotIntervalInHours+1:evalEndHourIndex/slotIntervalInHours),0));
                
                away_hourFlag = temp_occupancyData'==0;
                result_GTdata(away_hourFlag,day_idx)=1;
                result_GTdata(~away_hourFlag,day_idx)=2;
            else
                invalid_days(day_idx) = 1;
                continue;
            end
        end
        
        result_SMdata(:,invalid_days==1) = [];
        result_GTdata(:,invalid_days==1) = [];        
        
        temp_SMdata{idx} = result_SMdata;
        temp_GTdata{idx} = result_GTdata;
        temp_num_days(idx) = size(result_SMdata,2);
    end
    
    num_days = sum(temp_num_days);
    trainingSMdata = zeros(slot_num_in_day,num_days);
    trainingGTdata = zeros(slot_num_in_day,num_days);
    offset_day_idx = 0;
    for idx = 1:trainingHouses_num
        range = offset_day_idx+1:offset_day_idx+temp_num_days(idx);
        offset_day_idx = offset_day_idx + temp_num_days(idx);
        trainingSMdata(:,range) = temp_SMdata{idx};
        trainingGTdata(:,range) = temp_GTdata{idx};
    end   
    
    processedData = struct;
    processedData.testSMdata = testSMdata;
    processedData.testGTdata = testGTdata;
    processedData.testDateStrings = testDateStrings;
    processedData.trainingSMdata = trainingSMdata;
    processedData.trainingGTdata = trainingGTdata;
    
    save(filename,'processedData','smDataParams');
end
end

