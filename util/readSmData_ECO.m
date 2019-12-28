function [result] = readSmData_ECO(path_to_data,dataset, house_str, dates_strs, granularity_in_sec, option)

    % returns the smartmeter data of the specified household
    % granularity = data frequency (in seconds)
    % dates_strc = days
    % option = data type (power, current, ...)
    
    result = zeros(1, size(dates_strs,1)*(24*3600)/ granularity_in_sec);
    offset = 1;
    for day_idx=1:size(dates_strs,1)
        filename_sm = strcat(path_to_data, '/', dataset, '/smartmeter/', house_str, '/', dates_strs(day_idx), '.mat');
        % Smartmeter data
        if exist(filename_sm, 'file')
            vars = whos('-file',filename_sm);
            load(filename_sm);
            eval(['smartmeter_data=' vars.name ';']);
            eval(['clear ' vars.name ';']);
            if (granularity_in_sec > 1)
                 % powerallphases
                 eval(strcat('[mat_t,padded] = vec2mat(smartmeter_data.',option, ',', num2str(granularity_in_sec),');'));
                 assert(padded == 0, [num2str(granularity_in_sec), ' is not a permissable interval (does not divide into 24h)']);
                 result(1,offset:offset + (24 * 60 * 60)/granularity_in_sec -1) = mean(mat_t, 2);            
            else
                 eval(strcat('result(1,offset:offset + (24 * 60 * 60)/granularity_in_sec -1) = smartmeter_data.',option,';')); 
            end
        else
            result(1,offset:offset + (24 * 60 * 60)/granularity_in_sec -1) = -1;
        end
        offset = offset + (24 * 60 * 60) / granularity_in_sec;
    end
end
