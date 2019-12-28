function [controller_essParams] = getESSParams(controller_Params,controller_config)
cell_SOC_high = controller_config.cell_SOC_high;
cell_SOC_low = controller_config.cell_SOC_low;
cell_1C_power = controller_config.cell_1C_power; %in W

cellCapacityInAh = 1.8048;

soc_grid_boundaries_sim = linspace(cell_SOC_low,cell_SOC_high,9);
cell_pow_set_sim = linspace(-cell_1C_power,cell_1C_power,15);
slotIntervalInSeconds = controller_config.slotIntervalInSeconds;

deglifePartitions_num =  controller_config.deglifePartitions;
deglifePartitions = linspace(1,0.8,deglifePartitions_num+1);
allowedRelativeCapacityChange = (deglifePartitions(1)-deglifePartitions(2))*100;

cellSimData_all_partitions = cell(deglifePartitions_num,1);

for partition_idx = 1:deglifePartitions_num
    cellSimParams = struct;
    cellSimParams.soc_grid_boundaries = soc_grid_boundaries_sim;
    cellSimParams.cell_pow_set = cell_pow_set_sim;
    cellSimParams.initialRelCap = deglifePartitions(partition_idx)*100;
    cellSimParams.allowedRelativeCapacityChange = allowedRelativeCapacityChange; % 20
    cellSimParams.sample_num = controller_config.degSampleNum; % 20
    cellSimParams.slotIntervalInSeconds = slotIntervalInSeconds;
    cellSimParams.SOC_low = controller_config.cell_SOC_low;
    cellSimParams.SOC_high = controller_config.cell_SOC_high;
    cellSimParams.SOC_init = controller_config.cell_SOC_init;
    cellSimParams.cell_voltage_high = controller_config.cell_voltage_high;
    cellSimParams.cell_voltage_low = controller_config.cell_voltage_low;
    cellSimParams.deglifePartitions = controller_config.deglifePartitions;
    cellSimParams.timeAccelerationFactor = controller_config.timeAccelerationFactor;
    if(ispc)
        fileNamePrefix = strcat(controller_config.path_to_degradation_data,'\cellSimData');
    else
        fileNamePrefix = strcat(controller_config.path_to_degradation_data,'/cellSimData');
    end
    
    [filename,fileExists] = findFileName(cellSimParams,fileNamePrefix,'cellSimParams');
    if(fileExists)
        load(filename,'cellSimData');
        disp(strcat({'cellSimData loaded from '},filename,' .'));
    else
        error('Requires COMSOL simulation!');
    end
    cellSimData_all_partitions{partition_idx} = cellSimData;
end

%% Process cell degradation data
deglifePartitions_num = length(cellSimData_all_partitions);

[soc_num_sim,pow_num_sim,~] = size(cellSimData_all_partitions{1}.simTimeRatio_samples);

soc_grid_bin_mean_sim = zeros(soc_num_sim,1);
for soc_bin_idx_sim = 1:soc_num_sim
    soc_grid_bin_mean_sim(soc_bin_idx_sim) = (soc_grid_boundaries_sim(soc_bin_idx_sim) + soc_grid_boundaries_sim(soc_bin_idx_sim+1))/2;
end

energy_loss_rate_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in W
capacity_loss_rate_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in A
simTimeRatio_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
% voltage_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);
soc_kp1_mean_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num);

slotIntervalInHours = controller_config.slotIntervalInSeconds/3600;
paramsPrecisionDigits = controller_config.paramsPrecisionDigits;

for partition_idx = 1:deglifePartitions_num
    cellSimData = cellSimData_all_partitions{partition_idx};
    
    cell_energy_loss_samples = cellSimData.cell_energy_loss_samples; % in Wh
    cell_capacity_loss_rate_samples = cellSimData.capacity_loss_factor_samples*cellCapacityInAh/slotIntervalInHours; % in A
    simTimeRatio_samples = cellSimData.simTimeRatio_samples;    
%     cell_voltage_samples = cellSimData.cell_voltage_samples;
    z_kp1_idx_samples = cellSimData.z_kp1_idx_samples;
    
    simTimeThreshold_ene = 0;
    simTimeThreshold_cap = 0;
    
    for soc_bin_idx_sim = 1:soc_num_sim
        for pow_idx_sim = 1:pow_num_sim
            dataSamples = reshape(cell_energy_loss_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            simTimeRatioSamples = reshape(simTimeRatio_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            dataSamples = dataSamples./(simTimeRatioSamples*slotIntervalInHours); %in W
            dataSamples = dataSamples(~isinf(dataSamples)&~isnan(dataSamples)&simTimeRatioSamples>=simTimeThreshold_ene);
            if(~isempty(dataSamples))
                energy_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = abs(mean(dataSamples));
            end
            
            dataSamples = reshape(cell_capacity_loss_rate_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            dataSamples = dataSamples(~isnan(dataSamples)&simTimeRatioSamples>=simTimeThreshold_cap);
            if(~isempty(dataSamples))
                capacity_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = abs(mean(dataSamples));
            end 
            
%             dataSamples = reshape(cell_voltage_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
%             dataSamples = dataSamples(~isnan(dataSamples)&simTimeRatioSamples>=simTimeThreshold_cap);
%             if(~isempty(dataSamples))
%                 voltage_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = abs(mean(dataSamples));
%             end 
            
            dataSamples = reshape(z_kp1_idx_samples(soc_bin_idx_sim,pow_idx_sim,:),1,[]);
            dataSamples = dataSamples(~isnan(dataSamples)&simTimeRatioSamples>=simTimeThreshold_cap);
            if(~isempty(dataSamples))
                soc_kp1_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = mean(soc_grid_bin_mean_sim(dataSamples));
            end 
            
            simTimeRatioSamples(isinf(simTimeRatioSamples)|isnan(simTimeRatioSamples)) = 0;
            simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = abs(mean(simTimeRatioSamples));
        end
    end
    
    energy_loss_rate_mean_cell_partitions_sim(:,:,partition_idx) = round(max(energy_loss_rate_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
    capacity_loss_rate_mean_cell_partitions_sim(:,:,partition_idx) = round(max(capacity_loss_rate_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
%     voltage_mean_cell_partitions_sim(:,:,partition_idx) = round(max(voltage_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
    simTimeRatio_mean_cell_partitions_sim(:,:,partition_idx) = round(max(simTimeRatio_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
    soc_kp1_mean_cell_partitions_sim (:,:,partition_idx) = round(max(soc_kp1_mean_cell_partitions_sim(:,:,partition_idx),0),paramsPrecisionDigits);
end

res_energy_loss_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Wh
capacity_loss_cell_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Ah

simTimeRatioThreshold = 0.5;

for partition_idx = 1:deglifePartitions_num
    for soc_bin_idx_sim = 1:soc_num_sim
        for pow_idx_sim = 1:pow_num_sim
            soc_next = soc_kp1_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx);
            if(soc_next > soc_grid_boundaries_sim(end) || soc_next< soc_grid_boundaries_sim(1) || ...
                    simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) <=simTimeRatioThreshold)
                possible_state_transition = 0;
            else
                possible_state_transition = 1;
            end
            if(possible_state_transition)   
                res_energy_loss_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    energy_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;                
                capacity_loss_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    capacity_loss_rate_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;                
            end
        end
    end
end

%% Plot cell degradation maps

savefiles = 0;
ribbon_size = 1;
cell_cur_set_sim = linspace(-1,1,pow_num_sim);
plot_life_partitions_res_energy_loss_interp_percentage = 0;
if(plot_life_partitions_res_energy_loss_interp_percentage)
    filename = 'energyDegradation';
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    energy_loss_cell_paritions_percentage = zeros(size(res_energy_loss_cell_partitions_sim));
    for pow_idx_sim=1:pow_num_sim
        if(abs(cell_pow_set_sim(pow_idx_sim))>0)
            energy_loss_cell_paritions_percentage(:,pow_idx_sim,:) = res_energy_loss_cell_partitions_sim(:,pow_idx_sim,:)/abs(cell_pow_set_sim(pow_idx_sim))/slotIntervalInHours*100;
        end
    end
        
    y_axis_grid = cell_cur_set_sim;
    x_axis_grid = unique([0,soc_grid_boundaries_sim,1]);
    
    percentage_data = nan(length(y_axis_grid),length(x_axis_grid)-1,deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx_sim = 2:length(x_axis_grid)-2
            percentage_data(:,soc_bin_idx_sim,partition_idx) = energy_loss_cell_paritions_percentage(soc_bin_idx_sim-1,:,partition_idx);
        end
    end  
         
    energy_loss_plot = energy_loss_cell_paritions_percentage;
    
    energy_loss_min = min(energy_loss_plot(:));
    energy_loss_max = max(energy_loss_plot(:));
        
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        color_data_temp = cell(length(x_axis_grid)-1,1);
        for i = 1:length(x_axis_grid)-1
            color_data_temp{i} = [percentage_data(:,i,partition_idx),percentage_data(:,i,partition_idx)];
        end
        z_data = zeros(size(percentage_data(:,:,partition_idx)));
        
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');      
        h1 = ribbon(y_axis_grid,z_data,ribbon_size); 
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1, {'CData'}, color_data_temp,'FaceColor','interp','MeshStyle','column','FaceAlpha',1)
        xticks(0.5:2:10.5);
        xticklabels(0:0.2:1)
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0.5,10.5]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','YTick',[cell_cur_set_sim(1)  0  cell_cur_set_sim(end)],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
        colormap jet;
        caxis([energy_loss_min energy_loss_max]);
    end
    h = colorbar('Position',[0.888 0.068 0.0266666666666669 0.892],...
        'TickLabelInterpreter','latex','FontSize',fontSize);
    caxis([energy_loss_min energy_loss_max]);
    ylabel(h, 'Resistive energy loss ($\%$)','Interpreter','latex','FontSize',fontSize);
    colormap jet;
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

plot_life_partitions_capacity_loss_interp_fraction = 0;
if(plot_life_partitions_capacity_loss_interp_fraction)    
    filename = 'capacityDegradation';
    % plot figure--
    partition_idxs = 1:deglifePartitions_num;
    positions = {[0.0924678097197633 0.246165556113454 0.20936667244291 0.615405872457975],...
        [0.373264911169039 0.238165556113454 0.20936667244291 0.623405872457975],...
        [0.654062012618313 0.238165556113454 0.209366672442911 0.627405872457975]};    
    
    textpositions = {[0.232218011771345 1.10249990460461 0],...
        [0.181262597758605 1.10233123498622 0],...
        [0.24495686527453 1.09583826559003 0]...
        };
    
    capacity_loss_mean_cell_paritions_fraction = capacity_loss_cell_partitions_sim/cellCapacityInAh;
    
    capacity_loss_min = min(capacity_loss_mean_cell_paritions_fraction(:));
    capacity_loss_max = round(max(capacity_loss_mean_cell_paritions_fraction(:))*10^6)*10^-6+0.2*10^-6;
    
    y_axis_grid = cell_cur_set_sim;
    x_axis_grid = unique([0,soc_grid_boundaries_sim,1]);
    
    fraction_data = nan(length(y_axis_grid),length(x_axis_grid)-1,deglifePartitions_num);
    for partition_idx = 1:deglifePartitions_num
        for soc_bin_idx_sim = 2:length(x_axis_grid)-2
            fraction_data(:,soc_bin_idx_sim,partition_idx) = capacity_loss_mean_cell_paritions_fraction(soc_bin_idx_sim-1,:,partition_idx);
        end
    end  
        
    fontSize = 12;
    figure_xsize = 700;
    figure_ysize = 200;
    figure_position_offset = 50;
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    for partition_idx = partition_idxs
        color_data_temp = cell(length(x_axis_grid)-1,1);
        for i = 1:length(x_axis_grid)-1
            color_data_temp{i} = [fraction_data(:,i,partition_idx),fraction_data(:,i,partition_idx)];
        end
        z_data = zeros(size(fraction_data(:,:,partition_idx)));
        
        axes1 = axes('Parent',figure_1,...
            'Position',positions{partition_idx});
        box(axes1,'on');
        hold(axes1,'on');      
        h1 = ribbon(y_axis_grid,z_data,ribbon_size);     
        ylabel({'Current (A)'},'Interpreter','latex','Units', 'normalized','position',[-0.113092878602192,0.525806928450061,0]);
        xlabel({'SOC'},'Interpreter','latex','Units','normalized');
        set(h1, {'CData'}, color_data_temp,'FaceColor','interp','MeshStyle','column','FaceAlpha',1)
        xticks(0.5:2:10.5);
        xticklabels(0:0.2:1)
        axes1.YAxis.MinorTickValues = -1:0.5:1;
        xlim([0.5,10.5]);
        set(axes1,'FontSize',fontSize,'TickLabelInterpreter','latex',...
            'XMinorTick','on','YMinorTick','on','YTick',[cell_cur_set_sim(1)  0  cell_cur_set_sim(end)],...
            'Colormap',...
            [0 0 0.515625;0 0 0.552364864864865;0 0 0.58910472972973;0 0 0.625844594594595;0 0 0.662584459459459;0 0 0.699324324324324;0 0 0.736064189189189;0 0 0.772804054054054;0 0 0.809543918918919;0 0 0.846283783783784;0 0 0.883023648648649;0 0 0.919763513513513;0 0 0.956503378378378;0 0 0.993243243243243;0 0.0299831081081081 1;0 0.066722972972973 1;0 0.103462837837838 1;0 0.140202702702703 1;0 0.176942567567568 1;0 0.213682432432433 1;0 0.250422297297298 1;0 0.287162162162163 1;0 0.323902027027027 1;0 0.360641891891892 1;0 0.397381756756757 1;0 0.434121621621622 1;0 0.470861486486487 1;0 0.507601351351352 1;0 0.544341216216217 1;0 0.581081081081082 1;0 0.617820945945947 1;0 0.654560810810812 1;0 0.691300675675677 1;0 0.728040540540541 1;0 0.764780405405406 1;0 0.801520270270271 1;0 0.859375 1;0 0.871432648401826 1;0 0.883490296803653 1;0 0.895547945205479 1;0 0.907605593607306 1;0 0.919663242009132 1;0 0.931720890410959 1;0 0.943778538812785 1;0 0.955836187214611 1;0 0.967893835616438 1;0 0.979951484018264 1;0 0.992009132420091 1;0.00406678082191725 1 0.995933219178083;0.0161244292237437 1 0.983875570776256;0.0281820776255701 1 0.97181792237443;0.0402397260273966 1 0.959760273972603;0.052297374429223 1 0.947702625570777;0.0643550228310494 1 0.935644977168951;0.0764126712328759 1 0.923587328767124;0.0884703196347023 1 0.911529680365298;0.100527968036529 1 0.899472031963471;0.112585616438355 1 0.887414383561645;0.124643264840182 1 0.875356735159818;0.136700913242008 1 0.863299086757992;0.148758561643834 1 0.851241438356166;0.160816210045661 1 0.839183789954339;0.172873858447487 1 0.827126141552513;0.184931506849314 1 0.815068493150686;0.19698915525114 1 0.80301084474886;0.209046803652967 1 0.790953196347033;0.221104452054793 1 0.778895547945207;0.23316210045662 1 0.76683789954338;0.245219748858446 1 0.754780251141554;0.257277397260272 1 0.742722602739728;0.269335045662099 1 0.730664954337901;0.281392694063925 1 0.718607305936075;0.293450342465752 1 0.706549657534248;0.305507990867578 1 0.694492009132422;0.317565639269405 1 0.682434360730595;0.329623287671231 1 0.670376712328769;0.341680936073057 1 0.658319063926943;0.353738584474884 1 0.646261415525116;0.36579623287671 1 0.63420376712329;0.377853881278537 1 0.622146118721463;0.389911529680363 1 0.610088470319637;0.40196917808219 1 0.59803082191781;0.414026826484016 1 0.585973173515984;0.426084474885843 1 0.573915525114157;0.438142123287669 1 0.561857876712331;0.450199771689495 1 0.549800228310505;0.462257420091322 1 0.537742579908678;0.474315068493148 1 0.525684931506852;0.486372716894975 1 0.513627283105025;0.498430365296801 1 0.501569634703199;0.510488013698628 1 0.489511986301372;0.522545662100454 1 0.477454337899546;0.53460331050228 1 0.46539668949772;0.546660958904107 1 0.453339041095893;0.558718607305933 1 0.441281392694067;0.57077625570776 1 0.42922374429224;0.582833904109586 1 0.417166095890414;0.594891552511413 1 0.405108447488587;0.606949200913239 1 0.393050799086761;0.619006849315066 1 0.380993150684934;0.631064497716892 1 0.368935502283108;0.643122146118718 1 0.356877853881282;0.655179794520545 1 0.344820205479455;0.667237442922371 1 0.332762557077629;0.679295091324198 1 0.320704908675802;0.691352739726024 1 0.308647260273976;0.703410388127851 1 0.296589611872149;0.715468036529677 1 0.284531963470323;0.727525684931503 1 0.272474315068497;0.73958333333333 1 0.26041666666667;0.751640981735156 1 0.248359018264844;0.763698630136983 1 0.236301369863017;0.775756278538809 1 0.224243721461191;0.787813926940636 1 0.212186073059364;0.799871575342462 1 0.200128424657538;0.811929223744289 1 0.188070776255711;0.823986872146115 1 0.176013127853885;0.836044520547941 1 0.163955479452059;0.848102168949768 1 0.151897831050232;0.860159817351594 1 0.139840182648406;0.872217465753421 1 0.127782534246579;0.884275114155247 1 0.115724885844753;0.896332762557074 1 0.103667237442926;0.9083904109589 1 0.0916095890410999;0.920448059360726 1 0.0795519406392735;0.932505707762553 1 0.0674942922374471;0.944563356164379 1 0.0554366438356206;0.956621004566206 1 0.0433789954337942;0.968678652968032 1 0.0313213470319678;0.980736301369859 1 0.0192636986301413;0.992793949771685 1 0.00720605022831489;1 0.995148401826488 0;1 0.983090753424662 0;1 0.971033105022836 0;1 0.958975456621009 0;1 0.946917808219183 0;1 0.934860159817356 0;1 0.92280251141553 0;1 0.910744863013703 0;1 0.898687214611877 0;1 0.886629566210051 0;1 0.874571917808224 0;1 0.862514269406398 0;1 0.850456621004571 0;1 0.838398972602745 0;1 0.826341324200918 0;1 0.814283675799092 0;1 0.802226027397265 0;1 0.790168378995439 0;1 0.778110730593613 0;1 0.766053082191786 0;1 0.75399543378996 0;1 0.741937785388133 0;1 0.729880136986307 0;1 0.71782248858448 0;1 0.705764840182654 0;1 0.693707191780828 0;1 0.681649543379001 0;1 0.669591894977175 0;1 0.657534246575348 0;1 0.645476598173522 0;1 0.633418949771695 0;1 0.621361301369869 0;1 0.609303652968042 0;1 0.597246004566216 0;1 0.58518835616439 0;1 0.573130707762563 0;1 0.561073059360737 0;1 0.54901541095891 0;1 0.536957762557084 0;1 0.524900114155257 0;1 0.512842465753431 0;1 0.500784817351605 0;1 0.488727168949778 0;1 0.476669520547952 0;1 0.464611872146125 0;1 0.452554223744299 0;1 0.440496575342472 0;1 0.428438926940646 0;1 0.416381278538819 0;1 0.404323630136993 0;1 0.392265981735167 0;1 0.38020833333334 0;1 0.368150684931514 0;1 0.356093036529687 0;1 0.344035388127861 0;1 0.331977739726034 0;1 0.319920091324208 0;1 0.307862442922382 0;1 0.295804794520555 0;1 0.283747146118729 0;1 0.271689497716902 0;1 0.259631849315076 0;1 0.247574200913249 0;1 0.235516552511423 0;1 0.223458904109596 0;1 0.21140125570777 0;1 0.199343607305944 0;1 0.187285958904117 0;1 0.175228310502291 0;1 0.163170662100464 0;1 0.151113013698638 0;1 0.139055365296811 0;1 0.126997716894985 0;1 0.114940068493159 0;1 0.102882420091332 0;1 0.0908247716895056 0;1 0.0787671232876792 0;1 0.0667094748858528 0;1 0.0546518264840263 0;1 0.0425941780821999 0;1 0.0305365296803735 0;1 0.018478881278547 0;1 0.00642123287672058 0;0.994363584474894 0 0;0.982305936073068 0 0;0.970248287671241 0 0;0.958190639269415 0 0;0.946132990867588 0 0;0.934075342465762 0 0;0.922017694063936 0 0;0.909960045662109 0 0;0.897902397260283 0 0;0.885844748858456 0 0;0.87378710045663 0 0;0.861729452054803 0 0;0.849671803652977 0 0;0.83761415525115 0 0;0.825556506849324 0 0;0.813498858447498 0 0;0.801441210045671 0 0;0.789383561643845 0 0;0.777325913242018 0 0;0.765268264840192 0 0;0.753210616438365 0 0;0.741152968036539 0 0;0.729095319634713 0 0;0.717037671232886 0 0;0.70498002283106 0 0;0.692922374429233 0 0;0.680864726027407 0 0;0.66880707762558 0 0;0.656749429223754 0 0;0.644691780821927 0 0;0.632634132420101 0 0;0.620576484018275 0 0;0.608518835616448 0 0;0.596461187214622 0 0;0.584403538812795 0 0;0.572345890410969 0 0;0.560288242009142 0 0;0.548230593607316 0 0;0.53617294520549 0 0;0.524115296803663 0 0;0.512057648401837 0 0;0.50000000000001 0 0],...
            'YTickLabel',{'-1C','0','1C'});
        if(partition_idx == 1)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[0,1/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        elseif(partition_idx == deglifePartitions_num)
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(deglifePartitions_num-1),'/',num2str(deglifePartitions_num),',1]'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        else
            text('Units','normalized','BackgroundColor',[1 1 1],...
                'EdgeColor',[0 0 0],...
                'Interpreter','latex',...
                'String',strcat('age $\in$[',num2str(partition_idx-1),'/',num2str(deglifePartitions_num),',',num2str(partition_idx),'/',num2str(deglifePartitions_num),']'),...
                'Position',textpositions{partition_idx},...
                'Color',[0 0 0],'FontSize',fontSize);
        end
%         colormap jet;
        caxis([capacity_loss_min capacity_loss_max]);
    end
    h = colorbar('Position',[0.888 0.064 0.0266666666666666 0.816],...
        'TickLabelInterpreter','latex','Colormap',...
    [0 0 0.515625;0 0 0.552364864864865;0 0 0.58910472972973;0 0 0.625844594594595;0 0 0.662584459459459;0 0 0.699324324324324;0 0 0.736064189189189;0 0 0.772804054054054;0 0 0.809543918918919;0 0 0.846283783783784;0 0 0.883023648648649;0 0 0.919763513513513;0 0 0.956503378378378;0 0 0.993243243243243;0 0.0299831081081081 1;0 0.066722972972973 1;0 0.103462837837838 1;0 0.140202702702703 1;0 0.176942567567568 1;0 0.213682432432433 1;0 0.250422297297298 1;0 0.287162162162163 1;0 0.323902027027027 1;0 0.360641891891892 1;0 0.397381756756757 1;0 0.434121621621622 1;0 0.470861486486487 1;0 0.507601351351352 1;0 0.544341216216217 1;0 0.581081081081082 1;0 0.617820945945947 1;0 0.654560810810812 1;0 0.691300675675677 1;0 0.728040540540541 1;0 0.764780405405406 1;0 0.801520270270271 1;0 0.859375 1;0 0.871432648401826 1;0 0.883490296803653 1;0 0.895547945205479 1;0 0.907605593607306 1;0 0.919663242009132 1;0 0.931720890410959 1;0 0.943778538812785 1;0 0.955836187214611 1;0 0.967893835616438 1;0 0.979951484018264 1;0 0.992009132420091 1;0.00406678082191725 1 0.995933219178083;0.0161244292237437 1 0.983875570776256;0.0281820776255701 1 0.97181792237443;0.0402397260273966 1 0.959760273972603;0.052297374429223 1 0.947702625570777;0.0643550228310494 1 0.935644977168951;0.0764126712328759 1 0.923587328767124;0.0884703196347023 1 0.911529680365298;0.100527968036529 1 0.899472031963471;0.112585616438355 1 0.887414383561645;0.124643264840182 1 0.875356735159818;0.136700913242008 1 0.863299086757992;0.148758561643834 1 0.851241438356166;0.160816210045661 1 0.839183789954339;0.172873858447487 1 0.827126141552513;0.184931506849314 1 0.815068493150686;0.19698915525114 1 0.80301084474886;0.209046803652967 1 0.790953196347033;0.221104452054793 1 0.778895547945207;0.23316210045662 1 0.76683789954338;0.245219748858446 1 0.754780251141554;0.257277397260272 1 0.742722602739728;0.269335045662099 1 0.730664954337901;0.281392694063925 1 0.718607305936075;0.293450342465752 1 0.706549657534248;0.305507990867578 1 0.694492009132422;0.317565639269405 1 0.682434360730595;0.329623287671231 1 0.670376712328769;0.341680936073057 1 0.658319063926943;0.353738584474884 1 0.646261415525116;0.36579623287671 1 0.63420376712329;0.377853881278537 1 0.622146118721463;0.389911529680363 1 0.610088470319637;0.40196917808219 1 0.59803082191781;0.414026826484016 1 0.585973173515984;0.426084474885843 1 0.573915525114157;0.438142123287669 1 0.561857876712331;0.450199771689495 1 0.549800228310505;0.462257420091322 1 0.537742579908678;0.474315068493148 1 0.525684931506852;0.486372716894975 1 0.513627283105025;0.498430365296801 1 0.501569634703199;0.510488013698628 1 0.489511986301372;0.522545662100454 1 0.477454337899546;0.53460331050228 1 0.46539668949772;0.546660958904107 1 0.453339041095893;0.558718607305933 1 0.441281392694067;0.57077625570776 1 0.42922374429224;0.582833904109586 1 0.417166095890414;0.594891552511413 1 0.405108447488587;0.606949200913239 1 0.393050799086761;0.619006849315066 1 0.380993150684934;0.631064497716892 1 0.368935502283108;0.643122146118718 1 0.356877853881282;0.655179794520545 1 0.344820205479455;0.667237442922371 1 0.332762557077629;0.679295091324198 1 0.320704908675802;0.691352739726024 1 0.308647260273976;0.703410388127851 1 0.296589611872149;0.715468036529677 1 0.284531963470323;0.727525684931503 1 0.272474315068497;0.73958333333333 1 0.26041666666667;0.751640981735156 1 0.248359018264844;0.763698630136983 1 0.236301369863017;0.775756278538809 1 0.224243721461191;0.787813926940636 1 0.212186073059364;0.799871575342462 1 0.200128424657538;0.811929223744289 1 0.188070776255711;0.823986872146115 1 0.176013127853885;0.836044520547941 1 0.163955479452059;0.848102168949768 1 0.151897831050232;0.860159817351594 1 0.139840182648406;0.872217465753421 1 0.127782534246579;0.884275114155247 1 0.115724885844753;0.896332762557074 1 0.103667237442926;0.9083904109589 1 0.0916095890410999;0.920448059360726 1 0.0795519406392735;0.932505707762553 1 0.0674942922374471;0.944563356164379 1 0.0554366438356206;0.956621004566206 1 0.0433789954337942;0.968678652968032 1 0.0313213470319678;0.980736301369859 1 0.0192636986301413;0.992793949771685 1 0.00720605022831489;1 0.995148401826488 0;1 0.983090753424662 0;1 0.971033105022836 0;1 0.958975456621009 0;1 0.946917808219183 0;1 0.934860159817356 0;1 0.92280251141553 0;1 0.910744863013703 0;1 0.898687214611877 0;1 0.886629566210051 0;1 0.874571917808224 0;1 0.862514269406398 0;1 0.850456621004571 0;1 0.838398972602745 0;1 0.826341324200918 0;1 0.814283675799092 0;1 0.802226027397265 0;1 0.790168378995439 0;1 0.778110730593613 0;1 0.766053082191786 0;1 0.75399543378996 0;1 0.741937785388133 0;1 0.729880136986307 0;1 0.71782248858448 0;1 0.705764840182654 0;1 0.693707191780828 0;1 0.681649543379001 0;1 0.669591894977175 0;1 0.657534246575348 0;1 0.645476598173522 0;1 0.633418949771695 0;1 0.621361301369869 0;1 0.609303652968042 0;1 0.597246004566216 0;1 0.58518835616439 0;1 0.573130707762563 0;1 0.561073059360737 0;1 0.54901541095891 0;1 0.536957762557084 0;1 0.524900114155257 0;1 0.512842465753431 0;1 0.500784817351605 0;1 0.488727168949778 0;1 0.476669520547952 0;1 0.464611872146125 0;1 0.452554223744299 0;1 0.440496575342472 0;1 0.428438926940646 0;1 0.416381278538819 0;1 0.404323630136993 0;1 0.392265981735167 0;1 0.38020833333334 0;1 0.368150684931514 0;1 0.356093036529687 0;1 0.344035388127861 0;1 0.331977739726034 0;1 0.319920091324208 0;1 0.307862442922382 0;1 0.295804794520555 0;1 0.283747146118729 0;1 0.271689497716902 0;1 0.259631849315076 0;1 0.247574200913249 0;1 0.235516552511423 0;1 0.223458904109596 0;1 0.21140125570777 0;1 0.199343607305944 0;1 0.187285958904117 0;1 0.175228310502291 0;1 0.163170662100464 0;1 0.151113013698638 0;1 0.139055365296811 0;1 0.126997716894985 0;1 0.114940068493159 0;1 0.102882420091332 0;1 0.0908247716895056 0;1 0.0787671232876792 0;1 0.0667094748858528 0;1 0.0546518264840263 0;1 0.0425941780821999 0;1 0.0305365296803735 0;1 0.018478881278547 0;1 0.00642123287672058 0;0.994363584474894 0 0;0.982305936073068 0 0;0.970248287671241 0 0;0.958190639269415 0 0;0.946132990867588 0 0;0.934075342465762 0 0;0.922017694063936 0 0;0.909960045662109 0 0;0.897902397260283 0 0;0.885844748858456 0 0;0.87378710045663 0 0;0.861729452054803 0 0;0.849671803652977 0 0;0.83761415525115 0 0;0.825556506849324 0 0;0.813498858447498 0 0;0.801441210045671 0 0;0.789383561643845 0 0;0.777325913242018 0 0;0.765268264840192 0 0;0.753210616438365 0 0;0.741152968036539 0 0;0.729095319634713 0 0;0.717037671232886 0 0;0.70498002283106 0 0;0.692922374429233 0 0;0.680864726027407 0 0;0.66880707762558 0 0;0.656749429223754 0 0;0.644691780821927 0 0;0.632634132420101 0 0;0.620576484018275 0 0;0.608518835616448 0 0;0.596461187214622 0 0;0.584403538812795 0 0;0.572345890410969 0 0;0.560288242009142 0 0;0.548230593607316 0 0;0.53617294520549 0 0;0.524115296803663 0 0;0.512057648401837 0 0;0.50000000000001 0 0],...
    'FontSize',fontSize);
    caxis([capacity_loss_min capacity_loss_max]);
%     colormap jet;
    ylabel(h, 'Capacity loss ($\times$1C)','Interpreter','latex','FontSize',fontSize);
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    if(savefiles)
        saveas(gcf,strcat(filename,'.fig'));
        saveas(gcf,strcat(filename,'.eps'),'epsc');
    end
end

%% compute ESS losses  in time slot duration

cellsInSeries = controller_Params.cellsInSeries;
legsInParallel = controller_Params.legsInParallel;
slotIntervalInHours = controller_Params.slotIntervalInHours;
converterEfficiency = (controller_config.converterEfficiency)/100;     
out_pow_set = controller_Params.out_pow_set;  
paramsPrecisionDigits = controller_Params.paramsPrecisionDigits;
deglifePartitions_num =  controller_config.deglifePartitions;
d_num = controller_Params.d_num;
z_num = controller_Params.z_num;
d_offset = controller_Params.d_offset;
z_offset = controller_Params.z_offset;
p_pu = controller_Params.p_pu;
e_pu = controller_Params.e_pu;
batteryRatedCapacityInAh = controller_Params.batteryRatedCapacityInAh;
energyCostPer_Wh = controller_Params.energyCostPer_Wh;
capacityCostPer_Ah = controller_Params.capacityCostPer_Ah;
soc_grid_bin_mean = controller_Params.soc_grid_bin_mean;

bat_pow_set_sim = cell_pow_set_sim*cellsInSeries*legsInParallel;
out_pow_set_sim = zeros(1,pow_num_sim);
for pow_idx_sim = 1:pow_num_sim
    if(cell_pow_set_sim(pow_idx_sim)<0)
        out_pow_set_sim(pow_idx_sim) = bat_pow_set_sim(pow_idx_sim)*converterEfficiency;
    else
        out_pow_set_sim(pow_idx_sim) = bat_pow_set_sim(pow_idx_sim)/converterEfficiency;
    end
end

energy_loss_rate_mean_bat_partitions_sim = energy_loss_rate_mean_cell_partitions_sim*cellsInSeries*legsInParallel;
capacity_loss_rate_mean_ess_partitions_sim = capacity_loss_rate_mean_cell_partitions_sim*legsInParallel;

energy_loss_rate_mean_ess_partitions_sim = zeros(size(energy_loss_rate_mean_bat_partitions_sim));
energy_loss_ess_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Wh
capacity_loss_ess_partitions_sim = nan*zeros(soc_num_sim,pow_num_sim,deglifePartitions_num); %in Ah

for partition_idx = 1:deglifePartitions_num
    for soc_bin_idx_sim = 1:soc_num_sim
        for pow_idx_sim = 1:pow_num_sim
            if(simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) <=simTimeRatioThreshold)
                possible_state_transition = 0;
            else
                possible_state_transition = 1;
            end
            if(possible_state_transition)                
                out_pow = out_pow_set_sim(pow_idx_sim);
                bat_pow = bat_pow_set_sim(pow_idx_sim);                
                
                energy_loss_rate_mean_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    energy_loss_rate_mean_bat_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) + ...
                    abs(out_pow-bat_pow);
                
                energy_loss_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    energy_loss_rate_mean_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;
                
                capacity_loss_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx) = ...
                    capacity_loss_rate_mean_ess_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    simTimeRatio_mean_cell_partitions_sim(soc_bin_idx_sim,pow_idx_sim,partition_idx)*...
                    slotIntervalInHours;                
            end
        end
    end
end


%% interpolation
energy_loss_ess_partitions_interpolated_1 = zeros(soc_num_sim,d_num,deglifePartitions_num);
energy_loss_ess_partitions_interpolated = zeros(z_num,d_num,deglifePartitions_num);
capacity_loss_ess_partitions_interpolated_1 = zeros(soc_num_sim,d_num,deglifePartitions_num);
capacity_loss_ess_partitions_interpolated = zeros(z_num,d_num,deglifePartitions_num);

for partition_idx = 1:deglifePartitions_num
    for soc_bin_idx_sim = 1:soc_num_sim
        Yt = out_pow_set_sim';
        Zt = reshape(energy_loss_ess_partitions_sim(soc_bin_idx_sim,:,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3 && ~isempty(find(Zt>0, 1)))
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,pow_idx_sim_temp]=findpeaks(K);
            K(pow_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),out_pow_set','linear','extrap');
            energy_loss_ess_partitions_interpolated_1(soc_bin_idx_sim,:,partition_idx) = Zt2;
        end
        Yt = out_pow_set_sim';
        Zt = reshape(capacity_loss_ess_partitions_sim(soc_bin_idx_sim,:,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3)
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,pow_idx_sim_temp]=findpeaks(K);
            K(pow_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),out_pow_set','linear','extrap');
            capacity_loss_ess_partitions_interpolated_1(soc_bin_idx_sim,:,partition_idx) = Zt2;
        end
    end
    
    for pow_idx = 1:d_num
        Yt = soc_grid_bin_mean_sim;
        Zt = reshape(energy_loss_ess_partitions_interpolated_1(:,pow_idx,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3 && ~isempty(find(Zt>0, 1)))
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,soc_bin_idx_sim_temp]=findpeaks(K);
            K(soc_bin_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),soc_grid_bin_mean,'linear','extrap');
            energy_loss_ess_partitions_interpolated(:,pow_idx,partition_idx) = max(Zt2,0);
        end
        Yt = soc_grid_bin_mean_sim;
        Zt = reshape(capacity_loss_ess_partitions_interpolated_1(:,pow_idx,partition_idx),[],1);
        Yt(isnan(Zt)) = [];
        Zt(isnan(Zt)) = [];
        validcount = length(Yt);
        if(validcount>3)
            DT = delaunayTriangulation(Yt,Zt);
            [K,~] = convexHull(DT);
            [~,soc_bin_idx_sim_temp]=findpeaks(K);
            K(soc_bin_idx_sim_temp+1:end) = [];
            Zt2 = interp1(Yt(K),Zt(K),soc_grid_bin_mean,'linear','extrap');
            capacity_loss_ess_partitions_interpolated(:,pow_idx,partition_idx) = max(Zt2,0);
        end
    end
end

essUsageCost_ess_partitions_interpolated = energyCostPer_Wh*energy_loss_ess_partitions_interpolated + 5*capacityCostPer_Ah*capacity_loss_ess_partitions_interpolated;
essUsageCost_ess_partitions_interpolated = round(max(essUsageCost_ess_partitions_interpolated,0),paramsPrecisionDigits);
energy_loss_ess_partitions_interpolated = round(max(energy_loss_ess_partitions_interpolated,0),paramsPrecisionDigits);
capacity_loss_ess_partitions_interpolated = round(max(capacity_loss_ess_partitions_interpolated,0),paramsPrecisionDigits);

z_kp1_idxs_map = zeros(z_num,d_num,deglifePartitions_num);
for partition_idx = 1:deglifePartitions_num      
    for z_idx=1:z_num
        z_cur = (z_offset+z_idx)*e_pu;
        for d_idx=1:d_num
            out_pow = (d_idx+d_offset)*p_pu;
            en_loss_ess = energy_loss_ess_partitions_interpolated(z_idx,d_idx,partition_idx);
            z_next = z_cur + out_pow*slotIntervalInHours - en_loss_ess;
            z_next_idx = round(z_next/e_pu)-z_offset;
            if(z_next_idx>z_num)
                z_dif = z_next-(z_num+z_offset)*e_pu;
                if(z_dif <e_pu/2)
                    possible_state_transition = 1;
                    z_next_idx = z_num;
                else
                    possible_state_transition = 0;
                end
            elseif(z_next_idx<1)
                z_dif = (z_offset+1)*e_pu - z_next;
                if(z_dif <e_pu/2)
                    possible_state_transition = 1;
                    z_next_idx = 1;
                else
                    possible_state_transition = 0;
                end
            else
                possible_state_transition = 1;
            end
            if(possible_state_transition)
                z_kp1_idxs_map(z_idx,d_idx,partition_idx) = z_next_idx;
            else
                z_kp1_idxs_map(z_idx,d_idx,partition_idx) = nan;
                capacity_loss_ess_partitions_interpolated(z_idx,d_idx,partition_idx) = nan;
                energy_loss_ess_partitions_interpolated(z_idx,d_idx,partition_idx) = nan;
                essUsageCost_ess_partitions_interpolated(z_idx,d_idx,partition_idx) = nan;
            end            
        end
    end    
end

% Prepare params
controller_essParams{deglifePartitions_num} = struct;
for partition_idx = 1:deglifePartitions_num
    controller_essParams{partition_idx}.batteryRatedCapacityInAh = batteryRatedCapacityInAh;
    controller_essParams{partition_idx}.cellsInSeries = cellsInSeries;
    controller_essParams{partition_idx}.legsInParallel = legsInParallel;
    controller_essParams{partition_idx}.capacityLossInAh_map = capacity_loss_ess_partitions_interpolated(:,:,partition_idx);
    controller_essParams{partition_idx}.energyLossInWh_map = energy_loss_ess_partitions_interpolated(:,:,partition_idx);
    controller_essParams{partition_idx}.z_kp1_idxs_map = z_kp1_idxs_map(:,:,partition_idx);
    controller_essParams{partition_idx}.essUsageCost_map = essUsageCost_ess_partitions_interpolated(:,:,partition_idx);
end
end
