function [] = plotFigures__optimal_control(plotData)
controller_Params = plotData.controller_Params;
evaluationSMdata = plotData.evaluationSMdata;
soc_grid_bin_mean = controller_Params.soc_grid_bin_mean;

savefiles = 0;

plot_optimal_privacy_control = 1;
plot_optimal_demand_shaping = 1;

if (plot_optimal_privacy_control)
    filename = {'optimal_privacy_control'};
    simulatedControllerData = plotData.simulatedControllerData_privacy_control;
    numPlotDays = 10;
    
    k_num_in_day = controller_Params.k_num_in_day;
    
    unmodified_sm_data_vec = reshape(evaluationSMdata(:,1:numPlotDays), [], 1);        
    modified_sm_data = simulatedControllerData.modifiedSMdata(:,1:numPlotDays);    
    modified_sm_data_vec = reshape(modified_sm_data, [], 1);
        
    SOC_data = soc_grid_bin_mean(simulatedControllerData.z_k_idxs(:,1:numPlotDays));
    SOC_data_vec = reshape(SOC_data, [], 1);
    bayesRiskAveragedInHorizon_privacy_control = plotData.bayesRiskAveragedInHorizon_privacy_control;
    bayesRiskAveraged_vec = reshape(bayesRiskAveragedInHorizon_privacy_control(:,1:numPlotDays), [], 1);
    
    xlabel_strs = cell(1,numPlotDays);
    for day_idx = 1:numPlotDays
        xlabel_strs{day_idx} = num2str((day_idx));
    end
    
    figure_xsize = 800;
    figure_ysize = 600;
    fontSize = 16;
    
    % plot figure--
    x_data = 1:length(unmodified_sm_data_vec);
    figure_position_offset = 50;
    figure1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    axes1 = axes('Parent',figure1,...
        'Position',[0.13 0.735384615384615 0.775 0.217948717948718]);
    hold(axes1,'on');
    plot(x_data,unmodified_sm_data_vec,'DisplayName',' User demand','linewidth',1.5)
    plot(x_data,modified_sm_data_vec,'linewidth',2,'DisplayName',' SM readings')
    ylabel({'User demand (W)'},'Interpreter','latex','position',[-22,1044.3,-1.0000]);
    xlabel({'Days'},'Interpreter','latex');
    box(axes1,'on');
    grid(axes1,'on');
    ylim([0 2000])
    xlim([1 length(unmodified_sm_data_vec)])
    xticks(unique([1:k_num_in_day:length(unmodified_sm_data_vec),length(unmodified_sm_data_vec)]))
    xticklabels(xlabel_strs)
    legend_1 = legend('show');
    set(legend_1,'Interpreter','latex','FontSize',fontSize-2,'Location','northwest');
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 1:k_num_in_day/4:length(k_num_in_day);
    set(axes1,'TickLabelInterpreter','latex','FontSize',fontSize,'YMinorTick','on','YTick',...
        0:500:2000);
    
    axes2 = axes('Parent',figure1,...
        'Position',[0.13 0.426153846153846 0.775 0.187692307692307]);
    hold(axes2,'on');
    plot(x_data,SOC_data_vec,'linewidth',1.5)
    ylabel({'ESS SoC'}, 'position', [-22,0.5,-1],'Interpreter','latex');
    xlabel({'Days'},'Interpreter','latex');
    box(axes2,'on');
    grid(axes2,'on');
    xlim([1 length(unmodified_sm_data_vec)])
    xticks(unique([1:k_num_in_day:length(unmodified_sm_data_vec),length(unmodified_sm_data_vec)]))
    xticklabels(xlabel_strs)
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 1:k_num_in_day/4:length(unmodified_sm_data_vec);
    set(axes2,'TickLabelInterpreter','latex','FontSize',fontSize,'YMinorTick','on','YTick',...
        0:0.25:1);    
        
    axes3 = axes('Parent',figure1,...
        'Position',[0.13 0.109230769230769 0.775 0.201538461538461]);
    hold(axes3,'on');
    plot(x_data,bayesRiskAveraged_vec,'linewidth',1.5)
    ylabel({'Bayesian risk'}, 'position', [-22,0.5,-1],'Interpreter','latex');
    xlabel({'Days'},'Interpreter','latex');
    box(axes3,'on');
    grid(axes3,'on');
    ylim([0 1])
    xlim([1 length(unmodified_sm_data_vec)])
    xticks(unique([1:k_num_in_day:length(unmodified_sm_data_vec),length(unmodified_sm_data_vec)]))
    xticklabels(xlabel_strs)
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 1:k_num_in_day/4:length(unmodified_sm_data_vec);
    set(axes3,'TickLabelInterpreter','latex','FontSize',fontSize,'YMinorTick','on','YTick',...
        0:0.25:1);    
    
    fig_pos = figure1.PaperPosition;
    figure1.PaperSize = [fig_pos(3) fig_pos(4)];
    
    if(savefiles)
        saveas(gcf,strcat(filename{1},'.fig'));
        saveas(gcf,strcat(filename{1},'.eps'),'epsc');
    end
end

if (plot_optimal_demand_shaping)
    filename = {'optimal_demand_shaping'};
    simulatedControllerData = plotData.simulatedControllerData_demand_shaping;
    numPlotDays = 10;
    
    k_num_in_day = controller_Params.k_num_in_day;
    
    unmodified_sm_data_vec = reshape(evaluationSMdata(:,1:numPlotDays), [], 1);    
    
    modified_sm_data = simulatedControllerData.modifiedSMdata(:,1:numPlotDays);
    modified_sm_data_vec = reshape(modified_sm_data, [], 1);
        
    SOC_data = soc_grid_bin_mean(simulatedControllerData.z_k_idxs(:,1:numPlotDays));
    SOC_data_vec = reshape(SOC_data, [], 1);    
    
    bayesRiskAveragedInHorizon_demand_shaping = plotData.bayesRiskAveragedInHorizon_demand_shaping;
    bayesRiskAveraged_vec = reshape(bayesRiskAveragedInHorizon_demand_shaping(:,1:numPlotDays), [], 1);
    
    xlabel_strs = cell(1,numPlotDays);
    for day_idx = 1:numPlotDays
        xlabel_strs{day_idx} = num2str((day_idx));
    end
    
    figure_xsize = 800;
    figure_ysize = 600;
    fontSize = 16;
    
    % plot figure--
    x_data = 1:length(unmodified_sm_data_vec);
    figure_position_offset = 50;
    figure1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    
    axes1 = axes('Parent',figure1,...
        'Position',[0.13 0.735384615384615 0.775 0.217948717948718]);
    hold(axes1,'on');
    plot(x_data,unmodified_sm_data_vec,'DisplayName',' User demand','linewidth',1.5)
    plot(x_data,modified_sm_data_vec,'linewidth',2,'DisplayName',' SM readings')
    ylabel({'User demand (W)'},'Interpreter','latex','position',[-22,1044.3,-1.0000]);
    xlabel({'Days'},'Interpreter','latex');
    box(axes1,'on');
    grid(axes1,'on');
    ylim([0 2000])
    xlim([1 length(unmodified_sm_data_vec)])
    xticks(unique([1:k_num_in_day:length(unmodified_sm_data_vec),length(unmodified_sm_data_vec)]))
    xticklabels(xlabel_strs)
    legend_1 = legend('show');
    set(legend_1,'Interpreter','latex','FontSize',fontSize-2,'Location','northwest');
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 1:k_num_in_day/4:length(k_num_in_day);
    set(axes1,'TickLabelInterpreter','latex','FontSize',fontSize,'YMinorTick','on','YTick',...
        0:500:2000); 
    
    axes2 = axes('Parent',figure1,...
        'Position',[0.13 0.426153846153846 0.775 0.187692307692307]);
    hold(axes2,'on');
    plot(x_data,SOC_data_vec,'linewidth',1.5)
    ylabel({'ESS SoC'}, 'position', [-22,0.5,-1],'Interpreter','latex');
    xlabel({'Days'},'Interpreter','latex');
    box(axes2,'on');
    grid(axes2,'on');
    xlim([1 length(unmodified_sm_data_vec)])
    xticks(unique([1:k_num_in_day:length(unmodified_sm_data_vec),length(unmodified_sm_data_vec)]))
    xticklabels(xlabel_strs)
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 1:k_num_in_day/4:length(unmodified_sm_data_vec);
    set(axes2,'TickLabelInterpreter','latex','FontSize',fontSize,'YMinorTick','on','YTick',...
        0:0.25:1);    
    
    axes3 = axes('Parent',figure1,...
        'Position',[0.13 0.109230769230769 0.775 0.201538461538461]);
    hold(axes3,'on');
    plot(x_data,bayesRiskAveraged_vec,'linewidth',1.5)
    ylabel({'Bayesian risk'}, 'position', [-22,0.5,-1],'Interpreter','latex');
    xlabel({'Days'},'Interpreter','latex');
    box(axes3,'on');
    grid(axes3,'on');
    xlim([1 length(unmodified_sm_data_vec)])
    xticks(unique([1:k_num_in_day:length(unmodified_sm_data_vec),length(unmodified_sm_data_vec)]))
    xticklabels(xlabel_strs)
    ax = gca;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 1:k_num_in_day/4:length(unmodified_sm_data_vec);
    set(axes3,'TickLabelInterpreter','latex','FontSize',fontSize,'YMinorTick','on','YTick',...
        0:0.25:1);    
    
    fig_pos = figure1.PaperPosition;
    figure1.PaperSize = [fig_pos(3) fig_pos(4)];
    
    if(savefiles)
        saveas(gcf,strcat(filename{1},'.fig'));
        saveas(gcf,strcat(filename{1},'.eps'),'epsc');
    end
end

end

