function [] = plotFigures__deviation_privacy_tradeoff(plotData)
tradeOff_omega_withinUtility = plotData.tradeOff_omega_withinUtility;
batteryRatedCapacityInAh_values = plotData.batteryRatedCapacityInAh_values;
rmsd = plotData.rmsd;
risk = plotData.risk;
numSimulatedDays = plotData.numSimulatedDays;
cost = plotData.cost;

savefiles = 0;

plot_rmsd_vs_omega_vs_capacity = 1;
plot_risk_vs_omega_vs_capacity = 1;
plot_life_vs_omega_vs_capacity = 1;
plot_cost_vs_omega_vs_capacity = 1;

if (plot_rmsd_vs_omega_vs_capacity)
    filename = {'tradeoff_rmsd_omega_cap'};
    
    axis_position = [0.189230769230769 0.199999542236328 0.530619441804377 0.762500457763672];
    colorbar_position = [0.866769230769231,0.253055555555556,0.05,0.6775];
    
    colorbar_label_position = [-1.688334779305892,327.0963435861679,0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    z_data = rmsd;
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(z_data);
    xlabel({'Capacity (Ah)'},'Interpreter','latex','Units', 'normalized','rotation',36,'position',[0.5638,-0.195,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-36,'position',[0.372676811594203,-0.185380758017492,0]);
    zlabel({'RMSD (W)'},'Interpreter','latex','Units','normalized','position',[-0.2031,0.5228,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    zlim([0 600]);
    zticks(0:200:600)
    xlim([-0.5 6.5])
    ylim([-0.5 7.5])
    axes1 = gca;
    yticks([1,length(tradeOff_omega_withinUtility)]);
    yticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(batteryRatedCapacityInAh_values);
    xticks(1:2:(length(batteryRatedCapacityInAh_values)));
    xticklabels(batteryRatedCapacityInAh_values(1:2:end))
    h_cb = colorbar('peer',axes1,'Position',colorbar_position,...
        'TickLabelInterpreter','latex','FontSize',fontSize-6,'Limits',[0,600],'Ticks',0:200:600);
    ylabel(h_cb, 'RMSD (W)','Interpreter','latex','FontSize',fontSize,...
        'Position',colorbar_label_position);
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'Position',axis_position,'FontSize',fontSize,'TickLabelInterpreter','latex','XGrid','on',...
        'XMinorTick','on','YMinorTick','on','ZMinorTick','on','YGrid','on','ZGrid','on');
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    
    if(savefiles)
        saveas(gcf,strcat(filename{1},'.fig'));
        saveas(gcf,strcat(filename{1},'.eps'),'epsc');
    end
end

if (plot_risk_vs_omega_vs_capacity)
    filename = {'tradeoff_risk_omega_cap'};
    
    axis_position = [0.189230769230769 0.199999542236328 0.530619441804377 0.762500457763672];
    colorbar_position = [0.866769230769231,0.253055555555556,0.05,0.6775];
    
    colorbar_label_position = [-1.233789357272062,0.538966541212113,0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    z_data = risk;
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(z_data);
    xlabel({'Capacity (Ah)'},'Interpreter','latex','Units', 'normalized','rotation',36,'position',[0.5638,-0.195,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-36,'position',[0.372676811594203,-0.185380758017492,0]);
    zlabel({'Bayesian risk'},'Interpreter','latex','Units','normalized','position',[-0.222961359155138,0.513044807940734,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    zlim([0 1]);
    zticks(0:0.3:1)
    xlim([-0.5 6.5])
    ylim([-0.5 7.5])
    axes1 = gca;
    yticks([1,length(tradeOff_omega_withinUtility)]);
    yticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(batteryRatedCapacityInAh_values);
    xticks(1:2:(length(batteryRatedCapacityInAh_values)));
    xticklabels(batteryRatedCapacityInAh_values(1:2:end))
    h_cb = colorbar('peer',axes1,'Position',colorbar_position,...
        'TickLabelInterpreter','latex','FontSize',fontSize-6,'Limits',[0,1],'Ticks',0:0.25:1);
    ylabel(h_cb, 'Bayesian risk','Interpreter','latex','FontSize',fontSize,...
        'Position',colorbar_label_position);
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'Position',axis_position,'FontSize',fontSize,'TickLabelInterpreter','latex','XGrid','on',...
        'XMinorTick','on','ZMinorTick','on','YGrid','on','ZGrid','on');
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    
    if(savefiles)
        saveas(gcf,strcat(filename{1},'.fig'));
        saveas(gcf,strcat(filename{1},'.eps'),'epsc');
    end
end

if (plot_life_vs_omega_vs_capacity)
    filename = {'tradeoff_life_omega_cap'};
    
    axis_position = [0.2 0.199999542236328 0.519850211035146 0.762500457763672];
    colorbar_position = [0.849846153846154,0.264166666666667,0.05,0.6775];
    
    %     colorbar_label_position = [-1.476213603308707,1244.801350851528,0];
    colorbar_label_position = [-1.47621357440948 697.916084664767 0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    z_data = numSimulatedDays;
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(z_data);
    xlabel({'Capacity (Ah)'},'Interpreter','latex','Units', 'normalized','rotation',36,'position',[0.5638,-0.195,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-36,'position',[0.372676811594203,-0.185380758017492,0]);
    zlabel({'ESS life (in days)'},'Interpreter','latex','Units','normalized','position',[-0.243968228158067,0.622284036438086,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    zlimdata = [0 1200];
    zticksdata = 0:400:1200;
    zlim(zlimdata);
    zticks(zticksdata)
    
    xlim([-0.5 6.5])
    ylim([-0.5 7.5])
    axes1 = gca;
    yticks([1,length(tradeOff_omega_withinUtility)]);
    yticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(batteryRatedCapacityInAh_values);
    xticks(1:2:(length(batteryRatedCapacityInAh_values)));
    xticklabels(batteryRatedCapacityInAh_values(1:2:end))
    h_cb = colorbar('peer',axes1,'Position',colorbar_position,...
        'TickLabelInterpreter','latex','FontSize',fontSize-6,'Limits',zlimdata,'Ticks',zticksdata);
    ylabel(h_cb, 'ESS life (in days)','Interpreter','latex','FontSize',fontSize,...
        'Position',colorbar_label_position);
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'Position',axis_position,'FontSize',fontSize,'TickLabelInterpreter','latex','XGrid','on',...
        'YMinorTick','on','ZMinorTick','on','YGrid','on','ZGrid','on');
    
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    
    if(savefiles)
        saveas(gcf,strcat(filename{1},'.fig'));
        saveas(gcf,strcat(filename{1},'.eps'),'epsc');
    end
end

if (plot_cost_vs_omega_vs_capacity)
    filename = {'tradeoff_cost_omega_cap'};
    
    axis_position = [0.189230769230769 0.199999542236328 0.530619441804377 0.762500457763672];
    colorbar_position = [0.860615384615385,0.284166666666667,0.05,0.6775];
    
    colorbar_label_position = [-1.415607513803424,0.271753430171091,0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    z_data = cost;
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(z_data);
    xlabel({'Capacity (Ah)'},'Interpreter','latex','Units', 'normalized','rotation',36,'position',[0.5638,-0.195,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-36,'position',[0.372676811594203,-0.185380758017492,0]);
    zlabel({'ESS cost (euro/day)'},'Interpreter','latex','Units','normalized','position',[-0.2109,0.5814,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    
    zlimdata = [0 0.5];
    zticksdata = 0:0.25:0.5;
    
    zlim(zlimdata);
    zticks(zticksdata)
    xlim([-0.5 6.5])
    ylim([-0.5 7.5])
    axes1 = gca;
    yticks([1,length(tradeOff_omega_withinUtility)]);
    yticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(batteryRatedCapacityInAh_values);
    xticks(1:2:(length(batteryRatedCapacityInAh_values)));
    xticklabels(batteryRatedCapacityInAh_values(1:2:end))
    h_cb = colorbar('peer',axes1,'Position',colorbar_position,...
        'TickLabelInterpreter','latex','FontSize',fontSize-6,'Limits',zlimdata,'Ticks',zticksdata);
    ylabel(h_cb, 'ESS cost (euro/day)','Interpreter','latex','FontSize',fontSize,...
        'Position',colorbar_label_position);
    box(axes1,'on');
    grid(axes1,'on');
    set(axes1,'Position',axis_position,'FontSize',fontSize,'TickLabelInterpreter','latex','XGrid','on',...
        'YMinorTick','on','ZMinorTick','on','YGrid','on','ZGrid','on');
    
    fig_pos = figure_1.PaperPosition;
    figure_1.PaperSize = [fig_pos(3) fig_pos(4)];
    
    if(savefiles)
        saveas(gcf,strcat(filename{1},'.fig'));
        saveas(gcf,strcat(filename{1},'.eps'),'epsc');
    end
end

end

