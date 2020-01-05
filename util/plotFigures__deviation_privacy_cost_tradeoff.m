function [] = plotFigures__deviation_privacy_cost_tradeoff(plotData)
tradeOff_omega_withinUtility = plotData.tradeOff_omega_withinUtility;
tradeOff_sigma_forcost = plotData.tradeOff_sigma_forcost;
rmsd = plotData.rmsd;
risk = plotData.risk;
numSimulatedDays = plotData.numSimulatedDays;
cost = plotData.cost;

savefiles = 0;

plot_rmsd_vs_omega_vs_sigma = 1;
plot_risk_vs_omega_vs_sigma = 1;
plot_life_vs_omega_vs_sigma = 1;
plot_cost_vs_omega_vs_sigma = 1;

if (plot_rmsd_vs_omega_vs_sigma)
    filename = {'tradeoff_rmsd'};
    
    axis_position = [0.189230769230769 0.199999542236328 0.530619441804377 0.762500457763672];
    colorbar_position = [0.866769230769231,0.253055555555556,0.05,0.6775];
    
    colorbar_label_position = [-1.53681962779074,350.702900963217,0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    x_data = 1:length(tradeOff_sigma_forcost);
    
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    xlim([-0.5 7.5])
    ylim([-0.5 7.5])
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(rmsd);
    xlabel({'Trade-off ($\sigma$)'},'Interpreter','latex','Units', 'normalized','rotation',37,'position',[0.6015,-0.1422,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-34,'position',[0.372676811594203,-0.138772594752187,0]);
    zlabel({'RMSD (W)'},'Interpreter','latex','Units','normalized','position',[-0.2031,0.5228,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    zlim([0 600]);
    zticks(0:200:600)
    xticks([1,length(x_data)]);
    xticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(x_data);
    yticks([1,length(x_data)])
    yticklabels([0,1])
    axes1.YAxis.MinorTickValues = 1:1:length(x_data);
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

if (plot_risk_vs_omega_vs_sigma)
    filename = {'tradeoff_risk'};
    
    axis_position = [0.189230769230769 0.199999542236328 0.530619441804377 0.762500457763672];
    colorbar_position = [0.866769230769231,0.253055555555556,0.05,0.6775];
    
    colorbar_label_position = [-1.233789357272062,0.538966541212113,0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    x_data = 1:length(tradeOff_sigma_forcost);
    
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    xlim([-0.5 7.5])
    ylim([-0.5 7.5])
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(risk);
    xlabel({'Trade-off ($\sigma$)'},'Interpreter','latex','Units', 'normalized','rotation',37,'position',[0.6015,-0.1422,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-34,'position',[0.372676811594203,-0.138772594752187,0]);
    zlabel({'Bayesian risk'},'Interpreter','latex','Units','normalized','position',[-0.222961359155138,0.513044807940734,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    zlim([0 1]);
    zticks(0:0.3:1)
    xticks([1,length(x_data)]);
    xticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(x_data);
    yticks([1,length(x_data)])
    yticklabels([0,1])
    axes1.YAxis.MinorTickValues = 1:1:length(x_data);
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

if (plot_life_vs_omega_vs_sigma)
    filename = {'tradeoff_life'};
    max_life_in_plot = 2100;
    axis_position = [0.196923252986027,0.199999542236328,0.522926958049119,0.762500457763672];
    colorbar_position = [0.843692307692308,0.259722222222223,0.05,0.6775];
    
    colorbar_label_position = [-1.56710302829742 1281.21317402261 0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    x_data = 1:length(tradeOff_sigma_forcost);
    
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    xlim([-0.5 7.5])
    ylim([-0.5 7.5])
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(min(numSimulatedDays,max_life_in_plot));
    xlabel({'Trade-off ($\sigma$)'},'Interpreter','latex','Units', 'normalized','rotation',37,'position',[0.6015,-0.1422,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-34,'position',[0.372676811594203,-0.138772594752187,0]);
    zlabel({'ESS life (in days)'},'Interpreter','latex','Units','normalized','position',[-0.240810391830843,0.598960421277736,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    caxis([0 max_life_in_plot]);
    zlim([0 max_life_in_plot]);
    zticks(0:700:max_life_in_plot)
    xticks([1,length(x_data)]);
    xticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(x_data);
    yticks([1,length(x_data)])
    yticklabels([0,1])
    axes1.YAxis.MinorTickValues = 1:1:length(x_data);
    h_cb = colorbar('peer',axes1,'Position',colorbar_position,...
        'TickLabelInterpreter','latex','FontSize',fontSize-6,'Limits',[0,max_life_in_plot],'Ticks',0:700:max_life_in_plot);
    ylabel(h_cb, 'ESS life (in days)','Interpreter','latex','FontSize',fontSize,...
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

if (plot_cost_vs_omega_vs_sigma)
    filename = {'tradeoff_cost'};
    
    axis_position = [0.189230769230769 0.199999542236328 0.530619441804377 0.762500457763672];
    colorbar_position = [0.866769230769231,0.253055555555556,0.05,0.6775];
    
    colorbar_label_position = [-1.65803170204163 0.302573074180572 0];
    fontSize = 24;
    figure_xsize = 600;
    figure_ysize = 400;
    figure_position_offset = 50;
    
    x_data = 1:length(tradeOff_sigma_forcost);
    
    % plot figure--
    figure_1 = figure('Color','w','Renderer','painters');
    set(gcf, 'Position',  [figure_position_offset,figure_position_offset,figure_position_offset+figure_xsize, figure_position_offset+figure_ysize]);
    hold on
    axes1 = gca;
    xlim([-0.5 7.5])
    ylim([-0.5 7.5])
    set(axes1, 'Xdir', 'reverse')
    set(axes1, 'Ydir', 'reverse')
    view(axes1,[134 56]);
    h1 = bar3(cost);
    xlabel({'Trade-off ($\sigma$)'},'Interpreter','latex','Units', 'normalized','rotation',37,'position',[0.6015,-0.1422,0]);
    ylabel({'Trade-off ($\omega$)'},'Interpreter','latex','Units','normalized','rotation',-34,'position',[0.372676811594203,-0.138772594752187,0]);
    zlabel({'ESS cost (euro/day)'},'Interpreter','latex','Units','normalized','position',[-0.210913262324433,0.581448716735523,0]);
    set(h1, {'CData'}, get(h1,'ZData'), 'FaceColor','interp','MeshStyle','column')
    caxis([0 0.5])
    zlim([0 0.5]);
    zticks(0:0.25:0.6)
    xticks([1,length(x_data)]);
    xticklabels([0,1])
    axes1.XAxis.MinorTickValues = 1:1:length(x_data);
    yticks([1,length(x_data)])
    yticklabels([0,1])
    axes1.YAxis.MinorTickValues = 1:1:length(x_data);
    h_cb = colorbar('peer',axes1,'Position',colorbar_position,...
        'TickLabelInterpreter','latex','FontSize',fontSize-6,'Limits',[0,0.5],'Ticks',0:0.25:0.5);
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

