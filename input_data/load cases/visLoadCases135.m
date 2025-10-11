function visLoadCases135()
%Plot all load case time series on one plot

col = brewermap(5,'greys'); %colors
titlestr = {["High - Capacity"; "UUV Recharge"],"",["Next - generation Ocean"; "Observing and UUV Recharge"],"",["200 W"]};
%generate a dummy time vector
t1 = datetime(2010,1,1,0,0,0);
t2 = datetime(2016,1,1,0,0,0);
t1y = datetime(2010,3,1,0,0,0);
data.met.time = t1:hours(1):t2;
[load_case,load_series] = GenerateLoadCases_v4(data);
figure
set(gcf,'Units','inches','Position', [1, 1, 6.5, 4.5])
set(0,'DefaultTextFontname', 'times new roman')
set(0,'DefaultAxesFontName', 'times new roman')
tf = tiledlayout(3,1);
tf.Padding = 'compact';
tf.TileSpacing = 'tight';
fs = 10;
for i = 1:5
    if i == 1 || i == 3|| i == 5
        ax = nexttile;
        plot(data.met.time,load_series.L(i,:),'LineWidth',2,'Color',col(4,:))
        hold on
        mean(load_series.L(i,:))
        yline(mean(load_series.L(i,:)),'--','LineWidth',2,'Color',col(3,:))
        xlim([t1,t1y])
        ylim([0,1000])
        ylabel('[W]', 'FontSize',fs)
        %title(strcat("Load Case: ",string(i)))
        title(titlestr{i},'fontsize',fs)
        ax.FontSize = fs; 
        xtickformat('MMM-dd'); 
        ax.XAxis.SecondaryLabel.Visible = 'off';
        grid on
        box on
        if i == 1
            legend('Load','Average Load','fontsize',fs)
        end
    end
end
xlabel('Time [hr]','FontSize',fs)