function visLoadCases()
%Plot all load case time series on one plot

col = brewermap(5,'Set2'); %colors

%generate a dummy time vector
t1 = datetime(2010,1,1,0,0,0);
t2 = datetime(2016,1,1,0,0,0);
t1y = datetime(2010,3,1,0,0,0);
data.met.time = t1:hours(1):t2;
reltime = data.met.time - data.met.time(1);
[load_case,load_series] = GenerateLoadCases_v4(data);


fst = 10;
fs = 8;
figure
set(gcf,'Position', [100, 100, 400, 500])
tf = tiledlayout(3,1);
tf.Padding = 'compact';
tf.TileSpacing = 'tight';
load_case(5).name = 'Constant 322 W';
l = [1 3 5];
for j = 1:3
    i = l(j);
    ax(j) = nexttile;
    plot(reltime,load_series.L(i,:),'LineWidth',2,'Color','k')
    ylim([0,1000])
    ylh = ylabel('[W]','FontSize',fs,'Interpreter','latex');
    title(load_case(i).name,'FontSize',fst,'Interpreter','latex')
    xtickformat('h')
    xlim([reltime(1),reltime(1440)])
    grid on
    box on
    ax(j).TickLabelInterpreter = 'latex';
    set(gca,'FontSize',fs)
    set(ylh,'Rotation',0)
end
xlabel('Time [hr]','FontSize',fs,'Interpreter','latex')

