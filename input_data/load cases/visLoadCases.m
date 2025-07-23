function visLoadCases()
%Plot all load case time series on one plot

col = brewermap(5,'Set2'); %colors

%generate a dummy time vector
t1 = datetime(2010,1,1,0,0,0);
t2 = datetime(2016,1,1,0,0,0);
t1y = datetime(2010,3,1,0,0,0);
data.met.time = t1:hours(1):t2;
[load_case,load_series] = GenerateLoadCases_v4(data);
figure
set(gcf,'Position', [100, 100, 900, 500])
tf = tiledlayout(5,1);
tf.Padding = 'compact';
tf.TileSpacing = 'tight';
for i = 1:5
    nexttile
    plot(data.met.time,load_series.L(i,:),'LineWidth',2,'Color',col(i,:))
    xlim([t1,t1y])
    ylim([0,1000])
    ylabel('[W]')
    title(strcat("Load Case: ",string(i)))
end
xlabel('Time [hr]')

