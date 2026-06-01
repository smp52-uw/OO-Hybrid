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
fs = 10;
figure
set(gcf,'Units','inches')
set(gcf,'Position', [1, 1, 3.75, 4])
tf = tiledlayout(3,1);
tf.Padding = 'loose';
tf.TileSpacing = 'compact';
tf.InnerPosition = [0.17,0.13,0.75,0.8];
load_case(5).name = 'Constant 322 W';
load_case(3).name = {'Next-generation Ocean Observing and',' UUV recharge (OO + UUV Recharge)'};
load_case(1).name = {'High-capacity UUV Recharge',' (UUV Recharge)'};
tilelabel = {'(a)','(b)','(c)'};
l = [1 5 3];
for j = 1:3
    i = l(j);
    ax(j) = nexttile;
    plot(reltime,load_series.L(i,:),'LineWidth',2,'Color','k')
    ylim([0,1000])
    yl = ylabel({'$L_c$', '[W]'},'FontSize',fs,'Interpreter','latex','VerticalAlignment','middle','Rotation',0);
    set(yl,'units','normalized')
    ylp1 = yl.Position;
    ylp1(1) = ylp1(1)*1.5;
    set(yl,'position',ylp1)
    %title(load_case(i).name,'FontSize',fst,'Interpreter','latex')
    xtickformat('h')
    xlim([reltime(1),reltime(1440/2)])
    grid on
    box on
    ax(j).TickLabelInterpreter = 'latex';
    set(gca,'FontSize',fs)
    text(0.02, 0.85,tilelabel{j},'Units','Normalized', ...
    'VerticalAlignment','middle','FontSize',fs,'Interpreter','latex')
end
xlabel('Time [hr]','FontSize',fs,'Interpreter','latex')

