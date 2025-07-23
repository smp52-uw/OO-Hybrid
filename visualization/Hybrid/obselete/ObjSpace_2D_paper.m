%Make 2x2 2D vis plot
datacomp = figure;
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 8.5, 4.25])
set(gcf,'Color','w')
tiledlayout(2,2)
set(0,'defaulttextinterpreter','latex')
fs3 = 10; %annotation font size
%solar
nexttile
vis2DObjSpace(argBasin_inso,2,8)
ylabel('$G_{r,solar}$ [kW]')
text(0.85, 0.85,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex','Color','w')
%wind
nexttile
vis2DObjSpace(argBasin_wind,1,8)
ylabel('$G_{r,wind}$ [kW]')
text(0.85, 0.85,'(b)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex','Color','w')
%wave
nexttile
vis2DObjSpace(argBasin_wave,3,8)
ylabel('$G_{r,wave}$ [kW]')
text(0.85, 0.85,'(c)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex','Color','w')
%curr
nexttile
vis2DObjSpace(argBasin_curr15,5,15)
ylabel('$G_{r,current}$ [kW]')
text(0.85, 0.85,'(d)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex','Color','w')


print(datacomp,'C:\Users\smpal\Documents\GitHub\EWTEC_fig\2DResults_fig',  ...
    '-dpng','-r600')