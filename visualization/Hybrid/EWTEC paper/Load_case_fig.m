%Brian Polagye
%July 20, 2022
%Edits: Sarah Palmer 5/5/23

%Description: Generate time-resolved load cases
optInputs
%optInputs %load inputs
data = load(loc,loc);
data = data.(loc);
%load current data
curr = load(cloc);
data.curr = curr; %add current data to data structure

[data, opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);


[load_case,load_series] = GenerateLoadCases_v3(data);

%% Visualize load cases

L = load_series.L;
L_status = load_series.L_status ;
t = load_series.t ;
% 
% %Calculate maintenance interval
% dt = years(2);
% day_m1 = days(t(1)) + days(dt);
% day_m2 = days(t(1)) + 2*days(dt);
datacomp = figure;
set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultAxesFontName', 'cmr10')
set(gcf,'Color','w')
col = colormap(brewermap(5,'Pastel1')); %colors
fs = 6; %font size
set(gca,'FontSize',fs)
set(gcf,'Units','inches')
set(gcf,'Position', [1, 1, 5.5, 3.25])
f = tiledlayout(1,1);
ax(1) = nexttile(1,[1,1]);
i = 3;
tend = 500;
plot(days(t(1:tend)),L(i,1:tend),'-','linewidth',3,'Color',col(1,:))
hold on
yline(mean(L(i,1:tend)),'--','linewidth',3,'Color',col(2,:))
% 
% xline (day_m1,'-','linewidth',4)
% xline(day_m2,'-','linewidth',4)

%title(load_case(i).name,'fontweight','b')
yl = ylabel('$$L_c$$ [W]','Interpreter','latex');

grid on

xl = xlabel('Time [days]');
xlim([days(t(1)) days(t(tend))])
ax(1).TickLabelInterpreter = 'latex';

hl = legend('$$L_c(t)$$','$$L_{c,avg}$$','Interpreter','latex','location','northwest');
set(hl,'Interpreter','latex')
set(xl,'Interpreter','latex')
set(yl,'Interpreter','latex')

print(datacomp,'C:\Users\smpal\Documents\GitHub\EWTEC_fig\Load_fig_c',  ...
    '-dpng','-r600')

