function [] = visHybridSim(optStruct)

data = optStruct.data;
opt = optStruct.opt;
wave = optStruct.wave;
atmo = optStruct.atmo;
inso = optStruct.inso;
uc = optStruct.uc;
cturb = optStruct.cturb;
[data, opt] = prepHybrid(data,opt,uc,wave,atmo,inso,cturb);
output = optStruct.output;
length(data.met.time)
%extend time values
if ~isequal(length(output.min.Pdies),length(data.met.time))
    orig_l = length(data.met.time);
    vecMid = datevec(data.met.time(end));
    data.met.time = [data.met.time ; zeros(length(output.min.Pdies) ...
        - length(data.met.time),1)];
    for t = orig_l+1:length(data.met.time)
        vec = vecMid;
        vec(4) = vecMid(4) + t - orig_l;
        data.met.time(t) = datenum(vec);
    end
end

%I think this is redundant
% [data.met.wind_spd,data.met.time] = ...
%     extendToLifetime(data.met.wind_spd,data.met.time, ...
%     optStruct.uc.lifetime);

figure(1)
set(gcf,'Units','inches')
set(gcf,'Position', [1, 1, 4.2, 3.5])
set(gcf,'Color','w')
fs3 = 8;
col = colormap(brewermap(11,'Spectral')); %colors
col2 = colormap(brewermap(9,'Greys')); %colors
col3 = colormap(brewermap(9,'RdPu')); %colors
%set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultTextInterpreter', 'latex')
set(0,'DefaultAxesTickLabelInterpreter', 'latex')
set(0,'DefaultLegendInterpreter', 'latex')
%set(gcf, 'TickLabelInterpreter','latex')
%set(gcf, 'Position', [20, 20, 1300, 700])
tiledlayout(7,1,'TileSpacing','compact')

%Set length of time series
ts = 1:1:8760;
%ts = 1:1:length(data.met.time);
t1 = datetime(data.met.time(1),'ConvertFrom','datenum');
t2 = datetime(data.met.time(ts(end)),'ConvertFrom','datenum');
%STORAGE TIME SERIES
ax(1) = nexttile;
plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
    output.min.S1(1:ts(end))/1000,'Color',col(1,:), ... 
    'DisplayName','Battery 1','LineWidth',1)
%[255,69,0]/256
% hold on
% plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
%     output.min.S2(1:ts(end))/1000,'Color',col(3,:), ... 
%     'DisplayName','Battery 2','LineWidth',1)
% %[245,238,54]/256
% l1 = legend('show','location','southwest','Orientation','horizontal');
ylabel({'$S(t)$', '[kWh]'},'interpreter','latex')
% pos = get(l1,'Position');
% set(l1,'Position',pos + [0.347,0.037,0,0]);
drawnow
ylh = get(ax(1),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ylim([0 1.1*max(output.min.S2(1:end-1)/1000)])
yticks([0 round(max(output.min.S1(1:end-1)/1000),2)])
xlim([t1 t2])
set(ax(1),'XTick',[])
set(gca,'FontSize',6)
%set(gca,'LineWidth',2)
grid on
%add axis annotation
text(1.02, 0.5,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')

%POWER TIME SERIES
if opt.pm==5
    ax(2) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.Pdies(ts)/1000,'Color',col2(4,:), ... 
        'DisplayName','Diesel Power','LineWidth',1)
    %hold on
    %yline(output.min.kWd{end},':','Color',col2(4,:),'LineWidth',2)
    %[161,65,225]/256 
    %legend('show')
    ylabel({'$P_{\mbox{diesel}}$', '[kW]'},'interpreter','latex')
    ylh = get(ax(2),'ylabel');
    set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
        'VerticalAlignment','middle', ...
        'HorizontalAlignment','center')
    ylim([0 1.1*max(output.min.Pdies/1000)])
    xlim([t1 t2])
    set(ax(2),'XTick',[])
    set(gca,'FontSize',6)
    %set(gca,'LineWidth',2)
    grid on
    %add axis annotation
    text(1.02, 0.5,'(b)','Units','Normalized', ...
        'VerticalAlignment','middle','FontWeight','normal', ...
        'FontSize',fs3,'Interpreter','latex')
elseif opt.pm == 4
    ax(2) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.Pcurr(ts)/1000,'Color',col(11,:), ... 
        'DisplayName','Current Power','LineWidth',1)
    %hold on
    %yline(output.min.kWc{end},':','Color',col(11,:),'LineWidth',2)
    %[161,65,225]/256 
    %legend('show')
    ylabel({'$P_{\mbox{current}}$', '[kW]'},'interpreter','latex')
    ylh = get(ax(2),'ylabel');
    set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
        'VerticalAlignment','middle', ...
        'HorizontalAlignment','center')
    ylim([0 1.1*max(output.min.Pcurr/1000)])
    xlim([t1 t2])
    set(ax(2),'XTick',[])
    set(gca,'FontSize',6)
    %set(gca,'LineWidth',2)
    grid on
    %add axis annotation
    text(1.02, 0.5,'(b)','Units','Normalized', ...
        'VerticalAlignment','middle','FontWeight','normal', ...
        'FontSize',fs3,'Interpreter','latex')
end


%POWER TIME SERIES
ax(3) = nexttile;
plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
    output.min.Pinso(ts)/1000,'Color',col(5,:), ... 
    'DisplayName','Solar Power','LineWidth',1)
%hold on
%yline(output.min.kWi{end},':','Color',col(5,:),'LineWidth',2)
%[225,177,65]/256 
%legend('show')
ylabel({'$P_{\mbox{solar}}$', '[kW]'},'interpreter','latex')
ylh = get(ax(3),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ylim([0 1.1*max(output.min.Pinso/1000)])
xlim([t1 t2])
set(ax(3),'XTick',[])
set(gca,'FontSize',6)
%set(gca,'LineWidth',2)
grid on
%add axis annotation
text(1.02, 0.5,'(c)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')

%POWER TIME SERIES
ax(4) = nexttile;
plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
    output.min.Pwave(ts)/1000,'Color',col(10,:), ... 
    'DisplayName','Wave Power','LineWidth',1)
%hold on
%yline(output.min.kWwa{end},':','Color',col(10,:),'LineWidth',2)
%[65,105,225]/256 
%legend('show')
ylabel({'$P_{\mbox{wave}}$', '[kW]'},'interpreter','latex')
ylh = get(ax(4),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ylim([0 1.1*max(output.min.Pwave/1000)])
xlim([t1 t2])
set(ax(4),'XTick',[])
set(gca,'FontSize',6)
%set(gca,'LineWidth',2)
grid on
%add axis annotation
text(1.02, 0.5,'(d)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')

%POWER TIME SERIES
ax(5) = nexttile;
plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
    output.min.Pwind(ts)/1000,'Color',col(9,:), ... 
    'DisplayName','Wind Power','LineWidth',1)
%hold on
%yline(output.min.kWwi{end},':','Color',col(9,:),'LineWidth',2)
%[6,139,33]/256 
%%legend('show')
ylabel({'$P_{\mbox{wind}}$', '[kW]'},'interpreter','latex')
ylh = get(ax(5),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ylim([0 1.1*max(output.min.Pwind/1000)])
xlim([t1 t2])
set(ax(5),'XTick',[])
set(gca,'FontSize',6)
%set(gca,'LineWidth',2)
grid on
%add axis annotation
text(1.02, 0.5,'(e)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')

% %Current Power Time Series
% ax(6) = subplot(8,1,6);
% plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
%     output.min.Pcurr/1000,'Color',[6,139,33]/256, ... 
%     'DisplayName','Current Power','LineWidth',2)
% legend('show')
% ylabel('[kW]')
% xticks([])
% set(gca,'FontSize',10)
% set(gca,'LineWidth',2)
% grid on

%DUMPED POWER TIME SERIES
ax(6) = nexttile;
plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
    output.min.D(ts)/1000,'Color',col3(5,:), ... 
    'DisplayName','Power Dumped','LineWidth',1)
%[199,4,128]/256
%legend('show')
ylabel({'$P_{\mbox{discard}}$', '[kW]'},'interpreter','latex')
ylh = get(ax(6),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ylim([0 1.1*max(output.min.D/1000)])
xlim([t1 t2])
set(ax(6),'XTick',[])
set(gca,'FontSize',6)
%set(gca,'LineWidth',2)
grid on
%add axis annotation
text(1.02, 0.5,'(f)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')

%Load Time Series
ax(7) = nexttile;
plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
    output.min.L(ts)/1000,'Color',col2(7,:), ... 
    'DisplayName','Load','LineWidth',1)
%[212,92,176]/
%legend('show')
ylabel({'$L(t)$', '[kW]'},'interpreter','latex')
ylh = get(ax(7),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ylim([0 1.1*max(output.min.L/1000)])
yticks([0 1])
%xticks([])
t1 = datetime(data.met.time(1),'ConvertFrom','datenum');
t2 = datetime(data.met.time(ts(end)),'ConvertFrom','datenum');
xlim([t1 t2])
xlabel('Time','interpreter','latex')
set(gca,'FontSize',6)
%set(gca,'LineWidth',2)
grid on
%add axis annotation
text(1.02, 0.5,'(g)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')


% %Failure Time Series
% ax(8) = subplot(8,1,8);
% plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
%     output.min.F,'.','Color',[0,0,0], ... 
%     'DisplayName','Failures','LineWidth',2)
% legend('show')
% ylabel('[logical]')
% ylim([0,2])
% yticks([0 1])
% xlabel('Time')
% set(gca,'FontSize',10)
% set(gca,'LineWidth',2)
% grid on
% %EFFICIENCY TIME SERIES
% ax(4) = subplot(4,1,4);
% plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
%     output.min.eff_t,'Color',[255,20,147]/256, ... 
%     'DisplayName','Efficiency','LineWidth',2)
% legend('show')
% ylabel('[~]')
% xlabel('Time')
% ylim([0 1.25*max(output.min.eff_t)])
% set(gca,'FontSize',16)
% set(gca,'LineWidth',2)
% grid on
% 
% set(gcf, 'Position', [100, 100, 1400, 650])

linkaxes(ax,'x')
%linkaxes(ax(2:6),'y')

print(figure(1), '-dpng','-r600')
end

