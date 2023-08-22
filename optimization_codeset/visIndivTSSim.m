function [] = visIndivTSSim(optStruct,m)

% m
% 1:Wi 2:In 3:Wa 4:Di 5:Cu 6: dumped 7: load 8: S
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
if m == 6
    set(gcf,'Position', [1, 1, 4.31, 0.95])
else
    set(gcf,'Position', [1, 1, 4.31, 0.65])
end
%set(gcf,'Position', [1, 1, 4.2, 3.5])
set(gcf,'Color','w')
fs3 = 10;
fs1 = 10;
col = colormap(brewermap(11,'Spectral')); %colors
col2 = colormap(brewermap(9,'Greys')); %colors
col3 = colormap(brewermap(9,'RdPu')); %colors
%set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultTextInterpreter', 'latex')
set(0,'DefaultAxesTickLabelInterpreter', 'latex')
set(0,'DefaultLegendInterpreter', 'latex')
%set(gcf, 'TickLabelInterpreter','latex')
%set(gcf, 'Position', [20, 20, 1300, 700])
tiledlayout(1,1,'TileSpacing','compact')

%Set length of time series
ts = 1:1:(8760/2);
%ts = 1:1:length(data.met.time);
t1 = datetime(data.met.time(1),'ConvertFrom','datenum');
t2 = datetime(data.met.time(ts(end)),'ConvertFrom','datenum');
%STORAGE TIME SERIES
%ax(1) = nexttile;
if m == 8
    ax(1) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.S1(1:ts(end))/(1000*output.min.Smax{end}),'Color',col(1,:), ... 
        'DisplayName','Battery 1','LineWidth',1)
    yticks([0 1])
    %[255,69,0]/256
    % hold on
    % plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
    %     output.min.S2(1:ts(end))/1000,'Color',col(3,:), ... 
    %     'DisplayName','Battery 2','LineWidth',1)
    % %[245,238,54]/256
    % l1 = legend('show','location','southwest','Orientation','horizontal');
%     ylabel({'$S(t)$', '[kWh]'},'interpreter','latex')
    % pos = get(l1,'Position');
    % set(l1,'Position',pos + [0.347,0.037,0,0]);
    drawnow
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.S2(1:end-1)/1000)])
    %yticks([0 round(max(output.min.S1(1:end-1)/1000),2)])
    xlim([t1 t2])
    set(ax(1),'XTick',[])
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on
    %add axis annotation
%     text(1.02, 0.5,'(a)','Units','Normalized', ...
%         'VerticalAlignment','middle','FontWeight','normal', ...
%         'FontSize',fs3,'Interpreter','latex')

%POWER TIME SERIES
elseif m==4
    ax(1) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.Pdies(ts)/(1000*output.min.kWd{end}),'Color',col2(4,:), ... 
        'DisplayName','Diesel Power','LineWidth',1)
    %hold on
%     ylabel({'$P_{\mbox{diesel}}$', '[kW]'},'interpreter','latex')
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.Pdies/1000)])
    xlim([t1 t2])
    yticks([0 1])
    set(ax(1),'XTick',[])
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on
elseif m == 5
    ax(1) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.Pcurr(ts)/(1000*output.min.kWc{end}),'Color',col(11,:), ... 
        'DisplayName','Current Power','LineWidth',1)

%     ylabel({'$P_{\mbox{current}}$', '[kW]'},'interpreter','latex')
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.Pcurr/1000)])
    yticks([0 1])
    xlim([t1 t2])
    set(ax(1),'XTick',[])
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on


%POWER TIME SERIES
elseif m == 2
    ax(1) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.Pinso(ts)/(1000*output.min.kWi{end}),'Color',col(5,:), ... 
        'DisplayName','Solar Power','LineWidth',1)

%     ylabel({'$P_{\mbox{solar}}$', '[kW]'},'interpreter','latex')
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.Pinso/1000)])
    yticks([0 1])
    xlim([t1 t2])
    set(ax(1),'XTick',[])
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on

%POWER TIME SERIES
elseif m == 3
    ax(1) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.Pwave(ts)/(1000*output.min.kWwa{end}),'Color',col(10,:), ... 
        'DisplayName','Wave Power','LineWidth',1)

%     ylabel({'$P_{\mbox{wave}}$', '[kW]'},'interpreter','latex')
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.Pwave/1000)])
    yticks([0 1])
    xlim([t1 t2])
    set(ax(1),'XTick',[])
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on


%POWER TIME SERIES
elseif m == 1
    ax(1) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.Pwind(ts)/(1000*output.min.kWwi{end}),'Color',col(9,:), ... 
        'DisplayName','Wind Power','LineWidth',1)

%     ylabel({'$P_{\mbox{wind}}$', '[kW]'},'interpreter','latex')
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.Pwind/1000)])
    yticks([0 1])
    xlim([t1 t2])
    set(ax(1),'XTick',[])
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on

%DUMPED POWER TIME SERIES
elseif m == 6
    ax(1) = nexttile;
    sum_rated = output.min.kWi{end} + output.min.kWwi{end} + output.min.kWd{end}...
        + output.min.kWwa{end} + output.min.kWc{end};
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.D(ts)/(1000*sum_rated),'Color',col3(5,:), ... 
        'DisplayName','Power Dumped','LineWidth',1)
    yticks([0 1])
    %[199,4,128]/256
    %legend('show')
%     ylabel({'$P_{\mbox{discard}}$', '[kW]'},'interpreter','latex')
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.D/1000)])
    xlim([t1 t2])
    %set(ax(1),'XTick',[])

    t1 = datetime(data.met.time(1),'ConvertFrom','datenum');
    t2 = datetime(data.met.time(ts(end)),'ConvertFrom','datenum');
    xlim([t1 t2])
    xlabel('Time','interpreter','latex')
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on

%Load Time Series
elseif m == 7
    ax(1) = nexttile;
    plot(datetime(data.met.time(ts),'ConvertFrom','datenum'), ...
        output.min.L(ts)/1000,'Color',col2(7,:), ... 
        'DisplayName','Load','LineWidth',1)
    yticks([0 1])
    %[212,92,176]/
    %legend('show')
%     ylabel({'$L(t)$', '[kW]'},'interpreter','latex')
%     ylh = get(ax(1),'ylabel');
%     set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%         'VerticalAlignment','middle', ...
%         'HorizontalAlignment','center')
    %ylim([0 1.1*max(output.min.L/1000)])
    %yticks([0 1])
    xlim([t1 t2])
    xticks([])
    set(gca,'FontSize',fs1)
    %set(gca,'LineWidth',2)
    grid on
    %add axis annotation
%     text(1.02, 0.5,'(g)','Units','Normalized', ...
%         'VerticalAlignment','middle','FontWeight','normal', ...
%         'FontSize',fs3,'Interpreter','latex')
end

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

% linkaxes(ax,'x')
%linkaxes(ax(2:6),'y')

%print(figure(1), '-dpng','-r600')
end

