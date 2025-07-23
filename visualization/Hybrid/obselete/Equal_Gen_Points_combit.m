% optStruct = load('C:\Users\smpal\Documents\GitHub\EWTEC_fig\Res_62\argBasin_5_5_2_p2t_v16.mat');
% output = optStruct.argBasin_5_5_2_p2t_v16.output;
% opt = optStruct.argBasin_5_5_2_p2t_v16.opt;

optStruct = load('C:\Users\smpal\Documents\GitHub\EWTEC_fig\Res_62\argBasin_4_5_2_p2t_v16.mat');
output = optStruct.argBasin_4_5_2_p2t_v16.output;
opt = optStruct.argBasin_4_5_2_p2t_v16.opt;

surv = [output.surv{1,1}];
cost_a = [output.cost{1,1}];
Ki = [output.Ki_run{1,1}];
Kwi = [output.Kwi_run{1,1}];
Kd = [output.Kd_run{1,1}];
Kwa = [output.Kwa_run{1,1}];
Kc = [output.Kc_run{1,1}];
Smax = [output.S_run{1,1}];
total_G = [];

for i =2:5
    surv = [surv; output.surv{1,i}];
    cost_a = [cost_a; output.cost{1,i}];
    %cost_a{i}(cost_a{i}<1.2*min_cost{i}) = inf;

    Ki = [Ki; output.Ki_run{1,i}'];
    Kwi = [Kwi; output.Kwi_run{1,i}'];
    Kd = [Kd; output.Kd_run{1,i}'];
    Kwa = [Kwa; output.Kwa_run{1,i}'];
    Kc = [Kc; output.Kc_run{1,i}'];
    Smax = [Smax; output.S_run{1,i}'];

end
tot_G = Ki + Kwi + Kd + Kwa + Kc;
cost_a(surv<0.99) = inf;
cost_a(surv>1) = inf;
min_cost = min(cost_a);
Imin = find(cost_a == min(cost_a));
%cost_a(cost_a>1.2*min_cost) = inf;

%remove points that don't survive
Ki(cost_a==inf) = [];
Kwi(cost_a==inf) = [];
Kd(cost_a==inf) = [];
Kc(cost_a==inf) = [];
Kwa(cost_a==inf) = [];
Smax(cost_a==inf) = [];
tot_G(cost_a==inf) = [];
surv(cost_a == inf) = [];
cost_a(cost_a==inf) = [];
[sort_cost,IS] = sort(cost_a);
%[sort_surv,IS] = sort(surv);
%Make Plot Arrays

fs = 8;
fs3 = 8;
%col = colormap(brewermap(8,'Pastel1')); %colors
set(gcf,'Units','inches')
set(gcf,'Position', [1, 0.0, 4.25, 4.75])
col = colormap(brewermap(11,'Spectral')); %colors
col2 = colormap(brewermap(9,'Greys')); %colors
col3 = colormap(brewermap(9,'RdPu')); %colors
set(gcf,'Color','w')
figure(1)
tiledlayout(7,1,'TileSpacing','tight')
%tiledlayout(8,1,'TileSpacing','tight')
% p_max4 = max(cost_a)
% a_min4 = min(cost_a);
% p_max_mult4 = p_max4/a_min4;
% p_max5 = max(cost_a{5})
% a_min5 = min(cost_a{5});
% p_max_mult5 = p_max5/a_min5;
% lb = 0;
% IS_full = IS;
% AdvancedColormap('kkgw glww lww r rk k',length(IS_full), ...
%     [lb,lb+.025*(1-lb),lb+0.1*(1-lb),.35,.7,1]);
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
IS = IS(1:6000);
x4 = 1:1:length(IS);

%max1 = max([max(Ki{4}(IS{4})),max(Ki{5}(IS{5}))]);
ax7 = nexttile([1,1]);
scatter(x4,cost_a(IS),5,col2(7,:),'filled')
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$J_{\mbox{total}}$'},'FontSize',fs,'interpreter','latex');
ylh = get(ax7,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax7.TickLabelInterpreter = 'latex';
text(1.02, 0.5,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
set(gca,'FontSize',6)
grid on
box on

ax1 = nexttile;
scatter(x4,Ki(IS),5,col(5,:),'filled')
hold on
%find zeros
x4_zi = x4(Ki(IS) == 0);
cost_zi = sort_cost(x4_zi); 
Ki_z = zeros(length(x4_zi),1);
scatter(x4_zi,Ki_z,20,col2(9,:),'x')
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$G_{r\mbox{,solar}}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax1,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax1.TickLabelInterpreter = 'latex';
text(1.02, 0.5,'(b)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
set(gca,'FontSize',6)
%ylim([0 1])
grid on
box on

%max2 = max([max(Kwi{4}(IS{4})), max(Kwi{5}(IS{5}))]);
ax2 = nexttile;
scatter(x4,Kwi(IS),5,col(9,:),'filled')
hold on
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
%find zeros
x4_zwi = x4(Kwi(IS) == 0); 
Kwi_z = zeros(length(x4_zwi),1);
scatter(x4_zwi,Kwi_z,20,col2(9,:),'x')

ylabel({'$G_{r\mbox{,wind}}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax2,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax2.TickLabelInterpreter = 'latex';
text(1.02, 0.5,'(c)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
set(gca,'FontSize',6)
%ylim([0 max2])
grid on
box on

%max3 = max([max(Kwa{4}(IS{4})), max(Kwa{5}(IS{5}))]);
ax3 = nexttile;
scatter(x4,Kwa(IS),5,col(10,:),'filled')
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
hold on
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
%find zeros
x4_zwa = x4(Kwa(IS) == 0); 
Kwa_z = zeros(length(x4_zwa),1);
scatter(x4_zwa,Kwa_z,20,col2(9,:),'x')
ylabel({'$G_{r\mbox{,wave}}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax3,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax3.TickLabelInterpreter = 'latex';
text(1.02, 0.5,'(d)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
set(gca,'FontSize',6)
%ylim([0 max3])
grid on
box on

%max4 = max([max(Kd(IS)),max(Kd(IS))]);
if opt.pm == 5
    ax4 = nexttile;
    scatter(x4,Kd(IS),5,col2(4,:),'filled')
    hold on
    %caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
    %find zeros
    x4_zd = x4(Kd(IS) == 0); 
    Kd_z = zeros(length(x4_zd),1);
    scatter(x4_zd,Kd_z,20,col2(9,:),'x')
    %caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
    ylabel({'$G_{r\mbox{,diesel}}$','[kW]'},'FontSize',fs,'interpreter','latex');
    ylh = get(ax4,'ylabel');
    set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
        'VerticalAlignment','middle', ...
        'HorizontalAlignment','center')
    ax4.TickLabelInterpreter = 'latex';
    text(1.02, 0.5,'(e)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
    set(gca,'FontSize',6)
    %ylim([0 0.3])
    grid on
    box on

elseif opt.pm == 4
    ax4 = nexttile;
    scatter(x4,Kc(IS),5,col(11,:),'filled')
    hold on
    %caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
    %find zeros
    x4_zc = x4(Kc(IS) == 0); 
    Kc_z = zeros(length(x4_zc),1);
    scatter(x4_zc,Kc_z,20,col2(9,:),'x')
    %caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
    ylabel({'$G_{r\mbox{,current}}$','[kW]'},'FontSize',fs,'interpreter','latex');
    ylh = get(ax4,'ylabel');
    set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
        'VerticalAlignment','middle', ...
        'HorizontalAlignment','center')
    ax4.TickLabelInterpreter = 'latex';
    text(1.02, 0.5,'(e)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
    set(gca,'FontSize',6)
    %ylim([0 0.3])
    grid on
    box on
end

%max5 = max([max(tot_G(IS)),max(tot_G{5}(IS{5}))]);
ax5 = nexttile;
scatter(x4,tot_G(IS),5,col2(7,:),'filled')
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$G_{r\mbox{,total}}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax5,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax5.TickLabelInterpreter = 'latex';
text(1.02, 0.5,'(f)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
set(gca,'FontSize',6)
%ylim([0 max5])
grid on
box on

ax6 = nexttile;
%max6 = max([max(Smax{4}(IS{4})),max(Smax{5}(IS{5}))]);
scatter(x4,Smax(IS),5,col(1,:),'filled')
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$S_{m}$','[kWh]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax6,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax6.TickLabelInterpreter = 'latex';
text(1.02, 0.5,'(g)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
set(gca,'FontSize',6)
xlabel('Case Index','FontSize',fs,'interpreter','latex')
%ylim([0 max6])
grid on
box on

%ylim([min(cost_a{5}(IS{5})) max(cost_a{4}(IS{4}))])
% c4 = colorbar;
% set(c4.XLabel,{'String','Rotation','Position','Interpreter'},{'[kg]',0,[0.5 -0.01],'latex'})
% c4.Layout.Tile = 'south';
% c4.TickLabelInterpreter = 'latex';
% c4.FontSize = 6;
% pos4 = get(c4,'Position');
% set(c4,'Position',pos4 + [0,-0.12,0,-0.009]);
% ax8 = nexttile;
% scatter(x4,surv(IS),5,cost_a(IS),'filled')
% caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
% ylabel({'$a_{sim}$'},'FontSize',fs,'interpreter','latex');
% xlabel('Index in Grid','interpreter','latex')
% ylh = get(ax8,'ylabel');
% set(ylh,'Rotation',0,'Units','Normalized','Position',[-.1 .75 -1], ...
%     'VerticalAlignment','middle', ...
%     'HorizontalAlignment','center')
% ax8.TickLabelInterpreter = 'latex';
% %ylim([min(cost_a{5}(IS{5})) max(cost_a{4}(IS{4}))])
% grid on
% box on

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7],'x')

%print(datacomp, '-dpng','-r600')
print(figure(1), '-dpng','-r600')