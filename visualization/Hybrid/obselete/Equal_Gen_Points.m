optStruct = load('C:\Users\smpal\Documents\GitHub\EWTEC_fig\argBasin_5_5_2_p2t_v16.mat');
output = optStruct.argBasin_5_5_2_p2t_v16.output;

for i =1:5
    surv{i} = output.surv{1,i};
    cost_a{i} = output.cost{1,i};
    cost_a{i}(surv{i}<0.99) = inf;
    cost_a{i}(surv{i}>0.995) = inf;
    min_cost{i} = min(cost_a{i});
    Imin{i} = find(cost_a{i} == min(cost_a{i}));
    %cost_a{i}(cost_a{i}<1.2*min_cost{i}) = inf;

    Ki{i} = output.Ki_run{1,i};
    Kwi{i} = output.Kwi_run{1,i};
    Kd{i} = output.Kd_run{1,i};
    Kwa{i} = output.Kwa_run{1,i};
    Kc{i} = output.Kc_run{1,i};
    Smax{i} = output.S_run{1,i};
    tot_G{i} = Ki{i} + Kwi{i} + Kd{i} + Kwa{i} + Kc{i};

    %remove points that don't survive
    Ki{i}(cost_a{i}==inf) = [];
    Kwi{i}(cost_a{i}==inf) = [];
    Kd{i}(cost_a{i}==inf) = [];
    Kc{i}(cost_a{i}==inf) = [];
    Kwa{i}(cost_a{i}==inf) = [];
    Smax{i}(cost_a{i}==inf) = [];
    tot_G{i}(cost_a{i}==inf) = [];
    cost_a{i}(cost_a{i}==inf) = [];

    [sort_cost{i},IS{i}] = sort(cost_a{i});
end


fs = 8;
%col = colormap(brewermap(8,'Pastel1')); %colors
set(gcf,'Units','inches')
set(gcf,'Position', [1, 1, 8.0, 5.75])
set(gcf,'Color','w')
figure(1)
tiledlayout(7,2,'TileSpacing','tight')

p_max4 = max(cost_a{4})
a_min4 = min(cost_a{4});
p_max_mult4 = p_max4/a_min4;
p_max5 = max(cost_a{5})
a_min5 = min(cost_a{5});
p_max_mult5 = p_max5/a_min5;
lb = 0;
AdvancedColormap('kkgw glww lww r rk k',length(IS{4}), ...
    [lb,lb+.025*(1-lb),lb+0.1*(1-lb),.35,.7,1]);
%caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
x4 = 1:1:length(IS{4});
x5 = 1:1:length(IS{5});

max1 = max([max(Ki{4}(IS{4})),max(Ki{5}(IS{5}))]);
ax1 = nexttile;
scatter(x4,Ki{4}(IS{4}),2,cost_a{4}(IS{4}))
caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$G_{r,solar}$','[kW]'},'FontSize',fs,'interpreter','latex');
title('Iteration 4','FontSize',fs,'interpreter','latex');
ylh = get(ax1,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax1.TickLabelInterpreter = 'latex';
ylim([0 1])

ax8 = nexttile;
scatter(x5,Ki{5}(IS{5}),2,cost_a{5}(IS{5}))
caxis([a_min5 ceil(max(p_max_mult5)*a_min5)])
%ylabel('Solar')
title('Iteration 5','FontSize',fs,'interpreter','latex');
ax8.TickLabelInterpreter = 'latex';
ylim([0 1])

max2 = max([max(Kwi{4}(IS{4})), max(Kwi{5}(IS{5}))]);
ax2 = nexttile;
scatter(x4,Kwi{4}(IS{4}),2,cost_a{4}(IS{4}))
caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$G_{r,wind}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax2,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax2.TickLabelInterpreter = 'latex';
ylim([0 max2])

ax9 = nexttile;
scatter(x5,Kwi{5}(IS{5}),2,cost_a{5}(IS{5}))
caxis([a_min5 ceil(max(p_max_mult5)*a_min5)])
ax9.TickLabelInterpreter = 'latex';
ylim([0 max2])
%ylabel('Wind')

max3 = max([max(Kwa{4}(IS{4})), max(Kwa{5}(IS{5}))]);
ax3 = nexttile;
scatter(x4,Kwa{4}(IS{4}),2,cost_a{4}(IS{4}))
caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$G_{r,wave}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax3,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax3.TickLabelInterpreter = 'latex';
ylim([0 max3])

ax10 = nexttile;
scatter(x5,Kwa{5}(IS{5}),2,cost_a{5}(IS{5}))
caxis([a_min5 ceil(max(p_max_mult5)*a_min5)])
ax10.TickLabelInterpreter = 'latex';
ylim([0 max3])

max4 = max([max(Kd{4}(IS{4})),max(Kd{5}(IS{5}))]);
ax4 = nexttile;
scatter(x4,Kd{4}(IS{4}),2,cost_a{4}(IS{4}))
caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$G_{r,diesel}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax4,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax4.TickLabelInterpreter = 'latex';
ylim([0 0.3])

ax11 = nexttile;
scatter(x5,Kd{5}(IS{5}),2,cost_a{5}(IS{5}))
caxis([a_min5 ceil(max(p_max_mult5)*a_min5)])
ax11.TickLabelInterpreter = 'latex';
ylim([0 0.3])
%ylabel('Diesel')
% nexttile
% for i=4:5
%     plot(Kc{i},symbol{i})
%     ylabel('Current')
% end

max5 = max([max(tot_G{4}(IS{4})),max(tot_G{5}(IS{5}))]);
ax5 = nexttile;
scatter(x4,tot_G{4}(IS{4}),2,cost_a{4}(IS{4}))
caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$G_{r,total}$','[kW]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax5,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax5.TickLabelInterpreter = 'latex';
ylim([0 max5])

ax12 = nexttile;
scatter(x5,tot_G{5}(IS{5}),2,cost_a{5}(IS{5}))
caxis([a_min5 ceil(max(p_max_mult5)*a_min5)])
ax12.TickLabelInterpreter = 'latex';
ylim([0 max5])

ax6 = nexttile;
max6 = max([max(Smax{4}(IS{4})),max(Smax{5}(IS{5}))]);
scatter(x4,Smax{4}(IS{4}),2,cost_a{4}(IS{4}))
caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$S_{max}$','[kWh]'},'FontSize',fs,'interpreter','latex');
ylh = get(ax6,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax6.TickLabelInterpreter = 'latex';
ylim([0 max6])

ax13 = nexttile;
scatter(x5,Smax{5}(IS{5}),2,cost_a{5}(IS{5}))
caxis([a_min5 ceil(max(p_max_mult5)*a_min5)])
ax13.TickLabelInterpreter = 'latex';
ylim([0 max6])

ax7 = nexttile;
scatter(x4,cost_a{4}(IS{4}),2,cost_a{4}(IS{4}))
caxis([a_min4 ceil(max(p_max_mult4)*a_min4)])
ylabel({'$M_{total}$','[kg]'},'FontSize',fs,'interpreter','latex');
xlabel('Index in Grid','interpreter','latex')
ylh = get(ax7,'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
ax7.TickLabelInterpreter = 'latex';
ylim([min(cost_a{5}(IS{5})) max(cost_a{4}(IS{4}))])
c4 = colorbar;
%c4.Label.String = '[kg]';
c4.Location = 'south'
c4.TickLabelInterpreter = 'latex'
pos4 = get(c4,'Position');
set(c4,'Position',pos4 + [0,-0.12,0,-0.007]);

ax14 = nexttile;
scatter(x5,cost_a{5}(IS{5}),2,cost_a{5}(IS{5}))
caxis([a_min5 ceil(max(p_max_mult5)*a_min5)])
xlabel('Index in Grid','interpreter','latex')
ax14.TickLabelInterpreter = 'latex';
ylim([min(cost_a{5}(IS{5})) max(cost_a{4}(IS{4}))])
c5 = colorbar;
%c5.Label.String = '[kg]';
c5.Location = 'south'
c5.TickLabelInterpreter = 'latex'
pos5 = get(c5,'Position');
set(c5,'Position',pos5 + [0,-0.12,0,-0.007]);

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6 ax7],'x')
  linkaxes([ax8 ax9 ax10 ax11 ax12 ax13 ax14],'x')  
