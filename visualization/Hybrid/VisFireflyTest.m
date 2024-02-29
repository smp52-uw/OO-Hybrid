%Visualize Firefly Testing Data
clear
clc
%load csv of output data
data = readmatrix('FFA_TEST_916.csv');
numpop = data(:,3);
numit = data(:,4);
time = data(:,5);
cost = data(:,6);
itconv = data(:,7);
pos = data(:,8:13);
Kd = data(:,8);
Ki = data(:,9);
Kwi = data(:,10);
Kwa = data(:,11);
Kc = data(:,12);
Smax = data(:,13);

%calculate cost variation
testID = categorical({'50x10', '75x10', '100x10', '25x10', '25x12', '25x15', '25x20'});
testID = reordercats(testID,{'25x10', '25x12', '25x15', '25x20','50x10', '75x10', '100x10'});
cost_bin{1} = cost(1:5);
cost_bin{2} = cost(6:8);
cost_bin{3} = cost(9:11);
cost_bin{4} = cost(12:14);
cost_bin{5} = cost(15:17);
cost_bin{6} = cost(18:20);
cost_bin{7} = cost(20:23);

pos_bin{1} = pos(1:5);
pos_bin{2} = pos(6:8);
pos_bin{3} = pos(9:11);
pos_bin{4} = pos(12:14);
pos_bin{5} = pos(15:17);
pos_bin{6} = pos(18:20);
pos_bin{7} = pos(20:23);

time_bin{1} = time(1:5);
time_bin{2} = time(6:8);
time_bin{3} = time(9:11);
time_bin{4} = time(12:14);
time_bin{5} = time(15:17);
time_bin{6} = time(18:20);
time_bin{7} = time(20:23);

for i=1:7
    cost_std(i) = std(cost_bin{i});
    for j = 2:length(pos_bin{i})
        cost_diff{i}(j) = 100*abs(cost_bin{i}(j)-cost_bin{i}(1))/cost_bin{i}(1);
        pos_norm{i}(j) = norm((pos_bin{i}(1)-pos_bin{i}(j)));
    end
    cost_diff_avg(i) = mean(cost_diff{i});
    pos_norm_avg(i) = mean(pos_norm{i});
    time_avg(i) = mean((time_bin{i}./60));
    mean_cost(i) = mean(cost_bin{i});
end
cost_std(3) = 1;


%plot test ID vs. position variation
datacomp = figure;
fs = 8;
col = colormap(brewermap(9,'RdYlBu')); %colors
set(gcf,'Units','inches')
set(gcf,'Position', [0, 0, 18, 7.0])
set(gcf,'Color','w')

%t = tiledlayout(1,3,'TileSpacing','loose','Padding','tight')

%t.XLabel.String = 'Population Size - Number of Iterations';
%t.XLabel.Interpreter = 'Latex';
%t.XLabel.FontSize = 20;

ax1 = subplot(1,3,1);
set(gca,'Units','Inches')
set(gca, 'LineWidth', 1.5)
b = bar(testID,cost_diff_avg,'FaceColor',col(9,:));
yl = ylabel(ax1, {'$\overline{J_{diff}}$', '[\%]'},'Interpreter','Latex','rotation',0);
yl.Position(1) = -1.9;

%xlabel('Population Size - Number of Iterations','Interpreter','Latex')
ax1.TickLabelInterpreter = 'latex';
ax1.FontSize = 20;
b.FaceColor = 'flat';
b.CData(1,:) = col(3,:);
b.CData(2,:) = col(3,:);
b.CData(3,:) = col(3,:);
b.CData(4,:) = col(3,:);
pos1 = get(gca,'Position');
pos1(1) = 1.25;
pos1(2) = 0.85;
pos1(3) = 3.5;
set(gca,'Position',pos1)
set(gca, 'LineWidth', 1.5)
ax1.XAxis.TickLength = [0,0];

ax2 = subplot(1,3,2);
set(gca,'Units','Inches')
b = bar(testID,pos_norm_avg,'FaceColor',col(9,:));
yl = ylabel(ax2, {'$\overline{X_{diff}}$', '[-]'},'Interpreter','Latex','rotation',0);
yl.Position(1) = -2.1;
%xlabel('Population Size - Number of Iterations','Interpreter','Latex')
ax2.TickLabelInterpreter = 'latex';
ax2.FontSize = 20;
b.FaceColor = 'flat';
b.CData(1,:) = col(3,:);
b.CData(2,:) = col(3,:);
b.CData(3,:) = col(3,:);
b.CData(4,:) = col(3,:);
pos2 = get(gca,'Position');
pos2(1) = 6.25;
pos2(2) = 0.85;
pos2(3) = 3.5;
set(gca,'Position',pos2)
set(gca, 'LineWidth', 1.5)
ax2.XAxis.TickLength = [0,0];

ax3 = subplot(1,3,3);
set(gca,'Units','Inches')
b = bar(testID,time_avg,'FaceColor',col(9,:));

ylim([0 320])
yl = ylabel({'$\overline{T}$', '[min]'},'Interpreter','Latex','rotation',0);
yl.Position(1) = -1.9;
%xlabel('Population Size - Number of Iterations','Interpreter','Latex')
ax3.TickLabelInterpreter = 'latex';
ax3.FontSize = 20;
b.FaceColor = 'flat';
b.CData(1,:) = col(3,:);
% b.DisplayName = 'Local'
b.CData(2,:) = col(3,:);
b.CData(3,:) = col(3,:);
b.CData(4,:) = col(3,:);
pos3 = get(gca,'Position');
pos3(1) = 11.25;
pos3(2) = 0.85;
pos3(3) = 3.5;
set(gca,'Position',pos3)
set(gca, 'LineWidth', 1.5)
ax3.XAxis.TickLength = [0,0];
% b.DisplayName(5) = 'HPC'
%legend show

% datacomp2 = figure;
% col = colormap(brewermap(9,'RdYlBu')); %colors
% set(gcf,'Units','inches')
% set(gcf,'Position', [0, 0, 6, 7.0])
% set(gcf,'Color','w')
% 
% ax4 = gca;
% set(gca,'Units','Inches')
% b = bar(testID,mean_cost,'FaceColor',col(9,:));
% yl4 = ylabel(ax4, {"$\overline{J_{total}}$",'[-]'},'Interpreter','Latex','rotation',0);
% %yl4.Position(1)
% %yl4.Position(1) = 1
% ylim([600 640])
% ax4.TickLabelInterpreter = 'latex';
% ax4.FontSize = 20;
% b.FaceColor = 'flat';
% b.CData(1,:) = col(3,:);
% b.CData(2,:) = col(3,:);
% b.CData(3,:) = col(3,:);
% b.CData(4,:) = col(3,:);
% pos4 = get(gca,'Position');
% pos4(1) = 1.5;
% pos4(3) = 3.5;
% set(gca,'Position',pos4)
% set(gca, 'LineWidth', 1.5)
% ax4.XAxis.TickLength = [0,0];

% leg = legend('Orientation', 'Horizontal');
% leg.Layout.Tile = 'north';
print(datacomp,'C:\Users\smpal\Desktop\UW\Research\Hybrid\FFA_Test',  ...
    '-dpng','-r600')