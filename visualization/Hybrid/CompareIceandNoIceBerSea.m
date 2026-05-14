%Plot Bar Chart comparison of Bering Sea results with and without ice
%modeling

%No Ice data
noIce2D = 'C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\NoIce_BerSea\2DBF';
fileList = dir(fullfile(noIce2D,'*.mat'));
for i = 1:length(fileList)
    tmp = load(fullfile(noIce2D,fileList(i).name));
    nm = split(fileList(i).name,'.');
    nm = nm(1);
    optStruct(i) = tmp.(nm{1});
    costmin(i) = optStruct(i).output.min.cost;
end
[~,ind2D] = min(costmin);
optStruct2DNI = optStruct(ind2D);

noIce6D = 'C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\NoIce_BerSea\';
optStruct6DNI = load(fullfile(noIce6D,'SimplifiedResults.mat'));

%With Ice data
clear optStruct costmin fileList tmp nm loc
Ice2D = 'C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\AllLUP';
fileList = dir(fullfile(Ice2D,'*.mat'));
for i = 1:length(fileList)
    tmp = load(fullfile(Ice2D,fileList(i).name));
    nm = split(fileList(i).name,'.');
    nm = nm(1);
    loc = tmp.(nm{1}).loc;
    lc = tmp.(nm{1}).uc.loadcase;
    if strcmp(loc,'BerSea')
        if lc == 1
            optStruct(i) = tmp.(nm{1});
            costmin(i) = optStruct(i).output.min.cost;
        else
            costmin(i) = Inf;
        end
    else
        costmin(i) = Inf;
    end
end
[~,ind2DI] = min(costmin);
optStruct2DI = optStruct(ind2DI);

Ice6D = 'C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\FFA6DSensitivity_SW0\ffa6D_BerSeaLC1_SW0';
optStruct6DI = load(fullfile(Ice6D,'SimplifiedResults.mat'));

%% orgaize data for plotting
% X = categorical({'$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Single}\\\mathrm{Gen}\end{array}$',...
%     '$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Hybrid}\end{array}$',...
%     '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{Gen}\\\mathrm{Wind}\end{array}$',...
%     '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{Gen}\\\mathrm{Wave}\end{array}$',...
%     '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Hybrid}\end{array}$'});
% X = reordercats(X,{'$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Single}\\\mathrm{Gen}\end{array}$',...
%     '$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Hybrid}\end{array}$',...
%     '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{Gen}\\\mathrm{Wind}\end{array}$',...
%     '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{Gen}\\\mathrm{Wave}\end{array}$',...
%     '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Hybrid}\end{array}$'});


X = categorical({'$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Single}\\\mathrm{Gen}\end{array}$',...
    '$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Hybrid}\end{array}$',...
    '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{GenWi}\end{array}$',...
    '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{GenWa}\end{array}$',...
    '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Hybrid}\end{array}$'});
X = reordercats(X,{'$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Single}\\\mathrm{Gen}\end{array}$',...
    '$\begin{array}{c}\mathrm{No\:Ice}\\\mathrm{Hybrid}\end{array}$',...
    '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{GenWi}\end{array}$',...
    '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Single}\\\mathrm{GenWa}\end{array}$',...
    '$\begin{array}{c}\mathrm{Ice}\\\mathrm{Hybrid}\end{array}$'});
[~,minNI] = min(optStruct6DNI.cost);
[~,minI] = min(optStruct6DI.cost);

gen2DNI = [optStruct2DNI.output.min.kWwi{1}, optStruct2DNI.output.min.kWi{1}, optStruct2DNI.output.min.kWwa{1}, optStruct2DNI.output.min.kWd{1}, optStruct2DNI.output.min.kWc{1}];
gen2DI = [optStruct2DI.output.min.kWwi{1}, optStruct2DI.output.min.kWi{1}, optStruct2DI.output.min.kWwa{1}, optStruct2DI.output.min.kWd{1}, optStruct2DI.output.min.kWc{1}];
wavegen = [0,0,0.2244,0,0];

gen = [gen2DNI; optStruct6DNI.gen(minNI,:); gen2DI; wavegen; optStruct6DI.gen(minI,:)];
smax = [optStruct2DNI.output.min.Smax{1}, optStruct6DNI.smax(minNI), optStruct2DI.output.min.Smax{1},13, optStruct6DI.smax(minI)];
cost = [optStruct2DNI.output.min.cost, optStruct6DNI.cost(minNI), optStruct2DI.output.min.cost,2.75E5, optStruct6DI.cost(minI)];

%% plotting
colors = brewermap(11, 'Set3');
colpm(1,:) = colors(1,:);
colpm(2,:) = colors(6,:);
colpm(3,:) = colors(5,:);
colpm(4,:) = colors(9,:);
colpm(5,:) = colors(3,:);

colors = brewermap(11, 'Set1');
colc(1,:) = colors(9,:);
cols(1,:) = colors(8,:);

bw = 0.5;
fs = 10;

figure
set(gcf,'Units','Inches','Position', [0.05,0.05,4.5,4])
tf = tiledlayout(3,3);
tf.Padding = 'compact';
tf.TileSpacing = 'compact';

ax(1) = nexttile(2,[1 2]);
c = bar(X,cost./1000,bw);

ylh = ylabel({'Cost', '[\$1000]'},'Interpreter','latex');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.23 .5 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center','FontSize',fs)
set(c, 'FaceColor', 'Flat')
c(1).CData = colc(1,:);

grid on
box on
set(gca, 'XTickLabel', []);
ax(1).TickLabelInterpreter = 'latex';
ax(1).FontSize = fs;
ylim([150,300])

ax(2) = nexttile(5,[1 2]);
hold on
g = bar(X,gen,bw,'stacked');


g(1).CData = colpm(1,:);
g(2).CData = colpm(2,:);
g(3).CData = colpm(3,:);
g(4).CData = colpm(4,:);
g(5).CData = colpm(5,:);

ylh = ylabel({'Rated','Power',' [kW]'},'Interpreter','latex');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.23 .5 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center','FontSize',fs)
set(g, 'FaceColor', 'Flat')
grid on
box on
set(gca, 'XTickLabel', []);

ax(2).TickLabelInterpreter = 'latex';
ax(2).FontSize = fs;
lh = legend('Wind','Solar','Wave','','','Interpreter','Latex','FontSize',fs,'NumColumns',3,'Location','NorthEast');
lh.ItemTokenSize(1) = 10;
ylim([0,0.8])

ax(3) = nexttile(8,[1 2]);
s = bar(X,smax,bw);
xtickangle(0)
ylh = ylabel({'Battery','Max','Capacity','[kWh]'},'Interpreter','latex','Rotation',0,'HorizontalAlignment', 'center');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.23 .5 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center','FontSize',fs)
set(s, 'FaceColor', 'Flat')
s(1).CData = cols(1,:);
grid on
box on
ax(3).TickLabelInterpreter = 'latex';
ax(3).FontSize = fs;
ax(3).XAxis.TickLabelChild.VerticalAlignment = 'middle'
ylim([0,30])