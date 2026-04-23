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

Ice6D = 'C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\ffa6D_BerSeaLC1_SW0';
optStruct6DI = load(fullfile(Ice6D,'SimplifiedResults.mat'));

%% orgaize data for plotting
'$\begin{array}{c}Ice\\Single Generator\end{array}$';
X = categorical({'$\begin{array}{c}No Ice\\Single Generator\end{array}$','$\begin{array}{c}No Ice\\Hybrid\end{array}$','$\begin{array}{c}Ice\\Single Generator\end{array}$','$\begin{array}{c}Ice\\Hybrid\end{array}$'});
X = reordercats(X,{'$\begin{array}{c}No Ice\\Single Generator\end{array}$','$\begin{array}{c}No Ice\\Hybrid\end{array}$','$\begin{array}{c}Ice\\Single Generator\end{array}$','$\begin{array}{c}Ice\\Hybrid\end{array}$'});

[~,minNI] = min(optStruct6DNI.cost);
[~,minI] = min(optStruct6DI.cost);

gen2DNI = [optStruct2DNI.output.min.kWwi{1}, optStruct2DNI.output.min.kWi{1}, optStruct2DNI.output.min.kWwa{1}, optStruct2DNI.output.min.kWd{1}, optStruct2DNI.output.min.kWc{1}];
gen2DI = [optStruct2DI.output.min.kWwi{1}, optStruct2DI.output.min.kWi{1}, optStruct2DI.output.min.kWwa{1}, optStruct2DI.output.min.kWd{1}, optStruct2DI.output.min.kWc{1}];

gen = [gen2DNI; optStruct6DNI.gen(minNI,:); gen2DI; optStruct6DI.gen(minI,:)];
smax = [optStruct2DNI.output.min.Smax{1}, optStruct6DNI.smax(minNI), optStruct2DI.output.min.Smax{1},optStruct6DI.smax(minI)];
cost = [optStruct2DNI.output.min.cost, optStruct6DNI.cost(minNI), optStruct2DI.output.min.cost,optStruct6DI.cost(minI)];

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

figure
set(gcf,'Units','Inches','Position', [0.05,0.05,4,6])
tf = tiledlayout(3,1);
tf.Padding = 'loose';
tf.TileSpacing = 'compact';

ax(1) = nexttile;
c = bar(X,cost./1000);
ylabel({'Cost', '[\$1000]'},'Interpreter','latex','Rotation',0,'HorizontalAlignment', 'center')
set(c, 'FaceColor', 'Flat')
c.CData = colc(1,:);
grid on
box on
set(gca,'xtick',[])
ax(1).TickLabelInterpreter = 'latex';

ax(2) = nexttile;
g = bar(X,gen,'stacked');
ylabel({'Generation',' [kW]'},'Interpreter','latex','Rotation',0,'HorizontalAlignment', 'center')
set(g, 'FaceColor', 'Flat')
grid on
box on
set(gca,'xtick',[])
g(1).CData = colpm(1,:);
g(2).CData = colpm(2,:);
g(3).CData = colpm(3,:);
g(4).CData = colpm(4,:);
g(5).CData = colpm(5,:);
ax(2).TickLabelInterpreter = 'latex';
legend('Wind','Solar','Wave','Diesel','Current','Interpreter','Latex')

ax(3) = nexttile;
s = bar(X,smax);
ylabel({'Battery Max',' Capacity [kWh]'},'Interpreter','latex','Rotation',0,'HorizontalAlignment', 'center')
set(s, 'FaceColor', 'Flat')
s.CData = cols(1,:);
grid on
box on
ax(3).TickLabelInterpreter = 'latex';