clearvars -except allData
close all
set(0,'defaulttextinterpreter','none')
%set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultAxesFontName', 'cmr10')

optInputs
if ~exist('allData','var')
    %load data into structure

    load('BerSea')
    allData.BerSea = data;
    [allData.BerSea.data, allData.BerSea.opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);

    load('MidAtlSB')
    allData.MidAtlSB = data;
    [allData.MidAtlSB.data, allData.MidAtlSB.opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);

    load('PacWave')
    allData.PacWave = data;
    [allData.PacWave.data, allData.PacWave.opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);

    load('altPISCES')
    allData.altPISCES = data;
    [allData.altPISCES.data, allData.altPISCES.opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);

    load('altWETS')
    allData.altWETS = data;
    [allData.altWETS.data, allData.altWETS.opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);
end

%% retrieve fieldnames, number of locations and
fn = fieldnames(allData); %fieldnames
l = numel(fn); %locations

%initialize power density
Kwave = cell(1,l);
Kwind = cell(1,1);
Ksolar = cell(1,1);
Kcurr = cell(1,l);

%get monthly power densities and time signals
for i = 1:l
    [monthlyArray,Kwave{i}] = getMonthlyK(allData.(fn{i}),'wave');
    [~,Kwind{i}] = getMonthlyK(allData.(fn{i}),'wind');
    [~,Ksolar{i}] = getMonthlyK(allData.(fn{i}), 'inso');
    [~,Kcurr{i}] = getMonthlyK(allData.(fn{i}), 'curr');
end

%% unpack data structure
lats = zeros(1,l);
lons = zeros(1,l);
maptitles = {};
labels = categorical;
dists = zeros(1,l);
depths = zeros(1,l);

for i = l:-1:1
    lats(i) = allData.(fn{i}).lat;
    lons(i) = allData.(fn{i}).lon;
    maptitles{i} = allData.(fn{i}).title;
    labels(i) = allData.(fn{i}).title;
    %dists(i) = allData.(fn{i}).dist;    
    depths(i) = allData.(fn{i}).depth;
end
% %adjustments for plotting on map
maptitles{1} = {'Bering Sea'};
maptitles{2} = {'Mid-Atlantic','Shelf Break'};
maptitles{3} = {'Oregon'};
maptitles{4} = {'Washington'};
maptitles{5} = {'O''ahu'};

%plot settings
colb = brewermap(9,'YlGnBu');
colm = brewermap(9,'BuPu');
coln = brewermap(9,'RdPu');
colw = brewermap(9,'YlGn');
colo = brewermap(9,'Oranges');

col = [colb(7,:); colm(7,:); coln(6,:); colw(7,:); colo(7,:)];
%col = colormap(brewermap(l,'Set2')); %colors
ms = 9; %marker size
mlw = .5; %marker line width
fs = 6; %font size
fs2 = 8; %map font size
fs3 = 8; %annotation font size
lw = 1.2; %line width
bw_adj = .7; %bar width adjustment
%ah_adj = .8; %axis height adjustment
ah_adj = .6; %axis height adjustment
xl_adj = 1.05; %xlabel height adjustment

%% Location Figure
datacomp = figure(1);
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 6.5, 5])
set(gcf,'Color','w')
t = tiledlayout(4,3);
t.TileSpacing = 'loose';
t.Padding = 'compact';
%MAP
ax(1) = nexttile(1,[4,2]);
set(gca,'Units','Normalized')
wm = worldmap({'South America','Canada','USA'});
hold on
geoshow(ax(1),'landareas.shp','FaceColor',[0.93 0.93 0.93]);
set(findall(wm,'Tag','PLabel'),'visible','off')
set(findall(wm,'Tag','MLabel'),'visible','off')
set(ax(1),'LineWidth',10)
framem off %remove frame
gridm off %remove grid
%add point locations and text
posmod = [1.5 1.1 2.6 2.9 1.35 ; ...
          1 1 1 1 1]; %modify text position placement
for i = 1:l
    %add point locations
    pt = geoshow(lats(i),lons(i),'DisplayType','Point', 'Marker','d', ...
        'MarkerFaceColor',col(i,:),'MarkerEdgeColor','k', ... 
        'MarkerSize',ms,'LineWidth',mlw);
    %add text
    tx = text(pt.XData*posmod(1,i),pt.YData*posmod(2,i), ...
        maptitles{i},'Interpreter','latex');
    tx.FontSize = fs2;
end
text(0.05,0.95,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')

% MONTHLY K - WAVE
ax(1) = nexttile(3);
for i = 1:l
    plot(monthlyArray, Kwave{i},'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title])
    hold on
end
%xl = xlabel('Time');
xtickformat('yyyy');
ylabel({'Wave','Power','Density'},'FontSize',fs,'interpreter','latex');
%title('Wave','FontSize',fs,'interpreter','latex')
ylh = get(ax(1),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
text(-.175,.35,'$$\mathrm{\bigg[\frac{kW}{m}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
%hL = legend('show','location','eastoutside','box','off','Units','Inches');
drawnow
ax(1).TickLabelInterpreter = 'latex';
set(gca,'FontSize',fs)
%add axis annotation
text(1.05, 0.5,'(b)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - WIND
ax(2) = nexttile(6);
for i = 1:l
    plot(monthlyArray, ...
        Kwind{i}/1000,'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title])
    hold on
end
%xl = xlabel('Time');
xtickformat('yyyy');
ylabel({'Wind','Power','Density'},'FontSize',fs,'interpreter','latex');
%title('Wind','FontSize',fs,'interpreter','latex')
ylh = get(ax(2),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
text(-.175,.35,'$$\mathrm{\bigg[\frac{kW}{m^2}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
set(gca,'FontSize',fs)
ax(2).TickLabelInterpreter = 'latex';

%add axis annotation
text(1.05,.5,'(c)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - SOLAR
ax(3) = nexttile(9);
for i = 1:l
    plot(monthlyArray, ...
        Ksolar{i}/1000,'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title])
    hold on
end
%xl = xlabel('Time');
xtickformat('yyyy');
ylabel({'Solar','Power','Density'},'FontSize',fs,'interpreter','latex');
%title('Solar','FontSize',fs,'interpreter','latex')
ylh = get(ax(3),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
text(-.175,.35,'$$\mathrm{\bigg[\frac{kW}{m^2}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
set(gca,'FontSize',fs)
ax(3).TickLabelInterpreter = 'latex';

%add axis annotation
text(1.05,.5,'(d)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - CURRENT
ax(4) = nexttile(12);
for i = 1:l
    plot(monthlyArray, ...
        Kcurr{i}/1000,'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title]) %PLOTTING ONLY AT FIRST DEPTH (SHOULD IS BE THE 2ND DEPTH)
    hold on
end
xl = xlabel('Time','Interpreter','latex');
xtickformat('yyyy');
ylabel({'Current','Power','Density'},'FontSize',fs,'interpreter','latex');
%title('Solar','FontSize',fs,'interpreter','latex')
ylh = get(ax(4),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
text(-.175,.35,'$$\mathrm{\bigg[\frac{kW}{m^2}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
set(gca,'FontSize',fs)
ax(4).TickLabelInterpreter = 'latex';

%add axis annotation
text(1.05,.5,'(e)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

print(datacomp,['G:\My Drive\Hybrid Ocean Observation\OO-Hybrid Paper\DRAFT_ResourceHybrid'],  ...
    '-dpng','-r600')
