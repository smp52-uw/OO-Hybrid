%Arg Basin Cos Endurance IrmSea
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
    [data, opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);
    allData.BerSea.wave = opt.wave;

    load('MidAtlSB')
    allData.MidAtlSB = data;
    load('PacWave')
    allData.PacWave = data;
    load('PISCES')
    allData.PISCES = data;
    load('WETS')
    allData.WETS = data;
end

%retrieve fieldnames, number of locations and
fn = fieldnames(allData); %fieldnames
l = numel(fn); %locations

%initialize power density
Kwave = cell(1,l);
Kwind = cell(1,1);
Ksolar = cell(1,1);
Kcurr = cell(1,l);
Kwave_ts = cell(1,l);
Kwind_ts = cell(1,l);
Ksolar_ts = cell(1,l);
Kcurr_ts = cell(1,l);

%get monthly power densities and time signals
for i = 1:l
    Kwave{i} = getMonthlyK(allData.(fn{i}),'wave');
    Kwind{i} = getMonthlyK(allData.(fn{i}),'wind');
    Ksolar{i} = getMonthlyK(allData.(fn{i}), 'inso');
    Kcurr{i} = getMonthlyK(allData.(fn{i}), 'curr');

    Kwave_ts{i} = getPowerDensity(allData.(fn{i}),'wave');
    Kwind_ts{i} = getPowerDensity(allData.(fn{i}),'wind');
    Ksolar_ts{i} = getPowerDensity(allData.(fn{i}), 'inso');
    Kcurr_ts{i} = getPowerDensity(allData.(fn{i}), 'curr');
end

%unpack data structure
lats = zeros(1,l);
lons = zeros(1,l);
maptitles = {};
labels = categorical;
dists = zeros(1,l);
depths = zeros(1,l);
Kmean_wave = zeros(l,3);
Kmean_wind = zeros(1,3);
Kmean_solar = zeros(1,3);
for i = l:-1:1
    lats(i) = allData.(fn{i}).lat;
    lons(i) = allData.(fn{i}).lon;
    maptitles{i} = allData.(fn{i}).title;
    labels(i) = allData.(fn{i}).title;
    %dists(i) = allData.(fn{i}).dist;    
    depths(i) = allData.(fn{i}).depth;
    Kmean_wave(i,1) = mean(Kwave_ts{i}(:,2))/1000;
    Kmean_wave(i,2) = Kmean_wave(i,1) - prctile(Kwave_ts{i}(:,2),25)/1000;
    Kmean_wave(i,3) = prctile(Kwave_ts{i}(:,2),75)/1000 - Kmean_wave(i,1);
    Kmean_wind(i,1) = mean(Kwind_ts{i}(:,2))/1000;
    Kmean_wind(i,2) = Kmean_wind(i,1) - prctile(Kwind_ts{i}(:,2),25)/1000;
    Kmean_wind(i,3) = prctile(Kwind_ts{i}(:,2),75)/1000 - Kmean_wind(i,1);
    Kmean_solar(i,1) = mean(Ksolar_ts{i}(:,2),'omitnan')/1000;
    Kmean_solar(i,2) = Kmean_solar(i,1) - prctile(Ksolar_ts{i}(:,2),25)/1000;
    Kmean_solar(i,3) = prctile(Ksolar_ts{i}(:,2),75)/1000 - Kmean_solar(i,1);

    Kmean_curr(i,1) = mean(Kcurr_ts{i}(:,2),'omitnan')/1000; %THIS IS FOR THE FIRST DEPTH!!!
    Kmean_curr(i,2) = Kmean_curr(i,1) - prctile(Kcurr_ts{i}(:,2),25)/1000;
    Kmean_curr(i,3) = prctile(Kcurr_ts{i}(:,2),75)/1000 - Kmean_curr(i,1);
end
% %adjustments for plotting on map
maptitles{1} = {'Bering Sea'};
maptitles{2} = {'Mid-Atlantic','Shelf Break'};
maptitles{3} = {'PacWave'};
maptitles{4} = {'PISCES'};
maptitles{5} = {'WETS'};

%plot settings
col = colormap(brewermap(l,'Set2')); %colors
ms = 5; %marker size
mlw = .5; %marker line width
fs = 6; %font size
fs2 = 8; %map font size
fs3 = 8; %annotation font size
lw = 1.5; %line width
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
posmod = [1.5 1.1 2.6 2.4 1.4 ; ...
          1 1 1 1 1]; %modify text position placement
for i = 1:l
    %add point locations
    pt = geoshow(lats(i),lons(i),'DisplayType','Point', 'Marker','o', ...
        'MarkerFaceColor',col(i,:),'MarkerEdgeColor','k', ... 
        'MarkerSize',ms,'LineWidth',mlw);
    %add text
    tx = text(pt.XData*posmod(1,i),pt.YData*posmod(2,i), ...
        maptitles{i},'Interpreter','latex');
    tx.FontSize = fs2;
end


% MONTHLY K - WAVE
ax(1) = nexttile(3);
xt = [];
for i = 1:l
    plot(datetime(Kwave{i}(:,1),'ConvertFrom','datenum'), ...
        Kwave{i}(:,2)/1000,'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title])
    xt = [xt ; datetime(Kwave{i}(:,1),'ConvertFrom','datenum')];
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
text(-.175,.25,'$$\mathrm{\bigg[\frac{kW}{m}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
%hL = legend('show','location','eastoutside','box','off','Units','Inches');
drawnow
ax(1).TickLabelInterpreter = 'latex';
set(gca,'FontSize',fs)
%add axis annotation
text(1.05, 0.5,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - WIND
ax(2) = nexttile(6);
xt = [];
for i = 1:l
    plot(datetime(Kwind{i}(:,1),'ConvertFrom','datenum'), ...
        Kwind{i}(:,2)/1000,'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title])
    xt = [xt ; datetime(Kwind{i}(:,1),'ConvertFrom','datenum')];
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
text(-.175,.25,'$$\mathrm{\bigg[\frac{kW}{m^2}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
set(gca,'FontSize',fs)
ax(2).TickLabelInterpreter = 'latex';

%add axis annotation
text(1.05,.5,'(b)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - SOLAR
ax(3) = nexttile(9);
xt = [];
for i = 1:l
    plot(datetime(Ksolar{i}(:,1),'ConvertFrom','datenum'), ...
        Ksolar{i}(:,2)/1000,'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title])
    xt = [xt ; datetime(Ksolar{i}(:,1),'ConvertFrom','datenum')];
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
text(-.175,.25,'$$\mathrm{\bigg[\frac{kW}{m^2}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
set(gca,'FontSize',fs)
ax(3).TickLabelInterpreter = 'latex';

%add axis annotation
text(1.05,.5,'(c)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - CURRENT
ax(4) = nexttile(12);
xt = [];
for i = 1:l
    plot(datetime(Kcurr{i}(:,1),'ConvertFrom','datenum'), ...
        Kcurr{i}(:,2)/1000,'Color',col(i,:), ...
        'LineWidth',lw,'DisplayName',[allData.(fn{i}).title]) %PLOTTING ONLY AT FIRST DEPTH (SHOULD IS BE THE 2ND DEPTH)
    xt = [xt ; datetime(Kcurr{i}(:,1),'ConvertFrom','datenum')];
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
text(-.175,.25,'$$\mathrm{\bigg[\frac{kW}{m^2}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
set(gca,'FontSize',fs)
ax(4).TickLabelInterpreter = 'latex';

%add axis annotation
text(1.05,.5,'(d)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

print(datacomp,['C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\NAVFACReport\DRAFT_ResourceHybrid'],  ...
    '-dpng','-r600')
