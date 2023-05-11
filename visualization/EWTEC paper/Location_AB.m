%Arg Basin Cos Endurance IrmSea
clearvars -except allData
close all
set(0,'defaulttextinterpreter','none')
%set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultAxesFontName', 'cmr10')

if ~exist('allData','var')
    %load data into structure
    load('argBasin')

    allData.argBasin = argBasin;
end
allData.argBasin.curr = load('HYCOM_AB_mod_2018.mat');
%retrieve fieldnames, number of locations and
fn = fieldnames(allData); %fieldnames
l = 1; %locations

%initialize power densit
Kwave = cell(1,1);
kwind = cell(1,1);
Ksolar = cell(1,1);
Kwave_ts = cell(1,1);
Kcurr = cell(1,1);

%get monthly power densities

Kwave{1} = getMonthlyK(allData.(fn{1}),'wave');
Kwind{1} = getMonthlyK(allData.(fn{1}),'wind');
Ksolar{1} = getMonthlyK(allData.(fn{1}), 'inso');
Kcurr{1} = getMonthlyK(allData.(fn{1}), 'curr');
% Kwave_ts{1} = getPowerDensity(allData.(fn{1}),'wave');
% Kwind_ts{1} = getPowerDensity(allData.(fn{1}),'wind');
% Ksolar_ts{1} = getPowerDensity(allData.(fn{1}), 'inso');
% Kcurr{1} = get

%unpack data structure
lats = zeros(1,1);
lons = zeros(1,1);
maptitles = {};
labels = categorical;
for i = l:-1:1
    lats(i) = allData.(fn{i}).lat;
    lons(i) = allData.(fn{i}).lon;
    maptitles{i} = allData.(fn{i}).title;
    labels(i) = allData.(fn{i}).title;
end
%adjustments for plotting on map
maptitles{1} = {'Argentine','Basin'};

%plot settings
datacomp = figure;
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 5.5, 3.25])
set(gcf,'Color','w')
col = colormap(brewermap(8,'Pastel1')); %colors

%col = [col10(2,:); col10(4,:); col10(6,:); col10(8,:); col10(10,:);];
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


%MAP
t = tiledlayout(4,2)
t.TileSpacing = 'loose';
t.Padding = 'compact'
ax(1) = nexttile(1,[4,1]);
ax(1).Layout.Tile = 1;
% ax(1).Layout.TileSpan = [2 4];
%set(gca,'Units','Normalized')
wm = worldmap({'South America'});
hold on
geoshow(ax(1),'landareas.shp','FaceColor',[0.93 0.93 0.93]);
set(findall(wm,'Tag','PLabel'),'visible','off')
set(findall(wm,'Tag','MLabel'),'visible','off')
set(ax(1),'LineWidth',10)

framem off %remove frame
gridm off %remove grid
%add point locations and text
pt = geoshow(lats(1),lons(1),'DisplayType','Point', 'Marker','o', ...
    'MarkerFaceColor',col(1,:),'MarkerEdgeColor','k', ... 
    'MarkerSize',ms,'LineWidth',mlw);
%add text
tx = text(pt.XData*1.20,pt.YData*1.1, ...
    maptitles{1},'Interpreter','latex');
tx.FontSize = fs2;

ab_mod_lat = -41.25;
ab_mod_lon = -52.03;
pt_m = geoshow(ab_mod_lat,ab_mod_lon,'DisplayType','Point', 'Marker','o', ...
    'MarkerFaceColor',col(8,:),'MarkerEdgeColor','k', ... 
    'MarkerSize',ms,'LineWidth',mlw);
%add text
tx_m = text(pt_m.XData*1.15,pt_m.YData*0.87, ...
    {'Modified','Argentine Basin'},'Interpreter','latex');
tx_m.FontSize = fs2;

% % MONTHLY K - WAVE
ax(2) = nexttile(2,[1 1]);
xt = [];
plot(datetime(Kwave{1}(:,1),'ConvertFrom','datenum'), ...
    Kwave{1}(:,2)/1000,'Color',col(2,:), ...
    'LineWidth',lw)
xt = [xt ; datetime(Kwave{1}(:,1),'ConvertFrom','datenum')];
hold on
%xl = xlabel('Time');
xtickformat('yyyy');
ylabel({'Wave','Power','Density'},'FontSize',fs,'interpreter','latex');
ylh = get(ax(2),'ylabel');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.175 .75 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center')
text(-.175,.25,'$$\mathrm{\bigg[\frac{kW}{m}\bigg]}$$', ...
    'Units','Normalized','Interpreter','latex', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'HorizontalAlignment','center','FontSize',fs);
ylim([0 inf])
%drawnow

set(gca,'FontSize',fs)
ax(2).TickLabelInterpreter = 'latex';

%add axis annotation
text(1.02, 0.5,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - WIND
ax(3) = nexttile(4,[1 1]);
xt = [];
plot(datetime(Kwind{1}(:,1),'ConvertFrom','datenum'), ...
    Kwind{1}(:,2)/1000,'Color',col(3,:), ...
    'LineWidth',lw)
xt = [xt ; datetime(Kwind{1}(:,1),'ConvertFrom','datenum')];
hold on
%xl = xlabel('Time');
xtickformat('yyyy');
ylabel({'Wind','Power','Density'},'FontSize',fs,'interpreter','latex');
%title('Wind','FontSize',fs,'interpreter','latex')
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
text(1.02,.5,'(b)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on

% %MONTHLY K - SOLAR
ax(4) = nexttile(6,[1 1]);
xt = [];
plot(datetime(Ksolar{1}(:,1),'ConvertFrom','datenum'), ...
    Ksolar{1}(:,2)/1000,'Color',col(5,:), ...
    'LineWidth',lw)
xt = [xt ; datetime(Ksolar{1}(:,1),'ConvertFrom','datenum')];
hold on
%xl = xlabel('Time');
xtickformat('yyyy');
ylabel({'Solar','Power','Density'},'FontSize',fs,'interpreter','latex');
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
text(1.02,.5,'(c)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')
grid on
% % 

% %MONTHLY K - CURRENT
ax(4) = nexttile(8,[1 1]);
xt = [];
K_sz = size(Kcurr{1,1});
depth = [3]
col = colormap(brewermap(10,'Purples')); %colors
leg_txt = [];
for i = depth
    plot(datetime(Kcurr{1}(:,1),'ConvertFrom','datenum'), ...
        Kcurr{1}(:,i+1)/1000, 'Color',col(i+2,:),...
        'LineWidth',lw,'DisplayName',[strcat(num2str(allData.(fn{1}).curr.depth(i)),' m depth')])
    %leg_txt = [leg_txt, strcat(num2str(allData.(fn{1}).curr.depth(i)),'m')]
    hold on
end
% hL = legend('show','Location','northeast','box','off','Units','Inches')
% set(hL,'Interpreter','latex')
% hLPos = get(hL,'Position');
% set(hL,'Position',[hLPos(1)+0.15 hLPos(2)+0.1 hLPos(3) hLPos(4)+0.025])
xl = xlabel('Time');
set(xl,'Interpreter','latex')
xtickformat('yyyy');
st = datetime(Kcurr{1}(1,1),'ConvertFrom','datenum');
ed = datetime(Kcurr{1}(end,1),'ConvertFrom','datenum');
set(gca, 'XTick', [st:calmonths(2):ed]);
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
text(1.02,.5,'(d)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs3,'Interpreter','latex')

grid on
% % 


print(datacomp,'C:\Users\smpal\Documents\GitHub\EWTEC_fig\Loc_fig',  ...
    '-dpng','-r600')