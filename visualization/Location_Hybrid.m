%Arg Basin Cos Endurance IrmSea
clearvars -except allData
close all
set(0,'defaulttextinterpreter','none')
%set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultAxesFontName', 'cmr10')

if ~exist('allData','var')
    load('argBasin')
    load('cosEndurance_wa')
    load('irmSea')
    allData.argBasin = argBasin;
    allData.cosEndurance = cosEndurance_wa;
    allData.irmSea = irmSea;
end

%retrieve fieldnames, number of locations and
fn = fieldnames(allData); %fieldnames
l = numel(fn); %locations
allData.cosEndurance.title = 'Coastal Endurance'; %adjust CE title

%unpack data structure
lats = zeros(1,l);
lons = zeros(1,l);
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
maptitles{2} = {'CE'};
maptitles{3} = {'Irminger','Sea'};

%New Locations
lats(l+1) = 44.557; %N
lons(l+1) = -124.229; %W
maptitles{l+1} = {'PacWave'};
lats(l+2) = 21.474; %N
lons(l+2) = -157.754; %W
maptitles{l+2} = {'WETS'};
lats(l+3) = 35.545; %N
lons(l+3) = -74.817; %W
maptitles{l+3} = {'Mid-Atlantic','Shelf Break'};
lats(l+4) = 48.493; %N
lons(l+4) = -124.726; %W
maptitles{l+4} = {'PISCES'};
lats(l+5) = 34.129; %N
lons(l+5) = -119.223; %W LOC DIFFERENT THAN TABLE
maptitles{l+5} = {'Port','Hueneme'};
lats(l+6) = -58.783; %S
lons(l+6) = -64.563; %W
maptitles{l+6} = {'Antarctic','Circulation','Current'};
lats(l+7) = 57.511; %N
lons(l+7) = -163.789; %W
maptitles{l+7} = {'Bering', 'Sea'};

%plot settings
datacomp = figure;
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 6.5, 4])
set(gcf,'Color','w')
%col = colormap(brewermap(length(lats),'Pastel1')); %colors
col1 = [10 249 217]/255;
col2 = [255 0 247]/255;
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
%ax(1) = subplot(4,5,[1:2;6:7;11:12;16:17]);
%ax(1) = subplot(5,6,[1:2;7:8;13:14]);
ax(1) = subplot(1,1,1);
set(gca,'Units','Normalized')
wm = worldmap({'South America','Canada','Iceland','Pacific'});
%wm = worldmap('World');
hold on
geoshow(ax(1),'landareas.shp','FaceColor',[0.93 0.93 0.93]);
set(findall(wm,'Tag','PLabel'),'visible','off')
set(findall(wm,'Tag','MLabel'),'visible','off')
set(ax(1),'LineWidth',10)
% opos(1,:) = get(ax(1),'OuterPosition');
% opos(1,1) = 0.05;
% opos(1,3) = 0.29;
% opos(1,2) = 0.39;
% set(gca,'OuterPosition',opos(1,:))
framem off %remove frame
gridm off %remove grid
%add point locations and text
%AB, CE IS pacwave wets shelf-break pisces port antarctic bering
posmod = [1.10 -2 1.1 -4 1.85 1.1 -4.5 -0.7 1.1 1.65; ...
    1.1 0.95 0.9 0.85 1 0.9 1.05 .7 1.1 1]; %modify text position placement
for i = 1:length(lats)
    %add point locations
    if i <4
        col = col1;
    else
        col = col2;
    end
    pt = geoshow(lats(i),lons(i),'DisplayType','Point', 'Marker','o', ...
        'MarkerFaceColor',col,'MarkerEdgeColor','k', ... 
        'MarkerSize',ms,'LineWidth',mlw);
    %add text
    %tx = text(pt.XData, pt.YData, maptitles{i},'Interpreter','latex');
    tx = text(pt.XData*posmod(1,i),pt.YData*posmod(2,i), ...
        maptitles{i},'Interpreter','latex');
    tx.FontSize = fs2;
end