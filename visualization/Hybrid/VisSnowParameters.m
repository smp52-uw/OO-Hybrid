%Plot snow depth and temperature at BeringSea
close all
load('Site_6_Bering_Sea.mat')

%read in nc data for snow depth
% addpath(genpath("C:\Users\smpal\Downloads\cf3b5c5775e19dee5888e3531ea3d026"));
% file_name = "data_stream-oper_stepType-instant.nc";
addpath(genpath("C:\Users\smpal\Downloads"))
% file_name = '54856e723fe3243e893bd3796519aa8.nc'; %Land Snow Depth
file_name = 'f72172fc6ba5fb2f8923f4dc861906d5.nc'; %sea ice cover

%temp
file_name = '45c464209f132aa5578693103ac791ee.nc'; %temp

tmp = ncread(file_name,'valid_time'); %%THIS IS IN SECONDS SINCE 1970
date_ref = datenum(1970,1,1,0,0,0); %%WHY IS THIS REFERENCE DIFFERENT IN THIS FILE THAN THE DATA BRIAN DOWNLOADED
tmp = (double(tmp)/24/3600)+date_ref;
snow.time = datetime(tmp,'ConvertFrom','datenum','TimeZone','UTC');

snow.lat = ncread(file_name,'latitude');
snow.lon = ncread(file_name,'longitude');

temp = squeeze(ncread(file_name,'t2m'));
temp = squeeze(temp(2,2,:));
temp = temp - 273.15; %convert to C

figure
plot(snow.time, temp)

%ice = squeeze(ncread(file_name,'siconc'));
%ice = squeeze(ice(2,2,:));
%tmpsd = squeeze(ncread(file_name,'sd'));   %snow depth (SWE) 2,2
%snow.sd = squeeze(tmpsd(2,2,:));

%monthly average, and monthly hours with non-zero snow
for m = 1:12
    indm = find(month(solar.time) == m);
    indsd = find(month(snow.time) == m);
    meanSR(m) = mean(solar.mean_snowfall_rate(indm));
    meanS(m) = mean(solar.snowfall(indm));
    meansd(m) = mean(snow.sd(indsd));

    nonzerocount(m) = length(find(solar.snowfall(indm)>0));
    nonzeroPerc(m) = nonzerocount(m)./length(indm);
end

% figure %histograms 
% tiledlayout(3,4)
% for m = 1:12
%     clear edges N meanEdge
%     indm = find(month(solar.time) == m);
%     edges = linspace(0,max(solar.snowfall(indm)),15);
%     edges(1) = 0;
%     edges(2) = 1E-15;
%     if edges(3) < edges(2)
%         edges = edges(1:2);
%     end
%     [N,edges] = histcounts(solar.snowfall(indm),edges,'Normalization','probability');
%     for i = 1:length(edges)-1
%         meanEdge(i) = (edges(i) + edges(i+1))/2;
%     end
%     indLim = find(N < 0.01,1,'first');
%     nexttile
%     bar(meanEdge,N)
%     title(strcat("Month: ",string(m)))
%     ylim([0,1])
%     %xlim([0,meanEdge(indLim)])
% end

% figure %plot bar chart with hours of nonzero snowfall per month
% X = categorical({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% X = reordercats(X,{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
% h = bar(X, nonzeroPerc, 'FaceColor','flat');
% ylabel("Percent hours with nonzero snowfall")
% 

figure
tiledlayout(3,2)

ax(1) = nexttile;
plot(solar.time,solar.mean_snowfall_rate,'linewidth',1.5)
xlabel('Time')
ylabel('Mean Snowfall Rate [mm/s]') %the unit is kgm-2s-1 but they say it's equivalent to mm/2 of liquid water

nexttile
X = categorical({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
X = reordercats(X,{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
h = bar(X, meanSR, 'FaceColor','flat');
ylabel([{'Monthly Avg Mean'},{'Snowfall Rate [mm/s]'}])

ax(2) = nexttile;
plot(solar.time,solar.snowfall,'linewidth',1.5,'DisplayName','SWE')
legend
xlabel('Time')
ylabel('Snowfall [m]') %unit is m of equivalent water - the depth if it was spread evenly over the grid box

nexttile
X = categorical({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
X = reordercats(X,{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
h = bar(X, meanS, 'FaceColor','flat');
ylabel('Monthly Avg Snowfall [m]')

ax(3) = nexttile;
plot(snow.time,snow.sd,'linewidth',1.5,'DisplayName','SWE')
legend
xlabel('Time')
ylabel('Snowdepth [m]') %unit is m of equivalent water - the depth if it was spread evenly over the grid box

nexttile
X = categorical({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
X = reordercats(X,{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
h = bar(X, meansd, 'FaceColor','flat');
ylabel('Monthly Avg Snowdepth [m]')



linkaxes(ax,'x')




