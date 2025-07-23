clearvars -except allData
close all
set(0,'defaulttextinterpreter','none')
%set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultAxesFontName', 'cmr10')

load('argBasin')
allData.argBasin = argBasin;

%retrieve fieldnames, number of locations and
fn = fieldnames(allData); %fieldnames
l = numel(fn); %locations

%initialize power densit
Kwave = cell(1,l);
kwind = cell(1,1);
Ksolar = cell(1,1);
Kwave_ts = cell(1,l);

%get monthly power densities
for i = 1:l
    Kwave{i} = getMonthlyK(allData.(fn{i}),'wave');
    Kwind{i} = getMonthlyK(allData.(fn{i}),'wind');
    Ksolar{i} = getMonthlyK(allData.(fn{i}), 'inso');
    Kwave_ts{i} = getPowerDensity(allData.(fn{i}),'wave');
    Kwind_ts{i} = getPowerDensity(allData.(fn{i}),'wind');
    Ksolar_ts{i} = getPowerDensity(allData.(fn{i}), 'inso');
end

Kmean_wave = zeros(l,3);
Kmean_wind = zeros(1,3);
Kmean_solar = zeros(1,3);
for i = l:-1:1
    Kmean_wave(i,1) = mean(Kwave_ts{i}(:,2))/1000;
    Kmean_wave(i,2) = Kmean_wave(i,1) - prctile(Kwave_ts{i}(:,2),25)/1000;
    Kmean_wave(i,3) = prctile(Kwave_ts{i}(:,2),75)/1000 - Kmean_wave(i,1);
    Kmean_wind(i,1) = mean(Kwind_ts{i}(:,2))/1000;
    Kmean_wind(i,2) = Kmean_wind(i,1) - prctile(Kwind_ts{i}(:,2),25)/1000;
    Kmean_wind(i,3) = prctile(Kwind_ts{i}(:,2),75)/1000 - Kmean_wind(i,1);
    Kmean_solar(i,1) = mean(Ksolar_ts{i}(:,2),'omitnan')/1000;
    Kmean_solar(i,2) = Kmean_solar(i,1) - prctile(Ksolar_ts{i}(:,2),25)/1000;
    Kmean_solar(i,3) = prctile(Ksolar_ts{i}(:,2),75)/1000 - Kmean_solar(i,1);
end
xt = [];
col = colormap(brewermap(5,'Pastel1')); %colors
subplot(3,1,1)
for i = 1:l
    plot(datetime(Kwave{i}(:,1),'ConvertFrom','datenum'), ...
        Kwave{i}(:,2)/1000,'Color',col(2,:), ...
        'LineWidth',2,'DisplayName',[allData.(fn{i}).title])
    xt = [xt ; datetime(Kwave{i}(:,1),'ConvertFrom','datenum')];
    str = 'Wave $$ W/m $$';
    ylabel(str,'Interpreter','latex','FontWeight','bold')
end
subplot(3,1,2)
xt = [];
for i = 1:l
    plot(datetime(Kwind{i}(:,1),'ConvertFrom','datenum'), ...
        Kwind{i}(:,2)/1000,'Color',col(3,:), ...
        'LineWidth',2,'DisplayName',[allData.(fn{i}).title])
    xt = [xt ; datetime(Kwind{i}(:,1),'ConvertFrom','datenum')];
    str = 'Wind $$ W/m^2 $$';
    ylabel(str,'Interpreter','latex','FontWeight','bold')
end
subplot(3,1,3)
xt = [];
for i = 1:l
    plot(datetime(Ksolar{i}(:,1),'ConvertFrom','datenum'), ...
        Ksolar{i}(:,2)/1000,'Color',col(5,:), ...
        'LineWidth',2,'DisplayName',[allData.(fn{i}).title])
    xt = [xt ; datetime(Ksolar{i}(:,1),'ConvertFrom','datenum')];
    str = 'Solar $$ W/m^2 $$';
    ylabel(str,'Interpreter','latex','FontWeight','bold')
end