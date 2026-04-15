%Time Domain Visualization of Hybrid results

%some diagnostic plots to figure out what's going on with weird hybrid
%optimization resutls

selectedfolder = uigetdir();

fileList = dir(fullfile(selectedfolder,'*.mat'));

load(fullfile(selectedfolder,'SimplifiedResults.mat'))

[mincost,ii] = min(cost);

tmp = load(fullfile(selectedfolder,fileList(ii).name));
nm = split(fileList(ii).name,'.');
nm = nm(1);
optStruct = tmp.(nm{1});

%check that this is the min cost file
cost = optStruct.output.min.cost;
if cost ~= mincost
    error('Error - this is the wrong file')
end

%% plot snow information
figure
tiledlayout(3,1)

nexttile
plot(optStruct.opt.ice_ts)
ylabel('No Ice Trigger')

nexttile
plot(optStruct.output.min.Pinso./(optStruct.output.min.kWi{1}*1000))
ylabel('Norm Solar Power')

nexttile
plot(optStruct.output.min.Pwind./(optStruct.output.min.kWwi{1}*1000))
ylabel('Norm Wind Power')
xlabel('Hour')

%% Diesel Checks
rundies = zeros(size(optStruct.output.min.Pdies));
fueluse = zeros(size(optStruct.output.min.Pdies));
lph = polyval(optStruct.opt.p_dev.d_burn,optStruct.output.min.kWd{1});
for i = 1:length(optStruct.output.min.Pdies)
    if optStruct.output.min.Pdies(i) > 0
        rundies(i) = 1;
        fueluse(i) = lph;
    end
end
rundies(1:8760*2) = cumsum(rundies(1:8760*2)); %tally over time for first maintenance int
rundies(8760*2+1:8760*4) = cumsum(rundies(8760*2+1:8760*4)); %tally over time for second maintenance int
rundies(8760*4+1:end) = cumsum(rundies(8760*4+1:end)); %tally over time for third maintenance int

fueluse(1:8760*2) = cumsum(fueluse(1:8760*2)); %tally over time for first maintenance int
fueluse(8760*2+1:8760*4) = cumsum(fueluse(8760*2+1:8760*4)); %tally over time for second maintenance int
fueluse(8760*4+1:end) = cumsum(fueluse(8760*4+1:end)); %tally over time for third maintenance int

figure
tiledlayout(3,1)

nexttile
plot(optStruct.output.min.Pdies/(optStruct.output.min.kWd{1}*1000))
ylabel('Norm Dies Power')

nexttile
hold on
plot(rundies)
yline(250,'r')
ylabel('Runtime accumulation')

nexttile
hold on
plot(fueluse)
yline(800,'r')
ylabel('Fuel accumulation')

%% Compare wind and wave and solar power for this location
data = optStruct.data;
opt = optStruct.opt;
uc = optStruct.uc;
wave = optStruct.wave;
atmo = optStruct.atmo;
inso = optStruct.inso;
cturb = optStruct.cturb;
turb = optStruct.turb;
econ = optStruct.econ;


[data, opt] = prepHybrid(data,opt,uc,wave,atmo,inso,cturb); %run prepHybrid to get data to calc power

% % %% Test small PA
% % %load wec sim results into structure
% % pa2 = load('C:\Users\smpal\Downloads\point_a_data_2m.mat');
% % pa2 = pa2.point_a_data;
% % %preallocate scatter arrays
% % H_scat = [];
% % T_scat = [];
% % CWR_scat = [];
% % 
% % n = size(pa2,1)*size(pa2,2);
% % for i = 1:size(pa2,1)
% %     for j = 1:size(pa2,2)
% %         if ~isempty(pa2{i,j})
% %             H = pa2{i,j}.wave_height.value;
% %             T = pa2{i,j}.peak_period.value;
% % 
% %             CWR = pa2{i,j}.capture_width_ratio.value; %find cwr
% %             %populate scatter arrays
% %             T_scat = [T_scat ; T];
% %             H_scat = [H_scat ; H];    
% %             CWR_scat = [CWR_scat ; CWR];
% %         end
% %     end
% % end
% % 
% % F =  scatteredInterpolant(T_scat,H_scat,CWR_scat);
% % Gr = (2*wave.eta_ct*F(wave.Tp_ra,wave.Hs_ra)*opt.wave.wavepower_ra)/((1)*(1000));

%% 
clear wind
kW_wind = 0.265;
kW_inso = 0.83;
kW_wave = 0.2144;

Pwind = zeros(size(swso));
Pinso = zeros(size(swso));
Pwave = zeros(size(swso));

swso = data.swso;
wind = data.met.wind_spd; %[m/s]
if atmo.dyn_h %use log law to adjust wind speed based on rotor height
    wind = adjustHeight(wind,data.met.wind_ht, ...
    turb.clearance + sqrt(1000*2*kW_wind/(turb.eta*atmo.rho_a_c*pi*turb.ura^3)) ...
    ,'log',(data.wind.z0.*1000)); 

end
Pwind_ci = (1/2)*atmo.rho_a_c*turb.uci^3; %cut in power flux [W]
Pwind_co = (1/2)*atmo.rho_a_c*turb.uco^3; %cut out power flux [W]
Pwind_ra = (1/2)*atmo.rho_a_c*turb.ura^3; %rated power flux [W]
%find power from wind turbine 
if kW_wind ~= 0 
    Pflux = (1/2).*data.wind.rho.*wind.^3;
    Pwind = kW_wind.*1000.*data.wind.rho.*wind.^3./(atmo.rho_a_c.*turb.ura.^3); %[W]
    ind_wiCI = find(Pflux < Pwind_ci);
    ind_wiCO = find(Pflux > Pwind_co);
    ind_wiRA = find(Pflux > Pwind_ra);

    Pwind(ind_wiRA) = kW_wind*1000; %[W] cap at rated power
    Pwind(ind_wiCI) = 0; %[W] no power below CI
    Pwind(ind_wiCO) = 0; %[W] no power above CO
else
    Pwind = zeros(1,length(time));
end


wavepower = opt.wave.wavepower_ts; %wavepower timeseries [kW]
if wave.method == 1 %divide by B methodology - OUTDATED    
    disp('ERROR- method set to old divide by B method')
elseif wave.method == 2 %3d interpolation methodology
    %extract data
    Hs = opt.wave.Hs; %Hs timeseries
    Tp = opt.wave.Tp; %Tp timeseries
    %find width through rated power conditions
    width = interp1(opt.wave.B_func(2,:),opt.wave.B_func(1,:),kW_wave); %[m], B
    if kW_wave == 0 %Set width = 0 for no wave gen (in case the curve cross (0,0))
        width = 0;
    end
    cw = width.*opt.wave.F(Tp,Hs,width*ones(length(Tp),1)); %[m] cw ts
end

%compute wave power timeseries
Pwave = wave.eta_ct*cw.*wavepower - kW_wave*wave.house; %[kW] 
Pwave_CI = wave.CIR*kW_wave; %cut in power
Pwave(Pwave<Pwave_CI) = 0; %no negative power, no power below cut in
Pwave(Pwave>kW_wave) = kW_wave; %no larger than rated power
Pwave = Pwave*1000; %convert to watts
if kW_wave < 0.2144
    Pwave = zeros(1,length(time));
    %disp('Zero wave power due to WECSIM min')
end
%set panel degradation
eff = (1-((inso.deg/100)/8760)*(1:1:length(swso)));
%rain = repmat(linspace(0.5,0,24*30),[1,ceil(length(swso)/(24*30))]);
if econ.inso.scen == 1 %automated cleaning
    d_soil_eff = 0;
    econ.inso.marinization = econ.inso.marinization*2;
elseif econ.inso.scen == 2 %human cleaning
    disp("ERROR DON'T RUN HUMAN CLEANING")
    %d_soil_eff = (atmo.soil/100)/8760; %change in soil deg per hour
end
soil_eff = 1; %starting soil efficiency
%find solar efficiency
eff_t = eff.*soil_eff.*inso.eff;
eff_t = eff_t'; %transpose for vector calc
%find power from panel
if kW_inso ~= 0
    Pinso = eff_t./inso.eff.*kW_inso.*1000.*(swso./(inso.rated.*1000)); %[W]
    ind_insoRa = find(Pinso >= (eff_t./inso.eff.*kW_inso.*1000)); %indices with power above rated
    Pinso(ind_insoRa) = eff_t(ind_insoRa)./inso.eff.*kW_inso.*1000;
else %no solar panel
    Pinso = zeros(1,length(time));
end

% average power
Pwimean = mean(Pwind);
Pimean = mean(Pinso);
Pwamean = mean(Pwave);
Pavg = [Pwimean,Pimean,Pwamean; kW_wind*1000,kW_inso*1000,kW_wave*1000];
X = categorical({'Wind','Solar','Wave'});
X = reordercats(X,{'Wind','Solar','Wave'});

figure
tiledlayout(2,1)
nexttile
hold on
plot(Pinso,'linewidth',1.5)
plot(Pwind,'linewidth',1.5)
plot(Pwave,'linewidth',1.5)
ylabel('Power [W]')
xlabel('Hours')
legend('Solar','Wind','Wave')
grid on
nexttile
h = bar(X,Pavg);
ylabel('Power [W]')
legend('Average Power','Rated Power')


%% Look at power time sieres
%Define color options for 2% space
colors = brewermap(11, 'Set3');
colpm(1,:) = colors(1,:);
colpm(2,:) = colors(6,:);
colpm(3,:) = colors(5,:);
colpm(4,:) = colors(9,:);
colpm(5,:) = colors(3,:);

% Define color options for optimal points
c1 = brewermap(9,'PuBuGn');
coptpm(1,:) = c1(end,:);

c2 = brewermap(9,'Oranges');
coptpm(2,:) = c2(end,:);

c3 = brewermap(9,'PuBu');
coptpm(3,:) = c3(end,:);

c4 = brewermap(9,'Greys');
coptpm(4,:) = c4(7,:);

c5 = brewermap(9,'BuPu');
coptpm(5,:) = c5(end,:);

ct = colors(8,:);

figure
tiledlayout(4,2)

nexttile(1)
plot(optStruct.output.min.Pinso,'linewidth',1.2,'color',colpm(2,:))
ylabel('Solar Power [W]')
grid on

nexttile(3)
plot(optStruct.output.min.Pwind,'linewidth',1.2,'color',colpm(1,:))
ylabel('Wind Power [W]')
grid on

nexttile(5)
plot(optStruct.output.min.Pdies,'linewidth',1.2,'color',colpm(4,:))
ylabel('Diesel Power [W]')
grid on

nexttile(7)
plot(optStruct.output.min.Ptot,'linewidth',1.2,'color',ct)
ylabel('Total Power [W]')
grid on

nexttile(2)
plot(optStruct.output.min.L,'linewidth',1.2,'color','k')
ylabel('Load [W]')
grid on

