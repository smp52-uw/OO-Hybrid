%Visualize the economic model 

%needs optInputs, libraries, and use case
optInputs
opt.p_dev.t = calcDeviceVal('turbine',[],econ.wind_n);
opt.p_dev.d_cost = calcDeviceVal('dieselcost',[],econ.diescost_n);
opt.p_dev.d_mass = calcDeviceVal('dieselmass',[],econ.diesmass_n);
opt.p_dev.d_size = calcDeviceVal('dieselsize',[],econ.diessize_n);
opt.p_dev.d_burn = calcDeviceVal('dieselburn',[],econ.diesburn_n);
opt.p_dev.d_vol = calcDeviceVal('dieselvol',[],econ.diesvol_n);
opt.p_dev.b_size = calcDeviceVal('lfp_vol',[],econ.battsize_n);
[opt.p_dev.b,~,opt.p_dev.kWhmax] = calcDeviceVal('agm',[],econ.batt_n);

%Generator and battery set up
kW_inso = linspace(0,8,500);
kW_wind = linspace(0,8,500);
kW_dies = linspace(0,8,500);
kW_wave = linspace(0,8,500);
kW_curr = linspace(0,8,500);
Smax = linspace(1,500,500);

%Set up
nvi = 3;
runtime = [250,250,250];
lph = polyval(opt.p_dev.d_burn,kW_dies); %[l/h], burn rate
lph(1) = 0; %when generator is 0 kW
runtime_tot = sum(runtime);
newbatt = 0;
econ.wave.costmult = econ.wave.costmult_con;

%% Economic Model copied 4/7/2026
%% Calculate Mass of components (for weight approx optimization or for mooring model)
mass_solar = inso.wf.*kW_inso./(inso.rated*inso.eff); %[kg] - 1 array
mass_solar_E = 0; %electrical - for 1 buoy
%mass_solar_S = 3.15*kW_inso/(inso.rated*inso.eff); %structural - for 1 buoy
mass_solar_S = 0;
mass_wind = kW_wind.*turb.wf; %updated to include a tower

mass_curr = 0; %assume the turbine assembly is ballasted to neutral buoyancy (Brian)

mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of 1 generator [kg]
al_plate = 35.24; %[kg/m2] - 1.2" 6061 AL
gen_vol = econ.dies.volmult.*(polyval(opt.p_dev.d_size,kW_dies).^3);
mass_diesencl = 6.*(gen_vol.^(2/3)).*al_plate; %[kg]
dies_dens = 0.85; %[g/cm3] = [kg/L] from a chevron report
mass_fuel = max(runtime).*lph.*dies_dens; %[kg]

%generator 0 kW
mass_diesencl(1) = 0;
mass_fuel(1) = 0;
mass_dies(1) = 0;

mass_batt = Smax./(batt.se*batt.V/1000); %[kg] - 1 battery
batt_vol = polyval(opt.p_dev.b_size,Smax).*econ.batt.volmult;
batt_len = (batt_vol).^(1/3); %assume cube
mass_battencl = 6.*(batt_len.^2).*al_plate; %[kg]

% Pmtrl = econ.platform.wf* ...
%     (mass_dies + mass_diesencl + mass_fuel + mass_solar + mass_wind + mass_batt + ...
%     mass_battencl + mass_solar_E + mass_solar_S); % 1 platform material [kg]

mass_wec = kW_wave.*895.78./econ.wave.mass_mult; %based on RM3 weight with no float
mass_wec(1) = 0; %0kW wec

% %add instrumentation mass
% comp_plat_mass = comp_plat_mass + uc.loaddata(uc.loadcase).mass;

%% Economic Optimization

%% Costs
%costs - solar
Mcost_inso = 2*econ.inso.module.*kW_inso.*econ.inso.marinization*econ.inso.pcm; %module
Icost_inso = 2*econ.inso.installation.*kW_inso.*econ.inso.pcm; %installation
Ecost_inso = 2*econ.inso.electrical.*kW_inso.*econ.inso.marinization ...
    *econ.inso.pcm; %electrical infrastructure
Strcost_inso = 2*econ.inso.structural.*kW_inso.*econ.inso.marinization ...
    *econ.inso.pcm; %structural infrastructure
%costs - wind
kWcost_wind = 2.*polyval(opt.p_dev.t,kW_wind).*econ.wind.marinization ...
    *econ.wind.tcm; %turbine
kWcost_wind(1) = 0;
Icost_wind(1) = 0;
for i = 2:length(kW_wind)
    Icost_wind(i) = 2*(econ.wind.installed - 0.5*kWcost_wind(i)/ ...
        (kW_wind(i)*econ.wind.marinization*econ.wind.tcm))*kW_wind(i); %installation
end
%costs - wave
kWcost_wave = 2*econ.wave.costmult.*polyval(opt.p_dev.t,kW_wave); %wec
kWcost_wave(1) = 0;
Icost_wave(1) = 0;
for i = 2:length(kW_wave)
    Icost_wave(i) = 2*(econ.wind.installed - 0.5*kWcost_wave(i)/(kW_wave(i)*econ.wave.costmult))*kW_wave(i)*econ.wave.costmult; %installation 
end

%costs - current
kWcost_curr = 2*econ.curr.costmult.*polyval(opt.p_dev.t,kW_curr); %current
kWcost_curr(1) = 0;
Icost_curr(1) = 0;
for i = 1:length(kW_curr)
    Icost_curr(i) = 2*(econ.wind.installed - (0.5*kWcost_curr(i))/ ...
            (kW_curr(i)*econ.curr.costmult))*kW_curr(i)*econ.curr.costmult; %installation
end

%dies- costs
al_cost = 103.92; %$/sqft for 1/2" 6061 AL (2025)
al_cost = al_cost/0.0929; %convert to $/m2



kWcost_dies = 2.*(polyval(opt.p_dev.d_cost,kW_dies).*econ.dies.gcm + econ.dies.autostart); %generator (with spares provisioning: 2)
genencl = 2*6.*(gen_vol.^(2/3)).*al_cost;
fuel = runtime_tot.*lph.*econ.dies.fcost; %cost of consumed fuel
Icost_dies = 0.1.*(kWcost_dies + genencl); %installation

kWcost_dies(1) = 0;
genencl(1) = 0;
fuel(1) = 0;
Icost_dies(1) = 0;



Scost = (2)*batt.cost.*Smax; %2x so there is a spare on sure
battencl = 2*6.*(batt_len.^2).*al_cost;
Icost_batt = 0.1.*(Scost+battencl); %installation cost (10% is a guess from Brian)

%%Old mooring model
% %%New Mooring Model
% tempPmooring = nan(length(kW_dies),3);
% tempmass = nan(length(kW_dies),3);
% tempPmtrl = nan(length(kW_dies),3);
% tempMTRLMOOR = nan(length(kW_dies),3);
% Buoy_Mass = nan(length(kW_dies),1);
% dp = nan; %platform diameter is not needed for Task 2 modeling
% for b = 1:3 %loop through buoy sizes
%     tempPmooring(:,b) = 2*interp1(econ.platform.payloadmass(b,:),econ.platform.cost(b,:),comp_plat_mass,'linear'); %mooring cost - this should not extrapolate
%     tempmass(:,b) = interp1(econ.platform.payloadmass(b,:),econ.platform.mass(b,:),comp_plat_mass,'linear');
%     tempPmtrl(:,b) = 2*(tempmass(:,b)-comp_plat_mass)*econ.platform.steel;
%     tempMTRLMOOR(:,b) = tempPmtrl(:,b) + tempPmooring(:,b);
% end
% %[Pmooring,indMoor] = min(tempPmooring,[],2); %old - find min mooring cost and index of that min
% [~,indPM] = min(tempMTRLMOOR,[],2); %find index of lowers material + mooring cost
% 
% Pmooring = nan(length(indPM),1);
% for m = 1:length(indPM) %this is only a loop from the Hybrid Econ code which worked with matrices (length of indPM should always be 1)
%     Buoy_Mass(m,1) = tempmass(m,indPM(m)) - comp_plat_mass(m);
%     Pmooring(m,1) = tempPmooring(m,indPM(m));
%     if length(indPM) > 2
%         disp('WARNING: somehow the mooring cost is a vector')
%     end
% end
% 
% Pmtrl = 2*Buoy_Mass*econ.platform.steel;  %platform material cost - UPDATED BASED ON INFO FROM DEREK
% t_i = interp1(econ.platform.d_i,econ.platform.t_i,depth, ...
%     'linear','extrap'); %installation time
% if depth < 500
%     t_i = econ.platform.t_i (1); %make min time = 6 hours
% end
% Pinst = econ.vessel.speccost*(t_i/24); %platform instllation

% %triptime = dist*kts2mps(econ.vessel.speed)^(-1)*(1/86400); %[d]
% triptime = 0; %cost of trip attributed to instrumentation
% C_v = econ.vessel.speccost; %assuming that all trips require the spec vessel since everything is being replaced
% t_os = 2*t_i/24; %assuming recovery is similar to installation (have to recover and re-deploy)
% vesselcost = C_v*(nvi*(2*triptime + t_os)); %vessel cost

%Maintenance Costs
solarrepair = 1/4*0.5.*(Strcost_inso + Icost_inso + Ecost_inso).*(nvi); %solar refurb repair
solarrepair(1) = 0;

turbrepair = 1/2*0.5.*(kWcost_wind+Icost_wind).*(nvi); %turbine repair cost
turbrepair(1) = 0; 

cturbrepair = 1/2*0.5.*(kWcost_curr + Icost_curr).*(nvi); %current turbine refurb cost
cturbrepair(1) = 0;

wecrepair = 1/2*(0.5).*(kWcost_wave+Icost_wave).*(nvi); %wec repair cost
wecrepair(1) = 0; %if nvi = 0, wec repair must be 0

genrepair = 1/4*(0.5).*(kWcost_dies+Icost_dies).*(nvi); %generator repair cost
genrepair(1) = 0;
battreplace = Scost.*newbatt/2; %number of battery replacements
batteryrepair = 1/4*0.5.*(Scost + Icost_batt).*nvi; %battery repair
%mooringrepair = 1/4*0.5*Pmooring*nvi; %eventually need a mooring refurb cost

%Component Totals
BattTot = battencl + Scost + Icost_batt + battreplace + batteryrepair;
DiesTot = genencl + kWcost_dies + Icost_dies + fuel + genrepair;
SolarTot = Icost_inso + Mcost_inso + Ecost_inso + Strcost_inso + solarrepair;
WindTot = Icost_wind + kWcost_wind + turbrepair;
WaveTot = Icost_wave + kWcost_wave + wecrepair;
CurrTot = kWcost_curr + Icost_curr + cturbrepair;

BattMass = mass_batt + mass_battencl;
DiesMass = mass_diesencl + mass_fuel + mass_dies;
SolarMass = mass_solar + mass_solar_E + mass_solar_S;
WindMass = mass_wind;
WaveMass = mass_wec;
CurrMass = zeros(size(mass_wec));

%%
%Define color options for 2% space
colors = brewermap(11, 'Set3');
colors2 = brewermap(11, 'Set2');
colors3 = brewermap(11, 'PRGn');
colpm(1,:) = colors(1,:);
colpm(2,:) = colors(6,:);
colpm(3,:) = colors(5,:);
colpm(4,:) = colors2(8,:);
colpm(5,:) = colors3(4,:);
colb(1,:) = colors2(4,:);

fs = 10; %font size
figure
set(gcf,'Unit','Inches','OuterPosition',[0.3,0.3,5,4])
tf = tiledlayout(1,2);
tf.Padding = 'compact';
nexttile
hold on
plot(Smax./500,BattTot./1000,'linewidth',3,'DisplayName','Battery','color',colb(1,:))
plot(kW_dies./8,DiesTot./1000,'linewidth',3,'DisplayName','Diesel','color',colpm(4,:))
plot(kW_inso./8,SolarTot./1000,'linewidth',3,'DisplayName','Solar','color',colpm(2,:))
plot(kW_wind./8,WindTot./1000,'linewidth',3,'DisplayName','Wind','color',colpm(1,:))
plot(kW_wave./8,WaveTot./1000,'linewidth',3,'DisplayName','Wave','color',colpm(3,:))
plot(kW_curr./8,CurrTot./1000,'linewidth',3,'DisplayName','Current','color',colpm(5,:))
grid on
box on
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
xlabel('Fraction of max capacity','Interpreter','latex','FontSize',fs)
ylabel('Component Cost [\$1000]','Interpreter','latex','FontSize',fs)
text(0.05, 0.93,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')
lh = legend('Location','NorthOutside','Orientation','Horizontal','Interpreter','Latex','FontSize',fs,'NumColumns',6);
lh.Layout.Tile = 'North';
lh.ItemTokenSize = [15, 18]; % Makes the lines shorter

nexttile
hold on
plot(Smax./500,BattMass,'linewidth',3,'color',colb(1,:))
plot(kW_dies./8,DiesMass,'linewidth',3,'color',colpm(4,:))
plot(kW_inso./8,SolarMass,'linewidth',3,'color',colpm(2,:))
plot(kW_wind./8,WindMass,'linewidth',3,'color',colpm(1,:))
plot(kW_wave./8,WaveMass,'linewidth',3,'color',colpm(3,:))
plot(kW_curr./8,CurrMass,'linewidth',3,'color',colpm(5,:))
grid on
box on
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
xlabel('Fraction of max capacity','Interpreter','latex','FontSize',fs)
ylabel('Component Mass [kg]','Interpreter','latex','FontSize',fs)
text(0.05, 0.93,'(b)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')