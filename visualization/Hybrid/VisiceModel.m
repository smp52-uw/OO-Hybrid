%Visualize Ice model

%load data
optInputs
loc = 'BerSea';
data = load(loc,'data');
data = data.('data');

[data, opt] = prepHybrid(data,opt,uc(2),wave,atmo,inso,cturb);

%nominal generators
kW_wind = 2;
kW_inso = 2;

%Calculate Power for wind and solar
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

%% Set colors
colP = colormap(brewermap(12,'Set3')); %colors
colW = colormap(brewermap(11,'BrBG')); %colors
colT = colormap(brewermap(11,'PiYG')); %colors
fs = 10;
%% plot
U = data.wind.U; %this is the same as data.met.wind_spd but not affected by ice
time = data.met.time;
T = data.temperature.T2mw;
%ice models
figure
set(gcf,'Position', [100, 100, 400, 600])
tf = tiledlayout(2,1);
tf.Padding = 'compact';
tf.TileSpacing = 'compact';

ax(1) = nexttile(1); %met data
colororder([colW(9,:);colT(3,:)])

yyaxis left
plot(datetime(time,'convertfrom','datenum'),U,'linewidth',1.2)
ylabel('Wind Speed [m/s]','Interpreter','Latex','FontSize',fs)
grid on
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
yyaxis right
plot(datetime(time,'convertfrom','datenum'),T,'linewidth',1.2)
hold on
ylabel('Air Temperature [C]','Interpreter','Latex','FontSize',fs)
text(0.05, 0.93,'(a)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')

ax(2) = nexttile(2); %power data
plot(datetime(time,'convertfrom','datenum'),Pwind./kW_wind/1000,'linewidth',1.2,'color',colP(1,:),'DisplayName','Wind')
hold on
plot(datetime(time,'convertfrom','datenum'),Pinso./kW_inso/1000,'linewidth',1.2,'color',colP(6,:),'DisplayName','Solar')
ylabel('Normalized Power [-]','Interpreter','Latex','FontSize',fs)
xlabel('Time','Interpreter','Latex','FontSize',fs)
legend('Location','NorthEast','Interpreter','Latex','FontSize',fs)
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
text(0.05, 0.93,'(b)','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')
grid on

linkaxes(ax,'x')
