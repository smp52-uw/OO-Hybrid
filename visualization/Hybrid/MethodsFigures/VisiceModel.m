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
windNI = data.wind.U; %this is the same as data.met.wind_spd but not affected by ice
swsoNI = data.swso_noice; %this is the same as data.met.wind_spd but not affected by ice
swso = data.swso;
wind = data.met.wind_spd; %[m/s]
if atmo.dyn_h %use log law to adjust wind speed based on rotor height
    wind = adjustHeight(wind,data.met.wind_ht, ...
    turb.clearance + sqrt(1000*2*kW_wind/(turb.eta*atmo.rho_a_c*pi*turb.ura^3)) ...
    ,'log',(data.wind.z0.*1000)); 

    windNI = adjustHeight(windNI,data.met.wind_ht, ...
    turb.clearance + sqrt(1000*2*kW_wind/(turb.eta*atmo.rho_a_c*pi*turb.ura^3)) ...
    ,'log',(data.wind.z0.*1000)); 

end
Pwind_ci = (1/2)*atmo.rho_a_c*turb.uci^3; %cut in power flux [W]
Pwind_co = (1/2)*atmo.rho_a_c*turb.uco^3; %cut out power flux [W]
Pwind_ra = (1/2)*atmo.rho_a_c*turb.ura^3; %rated power flux [W]
%find power from wind turbine 
Pflux = (1/2).*data.wind.rho.*wind.^3;
Pwind = kW_wind.*1000.*data.wind.rho.*wind.^3./(atmo.rho_a_c.*turb.ura.^3); %[W]
ind_wiCI = find(Pflux < Pwind_ci);
ind_wiCO = find(Pflux > Pwind_co);
ind_wiRA = find(Pflux > Pwind_ra);

Pwind(ind_wiRA) = kW_wind*1000; %[W] cap at rated power
Pwind(ind_wiCI) = 0; %[W] no power below CI
Pwind(ind_wiCO) = 0; %[W] no power above CO

PfluxNI = (1/2).*data.wind.rho.*windNI.^3;
PwindNI = kW_wind.*1000.*data.wind.rho.*windNI.^3./(atmo.rho_a_c.*turb.ura.^3); %[W]
ind_wiCI = find(PfluxNI < Pwind_ci);
ind_wiCO = find(PfluxNI > Pwind_co);
ind_wiRA = find(PfluxNI > Pwind_ra);

PwindNI(ind_wiRA) = kW_wind*1000; %[W] cap at rated power
PwindNI(ind_wiCI) = 0; %[W] no power below CI
PwindNI(ind_wiCO) = 0; %[W] no power above CO


%set panel degradation
eff = (1-((inso.deg/100)/8760)*(1:1:length(swso)));
%rain = repmat(linspace(0.5,0,24*30),[1,ceil(length(swso)/(24*30))]);
if econ.inso.scen == 1 %automated cleaning
    d_soil_eff = 0;
    econ.inso.marinization = econ.inso.marinization*2;
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

    PinsoNI = eff_t./inso.eff.*kW_inso.*1000.*(swsoNI./(inso.rated.*1000)); %[W]
    ind_insoRa = find(PinsoNI >= (eff_t./inso.eff.*kW_inso.*1000)); %indices with power above rated
    PinsoNI(ind_insoRa) = eff_t(ind_insoRa)./inso.eff.*kW_inso.*1000;
else %no solar panel
    Pinso = zeros(1,length(time));
end

%% Set colors
colP = colormap(brewermap(12,'Set3')); %colors
colW = colormap(brewermap(11,'BrBG')); %colors
colT = colormap(brewermap(11,'PiYG')); %colors
colG = colormap(brewermap(9,'Set1')); %grey for no ice
fs = 10;
%% plot
time = data.met.time;
T = data.temperature.T2mw;

sttime = datetime(2018,01,08);
edtime = datetime(2018,02,21);
%ice models
figure
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 2.25, 2.5])
tf = tiledlayout(2,1);
tf.Padding = 'loose';
tf.TileSpacing = 'compact';
tf.InnerPosition = [0.31,0.17685,0.6,0.65];

indwind = windNI > 12;
indtemp = T < -2;
indaccumulation = find(indwind & indtemp);
inddiff = diff(indaccumulation);

ax(1) = nexttile(1); %met data
colororder([colW(9,:);colT(3,:)])
plot(datetime(time,'convertfrom','datenum'),windNI,'linewidth',1.2,'color',colW(9,:),'DisplayName','Data')
hold on
%yline(12,'--k','linewidth',1.2,'DisplayName','Ice Threshold')
for xx = 1:length(indaccumulation) - 1
    if inddiff(xx) <= 1
        xregion(datetime(time(indaccumulation(xx)),'convertfrom','datenum'),datetime(time(indaccumulation(xx+1)),'convertfrom','datenum'),'FaceColor', [0 0.447 0.741],'FaceAlpha',0.3)
    end
end
xlim([sttime,edtime])
yl = ylabel({'Wind',' Speed',' [m/s]'},'Interpreter','Latex','FontSize',fs,'Rotation',0,'VerticalAlignment','middle');
set(yl,'units','normalized')
ylp1 = yl.Position;
ylp1(1) = ylp1(1)*2.7;
set(yl,'position',ylp1)
grid on; box on
lg = legend('','Ice accumulation','Interpreter','Latex','FontSize',fs);
lg.ItemTokenSize(1)= 10;
lg.Layout.Tile = 'north';
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
text(0.05, 0.85,'\textbf{(a)}','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')
xticks([])

clear yl
ax(2) = nexttile(2); %met data
plot(datetime(time,'convertfrom','datenum'),T,'linewidth',1.2,'color',colT(3,:),'DisplayName','Data')
hold on
%yline(-2,'--k','linewidth',1.2,'DisplayName','Ice Threshold')
for xx = 1:length(indaccumulation) - 1
    if inddiff(xx) <= 1
        xregion(datetime(time(indaccumulation(xx)),'convertfrom','datenum'),datetime(time(indaccumulation(xx+1)),'convertfrom','datenum'),'FaceColor', [0 0.447 0.741],'FaceAlpha',0.3)
    end
end
yl = ylabel({'Air',' Temp.','[$^\circ$C]'},'Interpreter','Latex','FontSize',fs,'Rotation',0,'VerticalAlignment','middle');
set(yl,'units','normalized')
ylp2 = yl.Position;
ylp2(1) = ylp1(1);
set(yl,'position',ylp2)
%lg=legend('Interpreter','Latex','FontSize',fs);
%lg.ItemTokenSize(1)= 10;
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
text(0.05, 0.85,'\textbf{(b)}','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')
xlabel('Time','Interpreter','Latex','FontSize',fs)
grid on; box on
xlim([sttime,edtime])
xtickformat('MMM dd');

figure
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 2.25, 2.5])
tf = tiledlayout(2,1);
tf.Padding = 'loose';
tf.TileSpacing = 'compact';
tf.InnerPosition = [0.31,0.17685,0.6,0.65];

clear yl
ax(3) = nexttile(1); %power data
hold on
plot(datetime(time,'convertfrom','datenum'),PwindNI/1000,'linewidth',1.2,'color',colG(9,:),'DisplayName','Wind')
plot(datetime(time,'convertfrom','datenum'),Pwind/1000,'linewidth',1.2,'color',colP(1,:),'DisplayName','Wind')
yl=ylabel({'Wind',' Power',' [kW]'},'Interpreter','Latex','FontSize',fs,'Rotation',0,'VerticalAlignment','middle');
set(yl,'units','normalized')
ylp3 = yl.Position;
ylp3(1) = ylp1(1);
set(yl,'position',ylp3)
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
text(0.05, 0.85,'\textbf{(c)}','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')
grid on; box on
lg = legend('Neglecting Ice','','Interpreter','Latex','FontSize',fs);
lg.ItemTokenSize(1)= 10;
lg.Layout.Tile = 'north';
xlim([sttime,edtime])
xticks([])

clear yl
ax(4) = nexttile(2); %power data
hold on
plot(datetime(time,'convertfrom','datenum'),PinsoNI/1000,'linewidth',1.2,'color',colG(9,:),'DisplayName','Solar')
plot(datetime(time,'convertfrom','datenum'),Pinso/1000,'linewidth',1.2,'color',colP(6,:),'DisplayName','Solar')
yl = ylabel({'Solar',' Power',' [kW]'},'Interpreter','Latex','FontSize',fs,'Rotation',0,'VerticalAlignment','middle');
set(yl,'units','normalized')
ylp4 = yl.Position;
ylp4(1) = ylp1(1);
set(yl,'position',ylp4)
xlabel('Time','Interpreter','Latex','FontSize',fs)
set(gca,"TickLabelInterpreter",'latex','FontSize',fs)
text(0.05, 0.85,'\textbf{(d)}','Units','Normalized', ...
    'VerticalAlignment','middle','FontWeight','normal', ...
    'FontSize',fs,'Interpreter','latex')
grid on; box on
xlim([sttime,edtime])
xtickformat('MMM dd');
%linkaxes(ax,'x')
