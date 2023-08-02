function [cost,surv,CapEx,OpEx,kWcost_dies, kWcost_wave,kWcost_curr, kWcost_wind,Mcost_inso,Ecost_inso,Icost_inso,...
    Strcost_inso, Icost_wave, Icost_wind, Scost,Pmtrl,Pinst,Pmooring, ...
    vesselcost,genrepair,turbrepair, wecrepair, battreplace,battencl,genencl,fuel, ...
    triptime,runtime,nvi,batt_L1,batt_L2, batt_lft1,batt_lft2, nfr,noc,nbr,dp,S1,S2,Pdies,Pinso,Pwind,Pwave,Pcurr,Ptot,width,cw,D,L,F,eff_t,pvci] =  ...
    simHybridsg(kW_dies, kW_inso, kW_wind, kW_wave,kW_curr, Smax,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb) %#codegen

%% Created by Sarah Palmer Aug 2023 - started from Trent's OO-Tech code
%Surrogate simulation that doesn't take into account the battery
%degradation

ID = [kW_dies kW_inso kW_wind kW_wave kW_curr Smax];

%if physically impossible, set S_temp and C_temp to failed values
if opt.fmin && Smax < 0 || min([kW_dies,kW_inso,kW_wind,kW_wave,kW_curr]) < 0
    surv = 0;
    cost = inf;
    return
end

%set capture width modifier
%cw_mod = wave.cw_mod;
%Number of vessel intervention
nvi = (uc.lifetime*12)/uc.SI;
%set burn rate
lph = polyval(opt.p_dev.d_burn,kW_dies); %[l/h], burn rate
if kW_dies == 0  %Set burn rate = 0 for no diesel generator (in case the polyfit doesn't cross (0,0))
    lph = 0; 
end
%extract data
time = datenum(data.met.time);
%time = data.met.time; %THIS IS IN DATETIME
T = length(time);
dt = 24*(time(2) - time(1)); %time in hours
%dt = 1; %time in hours

%dist = data.dist; %[m] distance to shore
%depth = data.depth; %[m] water depth
swso = data.swso;
wind = data.met.wind_spd; %[m/s]
if atmo.dyn_h %use log law to adjust wind speed based on rotor height
    wind = adjustHeight(wind,data.met.wind_ht, ...
            turb.clearance + ...
            sqrt(1000*2*kW_wind/(turb.eta*atmo.rho_a*pi*turb.ura^3)) ...
            ,'log',atmo.zo);
end
wavepower = opt.wave.wavepower_ts; %wavepower timeseries
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
Pwave(Pwave<0) = 0; %no negative power
Pwave(Pwave>kW_wave) = kW_wave; %no larger than rated power
Pwave = Pwave*1000; %convert to watts
if kW_wave < 0.2144
    Pwave = zeros(1,length(time));
    %disp('Zero wave power due to WECSIM min')
end
%current power
t_depth = cturb.clearance + sqrt(1000*2*kW_curr/(cturb.eta*atmo.rho_w*cturb.ura^3)); %current turb height
[~,cturb_depth] = min(abs(t_depth - data.curr.depth));

c_speed = data.curr.speed6a(:,cturb_depth); 


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
t2 = tic;

%initialize diagnostic variables
S1 = ones(1,length(time))*Smax; %battery level timeseries

Ptot = zeros(1,length(time)); %power produced timeseries
Pdies = zeros(1,length(time));

D = zeros(1,length(time)); %power dumped timeseries
L = uc.draw; %power to load time series
F = zeros(1,length(time)); %failure series

%surv = 1;
charging = false;
runtime = zeros(1,nvi); %[h], amount of time spent running
dies_vol = zeros(1,nvi);

newbatt = 0; %number of new batteries that have to be bought
clear clean_ind batt_lft1

wind_p = wind;
wind_p(wind_p<turb.uci) = 0;
wind_p(wind_p>turb.uco) = 0;
wind_p(wind_p>turb.ura) = turb.ura;
Pwind = kW_wind.*1000.*wind_p.^3./turb.ura^3; %[W]

c_speed_p = c_speed;
c_speed_p(c_speed_p<cturb.uci) = 0;
c_speed_p(c_speed_p>cturb.uco) = 0;
c_speed_p(c_speed_p>cturb.ura) = cturb.ura;
Pcurr = kW_curr.*1000.*c_speed_p.^3./cturb.ura.^3; %[W]

    %find solar efficiency
%     soil_eff = soil_eff - d_soil_eff;
%     if soil_eff < 0 
%         soil_eff = 0; %no negative efficiency
%     end

eff_t = eff.*soil_eff.*inso.eff;
swso_p = swso;
swso_p(swso_p > inso.rated*1000) = inso.rated*1000;
Pinso = eff_t./inso.eff.*kW_inso.*1000.*(swso_p./(inso.rated*1000)); 

Prenew = Pwave + Pwind + Pinso + Pcurr; %total renewable power

%% run simulation
for t = 1:(length(time))
    sd1 = S1(t)*(batt.sdr/100)*(1/(30*24))*dt; %[Wh] self discharge
    if ~charging %generator off
        S1(t+1) = S1(t) + (Prenew(t) - uc.draw(t)).*dt - sd1;
    else %generator on
        Pdies(t) = kW_dies*1000;
        S1(t+1) = S1(t) + (Prenew(t)+Pdies(t) - uc.draw(t)).*dt - sd1;
        if opt.drun == 1
            charging = false;
        end
        %track fuel vol and runtime
        int_time = T/nvi; %dividing hours by number of intervention giving number of hours between each intervention?
        int_id = ceil(t/int_time); %calculated what intervention number we are on
        runtime(1,int_id) = runtime(1,int_id) + dt; %[h] when diesel gen is on
        dies_vol(1,int_id) = dies_vol(1,int_id) + lph*dt; %[L] burn rate multiplied by hours running should give total consumption
        if runtime(1,int_id) >= 250 || dies_vol(1,int_id) >= 800 %Brian set the maximums (250,800)
                %runtime(1,int_id)
                %dies_vol(1,int_id)
                surv = 0;
                cost = inf;
                %disp('Option not viable due to fuel volume or runtime')
                return
        end
    end
    if S1(t+1) >= (Smax*1000)
        D(t) = S1(t+1) - ((Smax*1000)); %[Wh]
        S1(t+1) = (Smax*1000);
        if opt.drun == 2
            charging = false; %turn off gen when battery is full
        end
    end
    Ptot(t) = Prenew(t) + Pdies(t); %total power
    if S1(t+1) <= Smax*batt.dmax*1000 %bottomed out
        S1(t+1) = dt*Ptot(t) + S1(t) - sd1;
        L(t) = 0; %load drops to zero
        F(t) = 1; %failure tracker
    end
end


%% Cost Calculation - COMMENTED OUT FOR WEIGHT APPROX
% % % %costs - solar
% % % Mcost_inso = 2*econ.inso.module*kW_inso*econ.inso.marinization*econ.inso.pcm; %module
% % % Icost_inso = 2*econ.inso.installation*kW_inso*econ.inso.pcm; %installation
% % % Ecost_inso = 2*econ.inso.electrical*kW_inso*econ.inso.marinization ...
% % %     *econ.inso.pcm; %electrical infrastructure
% % % Strcost_inso = 2*econ.inso.structural*kW_inso*econ.inso.marinization ...
% % %     *econ.inso.pcm; %structural infrastructure
% % % %costs - wind
% % % kWcost_wind = 2*polyval(opt.p_dev.t,kW_wind)*econ.wind.marinization ...
% % %     *econ.wind.tcm; %turbine
% % % if kW_wind == 0, kWcost_wind = 0; end
% % % if kW_wind ~=0 %can't calculate Icost with zero kWcost_wind
% % %     Icost_wind = 2*(econ.wind.installed - 0.5*kWcost_wind/ ...
% % %         (kW_wind*econ.wind.marinization*econ.wind.tcm))*kW_wind; %installation
% % %     if Icost_wind < 0, Icost_wind = 0; end
% % % else
% % %     Icost_wind = 0;
% % % end
% % % %costs - wave
% % % kWcost_wave = 2*econ.wave.costmult*polyval(opt.p_dev.t,kW_wave); %wec
% % % if kW_wave == 0, kWcost_wave = 0; end
% % % if kW_wave ~= 0 %no Icost if kW_wave = 0
% % %     Icost_wave = 2*(econ.wind.installed - (0.5*kWcost_wave)/ ...
% % %         (kW_wave*econ.wave.costmult))*kW_wave; %installation
% % %     if Icost_wave < 0, Icost_wave = 0; end
% % % else
% % %     Icost_wave = 0;
% % % end
% % % 
% % % %dies- costs
% % % if kW_dies == 0 
% % %     kWcost_dies = 0;
% % %     genencl = 0;
% % %     fuel = 0;
% % %     mass_dies = 0;
% % % else
% % %     kWcost_dies = polyval(opt.p_dev.d_cost,kW_dies)*2*econ.dies.gcm + ...
% % %         econ.dies.autostart; %generator (with spares provisioning: 2)
% % %     genencl = 2*polyval(opt.p_dev.d_size,kW_dies)^3* ...
% % %         (econ.dies.enclcost/econ.dies.enclcap); %generator enclosure
% % %     fuel = runtime_tot*lph*econ.dies.fcost; %cost of consumed fuel
% % %     mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of generator
% % % end
% % % 
% % % if bc == 1 %lead acid
% % %     if Smax < opt.p_dev.kWhmax %less than linear region
% % %         Scost = polyval(opt.p_dev.b,Smax);
% % %     else %in linear region
% % %         Scost = polyval(opt.p_dev.b,opt.p_dev.kWhmax)* ...
% % %             (Smax/opt.p_dev.kWhmax);
% % %     end
% % % elseif bc == 2 %lithium phosphate
% % %     Scost = (2+newbatt)*batt.cost*Smax; %2x so there is a spare on sure
% % % end
% % % battencl = 2*econ.batt.enclmult*Scost; %battery enclosure cost
% % % Pmtrl = 2*(1/1000)*econ.platform.wf*econ.platform.steel* ...
% % %     (mass_dies+inso.wf*kW_inso/(inso.rated*inso.eff)+kW_wind*turb.wf); %platform material
% % % 
% % % %% Will need new way to calculate dp
% % % dp = polyval(opt.p_dev.d_size,max(ID(1:end-1))); %Max ID 1:end-1 will give max of the rated generation
% % % if dp < 2, dp = 2; end 
% % % if dp > 8 %apply installtion time multiplier if big mooring - WRONG INSO.BOUNDARY needs to be updated
% % %     inst_mult = interp1([8 econ.platform.boundary_di], ...
% % %         [1 econ.platform.boundary_mf],dp);
% % % else
% % %     inst_mult = 1;
% % % end
% % % t_i = interp1(econ.platform.d_i,econ.platform.t_i,depth, ...
% % %     'linear','extrap')*inst_mult; %installation time
% % % Pinst = econ.vessel.speccost*(t_i/24); %platform instllation
% % % 
% % % %% Will need new mooring model
% % % Pmooring = 2*interp2(econ.platform.diameter, ...
% % %     econ.platform.depth, ...
% % %     econ.platform.cost,dp,depth,'linear'); %mooring cost
% % % 
% % % 
% % % %triptime = dist*kts2mps(econ.vessel.speed)^(-1)*(1/86400); %[d]
% % % triptime = 0; %cost of trip attributed to instrumentation
% % % C_v = econ.vessel.speccost; %assuming that all trips require the spec vessel since everything is being replaced
% % % 
% % % vesselcost = C_v*(nvi*(2*triptime + t_i)); %vessel cost
% % % %Maintenance Costs
% % % turbrepair = 1/2*0.5*(kWcost_wind+Icost_wind)*(nvi); %turbine repair cost
% % % if turbrepair < 0, turbrepair = 0; end
% % % wecrepair = 1/2*(0.5)*(kWcost_wave+Icost_wave)*(nvi); %wec repair cost
% % % if wecrepair < 0, wecrepair = 0; end %if nvi = 0, wec repair must be 0
% % % genrepair = 1/2*(0.5)*(kWcost_dies+genencl)*(nvi); %generator repair cost
% % % battreplace = Scost*nvi; %number of battery replacements
% % % mooringrepair = 0; %eventually need a mooring refurb cost
% % % 
% % % %Total Cost Calculations
% % % CapEx = Pmooring + Pinst + Pmtrl + battencl + Scost + ...
% % %     genencl + kWcost_dies + Icost_inso + Mcost_inso + Ecost_inso + ...
% % %     Strcost_inso + Icost_wind + kWcost_wind + Icost_wave + kWcost_wave;
% % % OpEx = fuel + battreplace + genrepair + vesselcost + turbrepair + wecrepair;
% % % cost = CapEx + OpEx;

%% Weight Approximation
if opt.tar == 1
    mass_solar = inso.wf*kW_inso/(inso.rated*inso.eff); %[kg] - 1 array
    mass_solar_E = 0; %electrical - for 1 buoy
    %mass_solar_S = 3.15*kW_inso/(inso.rated*inso.eff); %structural - for 1 buoy
    mass_solar_S = 0;
    %mass_wind = kW_wind*turb.wf; %[kg] - 1 turbine
    mass_wind = kW_wind*turb.wf; %updated to include a tower
    
    mass_curr = kW_curr*cturb.wf; %total guess
    
    mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of 1 generator [kg]
    if kW_dies == 0
        mass_dies = 0;
    end
    al_plate = 35.24; %[kg/m2]
    mass_diesencl = 6*(polyval(opt.p_dev.d_size,kW_dies)^2)*al_plate; %[kg]
    if kW_dies == 0
        mass_diesencl = 0;
    end
    dies_dens = 0.85; %[g/cm3] = [kg/L] from a chevron report
    mass_fuel = runtime_tot*lph*dies_dens; %[kg]
    
    mass_batt = Smax/(batt.se*batt.V/1000); %[kg] - 1 battery
    batt_vol = polyval(opt.p_dev.b_size,Smax);
    batt_len = batt_vol^(1/3); %assume cube
    mass_battencl = 6*(batt_len^2)*al_plate; %[kg]
    %old assumption is that 1 encl $ = batt $
    %mass_battencl = mass_batt;
    % Pmtrl = econ.platform.wf* ...
    %     (mass_dies + mass_diesencl + mass_fuel + mass_solar + mass_wind + mass_batt + ...
    %     mass_battencl + mass_solar_E + mass_solar_S); % 1 platform material [kg]
    %mass_wec = Pmtrl * 0.376; %[kg] based on the percent of the RM3 that isn't float mass
    mass_wec = kW_wave*895.78;
    if kW_wave == 0
        mass_wec = 0;
    end
    buoy_m = mass_solar + mass_solar_E + mass_solar_S + mass_wind + mass_curr + mass_dies + ...
        mass_diesencl + mass_fuel + mass_wec + mass_batt + mass_battencl;
    newbatt_m = newbatt * mass_batt; %[g]
    cost = buoy_m*2 + newbatt_m; %[kg] - 2 buoy weight + weight of extra batteries
    % if kW_wave < 0.2144 && kW_wave ~= 0
    %     cost = inf; %make cost infinity for bad wave power cases - since I'm only adjusting cost it shouldn't affect per_opt
    % end
else %Lowest gen capacity target
    mass_solar = turb.wf*kW_inso; %[kg] - 1 array
    mass_solar_E = 0; %electrical - for 1 buoy
    mass_solar_S = 0;
    mass_wind = kW_wind*turb.wf; %updated to include a tower
    
    mass_curr = kW_curr*turb.wf; %total guess
    
    mass_dies = turb.wf*kW_dies; %mass of 1 generator [kg]
    if kW_dies == 0
        mass_dies = 0;
    end
    mass_diesencl = 0;
    mass_fuel = 0; %[kg]
    al_plate = 35.24; %[kg/m2]
    mass_batt = Smax/(batt.se*batt.V/1000); %[kg] - 1 battery
    batt_vol = polyval(opt.p_dev.b_size,Smax);
    batt_len = batt_vol^(1/3); %assume cube
    mass_battencl = 6*(batt_len^2)*al_plate; %[kg]
    
    mass_wec = kW_wave*turb.wf;
    if kW_wave == 0
        mass_wec = 0;
    end
    buoy_m = mass_solar + mass_solar_E + mass_solar_S + mass_wind + mass_curr + mass_dies + ...
        mass_diesencl + mass_fuel + mass_wec + mass_batt + mass_battencl;
    newbatt_m = newbatt * mass_batt; %[g]
    cost = buoy_m*2 + newbatt_m; %[kg] - 2 buoy weight + weight of extra batteries
end
%% dummy output variables that aren't assigned without economics
CapEx = nan;
OpEx = nan;
kWcost_dies = mass_dies;
kWcost_wave = mass_wec;
kWcost_curr = mass_curr;
%kWcost_curr = nan;
kWcost_wind = mass_wind;
Mcost_inso = mass_solar;
Ecost_inso = mass_solar_E;
Icost_inso = nan;
Strcost_inso = mass_solar_S;
Icost_wave = nan;
Icost_wind = nan;
Scost = mass_batt;
Pinst = nan;
Pmooring = nan;
Pmtrl = nan;
vesselcost = nan;
genrepair = nan;
turbrepair = nan;
wecrepair = nan;
battreplace = newbatt_m;
battencl = mass_battencl;
genencl = mass_diesencl;
fuel = mass_fuel;
triptime = nan;
dp = nan;
pvci = nan;

%nan values due to surrogate sim
nvi = nan;
batt_lft1 = nan;
batt_lft2 = nan;
batt_L1 = nan;
batt_L2 = nan;
noc = nan;
nbr = nan;
nfr = nan;
S2 = nan;

%determine if desired uptime was met
surv = sum(L == uc.draw)/(length(L));
% if sum(L == uc.draw)/(length(L)) < uc.uptime
%     %surv = 0;
%     if opt.fmin
%         cost = inf;
%     end
% end

end
