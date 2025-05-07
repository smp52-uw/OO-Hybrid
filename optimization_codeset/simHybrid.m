function [cost,surv,CapEx,OpEx,kWcost_dies, kWcost_wave,kWcost_curr, kWcost_wind,Mcost_inso,Ecost_inso,Icost_inso,...
    Strcost_inso, Icost_wave, Icost_wind, Scost,Pmtrl,Pinst,Pmooring, ...
    vesselcost,genrepair,turbrepair, wecrepair, battreplace,battencl,genencl,fuel, ...
    triptime,runtime,nvi,batt_L1,batt_L2, batt_lft1,batt_lft2, nfr,noc,nbr,dp,S1,S2,Pdies,Pinso,Pwind,Pwave,Pcurr,Ptot,width,cw,D,L,F,eff_t,pvci] =  ...
    simHybrid(kW_dies, kW_inso, kW_wind, kW_wave,kW_curr, Smax,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb) %#codegen

%% Created by Sarah Palmer Jan 2023 - started from Trent's OO-Tech code

ID = [kW_dies kW_inso kW_wind kW_wave kW_curr Smax];
%disp('made it into SimHybrid')
% if kW_wave < 0.2144
%     surv = 0;
%     cost = inf;
%     disp('Option not viable due to WECSIM min')
%     return
% end
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
depth = data.depth; %[m] water depth
swso = data.swso;
wind = data.met.wind_spd; %[m/s]
if atmo.dyn_h %use log law to adjust wind speed based on rotor height
    for i = 1:length(wind)
        wind(i) = adjustHeight(wind(i),data.met.wind_ht, ...
            turb.clearance + ...
            sqrt(1000*2*kW_wind/(turb.eta*atmo.rho_a*pi*turb.ura^3)) ...
            ,'log',atmo.zo);
    end
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
% % [Pcurr_m,kW_testc,c_vel] = Current_Power(cturb,atmo);
%current_Pmatrix = load('current_Pmatrix.mat');
% Pcurr_m = current_Pmatrix.Pcurr.P;
% kW_testc = current_Pmatrix.Pcurr.kW_curr;
% c_vel = current_Pmatrix.Pcurr.c_vel;
%F_curr = opt.curr.F_curr;

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
S1 = zeros(1,length(time)); %battery level timeseries
S2 = zeros(1,length(time)); %battery level timeseries
Ptot = zeros(1,length(time)); %power produced timeseries
Pdies = zeros(1,length(time));
Pinso = zeros(1,length(time));
Pwind = zeros(1,length(time));
Pcurr= zeros(1,length(time));
Prenew = zeros(1,length(time));
D = zeros(1,length(time)); %power dumped timeseries
L = uc.draw; %power to load time series
F = zeros(1,length(time)); %failure series
batt_L1 = zeros(1,length(time)); %battery L (degradation) timeseries
fbi1 = 1; %fresh battery index
batt_L2 = zeros(1,length(time)); %battery L (degradation) timeseries
fbi2 = 1; %fresh battery index
eff_t = zeros(1,length(swso)); %[~] efficiency
%surv = 1;
charging = false;
runtime = zeros(1,nvi); %[h], amount of time spent running
dies_vol = zeros(1,nvi);
buoy1 = 1; %buoy1 = 1 buoy1 in water, =0 buoy 2 in water
if buoy1 == 1
    S1(1) = Smax*1000; %assume battery begins fully charged
    S2(1) = 0.5*Smax*1000; %assume battery begins fully charged
else
    disp('Error - battery 1 should start at sea')
end
newbatt = 0; %number of new batteries that have to be bought
clear clean_ind batt_lft1
%clean_ind = zeros(length(swso),1);

%% run simulation
for t = 1:(length(time))
    if t < fbi1 + batt.bdi - 1 %less than first interval after fresh batt
        batt_L1(t) = 0;
    elseif rem(t,batt.bdi) == 0 %evaluate degradation on interval
        batt_L1(t:t+batt.bdi) = batDegModel(S1(fbi1:t)/(1000*Smax), ...
            batt.T,3600*(t-fbi1),batt.rf_os,ID);
            %3600*t should be 3600*(t-fbi1)
        if batt_L1(t) > batt.EoL %battery is past its life
            disp('error: battery died')
            if ~exist('batt_lft1','var')
                batt_lft1 = t*dt*(1/8760)*(12); %[mo] battery lifetime
            end
        end
    end
    if t < fbi2 + batt.bdi - 1 %less than first interval after fresh batt
        batt_L2(t) = 0;
    elseif rem(t,batt.bdi) == 0 %evaluate degradation on interval
        %disp('evaluating batt 2 degredation')
        if t < T/nvi %use only calendar aging for the first 2 years with no cycling
           [L_cal, d_cal] = Calendar_degradation(S2(fbi2:t)/(1000*Smax), ...
                batt.T,3600*(t-fbi2));
            batt_L2(t:t+batt.bdi) = L_cal;
            %Calendar_degradation(SoC, T, t)
        else
            batt_L2(t:t+batt.bdi) = batDegModel(S2(fbi2:t)/(1000*Smax), ...
                batt.T,3600*(t-fbi2),batt.rf_os,ID);
        end
        if batt_L2(t) > batt.EoL %battery is past its life
            disp('error: battery died')
            if ~exist('batt_lft2','var')
                batt_lft2 = t*dt*(1/8760)*(12); %[mo] battery lifetime
            end
        end
    end

    cf1 = batt_L1(t)*Smax*1000; %[Wh] capacity fading
    sd1 = S1(t)*(batt.sdr/100)*(1/(30*24))*dt; %[Wh] self discharge
    cf2 = batt_L2(t)*Smax*1000; %[Wh] capacity fading
    sd2 = S2(t)*(batt.sdr/100)*(1/(30*24))*dt; %[Wh] self discharge

    %Calculate state of charge for buoy on shore
    if buoy1 == 1 %buoy 1 = 1 means buoy 2 is on shore
        S2(t+1) = S2(t)-sd2;
        if S2(t+1) > (Smax*1000) - cf2
            S2(t+1) = (Smax*1000) - cf2;
        end
    else
        S1(t+1) = S1(t)-sd1;
        if S1(t+1) > (Smax*1000) - cf1
            S1(t+1) = (Smax*1000) - cf1;
        end
    end

    %Calc Renewable Power
    %find solar efficiency
    soil_eff = soil_eff - d_soil_eff;
    if soil_eff < 0 
        soil_eff = 0; %no negative efficiency
    end
    %soil_eff = (1-(atmo.soil/100)/8760*rem(t,inso.pvci*(365/12)*24));
    %soil_eff = (1-soil_eff)*rain(t) + soil_eff; %rainfall clean
    eff_t(t) = eff(t)*soil_eff*inso.eff;
    %find power from panel
    if kW_inso ~= 0
        if swso(t) > inso.rated*1000 %rated irradiance
            Pinso(t) = eff_t(t)/inso.eff*kW_inso*1000; %[W]
        else %sub rated irradiance
            Pinso(t) = eff_t(t)/ ...
                inso.eff*kW_inso*1000*(swso(t)/(inso.rated*1000)); %[W]
        end
    end
    %find power from wind turbine
    if kW_wind ~= 0 
        if wind(t) < turb.uci %below cut out
            Pwind(t) = 0; %[W]
        elseif turb.uci < wind(t) && wind(t) <= turb.ura %below rated
            Pwind(t) = kW_wind*1000*wind(t)^3/turb.ura^3; %[W]
        elseif turb.ura < wind(t) && wind(t) <= turb.uco %above rated
            Pwind(t) = kW_wind*1000; %[W]
        else %above cut out
            Pwind(t) = 0; %[W]
        end
    end
    %Pcurr,kW_curr,c_vel
    if kW_curr ~= 0 
        %Pcurr(t) = interp2(kW_testc,c_vel,Pcurr_m,kW_curr,c_speed(t));
        %Pcurr(t) = F_curr(kW_curr,c_speed(t));
        if c_speed(t) < cturb.uci %below cut out
            Pcurr(t) = 0; %[W]
        elseif cturb.uci < c_speed(t) && c_speed(t) <= cturb.ura %below rated
            Pcurr(t) = kW_curr.*1000.*c_speed(t).^3./cturb.ura.^3; %[W]
        elseif cturb.ura < c_speed(t) && c_speed(t) <= cturb.uco %above rated
            Pcurr(t) = kW_curr.*1000; %[W]
        else %above cut out
            Pcurr(t) = 0; %[W]
        end
    end

    Prenew(t) = Pwave(t) + Pwind(t) + Pinso(t) + Pcurr(t); %total renewable power
    %Prenew(t) = 0; %test Configuration - only solar
    %Calculate state of charge
    if ~charging %generator off
        if buoy1 == 1
            S1(t+1) = dt*(Prenew(t)-uc.draw(t)) + S1(t) - sd1;
            if t < length(time)
                if S1(t+1) < uc.draw(t+1)*dt && kW_dies ~= 0 %if not enough power to run for the next hour
                    charging = true; %turn diesel generator on for the next hour
                end
            end
            if S1(t+1) >= (Smax*1000) - cf1
                D(t) = S1(t+1) - ((Smax*1000) - cf1); %[Wh]
                S1(t+1) = (Smax*1000) - cf1;
            end
        else
            S2(t+1) = dt*(Prenew(t)-uc.draw(t)) + S2(t) - sd2;
            if t < length(time)
                if S2(t+1) < uc.draw(t+1)*dt && kW_dies ~= 0 %if not enough power to run for the next hour
                    charging = true; %turn diesel generator on for the next hour
                end
            end
            if S2(t+1) >= (Smax*1000) - cf2
                D(t) = S2(t+1) - ((Smax*1000) - cf2); %[Wh]
                S2(t+1) = (Smax*1000) - cf2;
            end
        end
    else %generator on
        Pdies(t) = kW_dies*1000;
        if buoy1 == 1
            S1(t+1) = dt*(Pdies(t)+Prenew(t) - uc.draw(t)) + S1(t) - sd1; %[Wh]
            if S1(t+1) >= (Smax*1000 - cf1)
                D(t) = S1(t+1) - (Smax*1000 - cf1); %[Wh]
                S1(t+1) = Smax*1000 - cf1;
                if opt.drun == 2
                    charging = false; %turn off gen when battery is full
                end
            end
        else
            S2(t+1) = dt*(Pdies(t)+Prenew(t) - uc.draw(t)) + S2(t) - sd2; %[Wh]
            if S2(t+1) >= (Smax*1000 - cf2)
                D(t) = S2(t+1) - (Smax*1000 - cf2); %[Wh]
                S2(t+1) = Smax*1000 - cf2;
                if opt.drun == 2
                    charging = false; %turn off gen when battery is full
                end
            end
        end
        if opt.drun == 1
            charging = false; %only run the generator for 1 hour at a time
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
    Ptot(t) = Prenew(t) + Pdies(t); %total power
    if buoy1 == 1
        if S1(t+1) <= Smax*batt.dmax*1000 %bottomed out
            S1(t+1) = dt*Ptot(t) + S1(t) - sd1;
            L(t) = 0; %load drops to zero
            F(t) = 1; %failure tracker
        end
    else
        if S2(t+1) <= Smax*batt.dmax*1000 %bottomed out
            S2(t+1) = dt*Ptot(t) + S2(t) - sd2;
            L(t) = 0; %load drops to zero
            F(t) = 1; %failure trackerEoL
        end
    end
    if t == T/nvi || t == 2*T/nvi %if maintenance interval
        if buoy1 == 1 %buoy 1 going to shore
            S2(t+1) = (Smax*1000) - cf2; %buoy 2 charged up
            %if (1-abs((S1(t)-Smax)*1000 - cf1)/(Smax*1000)) >= batt.EoL %need new batt1
            if batt_L1(t) > batt.EoL %need new batt1
                S1(t+1) = (Smax*1000)*0.5; %stored at half total capacity
                newbatt = newbatt + 1;
                fbi1 = t+1; %set new battery interval
            else
                S1(t+1) = ((Smax*1000)-cf1)*0.5; %stored at half capacity
            end
        else %buoy 2 going to shore
            S1(t+1) = (Smax*1000) - cf1; %buoy 1 charged up
            %if (1-abs((S2(t)-Smax)*1000 - cf2)/(Smax*1000)) >= batt.EoL %need new batt2
            if batt_L2(t) > batt.EoL %need new batt2
                S2(t+1) = Smax*1000*0.5; %stored at half capacity
                newbatt = newbatt + 1;
                fbi2 = t+1;
            else
                S2(t+1) = ((Smax*1000)-cf2)*0.5; %stored at half capacity
            end
        end
        buoy1 = 1 - buoy1; %switch buoy case
    end
end
    
% battery degradation model
if batt.lcm == 1 %bolun's model
%     [batt_L,batt_lft] =  irregularDegradation(S/(Smax*1000), ...
%         data.wave.time',uc.lifetime,batt); %retrospective modeling (old)
    if ~exist('batt_lft1','var') %battery never reached EoL
        batt_lft1 = batt.EoL/batt_L1(t)*t*12/(8760); %[mo]
    end
    if ~exist('batt_lft12','var') %battery never reached EoL
        batt_lft2 = batt.EoL/batt_L2(t)*t*12/(8760); %[mo]
    end
elseif batt.lcm == 2 %dyanmic (old) model
    opt.phi = Smax/(Smax - (min(S)/1000)); %extra depth
    batt_lft1 = batt.lc_nom*opt.phi^(batt.beta); %new lifetime
    batt_lft1(batt_lft1 > batt.lc_max) = batt.lc_max; %no larger than max
else %fixed (really old) model
    batt_lft1 = batt.lc_nom; %[m]
end

time2 = toc(t2);

nbr1 = ceil((12*uc.lifetime/batt_lft1-1)); %number of battery replacements
nbr2 = ceil((12*uc.lifetime/batt_lft2-1)); %number of battery replacements
nbr = nbr1 + nbr2;
%% compute number of vessel requirements
runtime_tot = sum(runtime);
nfr = ceil(runtime_tot*lph/dies.fmax-1); %number of fuel replacements
if uc.lifetime/nfr > dies.ftmax/12 %fuel will go bad
    nfr = ceil(12*uc.lifetime/dies.ftmax-1);
end
noc = ceil(runtime_tot/dies.oilint-1); %number of oil changes

%Changing to a constant of vessel/life * lifetime
if nvi < (nbr-newbatt)
    disp('Warning Battery will die')
elseif nvi < nfr
    disp('Warning Fuel will run out')
elseif nvi < noc
    disp('Warning oil will run out')
end

%% Calculate Mass of components (for weight approx optimization or for mooring model)
mass_solar = inso.wf*kW_inso/(inso.rated*inso.eff); %[kg] - 1 array
mass_solar_E = 0; %electrical - for 1 buoy
%mass_solar_S = 3.15*kW_inso/(inso.rated*inso.eff); %structural - for 1 buoy
mass_solar_S = 0;
mass_wind = kW_wind*turb.wf; %updated to include a tower

mass_curr = 0; %assume the turbine assembly is ballasted to neutral buoyancy (Brian)

mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of 1 generator [kg]
al_plate = 35.24; %[kg/m2] - 1.2" 6061 AL
gen_vol = econ.dies.volmult*polyval(opt.p_dev.d_size,kW_dies)^3;
mass_diesencl = 6*(gen_vol^2/3)*al_plate; %[kg]
dies_dens = 0.85; %[g/cm3] = [kg/L] from a chevron report
mass_fuel = runtime_tot*lph*dies_dens; %[kg]
if kW_dies == 0
    mass_diesencl = 0;
    mass_fuel = 0;
    mass_dies = 0;
end
mass_batt = Smax/(batt.se*batt.V/1000); %[kg] - 1 battery
batt_vol = polyval(opt.p_dev.b_size,Smax)*econ.batt.volmult;
batt_len = (batt_vol)^(1/3); %assume cube
mass_battencl = 6*(batt_len^2)*al_plate; %[kg]

% Pmtrl = econ.platform.wf* ...
%     (mass_dies + mass_diesencl + mass_fuel + mass_solar + mass_wind + mass_batt + ...
%     mass_battencl + mass_solar_E + mass_solar_S); % 1 platform material [kg]

mass_wec = kW_wave*895.78;
if kW_wave == 0
    mass_wec = 0;
end
comp_plat_mass = mass_solar + mass_solar_E + mass_solar_S + mass_wind + mass_dies + ...
        mass_diesencl + mass_fuel + mass_batt + mass_battencl + mass_wec; %total mass of componenets on the platform

%% Economic Optimization
if opt.tar ==3 %economic
    %costs - solar
    Mcost_inso = 2*econ.inso.module*kW_inso*econ.inso.marinization*econ.inso.pcm; %module
    Icost_inso = 2*econ.inso.installation*kW_inso*econ.inso.pcm; %installation
    Ecost_inso = 2*econ.inso.electrical*kW_inso*econ.inso.marinization ...
        *econ.inso.pcm; %electrical infrastructure
    Strcost_inso = 2*econ.inso.structural*kW_inso*econ.inso.marinization ...
        *econ.inso.pcm; %structural infrastructure
    %costs - wind
    kWcost_wind = 2*polyval(opt.p_dev.t,kW_wind)*econ.wind.marinization ...
        *econ.wind.tcm; %turbine
    if kW_wind == 0, kWcost_wind = 0; end
    if kW_wind ~=0 %can't calculate Icost with zero kWcost_wind
        Icost_wind = 2*(econ.wind.installed - 0.5*kWcost_wind/ ...
            (kW_wind*econ.wind.marinization*econ.wind.tcm))*kW_wind; %installation
        if Icost_wind < 0, Icost_wind = 0; end
    else
        Icost_wind = 0;
    end
    %costs - wave
    kWcost_wave = 2*econ.wave.costmult*polyval(opt.p_dev.t,kW_wave); %wec
    if kW_wave == 0, kWcost_wave = 0; end
    if kW_wave ~= 0 %no Icost if kW_wave = 0
        Icost_wave = 2*(econ.wind.installed - (0.5*kWcost_wave)/ ...
            (kW_wave*econ.wave.costmult))*kW_wave*econ.wave.costmult; %installation
        if Icost_wave < 0, Icost_wave = 0; end
    else
        Icost_wave = 0;
    end

    %costs - current
    kWcost_curr = 2*econ.curr.costmult*polyval(opt.p_dev.t,kW_curr); %current
    if kW_curr == 0, kWcost_curr = 0; end
    
    if kW_curr ~= 0
        Icost_curr = 2*(econ.wind.installed - (0.5*kWcost_curr)/ ...
                (kW_curr*econ.curr.costmult))*kW_curr*econ.curr.costmult; %installation
        if Icost_curr < 0, Icost_wave = 0; end
    else
        Icost_curr = 0;
    end
    
    %dies- costs
    al_cost = 103.92; %$/sqft for 1/2" 6061 AL
    al_cost = al_cost*0.0929; %convert to $/m2

    if kW_dies == 0 
        kWcost_dies = 0;
        genencl = 0;
        fuel = 0;
        mass_dies = 0;
        Icost_dies = 0;
    else
        kWcost_dies = 2*(polyval(opt.p_dev.d_cost,kW_dies)*econ.dies.gcm + econ.dies.autostart); %generator (with spares provisioning: 2)
        %genencl = 2*polyval(opt.p_dev.d_size,kW_dies)^3* ...
        %    (econ.dies.enclcost/econ.dies.enclcap); %generator enclosure
        genencl = 2*6*(gen_vol^2/3)*al_cost;
        fuel = runtime_tot*lph*econ.dies.fcost; %cost of consumed fuel
        mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of generator
        Icost_dies = 0.1*(kWcost_dies + genencl); %installation
    end
    
    if bc == 1 %lead acid
        if Smax < opt.p_dev.kWhmax %less than linear region
            Scost = polyval(opt.p_dev.b,Smax);
        else %in linear region
            Scost = polyval(opt.p_dev.b,opt.p_dev.kWhmax)* ...
                (Smax/opt.p_dev.kWhmax);
        end
    elseif bc == 2 %lithium phosphate
        Scost = (2)*batt.cost*Smax; %2x so there is a spare on sure
    end
    battencl = 2*6*(batt_len^2)*al_cost;
    Icost_batt = 0.1*(Scost+battencl); %installation cost (10% is a guess from Brian)

    %%Old mooring model
    if any(strcmp(data.title,"Argentine Basin"))
        Buoy_Mass = econ.platform.wf*comp_plat_mass; %Trent's assumption that the platform is 5x the on platform equipment mass
        if opt.pm == 3 %wave
            dp = width;
        elseif opt.pm == 2 %solar
            dp = getInsoDiameter(kW_inso,inso);
            if dp < 2, dp = 2; end
        elseif opt.pm == 1 %wind
            dp = 0.8;
        elseif opt.pm == 4 %diesel
            dp = polyval(opt.p_dev.d_size,kW);
            if dp < 2, dp = 2; end
        end
        Pmooring = interp2(econ.platform.diameter, ...
        econ.platform.depth, ...
        econ.platform.cost,dp,depth,'linear');
    else
        %%New Mooring Model
        tempPmooring = nan(3,1);
        tempmass = nan(3,1);
        for b = 1:3 %loop through buoy sizes
            tempPmooring(b) = 2*interp1(econ.platform.payloadmass(b,:),econ.platform.cost(b,:),comp_plat_mass,'linear'); %mooring cost - this should not extrapolate
            tempmass(b) = interp1(econ.platform.payloadmass(b,:),econ.platform.mass(b,:),comp_plat_mass,'linear');
        end
        [Pmooring,indMoor] = min(tempPmooring);
        Buoy_Mass = tempmass(indMoor);
        dp = nan; %platform diameter is not needed for Task 2 modeling
    end
    Pmtrl = 2*Buoy_Mass*econ.platform.steel;  %platform material cost - UPDATED BASED ON INFO FROM DEREK
    t_i = interp1(econ.platform.d_i,econ.platform.t_i,depth, ...
        'linear','extrap'); %installation time
    if depth < 500
        t_i = econ.platform.t_i (1); %make min time = 6 hours
    end
    Pinst = econ.vessel.speccost*(t_i/24); %platform instllation
    
    %triptime = dist*kts2mps(econ.vessel.speed)^(-1)*(1/86400); %[d]
    triptime = 0; %cost of trip attributed to instrumentation
    C_v = econ.vessel.speccost; %assuming that all trips require the spec vessel since everything is being replaced
    t_os = 2*t_i/24; %assuming recovery is similar to installation (have to recover and re-deploy)
    vesselcost = C_v*(nvi*(2*triptime + t_os)); %vessel cost
    %Maintenance Costs
    solarrepair = 1/4*0.5*(Strcost_inso + Icost_inso + Ecost_inso); %solar refurb repair
    if solarrepair < 0, solarrepair = 0; end
    
    turbrepair = 1/2*0.5*(kWcost_wind+Icost_wind)*(nvi); %turbine repair cost
    if turbrepair < 0, turbrepair = 0; end

    cturbrepair = 1/2*0.5*(kWcost_curr + Icost_curr)*(nvi); %current turbine refurb cost
    if cturbrepair <0, cturbrepair = 0; end

    wecrepair = 1/2*(0.5)*(kWcost_wave+Icost_wave)*(nvi); %wec repair cost
    if wecrepair < 0, wecrepair = 0; end %if nvi = 0, wec repair must be 0

    genrepair = 1/4*(0.5)*(kWcost_dies+genencl)*(nvi); %generator repair cost
    if genrepair < 0, genrepair = 0; end
    battreplace = Scost*newbatt/2; %number of battery replacements
    batteryrepair = 1/4*0.5*(Scost + Icost_batt)*nvi; %battery repair
    mooringrepair = 1/4*0.5*Pmooring; %eventually need a mooring refurb cost
    
    %Total Cost Calculations
    CapEx = Pmooring + Pinst + Pmtrl + battencl + Scost + Icost_batt + ...
        genencl + kWcost_dies + Icost_dies + Icost_inso + Mcost_inso + Ecost_inso + ...
        Strcost_inso + Icost_wind + kWcost_wind + Icost_wave + kWcost_wave + kWcost_curr + Icost_curr;
    OpEx = fuel + battreplace + batteryrepair + genrepair + vesselcost + turbrepair + cturbrepair + wecrepair + solarrepair + mooringrepair;
    cost = CapEx + OpEx;

%% Weight Approximation
elseif opt.tar == 1
    buoy_m = mass_solar + mass_solar_E + mass_solar_S + mass_wind + mass_curr + mass_dies + ...
        mass_diesencl + mass_fuel + mass_wec + mass_batt + mass_battencl;
    newbatt_m = newbatt * mass_batt; %[g]
    cost = buoy_m*2 + newbatt_m; %[kg] - 2 buoy weight + weight of extra batteries
    % if kW_wave < 0.2144 && kW_wave ~= 0
    %     cost = inf; %make cost infinity for bad wave power cases - since I'm only adjusting cost it shouldn't affect per_opt
    % end
%% Lowest gen capacity target
else 
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
if opt.tar == 1 || opt.tar ==2
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
end
pvci = nan; %This variable is not needed anymore
%determine if desired uptime was met
surv = sum(L == uc.draw)/(length(L));
% if sum(L == uc.draw)/(length(L)) < uc.uptime  %This is handled in the
% optimization code (ffa or optHybrid)
%     %surv = 0;
%     if opt.fmin
%         cost = inf;
%     end
% end

if isnan(Pmooring)
    disp('Ahh - the platform is too big')
    cost = inf;
end

end
