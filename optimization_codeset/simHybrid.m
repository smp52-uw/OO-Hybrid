function [cost,surv,CapEx,OpEx,kWcost_dies, kWcost_wave,kWcost_wind,Mcost_inso,Ecost_inso,Icost_inso,...
    Strcost_inso, Icost_wave, Icost_wind, Scost,Pmtrl,Pinst,Pmooring, ...
    vesselcost,genrepair,turbrepair, wecrepair, battreplace,battencl,genencl,fuel, ...
    triptime,runtime,nvi,batt_L1,batt_L2, batt_lft1,batt_lft2, nfr,noc,nbr,dp,S1,S2,Pdies,Pinso,Pwind,Pwave,Ptot,width,cw,D,L,F,eff_t,pvci] =  ...
    simHybrid(kW_dies, kW_inso, kW_wind, kW_wave,Smax,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb)

%% Created by Sarah Palmer Jan 2023 - started from Trent's OO-Tech code

ID = [kW_dies kW_inso kW_wind kW_wave Smax];
disp('made it into SimHybrid')

%if physically impossible, set S_temp and C_temp to failed values
if opt.fmin && Smax < 0 || min([kW_dies,kW_inso,kW_wind,kW_wave]) < 0
    surv = 0;
    cost = inf;
    return
end

%set capture width modifier
cw_mod = wave.cw_mod;
%Number of vessel intervention
nvi = (uc.lifetime*12)/uc.SI;
%set burn rate
lph = polyval(opt.p_dev.d_burn,kW_dies); %[l/h], burn rate
if kW_dies == 0  %Set burn rate = 0 for no diesel generator (in case the polyfit doesn't cross (0,0))
    lph = 0; 
end
%extract data
time = datenum(data.met.time);
T = length(time);
dt = 24*(time(2) - time(1)); %time in hours
dist = data.dist; %[m] distance to shore
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
if wave.method == 1 %divide by B methodology - OUTDATED I THINK      
    disp('ERROR- method set to old divide by B method')
elseif wave.method == 2 %3d interpolation methodology
    %extract data
    Hs = opt.wave.Hs; %Hs timeseries
    Tp = opt.wave.Tp; %Tp timeseries
    %find width through rated power conditions
    width = interp1(opt.wave.B_func(2,:),opt.wave.B_func(1,:),kW_wave); %[m], B
    if kW_wave == 0 %Set width = 0 for no wave gen (in case the polyfit doesn't cross (0,0))
        width = 0;
    end
    cw = width.*opt.wave.F(Tp,Hs,width*ones(length(Tp),1)); %[m] cw ts
end

%compute wave power timeseries
Pwave = wave.eta_ct*cw.*wavepower - kW_wave*wave.house; %[kW] 
Pwave(Pwave<0) = 0; %no negative power
Pwave(Pwave>kW_wave) = kW_wave; %no larger than rated power
Pwave = Pwave*1000; %convert to watts

% set the cleaning interval based on the use case service interval
if uc.SI > 6
    %if service interval is long-term, guess the panels will be cleaned
    %once every thirty months, run shooting scheme
    inso.pvci = 10; %[months] initial guess
    over = true; %over/under indicator
    dm = 10; %change in pvci
    tol = inso.shoottol; %tolerance
    if inso.cleanstrat == 2
        mult = 2;
    else
        mult = 1;
    end
else
    %if service interval is short-term, assume the panels will be cleaned
    %every six months
    inso.pvci = 6;
end

%set panel degradation
eff = (1-((inso.deg/100)/8760)*(1:1:length(swso)));
%rain = repmat(linspace(0.5,0,24*30),[1,ceil(length(swso)/(24*30))]);
if econ.inso.scen == 1 %automated cleaning
    d_soil_eff = 0;
    econ.inso.marinization = econ.inso.marinization*2;
elseif econ.inso.scen == 2 %human cleaning
    d_soil_eff = (atmo.soil/100)/8760; %change in soil deg per hour
end
soil_eff = 1; %starting soil efficiency

% if inso.debug
%     disp([num2str(kW) ' kW'])
%     disp([num2str(Smax) ' kWh'])
%     pause
% end

cont = 1;
%t1 = tic;
t2 = tic;

while cont
    %initialize diagnostic variables
    S1 = zeros(1,length(time)); %battery level timeseries
    S2 = zeros(1,length(time)); %battery level timeseries
    S1(1) = Smax*1000; %assume battery begins fully charged
    S2(1) = Smax*1000; %assume battery begins fully charged
    Ptot = zeros(1,length(time)); %power produced timeseries
    Pdies = zeros(1,length(time));
    Pinso = zeros(1,length(time));
    Pwind = zeros(1,length(time));
    Prenew = zeros(1,length(time));
    D = zeros(1,length(time)); %power dumped timeseries
    L = uc.draw; %power to load time series
    F = zeros(1,length(time)); %failure series
    batt_L1 = zeros(1,length(time)); %battery L (degradation) timeseries
    fbi1 = 1; %fresh battery index
    batt_L2 = zeros(1,length(time)); %battery L (degradation) timeseries
    fbi2 = 1; %fresh battery index
    eff_t = zeros(1,length(swso)); %[~] efficiency
    surv = 1;
    charging = false;
    runtime = 0; %[h], amount of time spent running
    buoy1 = 1; %buoy1 = 1 buoy1 in water, =0 buoy 2 in water
    newbatt = 0; %number of new batteries that have to be bought
    %set cleaning interval
    clear clean_ind batt_lft1
    clean_ind = zeros(length(swso),1);
    if inso.cleanstrat == 3 || inso.cleanstrat == 4 && ...
            uc.SI > 6 %winter cleaning
        clean_ind(data.wint_clean_ind) = 1;
    else
        clean_ind(1:ceil((inso.pvci/12)*8760):end) = 1; %interval cleaning
    end

    %% run simulation
    for t = 1:length(time)
        if t < fbi1 + batt.bdi - 1 %less than first interval after fresh batt
            batt_L1(t) = 0;
        elseif rem(t,batt.bdi) == 0 %evaluate degradation on interval
            batt_L1(t:t+batt.bdi) = batDegModel(S1(fbi1:t)/(1000*Smax), ...
                batt.T,3600*(t-fbi1+1),batt.rf_os,ID);
                %3600*t should be 3600*(t-fbi1) but with the var that
                %changes with new battery
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
                    batt.T,3600*(t-fbi2+1));
                batt_L2(t:t+batt.bdi) = L_cal; %pretty sure this should be L not d
                %Calendar_degradation(SoC, T, t)
            else
                batt_L2(t:t+batt.bdi) = batDegModel(S2(fbi2:t)/(1000*Smax), ...
                    batt.T,3600*(t-fbi2+1),batt.rf_os,ID);
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
        if clean_ind(t) == 1
            soil_eff = 1; %panels cleaned
        end
        %soil_eff = (1-(atmo.soil/100)/8760*rem(t,inso.pvci*(365/12)*24));
        %soil_eff = (1-soil_eff)*rain(t) + soil_eff; %rainfall clean
        eff_t(t) = eff(t)*soil_eff*inso.eff;
        %find power from panel
        if kW_inso ~= 0 %removing the possibility of divide by zero error
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
        Prenew(t) = Pwave(t) + Pwind(t) + Pinso(t); %total renewable power
        %Prenew(t) = Pinso(t); %test Configuration - only solar
        %Calculate state of charge
        if ~charging %generator off
            if buoy1 == 1
                S1(t+1) = dt*(Prenew(t)-uc.draw(t)) + S1(t) - sd1;
                %if S(t+1) < dies.genon*Smax*1000 %turn generator on
                if S1(t+1) < uc.draw(t)*dt && kW_dies ~= 0 %if not enough power to run for the next hour
                    charging = true; %turn diesel generator on for the next hour
                elseif S1(t+1) >= (Smax*1000) - cf1
                    D(t) = S1(t+1) - ((Smax*1000) - cf1); %[Wh]
                    S1(t+1) = (Smax*1000) - cf1;
                end
            else
                S2(t+1) = dt*(Prenew(t)-uc.draw(t)) + S2(t) - sd2;
                %if S(t+1) < dies.genon*Smax*1000 %turn generator on
                if S2(t+1) < uc.draw(t)*dt && kW_dies ~= 0 %if not enough power to run for the next hour
                    charging = true; %turn diesel generator on for the next hour
                elseif S2(t+1) >= (Smax*1000) - cf2
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
                    charging = false;
                end
            else
                S2(t+1) = dt*(Pdies(t)+Prenew(t) - uc.draw(t)) + S2(t) - sd2; %[Wh]
                if S2(t+1) >= (Smax*1000 - cf2)
                    D(t) = S2(t+1) - (Smax*1000 - cf2); %[Wh]
                    S2(t+1) = Smax*1000 - cf2;
                    charging = false;
                end
            end
            %nvi = uc.lifetime*(uc.SI/12);
            int_time = T/nvi; %dividing hours by number of intervention giving number of hours between each intervention?
            runtime = zeros(1,nvi);
            dies_vol = zeros(1,nvi);
            int_id = ceil(t/int_time); %calculated what intervention number we are on
            runtime(int_id) = runtime(1,int_id) + dt; %[h] when diesel gen is on
            dies_vol(int_id) = dies_vol(1,int_id) + lph*dt; %[L] burn rate multiplied by hours running should give total consumption
            if runtime(1,int_id) >= 250 || dies_vol(1,int_id) > 800 %Brian set the maximums
                    surv = 0;
                    cost = inf;
                    disp('Option not viable due to fuel volume battor runtime')
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
                if (1-abs((S1(t)-Smax)*1000 - cf1)/(Smax*1000)) >= batt.EoL
                    S1(t+1) = (Smax*1000)*0.5; %stored at half total capacity
                    newbatt = newbatt + 1;
                    fbi1 = t+1; %set new battery interval
                else
                    S1(t+1) = ((Smax*1000)-cf1)*0.5; %stored at half capacity
                end
            else %buoy 2 going to shore
                %disp('buoy switch')
                %disp(cf1)
                S1(t+1) = (Smax*1000) - cf1; %buoy 1 charged up
                if (1-abs((S2(t)-Smax)*1000 - cf2)/(Smax*1000)) >= batt.EoL
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

    if inso.shootdebug
        disp(['pvci = ' num2str(inso.pvci)])
        disp(['battlc = ' num2str(batt_lft1)])
        disp(['dm = ' num2str(dm)])
        pause
    end
    
    if uc.SI == 6 || abs(batt_lft1 - mult*inso.pvci) < tol || ...
        inso.cleanstrat == 3 || inso.cleanstrat == 4
        cont = 0;
    elseif batt_lft1 < mult*inso.pvci %cleaning interval > battery life
        inso.pvci = inso.pvci - dm;
        if inso.shootdebug
            disp('Decreasing pvci...')
            %pause
        end
        if ~over
            dm = dm/2;
            over = true;
        end
    elseif batt_lft1 > mult*inso.pvci %cleaning interval < battery life
        %if batt lifetime is too long for cleaning to occur then...
        if inso.pvci > inso.cleanlim
            inso.pvci = inso.cleanlim; %set cleaning interval to limit
            cont = 0;
        else
            inso.pvci = inso.pvci + dm;
            if inso.shootdebug
                disp('Increasing pvci...')
                %pause
            end
            if over
                dm = dm/2;
                over = false;
            end
        end
    end
    
    %adjust if time is running too long
    %time1 = toc(t1);
    time2 = toc(t2);
    %     if time1 > 1
    %         dm = 20; %reset dm if it's taking too long
    %         disp('Resetting dm')
    %         t1 = tic; %reset timer
    %         %I don't think this improves convergence...
    %     end
    if time2 > 10 && ~inso.shootdebug
        error([num2str(kW_inso) ' kW and ' num2str(Smax) ...
            ' kWh do not converge'])
    end
end
pvci = inso.pvci; %pv cleaning interval
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
%nvi = uc.lifetime*(uc.SI/12);
if nvi < nbr
    disp('Warning Battery will die')
elseif nvi < nfr
    disp('Warning Fuel will run out')
elseif nvi < noc
    disp('Warning oil will run out')
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
mass_solar = inso.wf*kW_inso/(inso.rated*inso.eff); %[kg] - 1 array
mass_solar_E = 0; %electrical - for 1 buoy
mass_solar_S = 0; %structural - for 1 buoy
mass_wind = kW_wind*turb.wf; %[kg] - 1 turbine

mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of 1 generator [kg]
steel_plate = 98.8; %[kg/m2]
mass_diesencl = 6*(polyval(opt.p_dev.d_size,kW_dies)^2)*steel_plate; %[kg]
dies_dens = 0.85; %[g/cm3] = [kg/L] from a chevron report
mass_fuel = runtime_tot*lph*dies_dens; %[kg]

mass_batt = Smax/(batt.se*batt.V/1000); %[kg] - 1 battery
old assumption is that 1 encl $ = batt $
mass_battencl = 6*(polyval(opt.p_dev.b_size,Smax)^2)*steel_plate; %[kg]
Pmtrl = econ.platform.wf* ...
    (mass_dies + mass_diesencl + mass_fuel + mass_solar + mass_wind + mass_batt + ...
    mass_battencl + mass_solar_E + mass_solar_S); % 1 platform material [kg]
mass_wec = Pmtrl * 0.699; %[kg] based on the percent of the RM3 that isn't float mass

buoy_m = mass_solar + mass_solar_E + mass_solar_S + mass_wind + mass_dies + ...
    mass_diesencl + mass_fuel + mass_wec + mass_batt + Pmtrl;
newbatt_m = newbatt * mass_batt; %[g]
cost = buoy_m*2 + newbatt_m; %[kg] - 2 buoy weight + weight of extra batteries

%determine if desired uptime was met
surv = sum(L == uc.draw)/(length(L));
% if sum(L == uc.draw)/(length(L)) < uc.uptime
%     %surv = 0;
%     if opt.fmin
%         cost = inf;
%     end
% end

end
