function [cost,surv,CapEx,OpEx,kWcost_dies, kWcost_wave,kWcost_wind,Mcost_inso,Ecost_inso,Icost_inso,...
    Strcost_inso, Icost_wave, Icost_wind, Scost,Pmtrl,Pinst,Pmooring, ...
    vesselcost,genrepair,turbrepair, wecrepair, battreplace,battencl,genencl,fuel, ...
    triptime,runtime,nvi,batt_L,batt_lft,nfr,noc,dp,S,Pdies,Pinso,Pwind,Pwave,Ptot,width,cw,D,L,eff_t,pvci] =  ...
    simHybrid(kW_dies, kW_inso, kW_wind, kW_wave,Smax,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb)

%% Created by Sarah Palmer Jan 2023 - started from Trent's OO-Tech code

ID = [kW_dies kW_inso kW_wind kW_wave Smax];
display('made it into SimHybrid')

%set capture width modifier
cw_mod = wave.cw_mod;
%set burn rate
lph = polyval(opt.p_dev.d_burn,kW_dies); %[l/h], burn rate
%extract data
time = data.met.time;
T = length(time);
dt = 24*(data.met.time(2) - data.met.time(1)); %time in hours
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
if wave.method == 1 %divide by B methodology       
    cwr_b = wave.cw_mod.*opt.wave.cwr_b_ts; %[m^-1] eta timeseries (cwr/b)    
    %find width through rated power conditions
    width = sqrt(1000*kW_wave*(1+wave.house)/(wave.eta_ct* ...
        cw_mod*opt.wave.cwr_b_ra* ...
        opt.wave.wavepower_ra)); %[m] physical width of wec
    cw = cwr_b.*width^2; %[m] capture width timeseries
elseif wave.method == 2 %3d interpolation methodology
    %extract data
    Hs = opt.wave.Hs; %Hs timeseries
    Tp = opt.wave.Tp; %Tp timeseries
    %find width through rated power conditions
    width = interp1(opt.wave.B_func(2,:),opt.wave.B_func(1,:),kW_wave); %[m], B
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
    S = zeros(1,length(time)); %battery level timeseries
    S(1) = Smax*1000; %assume battery begins fully charged
    Ptot = zeros(1,length(time)); %power produced timeseries
    Pdies = zeros(1,length(time));
    Pinso = zeros(1,length(time));
    Pwind = zeros(1,length(time));
    Prenew = zeros(1,length(time));
    D = zeros(1,length(time)); %power dumped timeseries
    L = ones(1,length(time))*uc.draw; %power put to sensing timeseries
    batt_L = zeros(1,T); %battery L (degradation) timeseries
    fbi = 1; %fresh battery index
    eff_t = zeros(1,length(swso)); %[~] efficiency
    surv = 1;
    charging = false;
    runtime = 0; %[h], amount of time spent running
    
    %set cleaning interval
    clear clean_ind batt_lft
    clean_ind = zeros(length(swso),1);
    if inso.cleanstrat == 3 || inso.cleanstrat == 4 && ...
            uc.SI > 6 %winter cleaning
        clean_ind(data.wint_clean_ind) = 1;
    else
        clean_ind(1:ceil((inso.pvci/12)*8760):end) = 1; %interval cleaning
    end

    %run simulation
    for t = 1:length(time)
        if t < fbi + batt.bdi - 1 %less than first interval after fresh batt
            batt_L(t) = 0;
        elseif rem(t,batt.bdi) == 0 %evaluate degradation on interval
            batt_L(t:t+batt.bdi) = batDegModel(S(fbi:t)/(1000*Smax), ...
                batt.T,3600*t,batt.rf_os,ID);
            if batt_L > batt.EoL %new battery
                fbi = t+1;
                S(t) = Smax*1000; 
                if ~exist('batt_lft','var')
                    batt_lft = t*dt*(1/8760)*(12); %[mo] battery lifetime
                end
            end
        end
        cf = batt_L(t)*Smax*1000; %[Wh] capacity fading
        sd = S(t)*(batt.sdr/100)*(1/(30*24))*dt; %[Wh] self discharge

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
        if swso(t) > inso.rated*1000 %rated irradiance
            Pinso(t) = eff_t(t)/inso.eff*kW_inso*1000; %[W]
        else %sub rated irradiance
            Pinso(t) = eff_t(t)/ ...
                inso.eff*kW_inso*1000*(swso(t)/(inso.rated*1000)); %[W]
        end
        %find power from wind turbine
        if wind(t) < turb.uci %below cut out
            Pwind(t) = 0; %[W]
        elseif turb.uci < wind(t) && wind(t) <= turb.ura %below rated
            Pwind(t) = kW_wind*1000*wind(t)^3/turb.ura^3; %[W]
        elseif turb.ura < wind(t) && wind(t) <= turb.uco %above rated
            Pwind(t) = kW_wind*1000; %[W]
        else %above cut out
            Pwind(t) = 0; %[W]
        end

        Prenew(t) = Pwave(t) + Pwind(t) + Pinso(t); %total renewable power

        %Calculate state of charge
        if ~charging %generator off
            S(t+1) = dt*(Prenew(t)-uc.draw) + S(t) - sd;
            %if S(t+1) < dies.genon*Smax*1000 %turn generator on
            if S(t+1) < uc.draw*dt %if not enough power to run for the next hour
                charging = true; %turn diesel generator on for the next hour
            end
        else %generator on
            Pdies(t) = kW_dies*1000;
            S(t+1) = dt*(Pdies(t)+Prenew(t) - uc.draw) + S(t) - sd; %[Wh]
            if S(t+1) >= (Smax*1000 - cf)
                D(t) = S(t+1) - (Smax*1000 - cf); %[Wh]
                S(t+1) = Smax*1000 - cf;
                charging = false;
            end
            int_time = T/(econ.vessel.int*uc.lifetime); %dividing hours by number of intervention giving number of hours between each intervention?
            runtime = zeros(1,econ.vessel.int*uc.lifetime);
            dies_vol = zeros(1,econ.vessel.int*uc.lifetime);
            int_id = ceil(t/int_time); %calculated what intervention number we are on
            runtime(int_id) = runtime(1,int_id) + dt; %[h] when diesel gen is on
            dies_vol(int_id) = dies_vol(1,int_id) + lph*dt; %[L] burn rate multiplied by hours running should give total consumption
            if runtime(1,int_id) >= 250 || dies_vol(1,int_id) > 800 %Brian set the maximums
                    surv = 0;
                    cost = inf;
                    disp('Option not viable due to fuel volume or runtime')
                    return
            end
        end
        Ptot(t) = Prenew(t) + Pdies(t); %total power
        if S(t+1) <= Smax*batt.dmax*1000 %bottomed out
            S(t+1) = dt*Ptot(t) + S(t) - sd;
            L(t) = 0; %load drops to zero
        end
    end
    
    % battery degradation model
    if batt.lcm == 1 %bolun's model
    %     [batt_L,batt_lft] =  irregularDegradation(S/(Smax*1000), ...
    %         data.wave.time',uc.lifetime,batt); %retrospective modeling (old)
        if ~exist('batt_lft','var') %battery never reached EoL
            batt_lft = batt.EoL/batt_L(t)*t*12/(8760); %[mo]
        end
    elseif batt.lcm == 2 %dyanmic (old) model
        opt.phi = Smax/(Smax - (min(S)/1000)); %extra depth
        batt_lft = batt.lc_nom*opt.phi^(batt.beta); %new lifetime
        batt_lft(batt_lft > batt.lc_max) = batt.lc_max; %no larger than max
    else %fixed (really old) model
        batt_lft = batt.lc_nom; %[m]
    end

    if inso.shootdebug
        disp(['pvci = ' num2str(inso.pvci)])
        disp(['battlc = ' num2str(batt_lft)])
        disp(['dm = ' num2str(dm)])
        pause
    end
    
    if uc.SI == 6 || abs(batt_lft - mult*inso.pvci) < tol || ...
        inso.cleanstrat == 3 || inso.cleanstrat == 4
        cont = 0;
    elseif batt_lft < mult*inso.pvci %cleaning interval > battery life
        inso.pvci = inso.pvci - dm;
        if inso.shootdebug
            disp('Decreasing pvci...')
            %pause
        end
        if ~over
            dm = dm/2;
            over = true;
        end
    elseif batt_lft > mult*inso.pvci %cleaning interval < battery life
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
nbr = ceil((12*uc.lifetime/batt_lft-1)); %number of battery replacements

%% compute number of vessel requirements
runtime_tot = sum(runtime);
nfr = ceil(runtime_tot*lph/dies.fmax-1); %number of fuel replacements
if uc.lifetime/nfr > dies.ftmax/12 %fuel will go bad
    nfr = ceil(12*uc.lifetime/dies.ftmax-1);
end
noc = ceil(runtime_tot/dies.oilint-1); %number of oil changes

%Changing to a constant of vessel/life * lifetime
nvi = econ.vessel.int*uc.lifetime;
if nvi < nbr
    disp('Warning Battery will die')
elseif nvi < nfr
    disp('Warning Fuel will run out')
elseif nvi < noc
    disp('Warning oil will run out')
end

%% Cost Calculation
%costs - solar
Mcost_inso = econ.inso.module*kW_inso*econ.inso.marinization*econ.inso.pcm; %module
Icost_inso = econ.inso.installation*kW_inso*econ.inso.pcm; %installation
Ecost_inso = econ.inso.electrical*kW_inso*econ.inso.marinization ...
    *econ.inso.pcm; %electrical infrastructure
Strcost_inso = econ.inso.structural*kW_inso*econ.inso.marinization ...
    *econ.inso.pcm; %structural infrastructure
%costs - wind
kWcost_wind = 2*polyval(opt.p_dev.t,kW_wind)*econ.wind.marinization ...
    *econ.wind.tcm; %turbine
Icost_wind = 2*(econ.wind.installed - 0.5*kWcost_wind/ ...
    (kW_wind*econ.wind.marinization*econ.wind.tcm))*kW_wind; %installation
if Icost_wind < 0, Icost_wind = 0; end
%costs - wave
kWcost_wave = 2*econ.wave.costmult*polyval(opt.p_dev.t,kW_wave); %wec
Icost_wave = 2*(econ.wind.installed - (0.5*kWcost_wave)/ ...
    (kW_wave*econ.wave.costmult))*kW_wave; %installation
if Icost_wave < 0, Icost_wave = 0; end

%dies- costs
kWcost_dies = polyval(opt.p_dev.d_cost,kW_dies)*2*econ.dies.gcm + ...
    econ.dies.autostart; %generator (with spares provisioning: 2)
genencl = polyval(opt.p_dev.d_size,kW_dies)^3* ...
    (econ.dies.enclcost/econ.dies.enclcap); %generator enclosure
fuel = runtime_tot*lph*econ.dies.fcost; %cost of consumed fuel
if bc == 1 %lead acid
    if Smax < opt.p_dev.kWhmax %less than linear region
        Scost = polyval(opt.p_dev.b,Smax);
    else %in linear region
        Scost = polyval(opt.p_dev.b,opt.p_dev.kWhmax)* ...
            (Smax/opt.p_dev.kWhmax);
    end
elseif bc == 2 %lithium phosphate
    Scost = batt.cost*Smax;
end
battencl = econ.batt.enclmult*Scost; %battery enclosure cost
Pmtrl = (1/1000)*econ.platform.wf*econ.platform.steel* ...
    (polyval(opt.p_dev.d_mass,kW_dies)+inso.wf*kW_inso/inso.rated+kW_wind*turb.wf); %platform material

%% Will need new way to calculate dp
dp = polyval(opt.p_dev.d_size,max(ID(1:end-1))); %Max ID 1:end-1 will give max of the rated generation
if dp < 2, dp = 2; end 
if dp > 8 %apply installtion time multiplier if big mooring - WRONG INSO.BOUNDARY needs to be updated
    inst_mult = interp1([8 econ.platform.boundary_di], ...
        [1 econ.platform.boundary_mf],dp);
else
    inst_mult = 1;
end
t_i = interp1(econ.platform.d_i,econ.platform.t_i,depth, ...
    'linear','extrap')*inst_mult; %installation time
Pinst = econ.vessel.speccost*(t_i/24); %platform instllation

%% Will need new mooring model
Pmooring = interp2(econ.platform.diameter, ...
    econ.platform.depth, ...
    econ.platform.cost,dp,depth,'linear'); %mooring cost

if uc.SI < 12 %short term instrumentation
    triptime = 0; %attributed to instrumentation
    t_os = econ.vessel.t_ms/24; %[d]
    C_v = econ.vessel.speccost;
else %long term instrumentation and infrastructure
    triptime = dist*kts2mps(econ.vessel.speed)^(-1)*(1/86400); %[d]
    t_os = econ.vessel.t_mosv/24; %[d]
    C_v = econ.vessel.osvcost;
end
vesselcost = C_v*(nvi*(2*triptime + t_os)); %vessel cost
turbrepair = 1/2*0.5*(kWcost_wind+Icost_wind)*(nvi); %turbine repair cost
if turbrepair < 0, turbrepair = 0; end
wecrepair = 1/2*(0.5)*(kWcost_wave+Icost_wave)*(nvi); %wec repair cost
if wecrepair < 0, wecrepair = 0; end %if nvi = 0, wec repair must be 0
genrepair = 1/2*(0.5)*(kWcost_dies+genencl)*(nvi); %generator repair cost
battreplace = Scost*nbr; %number of battery replacements
CapEx = Pmooring + Pinst + Pmtrl + battencl + Scost + ...
    genencl + kWcost_dies + Icost_inso + Mcost_inso + Ecost_inso + ...
    Strcost_inso + Icost_wind + kWcost_wind + Icost_wave + kWcost_wave;
OpEx = fuel + battreplace + genrepair + vesselcost + turbrepair + wecrepair;
cost = CapEx + OpEx;

%determine if desired uptime was met. if not, output infinite cost.
if sum(L == uc.draw)/(length(L)) < uc.uptime
    surv = 0;
    if opt.fmin
        cost = inf;
    end
end

end
