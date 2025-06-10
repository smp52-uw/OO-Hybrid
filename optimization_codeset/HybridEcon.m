function [CapEx,OpEx,cost,costcomp,indPM,mass] = HybridEcon(data,econ,opt,inso,batt,turb,bc,kW_inso,kW_wind,kW_wave,kW_curr,kW_dies,Smax,runtime_tot,nvi,newbatt)
    %Hybrid Economic Model
    depth = data.depth;
    lph = polyval(opt.p_dev.d_burn,kW_dies); %[l/h], burn rate
    %% Calculate Mass of components (for weight approx optimization or for mooring model)
    mass_solar = inso.wf.*kW_inso./(inso.rated.*inso.eff); %[kg] - 1 array
    mass_solar_E = 0; %electrical - for 1 buoy
    %mass_solar_S = 3.15*kW_inso/(inso.rated*inso.eff); %structural - for 1 buoy
    mass_solar_S = 0;
    mass_wind = kW_wind.*turb.wf; %updated to include a tower
    
    mass_curr = 0; %assume the turbine assembly is ballasted to neutral buoyancy (Brian)
    
    mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of 1 generator [kg]
    al_plate = 35.24; %[kg/m2] - 1.2" 6061 AL
    gen_vol = econ.dies.volmult.*polyval(opt.p_dev.d_size,kW_dies).^3;
    mass_diesencl = 6.*(gen_vol.^2/3).*al_plate; %[kg]
    dies_dens = 0.85; %[g/cm3] = [kg/L] from a chevron report
    mass_fuel = runtime_tot.*lph*dies_dens; %[kg]
    if kW_dies == 0
        mass_diesencl = 0;
        mass_fuel = 0;
        mass_dies = 0;
    end
    mass_batt = Smax./(batt.se.*batt.V./1000); %[kg] - 1 battery
    batt_vol = polyval(opt.p_dev.b_size,Smax).*econ.batt.volmult;
    batt_len = (batt_vol).^(1/3); %assume cube
    mass_battencl = 6.*(batt_len.^2).*al_plate; %[kg]
    
    % Pmtrl = econ.platform.wf* ...
    %     (mass_dies + mass_diesencl + mass_fuel + mass_solar + mass_wind + mass_batt + ...
    %     mass_battencl + mass_solar_E + mass_solar_S); % 1 platform material [kg]
    
    mass_wec = kW_wave.*895.78; %based on RM3 weight with no float
    if kW_wave == 0
        mass_wec = 0;
    end
    comp_plat_mass = mass_solar + mass_solar_E + mass_solar_S + mass_wind + mass_dies + ...
            mass_diesencl + mass_fuel + mass_batt + mass_battencl + mass_wec; %total mass of componenets on the platform

    mass.total = comp_plat_mass;
    mass.solar= mass_solar;
    mass.solarE = mass_solar_E;
    mass.solarS = mass_solar_S;
    mass.wind = mass_wind;
    mass.dies = mass_dies;
    mass.diesencl = mass_diesencl;
    mass.fuel = mass_fuel;
    mass.batt = mass_batt;
    mass.battencl = mass_battencl;
    mass.wec = mass_wec;

    %% Costs
    %costs - solar
    Mcost_inso = 2.*econ.inso.module.*kW_inso.*econ.inso.marinization.*econ.inso.pcm; %module
    Icost_inso = 2.*econ.inso.installation.*kW_inso.*econ.inso.pcm; %installation
    Ecost_inso = 2.*econ.inso.electrical.*kW_inso.*econ.inso.marinization ...
        .*econ.inso.pcm; %electrical infrastructure
    Strcost_inso = 2.*econ.inso.structural.*kW_inso.*econ.inso.marinization ...
        .*econ.inso.pcm; %structural infrastructure
    %costs - wind
    kWcost_wind = 2.*polyval(opt.p_dev.t,kW_wind).*econ.wind.marinization ...
        .*econ.wind.tcm; %turbine
    if kW_wind == 0, kWcost_wind = 0; end
    if kW_wind ~=0 %can't calculate Icost with zero kWcost_wind
        Icost_wind = 2.*(econ.wind.installed - 0.5.*kWcost_wind./ ...
            (kW_wind.*econ.wind.marinization.*econ.wind.tcm)).*kW_wind; %installation
        if Icost_wind < 0, Icost_wind = 0; end
    else
        Icost_wind = 0;
    end
    %costs - wave
    kWcost_wave = 2.*econ.wave.costmult.*polyval(opt.p_dev.t,kW_wave); %wec
    if kW_wave == 0, kWcost_wave = 0; end
    if kW_wave ~= 0 %no Icost if kW_wave = 0
        Icost_wave = 2.*(econ.wind.installed - (0.5.*kWcost_wave)./ ...
            (kW_wave.*econ.wave.costmult)).*kW_wave.*econ.wave.costmult; %installation
        if Icost_wave < 0, Icost_wave = 0; end
    else
        Icost_wave = 0;
    end

    %costs - current
    kWcost_curr = 2.*econ.curr.costmult.*polyval(opt.p_dev.t,kW_curr); %current
    if kW_curr == 0, kWcost_curr = 0; end
    
    if kW_curr ~= 0
        Icost_curr = 2.*(econ.wind.installed - (0.5.*kWcost_curr)./ ...
                (kW_curr.*econ.curr.costmult)).*kW_curr.*econ.curr.costmult; %installation
        if Icost_curr < 0, Icost_wave = 0; end
    else
        Icost_curr = 0;
    end
    
    %dies- costs
    al_cost = 103.92; %$/sqft for 1/2" 6061 AL
    al_cost = al_cost.*0.0929; %convert to $/m2

    if kW_dies == 0 
        kWcost_dies = 0;
        genencl = 0;
        fuel = 0;
        mass_dies = 0;
        Icost_dies = 0;
    else
        kWcost_dies = 2.*(polyval(opt.p_dev.d_cost,kW_dies).*econ.dies.gcm + econ.dies.autostart); %generator (with spares provisioning: 2)
        %genencl = 2*polyval(opt.p_dev.d_size,kW_dies)^3* ...
        %    (econ.dies.enclcost/econ.dies.enclcap); %generator enclosure
        genencl = 2.*6.*(gen_vol.^2/3).*al_cost;
        fuel = runtime_tot.*lph.*econ.dies.fcost; %cost of consumed fuel
        mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of generator
        Icost_dies = 0.1.*(kWcost_dies + genencl); %installation
    end
    
    if bc == 1 %lead acid
        if Smax < opt.p_dev.kWhmax %less than linear region
            Scost = polyval(opt.p_dev.b,Smax);
        else %in linear region
            Scost = polyval(opt.p_dev.b,opt.p_dev.kWhmax).* ...
                (Smax./opt.p_dev.kWhmax);
        end
    elseif bc == 2 %lithium phosphate
        Scost = (2).*batt.cost.*Smax; %2x so there is a spare on sure
    end
    battencl = 2.*6.*(batt_len.^2).*al_cost;
    Icost_batt = 0.1.*(Scost+battencl); %installation cost (10% is a guess from Brian)

    %%Old mooring model
    if any(strcmp(data.title,"Argentine Basin"))
        Buoy_Mass = econ.platform.wf.*comp_plat_mass; %Trent's assumption that the platform is 5x the on platform equipment mass
        if opt.pm == 3 %wave
            dp = width;
        elseif opt.pm == 2 %solar
            dp = getInsoDiameter(kW_inso,inso);
            if dp < 2, dp = 2; end
        elseif opt.pm == 1 %wind
            dp = 0.8;
        elseif opt.pm == 4 %diesel
            dp = polyval(opt.p_dev.d_size,kW_dies);
            if dp < 2, dp = 2; end
        end
        Pmooring = interp2(econ.platform.diameter, ...
        econ.platform.depth, ...
        econ.platform.cost,dp,depth,'linear');
    else
        %%New Mooring Model
        tempPmooring = nan(length(kW_dies),3);
        tempmass = nan(length(kW_dies),3);
        dp = nan; %platform diameter is not needed for Task 2 modeling
        for b = 1:3 %loop through buoy sizes
            tempPmooring(:,b) = 2.*interp1(econ.platform.payloadmass(b,:),econ.platform.cost(b,:),comp_plat_mass,'linear'); %mooring cost - this should not extrapolate
            tempmass(:,b) = interp1(econ.platform.payloadmass(b,:),econ.platform.mass(b,:),comp_plat_mass,'linear');
            tempPmtrl(:,b) = 2.*(tempmass(:,b)-comp_plat_mass).*econ.platform.steel;
            tempMTRLMOOR(:,b) = tempPmtrl(:,b) + tempPmooring(:,b);
        end
        [Pmooring,indMoor] = min(tempPmooring,[],2); %find min mooring cost and index of that min
        [~,indPM] = min(tempMTRLMOOR,[],2); %find index of lowers material + mooring cost
        % for m = 1:length(indMoor)
        %     Buoy_Mass(m,1) = tempmass(m,indMoor(m));
        % end
        for m = 1:length(indPM)
            Buoy_Mass(m,1) = tempmass(m,indPM(m)) - comp_plat_mass(m);
            Pmooring(m,1) = tempPmooring(m,indPM(m));
        end
    end
    Pmtrl = 2.*Buoy_Mass.*econ.platform.steel;  %platform material cost - UPDATED BASED ON INFO FROM DEREK
    t_i = interp1(econ.platform.d_i,econ.platform.t_i,depth, ...
        'linear','extrap'); %installation time
    if depth < 500
        t_i = econ.platform.t_i (1); %make min time = 6 hours
    end
    Pinst = econ.vessel.speccost.*(t_i./24); %platform instllation
    
    %triptime = dist*kts2mps(econ.vessel.speed)^(-1)*(1/86400); %[d]
    triptime = 0; %cost of trip attributed to instrumentation
    C_v = econ.vessel.speccost; %assuming that all trips require the spec vessel since everything is being replaced
    t_os = 2.*t_i./24; %assuming recovery is similar to installation (have to recover and re-deploy)
    vesselcost = C_v.*(nvi.*(2.*triptime + t_os)); %vessel cost
    %Maintenance Costs
    solarrepair = 1/4.*0.5.*(Strcost_inso + Icost_inso + Ecost_inso); %solar refurb repair
    if solarrepair < 0, solarrepair = 0; end
    
    turbrepair = 1/2.*0.5.*(kWcost_wind+Icost_wind).*(nvi); %turbine repair cost
    if turbrepair < 0, turbrepair = 0; end

    cturbrepair = 1/2.*0.5.*(kWcost_curr + Icost_curr).*(nvi); %current turbine refurb cost
    if cturbrepair <0, cturbrepair = 0; end

    wecrepair = 1/2.*(0.5).*(kWcost_wave+Icost_wave).*(nvi); %wec repair cost
    if wecrepair < 0, wecrepair = 0; end %if nvi = 0, wec repair must be 0

    genrepair = 1/4.*(0.5).*(kWcost_dies+genencl).*(nvi); %generator repair cost
    if genrepair < 0, genrepair = 0; end
    battreplace = Scost.*newbatt./2; %number of battery replacements
    batteryrepair = 1/4.*0.5.*(Scost + Icost_batt).*nvi; %battery repair
    mooringrepair = 1/4.*0.5.*Pmooring; %eventually need a mooring refurb cost
    
    %Total Cost Calculations
    CapEx = Pmooring + Pinst + Pmtrl + battencl + Scost + Icost_batt + ...
        genencl + kWcost_dies + Icost_dies + Icost_inso + Mcost_inso + Ecost_inso + ...
        Strcost_inso + Icost_wind + kWcost_wind + Icost_wave + kWcost_wave + kWcost_curr + Icost_curr;
    OpEx = fuel + battreplace + batteryrepair + genrepair + vesselcost + turbrepair + cturbrepair + wecrepair + solarrepair + mooringrepair;
    cost = CapEx + OpEx;


    %package cost info
    costcomp.Pmooring = Pmooring;
    costcomp.Pinst = Pinst;
    costcomp.Pmtrl = Pmtrl;
    costcomp.battencl = battencl;
    costcomp.Scost = Scost;
    costcomp.Icost_batt = Icost_batt;
    costcomp.genencl = genencl;
    costcomp.kWcost_dies = kWcost_dies;
    costcomp.Icost_dies = Icost_dies;
    costcomp.Icost_inso = Icost_inso;
    costcomp.Mcost_inso = Mcost_inso;
    costcomp.Ecost_inso = Ecost_inso;
    costcomp.Strcost_inso = Strcost_inso;
    costcomp.Icost_wind = Icost_wind;
    costcomp.kWcost_wind = kWcost_wind;
    costcomp.Icost_wave = Icost_wave;
    costcomp.kWcost_wave = kWcost_wave;
    costcomp.kWcost_curr = kWcost_curr;
    costcomp.Icost_curr = Icost_curr;
    costcomp.fuel = fuel;
    costcomp.battreplace = battreplace;
    costcomp.batteryrepair = batteryrepair;
    costcomp.genrepair = genrepair;
    costcomp.vesselcost = vesselcost;
    costcomp.turbrepair = turbrepair;
    costcomp.cturbrepair = cturbrepair;
    costcomp.wecrepair = wecrepair;
    costcomp.solarrepair = solarrepair;
    costcomp.mooringrepair = mooringrepair;
end