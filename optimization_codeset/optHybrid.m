function [output,opt] = optHybrid(opt,data,atmo,batt,econ,uc,bc,dies,inso,turb,wave)
%% OptHybrid is a combination of the optDies/Inso/Wind/Wave functions
%individual functions written by Trent Dillon
%Combination and modification by Sarah May Palmer

%% Set boundaries of Mesh
%set Smax mesh 
opt.Smax_1 = 1;
opt.Smax_n = opt.bf.N; %[kWh]
%opt.Smax_n = 220; %for test

%set kW mesh - Dies
opt.dies.kW_1 = 0; %min size zero for hybrid sim
opt.dies.kW_m = dies.kWmax; %max size
%opt.dies.kW_m = 0.01; %zero diesel power for test

%set kW mesh - Inso
opt.inso.kW_1 = 0.0; %min size zero for hybrid sim (used to be 0.5)
if econ.platform.boundary == 2
    opt.inso.kW_m = 9.0477; %corresponds to 8m
else
    if econ.platform.boundary_di == 12
        opt.inso.kW_m = 20.3575; %corresponds to 12m
    else
        error('update econ.platform.boundary_di')
    end
end

%set kW mesh - wind
opt.wind.kW_1 = 0.0; %min size zero for hybrid sim (used to be 0.1)
opt.wind.kW_m = opt.bf.M; %[kW] (up to 8 kW for wind)

%set kW mesh - wave
opt.wave.kW_1 = 0.0; %lower limit for wecsim is 0.2143 (used to be 0.215 - set to zero for hybrid)
if ~opt.highresobj
    opt.wave.kW_m = opt.bf.M; %[kW]
else
    opt.bf.loc_ind = find(contains(opt.locations,data.loc, ...
        'IgnoreCase',false));
    opt.wave.kW_m = opt.bf.M_hros(opt.bf.loc_ind);
    opt.Smax_n = opt.bf.N_hros(opt.bf.loc_ind);
end

%% Sensitivity and Modifiers
%set sensitivity modifiers to 1 if absent and to value if existing
if ~isfield(wave,'cw_mod')
    wave.cw_mod = 1; %capture width modifier
%set sensitivity modifiers to value if existing
if isfield(data,'depth_mod')
    data.depth = data.depth_mod; %depth modifier
end
if isfield(data,'dist_mod')
    data.dist = data.dist_mod; %dist to coast modifier
end
if isfield(econ.vessel,'tmt_enf') && ...
        (opt.sens || opt.tdsens || opt.senssm) && ...
        isequal(opt.tuned_parameter,'tmt')
    econ.vessel.t_mosv = econ.vessel.tmt_enf; %osv maintenance time
    econ.vessel.t_ms = econ.vessel.tmt_enf; %spec maintenance time
end

% Commenting out the lambda adjustments cause the maintenance shceudle
% strategy doesn't consider unplanned failures
% %set econ scenario
% switch econ.wind.scen
%     case 1 %optimistic durability
%         econ.wind.lambda = econ.wind.lowfail; %vessel interventions
%     case 2 %conservative
%         econ.wind.lambda = econ.wind.highfail; %vessel interventions
% end
% %if sensitivity analysis
% if isfield(econ.wind,'lambda_mod')
%     econ.wind.lambda = econ.wind.lambda_mod; %lamdba modifier
% end

%set econ scenario
switch econ.wave.scen
    case 1 %conservative
        econ.wave.costmult = econ.wave.costmult_con; %cost multiplier
        %econ.wave.lambda = econ.wave.highfail; %vessel interventions
    case 2 %optimistic cost
        econ.wave.costmult = econ.wave.costmult_opt; %cost multiplier
        %econ.wave.lambda = econ.wave.highfail; %vessel interventions
    case 3 %optimistic durability
        econ.wave.costmult = econ.wave.costmult_con; %cost multiplier
        %econ.wave.lambda = econ.wave.lowfail; %vessel interventions
end
%if sensitivity analysis
% if isfield(econ.wave,'lambda_mod')
%     econ.wave.lambda = econ.wave.lambda_mod; %lamdba modifier
% end
if isfield(econ.wave,'costmult_mod')
    econ.wave.costmult = econ.wave.costmult_mod; %cost modifier
end

%% Check Coarse Mesh
%check to make sure coarse mesh will work
opt.fmin = false;
check_s = 0;
while ~check_s
    [~,check_s] = simHybrid(opt.dies.kW_m/2,opt.inso.kW_m/2,opt.wind.kW_m/2,opt.wave.kW_m/2, opt.Smax_n/2,opt,data, ... 
        atmo,batt,econ,uc,bc,dies,inso,wave,turb);
%     simHybrid(kW_dies, kW_inso, kW_wind, kW_wave,Smax,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb)
    if ~check_s
        disp('Adjusting maximum values for coarse mesh')
        opt.Smax_n = 2*opt.Smax_n; 
        opt.wind.kW_m = 2*opt.wind.kW_m;
        opt.wave.kW_m = 2*opt.wave.kW_m;
        if opt.wave.kW_m > wave.kW_max %no larger than max wec-sim value
            opt.wave.kW_m = wave.wave.kW_max;
        end
    end
end

% %% Run Simulation - BRUTE FORCE
% %initialize inputs/outputs and set up for parallelization
% %IS THERE ANY REASON WE WOULD HAVE A NON-SYMMETRIC GRID? - like 10 squares
% %in dies and 20 in solar or something
% j = opt.bf.j;
% k = opt.bf.k;
% l = opt.bf.l;
% m = opt.bf.m;
% n = opt.bf.n;
% opt.dies.kW = linspace(opt.dies.kW_1,opt.dies.kW_m,j);              %[kW]
% opt.inso.kW = linspace(opt.inso.kW_1,opt.inso.kW_m,k);              %[kW]
% opt.wind.kW = linspace(opt.wind.kW_1,opt.wind.kW_m,l);              %[kW]
% opt.wave.kW = linspace(opt.wave.kW_1,opt.wave.kW_m,m);              %[kW]
% opt.Smax = linspace(opt.Smax_1,opt.Smax_n,n);    %[kWh]
% %display('Before grid')
% %display(strcat(string(length(opt.dies.kW)),string(length(opt.inso.kW)),string(length(opt.wind.kW)),string(length(opt.wave.kW))))
% %[K,S] = meshgrid(opt.kW,opt.Smax);
% [Kd,Ki,Kwi,Kwa,S] = ndgrid(opt.dies.kW,opt.inso.kW,opt.wind.kW, opt.wave.kW, opt.Smax);
% Kd = reshape(Kd,[j*k*l*m*n 1]);
% Ki = reshape(Ki,[j*k*l*m*n 1]);
% Kwi = reshape(Kwi,[j*k*l*m*n 1]);
% Kwa = reshape(Kwa,[j*k*l*m*n 1]);
% S = reshape(S,[j*k*l*m*n 1]);
% C_temp = zeros(j*k*l*m*n,1);
% S_temp = zeros(j*k*l*m*n,1);
% X = zeros(j*k*l*m*n,1);
% %set number of cores
% if isempty(gcp('nocreate')) %no parallel pool running
%     cores = feature('numcores'); %find number of cofes
%     if cores > 4 %only start if using HPC
%         parpool(cores);
%     end
% end
% %parallel computing via parfor
% tGrid = tic;
% fmin_temp = opt.fmin; %making a temp variable to make the parfor loop happier
% disp(['Populating grid values: j=' num2str(j) ', k=' num2str(k) ', l=' num2str(l) ', m=' num2str(m) ', n=' num2str(n)])
% % disp('after grid')
% % display(strcat(string(length(Kd)),string(length(Ki)),string(length(Kwi)),string(length(Kwa)),string(length(S))))
% parfor (i = 1:j*k*l*m*n,opt.bf.maxworkers)
%     display(strcat('for loop index:',string(i)))
%     %if fmin is suggesting a negative input (physically impossible), set
%     %S_temp and C_temp to failed values
%     if fmin_temp && S(i) < 0 || min([Kd(i),Ki(i),Kwi(i),Kwa(i)]) < 0
%         S_temp(i) = 0;
%         C_temp(i) = inf;
%     else
%         [C_temp(i),S_temp(i)] = ...
%             simHybrid(Kd(i), Ki(i), Kwi(i), Kwa(i),S(i),opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb);
%         %simHybrid(kW_dies, kW_inso, kW_wind, kW_wave,Smax,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb)
%     end
%     if S_temp(i) == 0 %update obj val X
%         X(i) = inf;
%     else
%         X(i) = C_temp(i);
%     end
% end
% %REMOVED TRANSPOSE BC HIGH ORDER MATRIX CAN'T USE TRANSPOSE - MIGHT MESS UP
% %INDEXING
% %output.cost = reshape(C_temp,[j k l m n]); %return cost to matrix and structure
% %output.surv = reshape(S_temp,[j k l m n]); %return surv to matrix and structure
% %X = reshape(X,[j k l m n]); %return objval X to matrix
% output.tGrid = toc(tGrid);
% 
% disp('Brute forcing global minimum...')
% %[I(1),I(2),I(3),I(4),I(5)] = find(X == min(X(:)),1,'first');
% %hybrid approach
% I_min = find(X==min(X)); %find index of minimum cost - doesn't work if there's multiple minimums
% output.min.kWd = Kd(I_min);
% output.min.kWi = Ki(I_min);
% output.min.kWwi = Kwi(I_min);
% output.min.kWwa = Kwa(I_min);
% output.min.Smax = S(I_min);

%% Run Simulation - Telescope
%set number of cores
if isempty(gcp('nocreate')) %no parallel pool running
    cores = feature('numcores'); %find number of cofes
    if cores > 4 %only start if using HPC
        parpool(cores);
    end
end
tel_i = 1; %starting telescope index
tol = false; %logical variable for if telescoping tolerance has been met
while tol == false && tel_i <=opt.tel_max
    %initialize inputs/outputs and set up for parallelization
    if tel_i == 1 %can't send output to function before output is defined
        [opt] = telescope_opt(opt,tel_i,[]); %get arrays of kW and S for this grid
    else
        [opt] = telescope_opt(opt,tel_i,output); %get arrays of kW and S for this grid
    end
    j = opt.bf.j;
    k = opt.bf.k;
    l = opt.bf.l;
    m = opt.bf.m;
    n = opt.bf.n;
    [Kd,Ki,Kwi,Kwa,S] = ndgrid(opt.dies.kW{tel_i},opt.inso.kW{tel_i},opt.wind.kW{tel_i}, opt.wave.kW{tel_i}, opt.Smax{tel_i});
    Kd = reshape(Kd,[j*k*l*m*n 1]);
    Ki = reshape(Ki,[j*k*l*m*n 1]);
    Kwi = reshape(Kwi,[j*k*l*m*n 1]);
    Kwa = reshape(Kwa,[j*k*l*m*n 1]);
    S = reshape(S,[j*k*l*m*n 1]);
    C_temp = zeros(j*k*l*m*n,1);
    S_temp = zeros(j*k*l*m*n,1);
    X = zeros(j*k*l*m*n,1);

    %parallel computing via parfor
    tGrid = tic;
    fmin_temp = opt.fmin; %making a temp variable to make the parfor loop happier
    disp(['Telescope Iteration: ',num2str(tel_i)])
    disp(['Populating grid values: j=' num2str(j) ', k=' num2str(k) ', l=' num2str(l) ', m=' num2str(m) ', n=' num2str(n)])
    parfor (i = 1:j*k*l*m*n,opt.bf.maxworkers)
        %if physically impossible, set S_temp and C_temp to failed values
        if fmin_temp && S(i) < 0 || min([Kd(i),Ki(i),Kwi(i),Kwa(i)]) < 0
            S_temp(i) = 0;
            C_temp(i) = inf;
        else %Call Hybrid Sim function
            [C_temp(i),S_temp(i)] = ...
                simHybrid(Kd(i), Ki(i), Kwi(i), Kwa(i),S(i),opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb);
        end
        if S_temp(i) == 0 %update obj val X
            X(tel_i,i) = inf;
        else
            X(tel_i,i) = C_temp(i);
        end
    end
    %REMOVED TRANSPOSE BC HIGH ORDER MATRIX CAN'T USE TRANSPOSE - MIGHT MESS UP
    %INDEXING
    output.tGrid = toc(tGrid);
    disp('Brute forcing current minimum...')
    %hybrid approach
    I_min(tel_i) = find(X(tel_i,:)==min(X(tel_i,:))); %find index of minimum cost - doesn't work if there's multiple minimums
    output.min.kWd{tel_i} = Kd(I_min(tel_i));
    output.min.kWi{tel_i} = Ki(I_min(tel_i));
    output.min.kWwi{tel_i} = Kwi(I_min(tel_i));
    output.min.kWwa{tel_i} = Kwa(I_min(tel_i));
    output.min.Smax{tel_i} = S(I_min(tel_i));
    if tel_i > 1
        if abs(X(tel_i,I_min(tel_i)) - X((tel_i-1),I_min(tel_i-1)))/X(tel_i,I_min(tel_i)) <= opt.ctol
            tol = true;
            if abs(output.min.kWd{tel_i} - output.min.kWd{tel_i-1})/output.min.kWd{tel_i} > opt.kwtol
                tol = false;
            elseif abs(output.min.kWi{tel_i} - output.min.kWi{tel_i-1})/output.min.kWi{tel_i} > opt.kwtol
                tol = false;
            elseif abs(output.min.kWwi{tel_i} - output.min.kWwi{tel_i-1})/output.min.kWwi{tel_i} > opt.kwtol
                tol = false;
            elseif abs(output.min.kWwa{tel_i} - output.min.kWwa{tel_i-1})/output.min.kWwa{tel_i} > opt.kwtol
                tol = false;
            elseif abs(output.min.Smax{tel_i} - output.min.Smax{tel_i-1})/output.min.Smax{tel_i} > opt.kwtol
                tol = false;
            end
            if tol == true    
                disp('Telescoping Tolerance Met')
            end
        end
    end
    tel_i = tel_i +1; %update telescoping iteration
end
%% GET OUTPUT VALUES
[output.min.cost,output.min.surv,output.min.CapEx,output.min.OpEx,...
    output.min.kWcost_dies,output.min.kWcost_wave, output.min.kWcost_wind, output.min.Mcost_inso,...
    output.min.Ecost_inso, output.min.Icost_inso, output.min.Strcost_inso,...
    output.min.Icost_wave, output.min.Icost_wind, output.min.Scost,output.min.Pmtrl,...
    output.min.Pinst,output.min.Pmooring, ...
    output.min.vesselcost,output.min.genrepair, ...
    output.min.turbrepair, output.min.wecrepair,...
    output.min.battreplace,output.min.battencl,output.min.genencl, ...
    output.min.fuel,output.min.triptime,output.min.runtime, ... 
    output.min.nvi,output.min.batt_L,output.min.batt_lft, ...
    output.min.nfr,output.min.noc, output.min.nbr...
    output.min.dp,output.min.S,output.min.Pdies, output.min.Pinso, output.min.Pwind, output.min.Pwave,...
    output.min.Ptot, output.min.width,output.min.cw,output.min.D,output.min.L, output.min.F, output.min.eff_t, output.min.pvci] ...
    = simHybrid(output.min.kWd{end},output.min.kWi{end}, output.min.kWwi{end}, output.min.kWwa{end}, output.min.Smax{end}, ...
    opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb);

output.min.batt_dyn_lc = batt.lc_nom*(output.min.Smax{end}/ ...
    (output.min.Smax{end} - (min(output.min.S)/1000)))^batt.beta;
output.min.CFd = mean(output.min.Pdies)/(1000*output.min.kWd{end});
output.min.CFi = mean(output.min.Pinso)/(1000*output.min.kWi{end});
output.min.CFwi = mean(output.min.Pwind)/(1000*output.min.kWwi{end});
output.min.CFwa = mean(output.min.Pwave)/(1000*output.min.kWwa{end});
output.min.rotor_h = turb.clearance + ... 
    sqrt(1000*2*output.min.kWwi{end}/(atmo.rho_a*pi*turb.ura^3)); %rotor height
output.min.cw_avg = mean(output.min.cw); %average capture width
output.min.cwr_avg = mean(output.min.cw_avg/output.min.width); %average cwr

%cycles per year
% output.min.cyc60 = countCycles(output.min.S,output.min.Smax,60)/ ...
%     (length(data.wave.time)/8760);
% output.min.cyc80 = countCycles(output.min.S,output.min.Smax,80)/ ...
%     (length(data.wave.time)/8760);
% output.min.cyc100 = countCycles(output.min.S,output.min.Smax,100)/ ...
%     (length(data.wave.time)/8760);

end

