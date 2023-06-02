function [output,opt] = optHybrid(opt,data,atmo,batt,econ,uc,bc,dies,inso,turb,cturb, wave)
%% OptHybrid is a combination of the optDies/Inso/Wind/Wave functions
%individual functions written by Trent Dillon
%Combination and modification by Sarah May Palmer

%% Set boundaries of Mesh
%set Smax mesh 
opt.Smax_1 = 1;
opt.Smax_n = opt.bf.N; %[kWh]
%opt.Smax_n = 19; %for test

%set kW mesh - Dies
opt.dies.kW_1 = 0; %min size zero for hybrid sim
%opt.dies.kW_m = dies.kWmax; %max size
opt.dies.kW_m = opt.bf.M; %zero diesel power for test

%set kW mesh - Inso
opt.inso.kW_1 = 0.0; %min size zero for hybrid sim (used to be 0.5)
if econ.platform.boundary == 2
    %opt.inso.kW_m = 9.0477; %corresponds to 8m
    %opt.inso.kW_m = 8; %corresponds to 8m
    opt.inso.kW_m = opt.bf.M; %corresponds to 8m
else
    if econ.platform.boundary_di == 12
        disp('error do not use larger size')
        %opt.inso.kW_m = 20.3575; %corresponds to 12m
        opt.inso.kW_m = 8; %corresponds to 8m
    else
        error('update econ.platform.boundary_di')
    end
end

%set kW mesh - wind
opt.wind.kW_1 = 0.0; %min size zero for hybrid sim (used to be 0.1)
opt.wind.kW_m = opt.bf.M; %[kW] (up to 8 kW for wind)
%opt.wind.kW_m = 3;
%set kW mesh - current
opt.curr.kW_1 = 0.0; %min size zero for hybrid sim (used to be 0.1)
opt.curr.kW_m = opt.bf.M; %[kW] (up to 8 kW for wind)
%opt.curr.kW_m = 3;

%set kW mesh - wave
opt.wave.kW_1 = 0.0; %lower limit for wecsim is 0.2143 (used to be 0.215 - set to zero for hybrid)
if ~opt.highresobj
    opt.wave.kW_m = opt.bf.M; %[kW]
    %opt.wave.kW_m = 3; %[kW]
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
end
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
    [~,check_s] = simHybrid(opt.dies.kW_m/2,opt.inso.kW_m/2,opt.wind.kW_m/2,opt.wave.kW_m/2, opt.curr.kW_m/2, opt.Smax_n/2,opt,data, ... 
        atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);
%     simHybrid(kW_dies, kW_inso, kW_wind, kW_wave,Smax,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb)\
    %checks = surv -> if surv<uc.uptime
    if check_s < uc.uptime
        disp('Adjusting maximum values for coarse mesh')
        opt.Smax_n = 2*opt.Smax_n; 
        opt.wind.kW_m = 2*opt.wind.kW_m;
        opt.curr.kW_m = 2*opt.curr.kW_m;
        opt.wave.kW_m = 2*opt.wave.kW_m;
        if opt.wave.kW_m > wave.kW_max %no larger than max wec-sim value
            opt.wave.kW_m = wave.kW_max;
        end
    end
end


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
    if all(opt.alg == 'tel') || all(opt.alg == 'to2')
        j = opt.bf.j;
        k = opt.bf.k;
        l = opt.bf.l;
        m = opt.bf.m;
        n = opt.bf.n;
        o = opt.bf.o;
        if tel_i == 1 %can't send output to function before output is defined
            if all(opt.alg == 'tel')
                [opt] = telescope_opt(opt,tel_i,[]); %get arrays of kW and S for this grid
            else
                [opt] = telescope2_opt(opt,tel_i,[]); %get arrays of kW and S for this grid
            end
        else
            if all(opt.alg == 'tel')
                [opt] = telescope_opt(opt,tel_i,output); %get arrays of kW and S for this grid
            else
                [opt] = telescope2_opt(opt,tel_i,output); %get arrays of kW and S for this grid
            end
        end
        if opt.pd == 6
            [Kd,Ki,Kwi,Kwa,Kc,S] = ndgrid(opt.dies.kW{tel_i},opt.inso.kW{tel_i},opt.wind.kW{tel_i}, opt.wave.kW{tel_i}, opt.curr.kW{tel_i}, opt.Smax{tel_i});
            Kd = reshape(Kd,[j*k*l*m*n*o 1]);
            Ki = reshape(Ki,[j*k*l*m*n*o 1]);
            Kwi = reshape(Kwi,[j*k*l*m*n*o 1]);
            Kwa = reshape(Kwa,[j*k*l*m*n*o 1]);
            Kc = reshape(Kc,[j*k*l*m*n*o 1]);
            S = reshape(S,[j*k*l*m*n*o 1]);
            sim_run = ones(length(S),1);
            C_temp = zeros(j*k*l*m*n*o,1);
            S_temp = zeros(j*k*l*m*n*o,1);
            X(tel_i,:) = zeros(j*k*l*m*n*o,1);
        elseif opt.pd == 5
            if opt.pm == 4 %No Dies
                [Ki,Kwi,Kwa,Kc,S] = ndgrid(opt.inso.kW{tel_i},opt.wind.kW{tel_i}, opt.wave.kW{tel_i}, opt.curr.kW{tel_i}, opt.Smax{tel_i});
                j = 1;
                %Kd = reshape(Kd,[k*l*m*n*o 1]);
                Ki = reshape(Ki,[k*l*m*n*o 1]);
                Kwi = reshape(Kwi,[k*l*m*n*o 1]);
                Kwa = reshape(Kwa,[k*l*m*n*o 1]);
                Kc = reshape(Kc,[k*l*m*n*o 1]);
                S = reshape(S,[k*l*m*n*o 1]);
                sim_run = ones(length(S),1);
                C_temp = zeros(k*l*m*n*o,1);
                S_temp = zeros(k*l*m*n*o,1);
                X(tel_i,:) = zeros(k*l*m*n*o,1);
                Kd = zeros(length(S),1);
            elseif opt.pm == 5 %No Current
                [Kd,Ki,Kwi,Kwa,S] = ndgrid(opt.dies.kW{tel_i},opt.inso.kW{tel_i},opt.wind.kW{tel_i}, opt.wave.kW{tel_i}, opt.Smax{tel_i});
                o = 1;
                Kd = reshape(Kd,[j*k*l*m*n 1]);
                Ki = reshape(Ki,[j*k*l*m*n 1]);
                Kwi = reshape(Kwi,[j*k*l*m*n 1]);
                Kwa = reshape(Kwa,[j*k*l*m*n 1]);
                S = reshape(S,[j*k*l*m*n 1]);
                sim_run = ones(length(S),1);
                C_temp = zeros(j*k*l*m*n, 1);
                S_temp = zeros(j*k*l*m*n, 1);
                X(tel_i,:) = zeros(j*k*l*m*n, 1);
                Kc = zeros(length(S),1);
            end
            
        elseif opt.pd == 2
            %1:Wi 2:In 3:Wa 4:Di 5:Cu
            if opt.pm == 4
                opt.Gr{tel_i} = opt.dies.kW{tel_i};
                g = j;
                k = 1;
                l = 1;
                m = 1;
                o = 1;
            elseif opt.pm == 2
                opt.Gr{tel_i} = opt.inso.kW{tel_i};
                g = k;
                j = 1;
                l = 1;
                m = 1;
                o = 1;
            elseif opt.pm == 5
                opt.Gr{tel_i} = opt.curr.kW{tel_i};
                g = o;
                j = 1;
                k = 1;
                l = 1;
                m = 1;
            elseif opt.pm == 1
                opt.Gr{tel_i} = opt.wind.kW{tel_i};
                g = l;
                k = 1;
                j = 1;
                m = 1;
                o = 1;
            elseif opt.pm == 3
                opt.Gr{tel_i} = opt.wave.kW{tel_i};
                g = m;
                k = 1;
                l = 1;
                j = 1;
                o = 1;
            end
            [Gr,S] = ndgrid(opt.Gr{tel_i}, opt.Smax{tel_i});
            S = reshape(S,[g*n 1]);
            sim_run = ones(length(S),1);
            C_temp = zeros(g*n,1);
            S_temp = zeros(g*n,1);
            X(tel_i,:) = zeros(g*n,1);
            if opt.pm == 4
                Kd = reshape(Gr,[g*n 1]);
                Ki = zeros(length(S),1);
                Kwi = zeros(length(S),1);
                Kwa = zeros(length(S),1);
                Kc = zeros(length(S),1);
            elseif opt.pm == 2
                Ki = reshape(Gr,[g*n 1]);
                Kd = zeros(length(S),1);
                Kwi = zeros(length(S),1);
                Kwa = zeros(length(S),1);
                Kc = zeros(length(S),1);
            elseif opt.pm == 5
                Kc = reshape(Gr,[g*n 1]);
                Kd = zeros(length(S),1);
                Kwi = zeros(length(S),1);
                Kwa = zeros(length(S),1);
                Ki = zeros(length(S),1);
            elseif opt.pm == 1
                Kwi = reshape(Gr,[g*n 1]);
                Ki = zeros(length(S),1);
                Kd = zeros(length(S),1);
                Kwa = zeros(length(S),1);
                Kc = zeros(length(S),1);
            elseif opt.pm == 3
                Kwa = reshape(Gr,[g*n 1]);
                Ki = zeros(length(S),1);
                Kwi = zeros(length(S),1);
                Kd = zeros(length(S),1);
                Kc = zeros(length(S),1);
            end
        end
    elseif all(opt.alg == 'per')
        if tel_i == 1 %can't send output to function before output is defined
            %X = zeros(opt.tel_max, opt.bf.j*opt.bf.k*opt.bf.l*opt.bf.m*opt.bf.n*((2^(opt.tel_max-1))^opt.pd));
            [opt,Kd, Ki, Kwi, Kwa, Kc, S, C_temp, S_temp,X,sim_run] = persistence_opt_v2(opt,tel_i,[]); %get arrays of kW and S for this grid
        else
            %S_in = S_temp; %surv from previous iteration
            [opt,Kd, Ki, Kwi, Kwa, Kc, S, C_temp, S_temp,X,sim_run] = persistence_opt_v2(opt,tel_i,output); %get arrays of kW and S for this grid
        end
        j = opt.bf.j;
        k = opt.bf.k;
        l = opt.bf.l;
        m = opt.bf.m;
        n = opt.bf.n;
        o = opt.bf.o;
    elseif all(opt.alg == 'p2t')
        if tel_i == 1 %can't send output to function before output is defined
            %X = zeros(opt.tel_max, opt.bf.j*opt.bf.k*opt.bf.l*opt.bf.m*opt.bf.n*((2^(opt.tel_max-1))^opt.pd));
            [opt,Kd, Ki, Kwi, Kwa, Kc, S, C_temp, S_temp,X,sim_run] = per3_tel2(opt,tel_i,[],[]); %get arrays of kW and S for this grid
        else
            %S_in = S_temp; %surv from previous iteration
            [opt,Kd, Ki, Kwi, Kwa, Kc, S, C_temp, S_temp,X,sim_run] = per3_tel2(opt,tel_i,output,temp_min_cost); %get arrays of kW and S for this grid
        end
        j = opt.bf.j;
        k = opt.bf.k;
        l = opt.bf.l;
        m = opt.bf.m;
        n = opt.bf.n;
        o = opt.bf.o;
    end


    %parallel computing via parfor
    tGrid = tic;
    %fmin_temp = opt.fmin; %making a temp variable to make the parfor loop happier
    disp(['Optimization Iteration: ',num2str(tel_i)])
    disp(['Populating grid values: j=' num2str(j) ', k=' num2str(k) ', l=' num2str(l) ', m=' num2str(m) ', n=' num2str(n) ', o=' num2str(o)])
    sim_length = length(S);
    %sim_length = 10; %test case
    parfor (i = 1:sim_length,opt.bf.maxworkers)
        if sim_run(i) == 1 %only run points where the sim_run variable is true (used for persistence opt)
            %Call Hybrid Sim function
            %disp(i)
            [C_temp(i),S_temp(i)] = ...
                simHybrid(Kd(i), Ki(i), Kwi(i), Kwa(i), Kc(i), S(i),opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);
        %else %X only has dummy points for 5D persistence
            %disp("ERROR IN PERSISTENCE OPT!!!!!")
            
        end
        if S_temp(i) < uc.uptime %update obj val X - if sim_run skipped this i then the surv will fail and cost = inf
            X(tel_i,i) = inf;
        else
            X(tel_i,i) = C_temp(i);
        end


    end
    %REMOVED TRANSPOSE BC HIGH ORDER MATRIX CAN'T USE TRANSPOSE - MIGHT MESS UP
    output.cost{tel_i} = C_temp;
    output.surv{tel_i} = S_temp;
    %output.cost_a{tel_i} = X(tel_i,:);
    %saved outputs for persistence run
    output.surv_opt = S_temp;
    output.Kd_run{tel_i} = Kd;
    output.Ki_run{tel_i} = Ki;
    output.Kwi_run{tel_i} = Kwi;
    output.Kwa_run{tel_i} = Kwa;
    output.Kc_run{tel_i} = Kc;
    output.S_run{tel_i} = S;
    %INDEXING
    output.tGrid = toc(tGrid);
    disp('Brute forcing current minimum...')
    %hybrid approach
    %I_min(tel_i) = find(X(tel_i,:)==min(X(tel_i,:))); %find index of minimum cost - doesn't work if there's multiple minimums
    temp_I = find(X(tel_i,:)==min(X(tel_i,:)))
    if length(temp_I) == 1
        I_min(tel_i) = temp_I;
    elseif opt.tar == 1
        for p = 1:length(temp_I)
            C_tot(p) = Kd(temp_I(p)) + Ki(temp_I(p)) + Kwi(temp_I(p)) + Kwa(temp_I(p)) + Kc(temp_I(p)) + S(temp_I(p));
            gen_vec = [Kd(temp_I(p)) Ki(temp_I(p)) Kwi(temp_I(p)) Kwa(temp_I(p)) Kc(temp_I(p))];
            n_gen(p) = sum(gen_vec > 0);
        end
        min_n = find(n_gen == min(n_gen));
        if length(min_n)>1
            C_min_n = C_tot(min_n);
            min_C = find(C_tot == min(C_min_n));
            I_min(tel_i) = temp_I(min_C);
        else
            I_min(tel_i) = temp_I(min_n);
        end
    elseif opt.tar == 2
        for p = 1:length(temp_I)
            C_tot(p) = Kd(temp_I(p)) + Ki(temp_I(p)) + Kwi(temp_I(p)) + Kwa(temp_I(p)) + Kc(temp_I(p)) + S(temp_I(p));
            %gen_vec = [Kd(temp_I(p)) Ki(temp_I(p)) Kwi(temp_I(p)) Kwa(temp_I(p)) Kc(temp_I(p))];
            %n_gen(p) = sum(gen_vec > 0);
        end
        [min_g,min_gi] = min(C_tot); %minimum total generation
        I_min(tel_i) = temp_I(min_gi);
    end
    output.min.kWd{tel_i} = Kd(I_min(tel_i));
    output.min.kWi{tel_i} = Ki(I_min(tel_i));
    output.min.kWwi{tel_i} = Kwi(I_min(tel_i));
    output.min.kWwa{tel_i} = Kwa(I_min(tel_i));
    output.min.kWc{tel_i} = Kc(I_min(tel_i));
    output.min.Smax{tel_i} = S(I_min(tel_i));
    temp_min_cost(tel_i) = X(tel_i,I_min(tel_i));
    if tel_i > 1
        if any([all(opt.alg == 'per'), all(opt.alg == 'p2t')]) && X(tel_i,I_min(tel_i)) > temp_min_cost(tel_i-1)
            disp('Something Wrong - new it has higher min')
%             output.min.kWd{tel_i} = output.min.kWd{tel_i-1};
%             output.min.kWi{tel_i} = output.min.kWi{tel_i-1};
%             output.min.kWwi{tel_i} = output.min.kWwi{tel_i-1};
%             output.min.kWwa{tel_i} = output.min.kWwa{tel_i-1};
%             output.min.kWc{tel_i} = output.min.kWc{tel_i-1};
%             output.min.Smax{tel_i} = output.min.Smax{tel_i-1};
%             temp_min_cost(tel_i) = temp_min_cost(tel_i-1);
        end
    end
    if tel_i > 1 && any([all(opt.alg == 'tel'), all(opt.alg == 'to2')])
        %This check won't work for persistence right now because X(tel_i-1)
        %is not saved
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
            elseif abs(output.min.kWc{tel_i} - output.min.kWc{tel_i-1})/output.min.kWc{tel_i} > opt.kwtol
                tol = false;
            elseif abs(output.min.Smax{tel_i} - output.min.Smax{tel_i-1})/output.min.Smax{tel_i} > opt.kwtol
                tol = false;
            end
            if tol == true    
                disp('Optimization Tolerance Met')
            end
        end
    end

    tel_i = tel_i +1; %update optimization iteration
end
%% GET OUTPUT VALUES
[output.min.cost,output.min.surv,output.min.CapEx,output.min.OpEx,...
    output.min.kWcost_dies,output.min.kWcost_wave, output.min.kWcost_curr, output.min.kWcost_wind, output.min.Mcost_inso,...
    output.min.Ecost_inso, output.min.Icost_inso, output.min.Strcost_inso,...
    output.min.Icost_wave, output.min.Icost_wind, output.min.Scost,output.min.Pmtrl,...
    output.min.Pinst,output.min.Pmooring, ...
    output.min.vesselcost,output.min.genrepair, ...
    output.min.turbrepair, output.min.wecrepair,...
    output.min.battreplace,output.min.battencl,output.min.genencl, ...
    output.min.fuel,output.min.triptime,output.min.runtime, ... 
    output.min.nvi,output.min.batt_L1,output.min.batt_L2, output.min.batt_lft1, output.min.batt_lft2, ...
    output.min.nfr,output.min.noc, output.min.nbr,...
    output.min.dp,output.min.S1, output.min.S2, output.min.Pdies, output.min.Pinso, output.min.Pwind, output.min.Pwave,output.min.Pcurr,...
    output.min.Ptot, output.min.width,output.min.cw,output.min.D,output.min.L, output.min.F, output.min.eff_t, output.min.pvci] ...
    = simHybrid(output.min.kWd{end},output.min.kWi{end}, output.min.kWwi{end}, output.min.kWwa{end}, output.min.kWc{end}, output.min.Smax{end}, ...
    opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);
%nvi,batt_L1,batt_L2, batt_lft1,batt_lft2, nfr,noc,nbr,dp,S1,S2,Pdies
output.min.batt_dyn_lc1 = batt.lc_nom*(output.min.Smax{end}/ ...
    (output.min.Smax{end} - (min(output.min.S1)/1000)))^batt.beta;
output.min.batt_dyn_lc2 = batt.lc_nom*(output.min.Smax{end}/ ...
    (output.min.Smax{end} - (min(output.min.S2)/1000)))^batt.beta;
output.min.CFd = mean(output.min.Pdies)/(1000*output.min.kWd{end});
output.min.CFi = mean(output.min.Pinso)/(1000*output.min.kWi{end});
output.min.CFwi = mean(output.min.Pwind)/(1000*output.min.kWwi{end});
output.min.CFwa = mean(output.min.Pwave)/(1000*output.min.kWwa{end});
output.min.CFc = mean(output.min.Pcurr)/(1000*output.min.kWc{end});
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

