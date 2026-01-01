function [opt,Kd, Ki, Kwi, Kwa, Kc, S, C_temp, S_temp, X, sim_run] = per3_tel2(opt,tel_i, output,temp_min_cost)
%Persistence 3 Telescope 2 - Written 5-23-23
%  This function sets up the grid arrays for each grid
if tel_i == 1
    j = opt.bf.j;
    k = opt.bf.k;
    l = opt.bf.l;
    m = opt.bf.m;
    n = opt.bf.n;
    o = opt.bf.o;
    opt.dies.kW{tel_i} = linspace(opt.dies.kW_1,opt.dies.kW_m,j);              %[kW]
    opt.inso.kW{tel_i} = linspace(opt.inso.kW_1,opt.inso.kW_m,k);              %[kW]
    opt.wind.kW{tel_i} = linspace(opt.wind.kW_1,opt.wind.kW_m,l);              %[kW]
    opt.wave.kW{tel_i} = linspace(opt.wave.kW_1,opt.wave.kW_m,m);              %[kW]
    opt.curr.kW{tel_i} = linspace(opt.curr.kW_1,opt.curr.kW_m,o);              %[kW]
    opt.Smax{tel_i} = linspace(opt.Smax_1,opt.Smax_n,n);                       %[kWh]
    [Kd,Ki,Kwi,Kwa,Kc,S] = ndgrid(opt.dies.kW{tel_i},opt.inso.kW{tel_i},opt.wind.kW{tel_i}, opt.wave.kW{tel_i}, opt.curr.kW{tel_i}, opt.Smax{tel_i});
    Kd = reshape(Kd,[j*k*l*m*n*o 1]);
    Ki = reshape(Ki,[j*k*l*m*n*o 1]);
    Kwi = reshape(Kwi,[j*k*l*m*n*o 1]);
    Kwa = reshape(Kwa,[j*k*l*m*n*o 1]);
    Kc = reshape(Kc,[j*k*l*m*n*o 1]);
    S = reshape(S,[j*k*l*m*n*o 1]);
    C_temp = zeros(j*k*l*m*n*o,1);
    S_temp = zeros(j*k*l*m*n*o,1);
    %X(tel_i,((j*k*l*m*n)+1):end) = inf;
    sim_run = ones(length(S),1);
    X = zeros(tel_i,sum(sim_run));
elseif tel_i < 4 %Normal Persistence for 2&3
    %get surv data from previous run
    surv = output.surv_opt; 
    %zero out the not interesting cases
    surv(surv< opt.pl) = 0;
    surv(surv> opt.pr) = 0;
   
    %previous iteration reshaped grid style arrays
    Kd_old = output.Kd_run{tel_i-1};
    Ki_old = output.Ki_run{tel_i-1};
    Kwi_old = output.Kwi_run{tel_i-1};
    Kwa_old = output.Kwa_run{tel_i-1};
    Kc_old = output.Kc_run{tel_i-1};
    S_old = output.S_run{tel_i-1};

    %remove elements of grid style arrays where sruv = 0
    Kd_old(surv == 0) = [];
    Ki_old(surv == 0) = [];
    Kwi_old(surv == 0) = [];
    Kwa_old(surv == 0) = [];
    Kc_old(surv == 0) = [];
    S_old(surv == 0) = [];

    %Update the resolution (should be x2 each time but need to keep old
    %points!!!)
    opt.bf.j = (opt.bf.j*2)-1;
    opt.bf.k = (opt.bf.k*2)-1;
    opt.bf.l = (opt.bf.l*2)-1;
    opt.bf.m = (opt.bf.m*2)-1;
    opt.bf.n = (opt.bf.n*2)-1;
    opt.bf.o = (opt.bf.o*2)-1;
    j = opt.bf.j;
    k = opt.bf.k;
    l = opt.bf.l;
    m = opt.bf.m;
    n = opt.bf.n;
    o = opt.bf.o;
    %double the resolution on each dimension array
    opt.dies.kW{tel_i} = linspace(opt.dies.kW_1,opt.dies.kW_m,j);        %[kW]
    Kd_t = opt.dies.kW{tel_i};
    opt.inso.kW{tel_i} = linspace(opt.inso.kW_1,opt.inso.kW_m,k);        %[kW]
    Ki_t = opt.inso.kW{tel_i};
    opt.wind.kW{tel_i} = linspace(opt.wind.kW_1,opt.wind.kW_m,l);        %[kW]
    Kwi_t = opt.wind.kW{tel_i};
    opt.wave.kW{tel_i} = linspace(opt.wave.kW_1,opt.wave.kW_m,m);        %[kW]
    Kwa_t = opt.wave.kW{tel_i};
    opt.curr.kW{tel_i} = linspace(opt.curr.kW_1,opt.curr.kW_m,o);        %[kW]
    Kc_t = opt.curr.kW{tel_i};
    opt.Smax{tel_i} = linspace(opt.Smax_1,opt.Smax_n,n);                 %[kWh]
    S_t = opt.Smax{tel_i};
    
    %make the new grid style arrays for parfor loop
    %[Kd, Ki, Kwi, Kwa, S] = loopyloops(S_t, S_old, Kwa_t, Kwa_old, Kwi_t, Kwi_old, Ki_t, Ki_old, Kd_t, Kd_old,opt);
    [Kd, Ki, Kwi, Kwa, Kc, S, sim_run] = loopyloops_v2(S_t, S_old, Kc_t, Kc_old, Kwa_t, Kwa_old, Kwi_t, Kwi_old, Ki_t, Ki_old, Kd_t, Kd_old,opt);
    
    %initialize X, C_temp and S_temp
    X = zeros(tel_i,sum(sim_run));
    C_temp = zeros(length(S),1);
    S_temp = zeros(length(S),1);
elseif tel_i >= 4
    %Telescope
    if tel_i == 4
        %Smoosh together all the old data
        surv = [output.surv{1,1}; output.surv{1,2}; output.surv{1,3}];
        cost = [output.cost{1,1}; output.cost{1,2}; output.cost{1,3}]; %all cost from previous iteration
        min_cost = min([temp_min_cost(1), temp_min_cost(2), temp_min_cost(3)]); %min from all previous iterations
        %previous iteration reshaped grid style arrays
        Kd_old = [output.Kd_run{1,1}; output.Kd_run{1,2}'; output.Kd_run{1,3}'];
        Ki_old = [output.Ki_run{1,1}; output.Ki_run{1,2}'; output.Ki_run{1,3}'];
        Kwi_old = [output.Kwi_run{1,1}; output.Kwi_run{1,2}'; output.Kwi_run{1,3}'];
        Kwa_old = [output.Kwa_run{1,1}; output.Kwa_run{1,2}'; output.Kwa_run{1,3}'];
        Kc_old = [output.Kc_run{1,1}; output.Kc_run{1,2}'; output.Kc_run{1,3}'];
        S_old = [output.S_run{1,1}; output.S_run{1,2}'; output.S_run{1,3}'];
    else
        %get surv data from previous run
        %surv = output.surv_opt; 
        %zero out the not interesting cases

        surv = [output.surv{1,1}; output.surv{1,2}; output.surv{1,3}; output.surv{1,4}];
        cost = [output.cost{1,1}; output.cost{1,2}; output.cost{1,3}; output.cost{1,4}]; %all cost from previous iteration
        min_cost = min([temp_min_cost(1), temp_min_cost(2), temp_min_cost(3), temp_min_cost(4)]); %min from all previous iterations
        %previous iteration reshaped grid style arrays
        Kd_old = [output.Kd_run{1,1}; output.Kd_run{1,2}'; output.Kd_run{1,3}'; output.Kd_run{1,4}'];
        Ki_old = [output.Ki_run{1,1}; output.Ki_run{1,2}'; output.Ki_run{1,3}'; output.Ki_run{1,4}'];
        Kwi_old = [output.Kwi_run{1,1}; output.Kwi_run{1,2}'; output.Kwi_run{1,3}'; output.Kwi_run{1,4}'];
        Kwa_old = [output.Kwa_run{1,1}; output.Kwa_run{1,2}'; output.Kwa_run{1,3}'; output.Kwa_run{1,4}'];
        Kc_old = [output.Kc_run{1,1}; output.Kc_run{1,2}'; output.Kc_run{1,3}'; output.Kc_run{1,4}'];
        S_old = [output.S_run{1,1}; output.S_run{1,2}'; output.S_run{1,3}'; output.S_run{1,4}'];
    end

    surv(surv< opt.pl) = 0;
    surv(surv>opt.pr) = 0;
%     %% TEST FOR WA
%     %get surv data from previous run
%     surv = output.surv_opt; 
%     %zero out the not interesting cases
%     %surv(surv< 0.975) = 0;
%     %surv(surv>0.993) = 0;
%    
%     cost = output.cost{tel_i-1}; %all cost from previous iteration
%     min_cost = temp_min_cost(tel_i-1); %min from the previous iteration
%     %previous iteration reshaped grid style arrays
%     Kd_old = output.Kd_run{tel_i-1};
%     Ki_old = output.Ki_run{tel_i-1};
%     Kwi_old = output.Kwi_run{tel_i-1};
%     Kwa_old = output.Kwa_run{tel_i-1};
%     Kc_old = output.Kc_run{tel_i-1};
%     S_old = output.S_run{tel_i-1};
%     cost(surv == 0) = nan;
%     %% End TEST SECTION 

    %remove elements of grid style arrays where sruv = 0
    cost(surv == 0) = nan;
    %cost(cost > 1.2*min_cost) = nan; %remove all points that are too heavy
    cost(cost< (min_cost*opt.tl) | cost > (min_cost*opt.tr)) = nan; 
    Kd_old(isnan(cost)) = nan;
    Ki_old(isnan(cost)) = nan;
    Kwi_old(isnan(cost)) = nan;
    Kwa_old(isnan(cost)) = nan;
    Kc_old(isnan(cost)) = nan;
    S_old(isnan(cost)) = nan;

    %remove any nan
    Kd_old(isnan(Kd_old)) = [];
    Ki_old(isnan(Ki_old)) = [];
    Kwi_old(isnan(Kwi_old)) = [];
    Kwa_old(isnan(Kwa_old)) = [];
    Kc_old(isnan(Kc_old)) = [];
    S_old(isnan(S_old)) = [];

    %find new min and max points of next telescoped grid (accounts for a
    %min at the edge of the tele_i-1 grid
    if tel_i == 4
        j = 5; k = 5; l = 5; m= 5; n = 5; o = 5; %going to be 9 points - CAN ONLY BE AN ODD NUMBER
    else
        j = 3; k = 3; l = 3; m= 3; n = 3; o = 3; %going to be 9 points - CAN ONLY BE AN ODD NUMBER
    end
    %Might be funky for 5th it
    disp('new grid spacing...')
    if tel_i == 4
        %opt.bf.j %number of points
        dx = opt.bf.M/opt.bf.j; %number of points
        ds = opt.bf.N/opt.bf.n; %number of points
    elseif tel_i == 5
        dx = 2*(opt.bf.M/opt.bf.j)/j;
        ds = 2*(opt.bf.N/opt.bf.n)/n;
    end

    Big_Matrix = [];
    %disp('This will create: ',num2str((j^5)*length(Kd_old)),'points')
    %matrix of old points before adjustments
    old_dim = size(Kd_old)
    if old_dim(1) > 1
        old_points = [Kd_old Ki_old Kwi_old Kwa_old Kc_old S_old];
    else
        old_points = [Kd_old' Ki_old' Kwi_old' Kwa_old' Kc_old' S_old'];
    end
    parfor (gp = 1 : length(Kd_old),opt.bf.maxworkers)
    %for gp = 1:length(Kd_old)
        %define grid boundaries 1 dx away from point
        j_1 = Kd_old(gp)-dx;
        j_max = Kd_old(gp)+dx;

        k_1 = Ki_old(gp)-dx;
        k_max = Ki_old(gp)+dx;

        l_1 = Kwi_old(gp)-dx;
        l_max = Kwi_old(gp)+dx;

        m_1 = Kwa_old(gp)-dx;
        m_max = Kwa_old(gp)+dx;

        o_1 = Kc_old(gp)-dx;
        o_max = Kc_old(gp)+dx;

        n_1 = S_old(gp)-ds;
        n_max = S_old(gp)+ds;
    
        %Make new arrays
        kd_temp = linspace(j_1,j_max,j);             %[kW]

        ki_temp = linspace(k_1,k_max,k);             %[kW]

        kwi_temp = linspace(l_1,l_max,l);            %[kW]

        kwa_temp = linspace(m_1,m_max,m);   

        kc_temp = linspace(o_1,o_max,o);   

        smax_temp = linspace(n_1,n_max,n);
        
        [Kd,Ki,Kwi,Kwa,Kc,S] = ndgrid(kd_temp,ki_temp, kwi_temp, kwa_temp, kc_temp, smax_temp);
        Kd = reshape(Kd,[((j)^6) 1]);
        Ki = reshape(Ki,[((j)^6) 1]);
        Kwi = reshape(Kwi,[((j)^6) 1]);
        Kwa = reshape(Kwa,[((j)^6) 1]);
        Kc = reshape(Kc,[((j)^6) 1]);
        S = reshape(S,[((j)^6) 1]);

        %remove bad points
        max_x = opt.bf.M;
        max_s = opt.bf.N;
        ind_rm = Kd<0 | Ki<0 | Kwi <0 |Kwa<0 | Kc < 0 | S<0;
        Kd(ind_rm) = [];
        Ki(ind_rm)  = [];
        Kwi(ind_rm)  = [];
        Kwa(ind_rm) = [];
        Kc(ind_rm)  = [];
        S(ind_rm)  = [];
        
        ind_rm2 = Kd> max_x | Ki>max_x | Kwi > max_x |Kwa > max_x | Kc > max_x | S> max_s;
        Kd(ind_rm2) = [];
        Ki(ind_rm2)  = [];
        Kwi(ind_rm2)  = [];
        Kwa(ind_rm2) = [];
        Kc(ind_rm2)  = [];
        S(ind_rm2)  = [];

        %make any wave points below WECSIM be zero wave
        Kwa(Kwa<0.2144) = 0;

        gp_grid = [Kd'; Ki'; Kwi'; Kwa'; Kc'; S';];
        Big_Matrix = [Big_Matrix gp_grid];
    end
    Big_Matrix = round(Big_Matrix,4);
    new_points = unique(Big_Matrix','rows','stable');
    %new_points = new_points;
    old_points = round(old_points,4);
    [extra_points,in,io] = intersect(new_points, old_points,'stable','rows');
    new_points(in,:) = [];
    new_points = new_points';
    disp(size(Big_Matrix))
    disp(size(new_points))
    Kd = new_points(1,:);
    Ki = new_points(2,:);
    Kwi = new_points(3,:);
    Kwa = new_points(4,:);
    Kc = new_points(5,:);
    S = new_points(6,:);

    %initialize X, C_temp and S_temp
    disp(size(S))
    sim_run = ones(length(S),1);
    X = zeros(tel_i,sum(sim_run));
    C_temp = zeros(length(S),1);
    S_temp = zeros(length(S),1);
    
end

if opt.pd == 5
    if opt.pm == 4 %no dies
        sim_run(Kd>0) = 0; %dont run simHybrid is Kd is greater than zero
    elseif opt.pm == 5 %no current
        sim_run(Kc>0) = 0;
    end
elseif opt.pd == 3
    if opt.pm == 12 %wind+inso
        sim_run(Kd>0) = 0; %dont run simHybrid is Kd is greater than zero
        sim_run(Kwa>0) = 0; %dont run simHybrid is Kd is greater than zero
        sim_run(Kc>0) = 0; %dont run simHybrid is Kd is greater than zero
    end
end
end