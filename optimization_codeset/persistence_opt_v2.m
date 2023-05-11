function [opt,Kd, Ki, Kwi, Kwa, Kc, S, C_temp, S_temp, X, sim_run] = persistence_opt_v2(opt,tel_i, output)
%Persistence_opt written by Sarah Palmer Apr 2023
%  This function sets up the grid arrays for each telescoping grid
%X = zeros(opt.tel_max, opt.bf.j*opt.bf.k*opt.bf.l*opt.bf.m*opt.bf.n*(2^opt.tel_max));
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
else
    %get surv data from previous run
    surv = output.surv_opt; 
    %zero out the not interesting cases
%     surv(surv< 0.95) = 0;
%     surv(surv>0.998) = 0;
    surv(surv< 0.975) = 0;
    surv(surv>0.993) = 0;
   
    %previous iteration reshaped grid style arrays
    Kd_old = output.Kd_run;
    Ki_old = output.Ki_run;
    Kwi_old = output.Kwi_run;
    Kwa_old = output.Kwa_run;
    Kc_old = output.Kc_run;
    S_old = output.S_run;

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
    [Kd, Ki, Kwi, Kwa, Kc, S] = loopyloops_v2(S_t, S_old, Kc_t, Kc_old, Kwa_t, Kwa_old, Kwi_t, Kwi_old, Ki_t, Ki_old, Kd_t, Kd_old,opt);
    %sim_run is 1 for all points added to the new grid like arrays
    sim_run = ones(length(S),1);
    
    %initialize X, C_temp and S_temp
    X = zeros(tel_i,sum(sim_run));
    C_temp = zeros(length(S),1);
    S_temp = zeros(length(S),1);


end
end