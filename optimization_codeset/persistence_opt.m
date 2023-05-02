function [opt,Kd, Ki, Kwi, Kwa, S, C_temp, S_temp, X, sim_run] = persistence_opt(opt,tel_i, output)
%Telescope_opt written by Sarah Palmer Jan 2023
%  This function sets up the grid arrays for each telescoping grid
%X = zeros(opt.tel_max, opt.bf.j*opt.bf.k*opt.bf.l*opt.bf.m*opt.bf.n*(2^opt.tel_max));
if tel_i == 1
    j = opt.bf.j;
    k = opt.bf.k;
    l = opt.bf.l;
    m = opt.bf.m;
    n = opt.bf.n;
    opt.dies.kW{tel_i} = linspace(opt.dies.kW_1,opt.dies.kW_m,j);              %[kW]
    opt.inso.kW{tel_i} = linspace(opt.inso.kW_1,opt.inso.kW_m,k);              %[kW]
    opt.wind.kW{tel_i} = linspace(opt.wind.kW_1,opt.wind.kW_m,l);              %[kW]
    opt.wave.kW{tel_i} = linspace(opt.wave.kW_1,opt.wave.kW_m,m);              %[kW]
    opt.Smax{tel_i} = linspace(opt.Smax_1,opt.Smax_n,n);                       %[kWh]
    [Kd,Ki,Kwi,Kwa,S] = ndgrid(opt.dies.kW{tel_i},opt.inso.kW{tel_i},opt.wind.kW{tel_i}, opt.wave.kW{tel_i}, opt.Smax{tel_i});
    Kd = reshape(Kd,[j*k*l*m*n 1]);
    Ki = reshape(Ki,[j*k*l*m*n 1]);
    Kwi = reshape(Kwi,[j*k*l*m*n 1]);
    Kwa = reshape(Kwa,[j*k*l*m*n 1]);
    S = reshape(S,[j*k*l*m*n 1]);
    C_temp = zeros(j*k*l*m*n,1);
    S_temp = zeros(j*k*l*m*n,1);
    %X(tel_i,((j*k*l*m*n)+1):end) = inf;
    sim_run = ones(length(S),1);
    X = zeros(tel_i,sum(sim_run));
else
    %get surv data from previous run
    surv = output.surv_opt; 
    %zero out the not interesting cases
    surv(surv< 0.95) = 0;
    surv(surv>0.998) = 0;
    Kd_run = opt.Kd_run;
    [Kd_1,Ki_1,Kwi_1,Kwa_1,S_1] = ndgrid(opt.dies.kW{tel_i-1},opt.inso.kW{tel_i-1},opt.wind.kW{tel_i-1}, opt.wave.kW{tel_i-1}, opt.Smax{tel_i-1});
    Kd_old = reshape(Kd_1,[(j*k*l*m*n) 1]);
    %surv_g = reshape(surv, [length(Kd_1) length(Ki_1) length(Kwi_1) length(Kwa_1) length(S_1)]);
    %interp surv_g from run kd to full array kd from last it
    %if not in run kd make surv_g =0
    %reshape
    %Update the resolution (should be x2 each time tel_i increases by 1)
    opt.bf.j = opt.bf.j*2;
    opt.bf.k = opt.bf.k*2;
    opt.bf.l = opt.bf.l*2;
    opt.bf.m = opt.bf.m*2;
    opt.bf.n = opt.bf.n*2;
    j = opt.bf.j;
    k = opt.bf.k;
    l = opt.bf.l;
    m = opt.bf.m;
    n = opt.bf.n;
    %double the grid 
    %tel_i = 2, j*(2^1), tel_1 = 3, j*(2^2)
    opt.dies.kW{tel_i} = linspace(opt.dies.kW_1,opt.dies.kW_m,j);        %[kW]
    opt.inso.kW{tel_i} = linspace(opt.inso.kW_1,opt.inso.kW_m,k);        %[kW]
    opt.wind.kW{tel_i} = linspace(opt.wind.kW_1,opt.wind.kW_m,l);        %[kW]
    opt.wave.kW{tel_i} = linspace(opt.wave.kW_1,opt.wave.kW_m,m);        %[kW]
    opt.Smax{tel_i} = linspace(opt.Smax_1,opt.Smax_n,n);                 %[kWh]

    [Kd,Ki,Kwi,Kwa,S] = ndgrid(opt.dies.kW{tel_i},opt.inso.kW{tel_i},opt.wind.kW{tel_i}, opt.wave.kW{tel_i}, opt.Smax{tel_i});
    [Kd_1,Ki_1,Kwi_1,Kwa_1,S_1] = ndgrid(opt.dies.kW{tel_i-1},opt.inso.kW{tel_i-1},opt.wind.kW{tel_i-1}, opt.wave.kW{tel_i-1}, opt.Smax{tel_i-1});
    disp(size(Kd))
    disp((j*k*l*m*n))

    S_f = griddedInterpolant(Kd_1,Ki_1,Kwi_1,Kwa_1,S_1, surv_g);
    sim_run = S_f(Kd, Ki, Kwi, Kwa,S); %interpolate to the denser grid
    sim_run(sim_run>0) = 1; %any non-zero value should be 1

    Kd = reshape(Kd,[(j*k*l*m*n) 1]);
    Ki = reshape(Ki,[(j*k*l*m*n) 1]);
    Kwi = reshape(Kwi,[(j*k*l*m*n) 1]);
    Kwa = reshape(Kwa,[(j*k*l*m*n) 1]);
    S = reshape(S,[(j*k*l*m*n) 1]);
    C_temp = zeros(j*k*l*m*n,1);
    S_temp = zeros(j*k*l*m*n,1);
%     X(tel_i,((j*k*l*m*n)+1):end) = inf;
    sim_run = reshape(sim_run,[(j*k*l*m*n) 1]);
    X = zeros(tel_i,sum(sim_run));
    %remove all indices where sim run is zero
    Kd(sim_run == 0) = [];
    Ki(sim_run == 0) = [];
    Kwi(sim_run == 0) = [];
    Kwa(sim_run == 0) = [];
    S(sim_run == 0) = [];
    C_temp(sim_run == 0) = [];
    S_temp(sim_run == 0) = [];
    sim_run(sim_run == 0) = [];

    %set opt.simrun to 0 or 1 depending on data (linearly interpolate 0/1 from coarse grid and
    %then round up should give the right 0/1 in the new grid)


end
end