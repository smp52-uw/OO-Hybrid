function [opt,Kd, Ki, Kwi, Kwa, S, C_temp, S_temp, X, sim_run] = persistence_opt(opt,tel_i,output,S_in)
%Telescope_opt written by Sarah Palmer Jan 2023
%  This function sets up the grid arrays for each telescoping grid
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
    X = zeros(j*k*l*m*n,1);
    sim_run = ones(length(S),1);
else
    %get surv data from previous run
    surv = output.surv_opt; 
    %zero out the not interesting cases
    surv(surv< 0.95) = 0;
    surv(surv>0.998) = 0;
    surv_g = reshape(surv, [opt.bf.j opt.bf.k opt.bf.l opt.bf.m opt.bf.n]);
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
    disp(size(Kd))
    disp((j*k*l*m*n))

    S_f = griddedInterpolant(surv_g);
    sim_run = S_f(Kd, Ki, Kwi, Kwa,S); %interpolate to the denser grid
    sim_run(sim_run>0) = 1; %any non-zero value should be 1

    Kd = reshape(Kd,[(j*k*l*m*n) 1]);
    Ki = reshape(Ki,[(j*k*l*m*n) 1]);
    Kwi = reshape(Kwi,[(j*k*l*m*n) 1]);
    Kwa = reshape(Kwa,[(j*k*l*m*n) 1]);
    S = reshape(S,[(j*k*l*m*n) 1]);
    C_temp = zeros(j*k*l*m*n,1);
    S_temp = zeros(j*k*l*m*n,1);
    X = zeros(j*k*l*m*n,1);
    sim_run = reshape(sim_run,[(j*k*l*m*n) 1]);
    %set opt.simrun to 0 or 1 depending on data (linearly interpolate 0/1 from coarse grid and
    %then round up should give the right 0/1 in the new grid)


end
end