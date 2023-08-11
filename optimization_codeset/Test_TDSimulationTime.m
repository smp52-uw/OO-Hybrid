%%%Test simulation run times

%%Call Hybrid Opt functions needed
optInputs

data = load(loc,loc);
data = data.(loc);
%load current data
curr = load(cloc);
data.curr = curr; %add current data to data structure

%curve-fit device scatters, find polyvals
opt.p_dev.t = calcDeviceVal('turbine',[],econ.wind_n);
opt.p_dev.d_cost = calcDeviceVal('dieselcost',[],econ.diescost_n);
opt.p_dev.d_mass = calcDeviceVal('dieselmass',[],econ.diesmass_n);
opt.p_dev.d_size = calcDeviceVal('dieselsize',[],econ.diessize_n);
opt.p_dev.d_burn = calcDeviceVal('dieselburn',[],econ.diesburn_n);
opt.p_dev.d_vol = calcDeviceVal('dieselvol',[],econ.diesvol_n);
opt.p_dev.b_size = calcDeviceVal('lfp_vol',[],econ.battsize_n);
[opt.p_dev.b,~,opt.p_dev.kWhmax] = calcDeviceVal('agm',[],econ.batt_n);

%HYBRID Prep Function Calls
[data, opt] = prepHybrid(data,opt,uc(c),wave,atmo,inso,cturb);
%Hybrid load case call
[uc(c).loaddata, loadseries] = GenerateLoadCases_v4(data); %updated to even hour loads
uc(c).draw = loadseries.L(uc(c).loadcase,:);
opt.fmin = false;

%Set Test Point
Kd = 0.05;
Ki = 1.5;
Kwi = 2;
Kwa = 1.25;
Kc = 0;
S = 10;

%Run TD Sim
i = 1;
tTot = tic;
[C_temp(i),S_temp(i)] = ...
                simHybrid(Kd(i), Ki(i), Kwi(i), Kwa(i), Kc(i), S(i),opt,data,atmo,batt,econ,uc(c),bc,dies,inso,wave,turb,cturb);
disp(['Time Domain Simulation complete after ' ...
    num2str(round(toc(tTot),4)) ' seconds.'])
C_temp(i)
S_temp(i)
%Run Surrogate Sim
tTots = tic;
[C_temp(i),S_temp(i)] = ...
                simHybridsg(Kd(i), Ki(i), Kwi(i), Kwa(i), Kc(i), S(i),opt,data,atmo,batt,econ,uc(c),bc,dies,inso,wave,turb,cturb);
disp(['Surrogate Simulation after ' ...
    num2str(round(toc(tTots),4)) ' seconds.'])
C_temp(i)
S_temp(i)