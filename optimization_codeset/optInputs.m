%simulation settings
%interactive job
econ.wave.scen = 1; %scenario indicator 1:C, 2:OC, 3:OD
econ.inso.scen = 1; %scenario indicator 1:AU, 2:HU (don't use human)
econ.wind.scen = 2; %scenario indicator 1:OD, 2:C

opt.bf.j = 9;
opt.bf.k = 9;
opt.bf.l = 9;
opt.bf.m = 9;
opt.bf.n = 9;
opt.bf.o = 9;
opt.ffa.pop = 25; %population size 

%% Optimization Algorithm
opt.tel_max = 3; %maximum number of telescoping iterations
opt.ffa.max = 10; %max number of firefly iterations
opt.ctol = 1/100; %Tolerance on minimum cost [1% of cost]
opt.kwtol = 1/100; %Tolerance on kW or kWh of minimum system [1% of kW or kWh]
opt.alg = 'ffa'; %'tel' -Telescope, 'per' -persistence band, 'to2' -tel 2 box, 'p2t - per to tel, 'ffa'-firefly

opt.pd = 6; %6 = 6D hybrid sim, 2 = 1 gen + batt, 3 = 2 gen + batt
opt.pm = 4; %power module (for 2D sim), 1:Wi 2:In 3:Wa 4:Di 5:Cu 12:Wi+In
opt.tar = 2; %1 = mass, 2 = gen cap 
opt.pl = 0.975; %persistence left side
opt.pr = 0.995; %persistence right side
opt.tl = 0.97; %telescope left side
opt.tr = 1.03; %telescope right side
opt.drun = 1; %Diesel run method: 1=1 hour, 2=til batt full

opt.bf.M = 8; %[kW] max kW in grid
opt.bf.N = 500; %[kWh] max Smax in grid
%% Run Inputs
opt.allscenuses = 0;
opt.alllocuses = 0;
opt.sens = 0;
opt.tdsens = 0;
opt.senssm = 0;
opt.highresobj = 0;
%pm = 4; %power module, 1:Wi 2:In 3:Wa 4:Di
c = 2;  %use case 1:ST 2:LT (Only use LT)
loc = 'argBasin'; %location
%loc = 'cosEndurance_wa'
cloc = 'HYCOM_AB_mod_2018'; %current data location
%batch = false;
if ~exist('batchtype','var')
    batchtype = [];
    batchscen = [];
    batchloc = [];
    batchc = [];
end
if isequal(batchtype,'ssm')
    econ.wave.scen = batchscen;
    econ.inso.scen = batchscen;
    econ.wind.scen = batchscen;
    opt.bf.j = 11;
    opt.bf.k = 11;
    opt.bf.l = 11;
    opt.bf.m = 11;
    opt.bf.n = 11;
    opt.allscenuses = 0;
    opt.alllocuses = 0;
    opt.sens = 0;
    opt.tdsens = 0;
    opt.senssm = 1;
    opt.highresobj = 0;
    %pm = batchpm;
    c = batchc;
    loc = batchloc;
    %batch = true;
elseif isequal(batchtype,'alllocuses')
    econ.wave.scen = batchscen; 
    econ.inso.scen = batchscen;
    econ.wind.scen = batchscen;
    opt.bf.j = 11;
    opt.bf.k = 11;
    opt.bf.l = 11;
    opt.bf.m = 11;
    opt.bf.n = 11;
    opt.allscenuses = 0;
    opt.alllocuses = 1;
    opt.sens = 0;
    opt.tdsens = 0;
    opt.senssm = 0;
    opt.highresobj = 0;
    %pm = batchpm;
    c = [];
    loc = [];
    %batch = true;
elseif isequal(batchtype,'hros')
    econ.wave.scen = 1; %scenario indicator 1:C,2:OC,3:OD
    opt.bf.j = 11;
    opt.bf.k = 11;
    opt.bf.l = 11;
    opt.bf.m = 11;
    opt.bf.n = 11;
    opt.allscenuses = 0;
    opt.alllocuses = 1;
    opt.sens = 0;
    opt.tdsens = 0;
    opt.senssm = 0;
    opt.highresobj = 1;
    %pm = 3;
    c = batchc;
    loc = batchloc;
    %batch = true;
elseif isequal(batchtype,'sens')
    opt.tuning_array = linspace(0,2.25,10);
    opt.tuned_parameter = 'wiv';
    econ.wave.scen = batchscen; 
    opt.bf.j = 11;
    opt.bf.k = 11;
    opt.bf.l = 11;
    opt.bf.m = 11;
    opt.bf.n = 11;
    opt.allscenuses = 0;
    opt.alllocuses = 0;
    opt.sens = 1;
    opt.tdsens = 0;
    opt.senssm = 0;
    opt.highresobj = 0;
    %pm = 3;
    c = batchc;
    loc = batchloc;
    %batch = true;
end

%check to see if HPC
if feature('numcores') < 10
    opt.bf.j = 3;
    opt.bf.k = 3;
    opt.bf.l = 3;
    opt.bf.n = 3;
    opt.bf.m = 3;
    opt.bf.o = 3;
end

%% strings
opt.locations = {'argBasin';'cosEndurance_wa'; ...
    'cosPioneer';'irmSea';'souOcean'};
%opt.powermodules = {'wind';'inso';'wave';'dies'};
opt.usecases = {'short term';'long term'};
opt.wavescens = {'Conservative';'Optimistic Cost';'Optimistic Durability'};
% if pm == 1
%     opt.scens = {'optimistic durability','conservative'};
% elseif pm == 2
%     opt.scens = {'automated','human'};
% elseif pm == 3
%     opt.scens = opt.wavescens;
% elseif pm == 4
%     opt.scens = {'default'};
% end

%% ECONOMIC
%polynomial fits
econ.batt_n = 1;                    %[~]
econ.battsize_n = 1;                %[~]
econ.wind_n = 1;                    %[~]
econ.diescost_n = 1;                %[~]
econ.diesmass_n = 1;                %[~]   
econ.diessize_n = 1;                %[~]   
econ.diesburn_n = 1;                %[~]   
econ.diesvol_n = 1;                %[~]  

%HYBRID PLATFORM DEFINITION IS NEEDED - THIS IS A TEMPORARY FIX
load('mdd_output_wave.mat')
econ.platform.cost = cost;
econ.platform.depth = depth;
econ.platform.diameter = diameter;

econ.platform.boundary = 2; %1: multi-mooring, 2: 8m diameter limit
econ.platform.boundary_di = 12; %[m] for multi-mooring
econ.platform.boundary_mf = 3; %multi line factor
%End of temporary fix section

clear cost depth diameter
econ.platform.wf = 5;               %weight factor (of light ship)
econ.platform.steel = 2000;          %[$/metric ton], steelbenchmarker
econ.platform.t_i = [6 12];         %[h] added h for inst
econ.platform.d_i = [500 5000];     %[m] depth for inst cost

%econ.vessel.osvcost = 15000*1.15;  %[$/day] 2020->2022
%econ.vessel.speed = 10;             %[kts]
%econ.vessel.t_mosv = 6;             %[h] time on site for maint (osv)
econ.vessel.speccost = 50000*1.15;       %[$/day] 2020->2022
%econ.vessel.t_ms = 6;               %[h] time on site  (spec) - obselete
%econ.vessel.int = 1;                %[vessel int/year] -obselete

%battery 
econ.batt.enclmult = 1;             %multiplier on battery cost for encl
%wind
econ.wind.installed = 5120;         %[$/kW] installed cost (DWR, 2022)
econ.wind.tcm = 1;                  %turbine cost multiplier (sens var)
%econ.wind.mim = 137/49;             %marine installment multiplier (CoWR)
econ.wind.marinization = 1.8;       %[CoWR]
%solar
econ.inso.module = 480;             %[$/kW], all SCB, Q12022
econ.inso.installation = 160;       %[$/kW]
econ.inso.electrical = 310;         %[$/kW]
econ.inso.structural = 90;         %[$/kW]
econ.inso.marinization = 1.2;       %[~]
econ.inso.pcm = 1;                  %cost multiplier (sens var)
%wave costs
econ.wave.scenarios = 3;            %number of scenarios
econ.wave.costmult_con = 10;         %conservative cost multiplier
econ.wave.costmult_opt = 4;         %optimistic cost multiplier
%diesel costs
econ.dies.fcost = 1.4;              %[$/L] diesel fuel cost
econ.dies.enclcost = 5000*1.19;     %[$], 2018->2022
econ.dies.enclcap = 1.5;            %[m^3]
econ.dies.autostart = 3000*1.15;    %[$], 2020->2022
%econ.dies.fail = .2;                %failures per year
econ.dies.gcm = 1;                  %generator cost multiplier (sens var)

%ENERGY
%wind parameters
turb.uci = 3;               %[m/s] guess
turb.ura = 11;              %[m/s] awea
turb.uco = 30;              %[m/s] guess
turb.eta = 0.35;            %[~] guess
turb.clearance = 4;         %[m] surface to bottom of swept area clearance
turb.wf = 80;               %[kg/kW] - updated with tower
%turb.nu = 0.26;
% turb.spar_t = 0.04;         %[m] spar thickness
% turb.spar_ar = 6;           %aspect ratio 
% turb.spar_bm = 10;           %buoyancy multiplier
% current parameters

cturb.uci = 0.5;              %[m/s] 
cturb.ura = 2;                %[m/s] 
cturb.uco = 3;                %[m/s] 
cturb.eta = 0.4*0.7;          %[~] guess (Brian)
cturb.wf = 96;               %[kg/kW] RM2
cturb.clearance = 1.5;        % [m] under the water (from Brian)

%solar parameters
inso.rated = 1;             %[kW/m^2] from Brian
inso.eff = 0.18;            %[~] from Devin (may trail off when off of MPP)
inso.deg = 0.5;             %[%/year]
%inso.pvci = 24;             %[months] cleaning interval
inso.wf = 30;               %[kg/m^2] weight factor
inso.debug = false;          %toggle debugging kW/kWh combo for shooter
inso.shootdebug = false;    %toggle debugging pvci shooter
inso.shoottol = 5;          %months
%inso.ct_eval = false;       %evaluate/compare trips for cleaning
inso.cleanstrat = 4;        %panel cleaning strategy 1:NC, 2:CT, 3:CTW
inso.cleanlim = 20;         %[mo] maximum limit for cleaning
%inso.nu = 1.01;             %[m/kW]
%wave energy parameters
wave.method = 2;            %1: divide by B, 2: 3d interpolation, 3: Trevor PM
wave.B_func_n = 1000;       %number of points in B(Gr) function
wave.Hs_ra = 4;             %[m], rated wave height
wave.Tp_ra = 9;            %[s], rated peak period
wave.eta_ct = 0.6;          %[~] wec efficiency
wave.house = 0.10;          %percent of rated power as house load
wave.kW_max = 17;           %[kW] maximum limit for wec-sim output
% wave.wsr = 'struct3m_opt';  %wec sim run
% wave.wsHs = 3;              %[m] wec sim Hs
%diesel parameters
dies.fmax = 800;            %[liters] fuel capacity
dies.ftmax = 18;            %[m] fuel can sit idle before going "bad"
%dies.lph = 2;               %[l/h]
dies.oilint = 250;          %[hours] maintenance interval

%dies.genon = 0.1;           %battery level generator turns on at -
%obselete because the generator now turns on when storage is too low to
%satisfy load

dies.kWmax = 8;            %maximum power generation
dies.kWmin = 1;             %minimum power generation
dies.bm = 4;                %barge multiplier
% %AGM parameters
% agm.V = 12;                %[V] Voltage
% agm.se = 3.3;              %[Ah/kg] specific energy factor
% agm.lc_nom = 18;           %[months] nominal life cycle
% agm.beta = 6/10;           %decay exponential for life cycle
% agm.lc_max = 12*5;        %maximum months of operation
% agm.sdr = 5;               %[%/month] self discharge rate
% %agm.dyn_lc = true;         %toggle dynamic life cycle
% agm.dmax = .2;             %maximum depth of discharge
%LFP parameters
lfp.V = 12;                 %[V] Voltage
lfp.se = 8.75;              %[Ah/kg] specific energy factor
lfp.lc_nom = 18;            %[months] nominal life cycle
lfp.beta = 1;               %decay exponential for life cycle
lfp.lc_max = 12*5;          %maximum months of operation
lfp.sdr = 3;                %[%/month] self discharge rate
%lfp.dyn_lc = true;         %toggle dynamic life cycle
lfp.dmax = 0;              %maximum depth of discharge
lfp.cost = 466;             %[$/kWh] - irena2020electricty
lfp.lcm = 1;%battery life cycle model, 1:bolun 2:dyn_lc 3:fixed_lc
lfp.T = 15;                 %[C] temperature
lfp.EoL = 0.2;              %battery end of life
lfp.rf_os = true;           %toggle using open source rainflow
lfp.bdi = 2190;              %battery degradation evaluation interaval
bc = 2; %battery chemistry 1:AGM 2:LFP
if bc == 1 %agm chemistry
    batt = agm;
elseif bc == 2 %lfp chemistry
    batt = lfp;
end

%atmospheric parameters
atmo.rho_a = 1.225;         %[kg/m^3] density of air
atmo.rho_w = 1025;          %[kg/m^3] density of water
atmo.g = 9.81;              %[m/s^2]
%atmo.h = 4;                 %[m]
atmo.zo = 0.2;             %[mm]
atmo.dyn_h = true;          %toggle dynamic hub height
atmo.soil = 35;             %[%/year]
%atmo.clean = 0.5;           %heavy rain cleans X amt of soil

%% USE CASES
%short term instrumentation
%uc(1).draw = 200;               %[W] - load now defined in optRun
uc(1).lifetime = 6;             %[y]
uc(1).loadcase = 3;             %1=HCUUV, 2=HFUUV, 3=OOUUV, 4=HF radar, 5 = 200W
uc(1).SI = 24;                   %[months] service interval
uc(1).uptime = .99;             %[%] uptime
%long term instrumentation
%uc(2).draw = 200;               %[W] - secondary node
uc(2).lifetime = 6;             %[y]
uc(2).loadcase = 3;             %1=HCUUV, 2=HFUUV, 3=OOUUV, 4=HF radar, 5 = 200W
uc(2).SI = 24;                  %[months] service interval
uc(2).uptime = .99;             %[%] uptime

%sensitivity analaysis
if ~isfield(opt,'tuning_array') && ~isfield(opt,'tuned_parameter')
% opt.tuning_array = [100 95 90 85 80 75 70];
% opt.tuned_parameter = 'wcp'; %wave cutout percentile
% opt.tuning_array = [1 2 3 4 5 6 7 8 9 10];
% opt.tuned_parameter = 'wcm'; %wave cost multiplier
% opt.tuning_array = [45 50 55 60 65 70 75 80 85 90];
% opt.tuned_parameter = 'wrp'; %wave rated percentile
% opt.tuning_array = linspace(.80,1,10);
% opt.tuned_parameter = 'utp';
% opt.tuning_array = [10:10:200];
% opt.tuned_parameter = 'load';
% opt.tuning_array = [0.01,0.2,.5];
% opt.tuned_parameter = 'zo';
% opt.tuning_array = [0,1,2,3,4,5];
% opt.tuned_parameter = 'utf';
% opt.tuning_array = [0 .01 .025 .05 .075 .1 .15 .2 .25];
% opt.tuned_parameter = 'whl'; %wec house load
% opt.tuning_array = [1 1.2 1.4 1.6 1.8 2];
% opt.tuned_parameter = 'imf'; %inso marinization factor
% opt.tuning_array = linspace(0.1,10,10);
% opt.tuned_parameter = 'btm'; %battery time slope
% opt.tuning_array = [10 20 30 40 50 60 70 80 90 100];
% opt.tuned_parameter = 'mbt'; %minimum battery for time added
% opt.tuning_array = linspace(1/2,2,10);
% opt.tuned_parameter = 'cwm'; %capture width multiplier
% opt.tuning_array = 1:1:10;
% opt.tuned_parameter = 'wcm'; %wave cost multiplier
% opt.tuning_array = linspace(0,9,50);
% opt.tuned_parameter = 'wiv'; %wec interventions
% opt.tuning_array = linspace(1/2,2,10);
% opt.tuned_parameter = 'dep'; %depth modifier
% opt.tuning_array = linspace(uc(c).lifetime-3,uc(c).lifetime+3,10);
% opt.tuned_parameter = 'lft'; %lifetime
    opt.tuning_array = linspace(10,1400,10)*1000;
    opt.tuned_parameter = 'dtc'; %distance to coast [OPEX]
end

%opt 2D sens
% opt.tdsens_ta(1,:) = 0.1:0.04:1.7;
% opt.tdsens_ta(2,:) = 40:2:120;
% opt.tdsens_tp{1} = 'btm'; %battery time slope
% opt.tdsens_tp{2} = 'mbt'; %minimum battery for time added
opt.tdsens_ta(1,:) = 2:1:7;
opt.tdsens_ta(2,:) = 6:1:11;
opt.tdsens_tp{1} = 'hra'; %rated Hs
opt.tdsens_tp{2} = 'tra'; %rated Tp

%optimization parameters
opt.V = 2;
opt.bf.M_hros = [2 4 5 1 1.75]; %[kW], high res os
opt.bf.N_hros = [350 500 500 300 350]; %[kWh], high res os
opt.bf.maxworkers = 36; %maximum cores
