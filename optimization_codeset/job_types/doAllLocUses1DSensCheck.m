function [allLUP] = doAllLocUses1DSensCheck(prepath,name,batchtype,Aparams)
tTot = tic;
%%This should only be run with array job types

namepts = split(name,'_');
arrayID = namepts{4};
date = namepts{5};
opt.jobstore = strcat("/gscratch/scrubbed/local_cluster_jobs/Job",arrayID,date);
%parse inputs from array job
if ~isempty(Aparams)
    Aloc = Aparams(1); %location
    Aload = Aparams(2); %load case
    Asens = Aparams(3); %sens var
    Aval = Aparams(4); %sens val

    opt.ffa.sens{1} = Aparams(5);
    opt.ffa.sens{2} = Aparams(6);
    opt.ffa.sens{3} = Aparams(7);
    opt.ffa.sens{4} = Aparams(8);
    opt.ffa.sens{5} = Aparams(9);
    opt.ffa.sens{6} = Aparams(10);
end

locoptions = {'PacWave','MidAtlSB','BerSea','altWETS','altPISCES'};

%sensitivity variables to check
%'lfp.se','econ.wave.mass_mult','econ.platform.steel','lfp.cost','econ.refurb_mult','econ.dies.fcost','econ.wind.marinization','econ.wave.costmult_con','econ.curr.costmult'
sensmin = [1/4, 10, 1/10,1/10,1/2,1/5,1/2,1/10,1/3];
sensmax = [10, 1/10, 10,10,10,10,10,10,10];


%initialize outputs
clear allLUP
allLUP(1,1) = struct();

loc = locoptions{Aloc(ll)}; %reset the location
data = load(loc,'data');
data = data.('data');
optInputs %load inputs
uc(c).loadcase = Aload(uu); %reset the load case

%overwrite array parameters
opt.ffa.max = opt.ffa.sens{1};
opt.ffa.pop = opt.ffa.sens{2}; 
opt.ffa.gamma = opt.ffa.sens{3};
opt.ffa.beta0 = opt.ffa.sens{4};
opt.ffa.alpha = opt.ffa.sens{5};
opt.ffa.adamp = opt.ffa.sens{6};

if Aval == 1
    mult = sensmin;
else
    mult = sensmax;
end

switch Asens
    case 1
        disp('changing lfp.se')
        lfp.se = lfp.se * mult(Asens)
    case 2
        disp('changing econ.wave.mass_mult')
        econ.wave.mass_mult =  econ.wave.mass_mult * mult(Asens)
    case 3
        disp('changing econ.platform.steel')
        econ.platform.steel = econ.platform.steel * mult(Asens)
    case 4
        disp('changing lfp.cost')
        lfp.cost = lfp.cost * mult(Asens)
    case 5
        disp('changing econ.refurb_mult')
        econ.refurb_mult = econ.refurb_mult * mult(Asens)
    case 6
        disp('changing econ.dies.fcost')
        econ.dies.fcost = econ.dies.fcost * mult(Asens)
    case 7
        disp('changing econ.wind.marinization')
        econ.wind.marinization = econ.wind.marinization * mult(Asens)
    case 8
        disp('changing econ.wave.costmult_con')
        econ.wave.costmult_con = econ.wave.costmult_con * mult(Asens)
    case 9
        disp('changing econ.curr.costmult')
        econ.curr.costmult = econ.curr.costmult * mult(Asens)
end

disp(['Optimization at ' char(loc) ...
    ' for power module Hybrid beginning after ' ...
    num2str(round(toc(tTot),2)) ' seconds.'])

[allLUP(1).output, ...
    allLUP(1).opt] = optRun(opt,data,atmo,batt,econ,uc(c),bc,inso,turb,cturb, wave,dies);
allLUP(1).data = data;
allLUP(1).atmo = atmo;
allLUP(1).batt = batt;
allLUP(1).econ = econ;
allLUP(1).uc = uc(c);
allLUP(1).loc = loc;
allLUP(1).turb = turb;
allLUP(1).inso = inso;
allLUP(1).wave = wave;
allLUP(1).dies = dies;
allLUP(1).cturb = cturb;
allLUP(1).comptime = num2str(round(toc(tTot),2)); %seconds
    
%Periodic Save so in case the simulation is longer than the wall time
stru.(name) = allLUP;
nameext = strcat(prepath, name,'.mat');
save(nameext, '-struct','stru','-v7.3')

disp(['Optimization at ' char(char(loc)) ...
    ' for Hybrid application complete after ' ...
    num2str(round(toc(tTot),2)) ' seconds.'])

end
