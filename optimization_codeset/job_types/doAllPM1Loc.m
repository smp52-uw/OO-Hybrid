function [allPM1Loc] = doAllPM1Loc(prepath,name,batchtype)
tTot = tic;
optInputs %load inputs

data = load(loc,'data');
data = data.('data');

%PM
task2pm = 1:1:5;

%initialize outputs
clear allPM1Loc
allPM1Loc(1,length(task2pm)) = struct();

for PowMod = 1:length(task2pm)
    opt.pm = task2pm(PowMod); %Re-set load case

    disp(['Optimization at ' char(loc) ...
        ' for power module ' ...
        string(task2pm(PowMod)) ' beginning after ' ...
        num2str(round(toc(tTot),2)) ' seconds.'])
    [allPM1Loc(1,PowMod).output, ...
        allPM1Loc(1,PowMod).opt] = optRun(opt,data,atmo,batt,econ,uc(c),bc,inso,turb,cturb, wave,dies);
    allPM1Loc(1,PowMod).data = data;
    allPM1Loc(1,PowMod).atmo = atmo;
    allPM1Loc(1,PowMod).batt = batt;
    allPM1Loc(1,PowMod).econ = econ;
    allPM1Loc(1,PowMod).uc = uc(c);
    allPM1Loc(1,PowMod).loc = loc;
    allPM1Loc(1,PowMod).turb = turb;
    allPM1Loc(1,PowMod).inso = inso;
    allPM1Loc(1,PowMod).wave = wave;
    allPM1Loc(1,PowMod).dies = dies;
    allPM1Loc(1,PowMod).cturb = cturb;
    allPM1Loc(1,PowMod).comptime = num2str(round(toc(tTot),2)); %seconds

    %Periodic Save so in case the simulation is longer than the wall time
    stru.(name) = allPM1Loc;
    save([prepath name '.mat'], '-struct','stru','-v7.3')

    disp(['Optimization at ' char(opt.locations(1)) ...
        ' for ' ...
        string(task2pm(PowMod)) ' application complete after ' ...
        num2str(round(toc(tTot),2)) ' seconds.'])
end
end
