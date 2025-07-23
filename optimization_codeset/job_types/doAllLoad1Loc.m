function [allLoad1Loc] = doAllLoad1Loc(prepath,name,batchtype)
tTot = tic;
optInputs %load inputs

data = load(loc,'data');
data = data.('data');

%initialize outputs
clear allLoad1Loc
allLoad1Loc(1,5) = struct();

for LoadC = 1:5
    uc(c).loadcase = LoadC; %Re-set load case

    disp(['Optimization at ' char(loc) ...
        ' for load case ' ...
        char(LoadC) ' beginning after ' ...
        num2str(round(toc(tTot),2)) ' seconds.'])
    [allLoad1Loc(1,LoadC).output, ...
        allLoad1Loc(1,LoadC).opt] = optRun(opt,data,atmo,batt,econ,uc(c),bc,inso,turb,cturb, wave,dies);
    allLoad1Loc(1,LoadC).data = data;
    allLoad1Loc(1,LoadC).atmo = atmo;
    allLoad1Loc(1,LoadC).batt = batt;
    allLoad1Loc(1,LoadC).econ = econ;
    allLoad1Loc(1,LoadC).uc = uc(c);
    allLoad1Loc(1,LoadC).loc = loc;
    allLoad1Loc(1,LoadC).turb = turb;
    allLoad1Loc(1,LoadC).inso = inso;
    allLoad1Loc(1,LoadC).wave = wave;
    allLoad1Loc(1,LoadC).dies = dies;
    allLoad1Loc(1,LoadC).cturb = cturb;
    allLoad1Loc(1,LoadC).comptime = num2str(round(toc(tTot),2)); %seconds

    %Periodic Save so in case the simulation is longer than the wall time
    stru.(name) = allLoad1Loc;
    save([prepath name '.mat'], '-struct','stru','-v7.3')

    disp(['Optimization at ' char(opt.locations(1)) ...
        ' for ' ...
        char(LoadC) ' application complete after ' ...
        num2str(round(toc(tTot),2)) ' seconds.'])
end
end

