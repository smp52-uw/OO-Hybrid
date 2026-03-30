function [allPM1Loc] = doAllPM1Loc(prepath,name,batchtype,Aparams)
tTot = tic;

%parse inputs from array job
if ~isempty(Aparams)
    Aloc = Aparams(1);
    Apm = Aparams(2);
    Aload = Aparams(3);
end

locoptions = {'PacWave','MidAtlSB','BerSea','altWETS','altPISCES'};

%initialize outputs
clear allPM1Loc
allPM1Loc(1,length(task2pm)) = struct();

for ll = 1:length(Aloc)
    optInputs %load inputs
    loc = locoptions{Aloc(ll)}; %reset the location
    data = load(loc,'data');
    data = data.('data');

    for pp = 1:length(Apm)
        for uu = l:length(Aload)
            opt.pm = Apm(pp); %Re-set power module
            uc(c).loadcase = Aload(uu); %reset the load case

            disp(['Optimization at ' char(loc) ...
                ' for power module ' ...
                string(Apm(pp)) ' beginning after ' ...
                num2str(round(toc(tTot),2)) ' seconds.'])
            [allPM1Loc(1,pp).output, ...
                allPM1Loc(1,pp).opt] = optRun(opt,data,atmo,batt,econ,uc(c),bc,inso,turb,cturb, wave,dies);
            allPM1Loc(1,pp).data = data;
            allPM1Loc(1,pp).atmo = atmo;
            allPM1Loc(1,pp).batt = batt;
            allPM1Loc(1,pp).econ = econ;
            allPM1Loc(1,pp).uc = uc(c);
            allPM1Loc(1,pp).loc = loc;
            allPM1Loc(1,pp).turb = turb;
            allPM1Loc(1,pp).inso = inso;
            allPM1Loc(1,pp).wave = wave;
            allPM1Loc(1,pp).dies = dies;
            allPM1Loc(1,pp).cturb = cturb;
            allPM1Loc(1,pp).comptime = num2str(round(toc(tTot),2)); %seconds
        
            %Periodic Save so in case the simulation is longer than the wall time
            stru.(name) = allPM1Loc;
            save([prepath name '.mat'], '-struct','stru','-v7.3')
        
            disp(['Optimization at ' char(char(loc)) ...
                ' for ' ...
                string(Apm(pp)) ' application complete after ' ...
                num2str(round(toc(tTot),2)) ' seconds.'])
        end
    end
end
end
