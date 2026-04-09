function [allLUP] = doAllLocUsesPMs(prepath,name,batchtype,Aparams)
tTot = tic;

namepts = split(name,'_');
arrayID = namepts{4};
date = namepts{5};
opt.jobstore = strcat("/gscratch/scrubbed/local_cluster_jobs/Job",arrayID,date);
%parse inputs from array job
if ~isempty(Aparams)
    Aloc = Aparams(1);
    Apm = Aparams(3);
    Aload = Aparams(2);
else
    Aloc = 1:1:5;
    Apm = 1:1:5;
    Aload = [1 3 5];
end

locoptions = {'PacWave','MidAtlSB','BerSea','altWETS','altPISCES'};

%initialize outputs
clear allLUP
allLUP(length(Aloc),length(Apm),length(Aload)) = struct();

for ll = 1:length(Aloc)
    optInputs %load inputs
    loc = locoptions{Aloc(ll)}; %reset the location
    data = load(loc,'data');
    data = data.('data');
    for pp = 1:length(Apm)
        for uu = 1:length(Aload)
            opt.pm = Apm(pp); %Re-set power module
            uc(c).loadcase = Aload(uu); %reset the load case

            disp(['Optimization at ' char(loc) ...
                ' for power module ' ...
                string(Apm(pp)) ' beginning after ' ...
                num2str(round(toc(tTot),2)) ' seconds.'])

            [allLUP(ll,pp,uu).output, ...
                allLUP(ll,pp,uu).opt] = optRun(opt,data,atmo,batt,econ,uc(c),bc,inso,turb,cturb, wave,dies);
            allLUP(ll,pp,uu).data = data;
            allLUP(ll,pp,uu).atmo = atmo;
            allLUP(ll,pp,uu).batt = batt;
            allLUP(ll,pp,uu).econ = econ;
            allLUP(ll,pp,uu).uc = uc(c);
            allLUP(ll,pp,uu).loc = loc;
            allLUP(ll,pp,uu).turb = turb;
            allLUP(ll,pp,uu).inso = inso;
            allLUP(ll,pp,uu).wave = wave;
            allLUP(ll,pp,uu).dies = dies;
            allLUP(ll,pp,uu).cturb = cturb;
            allLUP(ll,pp,uu).comptime = num2str(round(toc(tTot),2)); %seconds
        
            %Periodic Save so in case the simulation is longer than the wall time
            stru.(name) = allLUP;
            nameext = strcat(prepath, name,'.mat');
            save(nameext, '-struct','stru','-v7.3')
        
            disp(['Optimization at ' char(char(loc)) ...
                ' for ' ...
                string(Apm(pp)) ' application complete after ' ...
                num2str(round(toc(tTot),2)) ' seconds.'])
        end
    end
end
end
