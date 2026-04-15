function [multStruct] = doFailSurvsens(prepath,name,batchtype,Aparams)
optInputs %load inputs
data = load(loc,'data');
data = data.('data');


namepts = split(name,'_');
arrayID = namepts{4};
date = namepts{5};
opt.jobstore = strcat("/gscratch/scrubbed/smp52/local_cluster_jobs/Job",arrayID,date);

failopts = [4E8 4E7 4E6 2];
if ~isempty(Aparams)
    opt.failsurv = failopts(Aparams(1));
end

%initalize outputs
clear multStruct
multStruct = struct();

tTot = tic;
disp('Next point in failsurv starting now')


[multStruct.output,multStruct.opt] = ...
optRun(opt,data,atmo,batt,econ,uc(c),bc,inso,turb,cturb, wave,dies);

multStruct.data = data;
multStruct.atmo = atmo;
multStruct.batt = batt;
multStruct.econ = econ;
multStruct.uc = uc(c);
multStruct.c = c;
multStruct.loc = loc;
multStruct.turb = turb;
multStruct.cturb = cturb;
multStruct.inso = inso;
multStruct.wave = wave;
multStruct.dies = dies;
multStruct.comptime = num2str(round(toc(tTot),2)); %seconds
    
stru.(name) = multStruct;
nameext = strcat(prepath, name,'.mat');
save(nameext, '-struct','stru','-v7.3')

end

