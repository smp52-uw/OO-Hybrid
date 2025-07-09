function [multStruct] = doFFAsens(prepath,name,batchtype)
optInputs %load inputs
data = load(loc,'data');
data = data.('data');

opt.S1 = length(opt.ffa.sens{1});
opt.S2 = length(opt.ffa.sens{2});
opt.S3 = length(opt.ffa.sens{3});
opt.S4 = length(opt.ffa.sens{4});
opt.S5 = length(opt.ffa.sens{5});
opt.S6 = length(opt.ffa.sens{6});

%initalize outputs
clear multStruct
multStruct(opt.S1,opt.S2,opt.S3,opt.S4,opt.S5,opt.S6) = struct();
for p = 1:opt.S1
    for q = 1:opt.S2
        for r = 1:opt.S3
            for s = 1:opt.S4
                for u = 1:opt.S5
                    for v = 1:opt.S6

                        opt.ffa.max = opt.ffa.sens{1}(p);
                        opt.ffa.pop = opt.ffa.sens{2}(q); 
                        opt.ffa.gamma = opt.ffa.sens{3}(r);
                        opt.ffa.beta0 = opt.ffa.sens{4}(s);
                        opt.ffa.alpha = opt.ffa.sens{5}(u);
                        opt.ffa.adamp = opt.ffa.sens{6}(v);
                
                        tTot = tic;
                        disp('Next point in FFA sensitivity starting now')
                
                        [multStruct(p,q,r,s,u,v).output,multStruct(p,q,r,s,u,v).opt] = ...
                        optRun(opt,data,atmo,batt,econ,uc(c),bc,inso,turb,cturb, wave,dies);
                
                        multStruct(p,q,r,s,u,v).data = data;
                        multStruct(p,q,r,s,u,v).atmo = atmo;
                        multStruct(p,q,r,s,u,v).batt = batt;
                        multStruct(p,q,r,s,u,v).econ = econ;
                        multStruct(p,q,r,s,u,v).uc = uc(c);
                        multStruct(p,q,r,s,u,v).c = c;
                        multStruct(p,q,r,s,u,v).loc = loc;
                        multStruct(p,q,r,s,u,v).turb = turb;
                        multStruct(p,q,r,s,u,v).cturb = cturb;
                        multStruct(p,q,r,s,u,v).inso = inso;
                        multStruct(p,q,r,s,u,v).wave = wave;
                        multStruct(p,q,r,s,u,v).dies = dies;
                        multStruct(p,q,r,s,u,v).comptime = num2str(round(toc(tTot),2)); %seconds

                        stru.(name) = multStruct;
                        save([prepath name '.mat'], '-struct','stru','-v7.3')
                    end
                end
            end
        end
    end
end
end

