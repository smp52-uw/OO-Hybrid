%organize the data for tables for the NAVFAC Report
clear
clc
close all

folder = "C:\Users\smpal\Documents\OO-Hybrid-Research\NAVFAC_ReportResults\";

task2loc = {'MidAtlSB','BerSea','PacWave','PISCES','WETS'};

pms = {'1','2','3','5'};

addpath('C:\Users\smpal\OneDrive\Documents\GitHub\WECFloatHydro\DrosteEffect-BrewerMap-b6a6efc\')

lc = 1;

for l = 1:max(size(task2loc))
    for p = 1:max(size(pms))
        filename = strcat(task2loc{l},'_PD2PM',pms{p},'_AllLC_*_*_*.mat');
        filename = dir(strcat(folder, filename));
        if max(size(filename)) == 2
            xn = split(filename(1).name,"_");
            temp = load(fullfile(folder,filename(1).name));
            if any(strcmp(xn,'FFAOpt'))
                disp('success')
                temp = load(fullfile(folder,filename(1).name));
                fn = fieldnames(temp);
                optStructFFA{l,p} = temp.(fn{1});

                temp = load(fullfile(folder,filename(2).name));
                fn = fieldnames(temp);
                optStructBF{l,p} = temp.(fn{1});

            elseif any(strcmp(xn,'BF500'))
                temp = load(fullfile(folder,filename(1).name));
                fn = fieldnames(temp);
                optStructBF{l,p} = temp.(fn{1});

                temp = load(fullfile(folder,filename(2).name));
                fn = fieldnames(temp);
                optStructFFA{l,p} = temp.(fn{1});

            end
        else
            xn = split(filename(1).name, "_");
            if sum(strcmp(xn,"FFAOpt")) == 1
                temp = load(fullfile(folder,filename(1).name));
                fn = fieldnames(temp);
                optStructFFA{l,p} = temp.(fn{1});
            else
                temp = load(fullfile(folder,filename(1).name));
                fn = fieldnames(temp);
                optStructBF{l,p} = temp.(fn{1});
            end

        end
    end
end
%%
format longg
costtab = cell(6,5);
gentab = cell(6,5);
stortab = cell(6,5);
gen = {"kWwi","kWi","kWwa","kWc"};
loadcases = [1, 3, 5];
c = 2;

for l = 1:max(size(task2loc))
    costtab{l+1,1} = task2loc{l};
    for p = 1:max(size(pms))
        costtab{1,p+1} = pms{p};
        if ~isempty(optStructFFA{l,p})
            if optStructFFA{l,p}(c).uc.loadcase == loadcases(c)
                costtab{l+1,p+1} = round(optStructFFA{l,p}(c).output.min.cost./1000,2);
                gentab{l+1,p+1} = round(optStructFFA{l,p}(c).output.min.(gen{p}){1},2);
                stortab{l+1,p+1} = round(optStructFFA{l,p}(c).output.min.Smax{1},2);
            else
                cc = loadcases(c)
                costtab{l+1,p+1} = round(optStructFFA{l,p}(cc).output.min.cost./1000,2);
                gentab{l+1,p+1} = round(optStructFFA{l,p}(cc).output.min.(gen{p}){1},2);
                stortab{l+1,p+1} = round(optStructFFA{l,p}(cc).output.min.Smax{1},2);
            end

        else
            if optStructBF{l,p}(c).uc.loadcase == loadcases(c)
                costtab{l+1,p+1} = round(optStructBF{l,p}(c).output.min.cost./1000,2);
                gentab{l+1,p+1} = round(optStructBF{l,p}(c).output.min.(gen{p}){1},2);
                stortab{l+1,p+1} = round(optStructBF{l,p}(c).output.min.Smax{1},2);
            else
                cc = loadcases(c)
                costtab{l+1,p+1} = round(optStructBF{l,p}(cc).output.min.cost./1000,2);
                gentab{l+1,p+1} = round(optStructBF{l,p}(cc).output.min.(gen{p}){1},2);
                stortab{l+1,p+1} = round(optStructBF{l,p}(cc).output.min.Smax{1},2);
            end

        end
    end
end