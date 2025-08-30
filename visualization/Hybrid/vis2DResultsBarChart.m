%plot all the 2D data in bar chart form
clear
clc
close all

folder = "C:\Users\smpal\Documents\OO-Hybrid-Research\NAVFAC_ReportResults\";

task2loc = {'WETS','PISCES','PacWave','MidAtlSB','BerSea'};

pms = {'1','2','3','5'};

addpath('C:\Users\smpal\OneDrive\Documents\GitHub\WECFloatHydro\DrosteEffect-BrewerMap-b6a6efc\')
col{1} = brewermap(9,'Set3'); 
cols = [col{1}(1,:); col{1}(2,:); col{1}(5,:); col{1}(3,:)];
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

%% plot the minimum cost
loadcase3 = [1,3,5];
figure
tiledlayout(max(size(task2loc)),3)
for l = 1:max(size(task2loc))
    for lc = 1:3
        LL = loadcase3(lc);
        nexttile
        y = nan(max(size(pms)),2);
        for p = 1:max(size(pms))
            x(p) = str2double(pms{p});
            if isempty(optStructFFA{l,p}) == 0
                y(p,1) = optStructFFA{l,p}(lc).output.min.cost/1000;
            else
                y(p,1) = nan;
            end
            if isempty(optStructBF{l,p}) == 0
                if max(size(optStructBF{l,p})) == 5
                    y(p,2) = optStructBF{l,p}(LL).output.min.cost/1000;
                else
                    y(p,2) = optStructBF{l,p}(lc).output.min.cost/1000;
                end
            else
                y(p,2) = nan;
            end
        end
        b = bar(x,y);
        b(1).FaceColor = 'flat'; b(2).FaceColor = 'flat';
        b(1).CData = cols; b(2).CData = cols;

        xtips2 = b(1).XEndPoints;
        ytips2 = b(1).YEndPoints;
        labels2 = string(round(b(1).YData,2));
        text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')

        xtips2 = b(2).XEndPoints;
        ytips2 = b(2).YEndPoints;
        labels2 = string(round(b(2).YData,2));
        text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')

        ylabel('Cost [$1000]')
        xlabel('Power Module')
        %title(strcat(task2loc{l}," Load Case: ",string(optStruct{l,p}(lc).uc.loadcase)))
        ylim([0,2000])
    end
end

%Plot min battery size
figure
tiledlayout(max(size(task2loc)),3)

for l = 1:max(size(task2loc))
    for lc = 1:3
        nexttile
        for p = 1:max(size(pms))
            x(p) = str2double(pms{p});
            y(p) = optStruct{l,p}(lc).output.min.Smax{1};
        end
        b = bar(x,y);
        b.FaceColor = 'flat';
        b.CData = cols;
        xtips2 = b.XEndPoints;
        ytips2 = b.YEndPoints;
        labels2 = string(round(b.YData,2));
        text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
        ylabel('S_max [kWh]')
        xlabel('Power Module')
        title(strcat(task2loc{l}," Load Case: ",string(optStruct{l,p}(lc).uc.loadcase)))
        ylim([0,500])
    end
end

%Plot min generator capacity
figure
tiledlayout(max(size(task2loc)),3)

for l = 1:max(size(task2loc))
    for lc = 1:3
        nexttile
        for p = 1:max(size(pms))
            x(p) = str2double(pms{p});
            if x(p) == 1
                y(p) = optStruct{l,p}(lc).output.min.kWwi{1};
            elseif x(p) == 2
                y(p) = optStruct{l,p}(lc).output.min.kWi{1};
            elseif x(p) == 3
                y(p) = optStruct{l,p}(lc).output.min.kWwa{1};
            elseif x(p) == 5
                y(p) = optStruct{l,p}(lc).output.min.kWc{1};
            end
        end
        b = bar(x,y);
        b.FaceColor = 'flat';
        b.CData = cols;
        xtips2 = b.XEndPoints;
        ytips2 = b.YEndPoints;
        labels2 = string(round(b.YData,2));
        text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
        ylabel('[kW]')
        xlabel('Power Module')
        title(strcat(task2loc{l}," Load Case: ",string(optStruct{l,p}(lc).uc.loadcase)))
        ylim([0,10])
    end
end
