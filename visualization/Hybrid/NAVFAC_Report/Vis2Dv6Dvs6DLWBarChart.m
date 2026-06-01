%Compare the 2D, 6D and 6D low WEC cost/mass results for Bering Sea
close all
clear

set(0,'defaulttextinterpreter','none')
%set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'calibri')
set(0,'DefaultAxesFontName', 'calibri')

%% merge structures 
folder = "C:\Users\smpal\Documents\OO-Hybrid-Research\NAVFAC_ReportResults\";
task2loc = {'WETS','PISCES','PacWave','MidAtlSB','BerSea'};
pmtext = {"Wind","Solar","Wave","Current"};
locationI = 5;
pms = {'1','2','3','5'};
%pull 2D data
for p = 1:max(size(pms))
    filename = strcat(task2loc{locationI},'_PD2PM',pms{p},'_AllLC_*_*_*.mat');
    filename = dir(strcat(folder, filename));
    if max(size(filename)) == 2 %if the pm was run twice then I only want the FFA Opt solution
        xn = split(filename(1).name,"_");
        temp = load(fullfile(folder,filename(1).name));
        if any(strcmp(xn,'FFAOpt'))
            temp = load(fullfile(folder,filename(1).name));
            fn = fieldnames(temp);
            optStructFFA{p} = temp.(fn{1});

        elseif any(strcmp(xn,'BF500'))
            temp = load(fullfile(folder,filename(2).name));
            fn = fieldnames(temp);
            optStructFFA{p} = temp.(fn{1});

        end

        if max(size(optStructFFA{p})) == 3
            mincost(1,p) = optStructFFA{p}(1).output.min.cost;
            mincost(2,p) = optStructFFA{p}(2).output.min.cost;
            mincost(3,p) = optStructFFA{p}(3).output.min.cost;
        else
            mincost(1,p) = optStructFFA{p}(1).output.min.cost;
            mincost(2,p) = optStructFFA{p}(3).output.min.cost;
            mincost(3,p) = optStructFFA{p}(5).output.min.cost;
        end

    else 
        xn = split(filename(1).name, "_");
        if sum(strcmp(xn,"BF500")) == 1 %if this is MidATLSB then only BF was run but I'm gonna get all of the pms
            temp = load(fullfile(folder,filename(1).name));
            fn = fieldnames(temp);
            optStructBF{p} = temp.(fn{1});
            if max(size(optStructBF{p})) == 3
                mincost(1,p) = optStructBF{p}(1).output.min.cost;
                mincost(2,p) = optStructBF{p}(2).output.min.cost;
                mincost(3,p) = optStructBF{p}(3).output.min.cost;
            else
                mincost(1,p) = optStructBF{p}(1).output.min.cost;
                mincost(2,p) = optStructBF{p}(3).output.min.cost;
                mincost(3,p) = optStructBF{p}(5).output.min.cost;
            end
        else
            temp = load(fullfile(folder,filename(1).name));
            fn = fieldnames(temp);
            optStructFFA{p} = temp.(fn{1});
            if max(size(optStructFFA{p})) == 3
                mincost(1,p) = optStructFFA{p}(1).output.min.cost;
                mincost(2,p) = optStructFFA{p}(2).output.min.cost;
                mincost(3,p) = optStructFFA{p}(3).output.min.cost;
            else
                mincost(1,p) = optStructFFA{p}(1).output.min.cost;
                mincost(2,p) = optStructFFA{p}(3).output.min.cost;
                mincost(3,p) = optStructFFA{p}(5).output.min.cost;
            end
        end

    end
end

%if we have the BF structure we need to find the opt2D
if exist('optStructBF','var') == 1
    [~,indp] = min(mincost,[],2);
    optStruct2D(1,1) = optStructBF{indp(1)}(1);
    optStruct2D(1,2) = optStructBF{indp(2)}(2);
    optStruct2D(1,3) = optStructBF{indp(3)}(3);
else
    [~,indp] = min(mincost,[],2);
    optStruct2D(1,1) = optStructFFA{indp(1)}(1);
    optStruct2D(1,2) = optStructFFA{indp(2)}(2);
    optStruct2D(1,3) = optStructFFA{indp(3)}(3);
end

%pull hybrid data
filename = strcat(task2loc{locationI},'_PD6','_*.mat');
%BerSea_PD6_AllLC_FFAOpt_Low1010MassCostWEC_09302025
filename = dir(strcat(folder, filename));
for j = 1:max(size(filename))
    xn = split(filename(j).name,"_");
    if sum(strcmp(xn,"Low1010MassCostWEC")) == 1
        temp = load(fullfile(folder,filename(j).name));
        fn = fieldnames(temp);
        optStruct6DLW = temp.(fn{1});
    else %This is the normal 6D Run
        temp = load(fullfile(folder,filename(j).name));
        fn = fieldnames(temp);

        if length(size(temp.(fn{1}))) == 3
            tempopt = temp.(fn{1})(2,1,1);
        else
            tempopt = temp.(fn{1})(2,1,1,2);
        end

        if sum(strcmp(xn,"LC1")) == 1
            optStruct6D(1,1) = tempopt;
        elseif sum(strcmp(xn,"LC3")) == 1
            optStruct6D(1,2) = tempopt;
        else
            optStruct6D(1,3) = tempopt;
        end

    end
end


%%

%% plotting setup
results = figure;
set(gcf,'Units','inches')
set(gcf, 'Position', [1, 1, 7, 5])

tfig = tiledlayout(1,3);
tfig.Padding = 'compact';
title(tfig,strcat("Cost Comparison: ",task2loc{locationI}))

%colors
addpath(genpath("C:\Users\smpal\Documents\cmasher"))
%wind = green
%solar = yellow
%wave = blue
%current = purple
%hybrid = pink
col{1} = cmasherImport('jungle', 50);
col{2} = cmasherImport('amber', 70);
col{3} = cmasherImport('arctic', 50);
col{4} = cmasherImport('amethyst', 50);
col{5} = cmasherImport('flamingo', 50);
col{6} = cmasherImport('gothic', 50);

col2D{1} = [col{1}(25,:) ; col{1}(30,:); col{1}(35,:); col{1}(40,:)];
col2D{2} = [col{2}(55,:) ; col{2}(60,:); col{2}(65,:); col{2}(70,:)];
col2D{3} = [col{3}(25,:) ; col{3}(30,:); col{3}(35,:); col{3}(40,:)];
col2D{4} = [col{4}(25,:) ; col{4}(30,:); col{4}(35,:); col{4}(40,:)];
col6D = [col{5}(30,:) ; col{5}(35,:); col{5}(40,:); col{5}(45,:)];
col6DLW = [col{6}(30,:) ; col{6}(35,:); col{6}(40,:); col{6}(45,:)];

%bar chart settings
MaxGroupWidth = 0.75;
groupOffset = MaxGroupWidth;

for lcind = 1:3
    %unpack cost variables into costdata
    %Cost Types:
    %Gen CapEx
    %Battery CapEx
    %Mooring Cost CapEx
    %Repair OpEx
    %Vessel Cost (for installation and repairs)
    %2D
    SG = optStruct2D(lcind).output.min;

    al_cost = 103.92; %$/sqft for 1/2" 6061 AL
    al_cost = al_cost*0.0929; %convert to $/m2
    batt_vol = polyval(optStruct2D(lcind).opt.p_dev.b_size,SG.Smax{1})*optStruct2D(lcind).econ.batt.volmult;
    batt_len = (batt_vol)^(1/3); %assume cube
    battencl = 2*6*(batt_len^2)*al_cost;
    SG.Icost_batt = 0.1*(SG.Scost+battencl); %installation cost (10% is a guess from Brian)

    costdata(1,1) = SG.CapEx - (SG.Pmooring + SG.Pinst + SG.Pmtrl + SG.battencl + SG.Scost + SG.Icost_batt);
    costdata(2,1) = SG.battencl + SG.Scost + SG.Icost_batt; %
    costdata(3,1) = SG.Pmooring + SG.Pmtrl;
    costdata(4,1) = SG.OpEx - SG.vesselcost;
    %costdata(5,1) = SG.vesselcost + SG.Pinst;
    
    %6D
    HG = optStruct6D(lcind).output.min;

    batt_vol = polyval(optStruct6D(lcind).opt.p_dev.b_size,HG.Smax{1})*optStruct6D(lcind).econ.batt.volmult;
    batt_len = (batt_vol)^(1/3); %assume cube
    battencl = 2*6*(batt_len^2)*al_cost;
    HG.Icost_batt = 0.1*(HG.Scost+battencl); %installation cost (10% is a guess from Brian)

    costdata(1,2) = HG.CapEx - (HG.Pmooring + HG.Pinst + HG.Pmtrl + HG.battencl + HG.Scost + HG.Icost_batt); %
    costdata(2,2) = HG.battencl + HG.Scost + HG.Icost_batt; %
    costdata(3,2) = HG.Pmooring + HG.Pmtrl;
    costdata(4,2) = HG.OpEx - HG.vesselcost;
    %costdata(5,2) = HG.vesselcost + HG.Pinst;

    %6D - Low WEC mass and cost
    HGLW = optStruct6DLW(lcind).output.min;

    batt_vol = polyval(optStruct6DLW(lcind).opt.p_dev.b_size,HGLW.Smax{1})*optStruct6DLW(lcind).econ.batt.volmult;
    batt_len = (batt_vol)^(1/3); %assume cube
    battencl = 2*6*(batt_len^2)*al_cost;
    HGLW.Icost_batt = 0.1*(HGLW.Scost+battencl); %installation cost (10% is a guess from Brian)

    costdata(1,3) = HGLW.CapEx - (HGLW.Pmooring + HGLW.Pinst + HGLW.Pmtrl + HGLW.battencl + HGLW.Scost + HGLW.Icost_batt); %
    costdata(2,3) = HGLW.battencl + HGLW.Scost + HGLW.Icost_batt; %
    costdata(3,3) = HGLW.Pmooring + HGLW.Pmtrl;
    costdata(4,3) = HGLW.OpEx - HGLW.vesselcost;
    %costdata(5,2) = HG.vesselcost + HG.Pinst;
    

    
    %convert all costs to 1000 $
    costdata = costdata ./1000;
    %plot
    
    ax(lcind) = nexttile();
    hold on
    
    X = categorical({'Generation CapEx','Storage CapEx','Mooring CapEx','Repairs OpEx'});
    X = reordercats(X,{'Generation CapEx','Storage CapEx','Mooring CapEx','Repairs OpEx'});
    h = bar(X, costdata, 'FaceColor','flat');
    set(h,'BarWidth',groupOffset);

    h(1).CData = col2D{optStruct2D(lcind).opt.pm};
    h(2).CData = col6D;
    h(3).CData = col6DLW;

        
    hold off;
    grid on
    title(strcat("Load Case: ",string(optStruct2D(lcind).uc.loadcase)))

    if lcind == 1
        ylabel('Cost [$1000]')
    end

    if lcind == 2
        sgleg = strcat(pmtext{optStruct2D(lcind).opt.pm}," + Battery");
        legend(sgleg, "Hybrid","Hybrid - Low WEC",'location','northoutside','NumColumns',3)
    end
end


