%load in hybrid simulations
%set the folder with the sensitivity results
%selectedfolder = "C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\AllLUP6D";
%selectedfolder = "C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\AllLUP6D\LC3_NoBatteryRule"; %Adjusted LC 3 with no battery rule results
selectedfolder = "C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\All_MidAtlSBE_UP\6D"';
fileList = dir(fullfile(selectedfolder,'*.mat'));
numfiles = max(size(fileList));

%% load 2D results
folder2D = "C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\AllLUP";
min2D= load(fullfile(folder2D,"AllLUP_minSys.mat"));
minSys = min2D.minSys;
clear min2D

%%Set up locations and load cases
locoptions = {'PacWave','MidAtlSB','BerSea','altWETS','altPISCES','MidAtlSB_E'};
locdisp =  {{'Coastal','Oregon'},{'Mid-Atlantic',' Shelf Break'},{'Bering',' Sea'},{'Coastal','O''ahu'},{'Coastal', 'Washington'}};
loadoptions = [1 3 5];

%% either load the simplified mat file or load all files
results = cell(6,3,5);
try
    load(fullfile(selectedfolder,'SimplifiedResults.mat'))
    savefile = 0;
    numfiles = numfiles - 1;
catch
    disp('No simplified results - loading all files')
    i = 1;
    for j = 1:numfiles
        try
            tmp = load(fullfile(selectedfolder,fileList(j).name));

        catch
            warning(strcat('Problem loading file:',string(fileList(j).name)));
            tmp = [];
        end
        if ~isempty(tmp)
            nm = split(fileList(j).name,'.');
            nm = nm(1);
            loc{i} = tmp.(nm{1}).loc;
            lc(i) = tmp.(nm{1}).uc.loadcase;

            %indexing
            li = find(strcmp(loc{i},locoptions));
            ui = find(lc(i) == loadoptions);
            pi = tmp.(nm{1}).opt.pm;

            results{li,ui,pi} = tmp.(nm{1});
        
            clear tmp
            savefile = 1;
            i = i + 1;
        end
    end
end

%%
c1 = brewermap(9,'Set3');
coptpm(1,:) = c1(1,:);

coptpm(2,:) = c1(6,:);
coptpm(3,:) = c1(5,:);

c4 = brewermap(9,'Set2');
coptpm(4,:) = c4(8,:);

c5 = brewermap(9,'PRGn');
coptpm(5,:) = c5(4,:);
coptpm(6,:) = c4(4,:);

locloop = unique(loc);
loadloop = unique(lc);

analysis.loadcases = {'UUV','OO','OO+UUV'};
auu = [1,3,2];
for xx = 1:length(locloop)
    ll = strcmp(locloop{xx}, locoptions);
    ll = find(ll);
    for uu  = 1:length(loadloop)
        emptyresults = [isempty(results{ll,uu,1}), isempty(results{ll,uu,2}), isempty(results{ll,uu,3}), isempty(results{ll,uu,4}), isempty(results{ll,uu,5})];
        for rr = 1:5
            if ~emptyresults(rr)
                tmpcost(rr) = results{ll,uu,rr}.output.min.cost;
                allgen(rr,:) = [results{ll,uu,rr}.output.min.kWwi{1},results{ll,uu,rr}.output.min.kWi{1},results{ll,uu,rr}.output.min.kWwa{1},results{ll,uu,rr}.output.min.kWd{1},results{ll,uu,rr}.output.min.kWc{1}];
                allsmax(rr) = results{ll,uu,rr}.output.min.Smax{1};
            else
                tmpcost(rr) = nan;
                allgen(rr,:) = [nan,nan,nan,nan,nan];
                allsmax(rr) = nan;
            end
            %visFireflies(results{ll,uu,rr},strcat(locoptions{ll}," LC: ",string(loadoptions(uu))," R: ",string(rr)))
        end
        [cost6D,ind] = min(tmpcost);
       %close all
        tmpgen = [results{ll,uu,ind}.output.min.kWwi{1},results{ll,uu,ind}.output.min.kWi{1},results{ll,uu,ind}.output.min.kWwa{1},results{ll,uu,ind}.output.min.kWd{1},results{ll,uu,ind}.output.min.kWc{1}];

        % Cost Variation in repeat solutions
        analysis.mincostvar(ll,auu(uu)) = std(tmpcost)/mean(tmpcost,'omitnan');

        % Hybridization Benefit
        % cost2D = min(minSys{ll,uu}.cost);
        % analysis.hybridbenefit(ll,auu(uu)) = abs(cost6D-cost2D)/(cost2D);
        % tmpgen(tmpgen >1E-3) = 0; %zeroing out any generators smaller than 1 W
        % if sum(tmpgen>0) > 1
        %     analysis.hybridbenefit(ll,auu(uu)) = 0;
        % end

        %plot results
        figure
        tf = tiledlayout(3,1);
        title(tf,strcat(locoptions{ll}," : ",string(loadoptions(uu))))
        nexttile
        bar(tmpcost,'FaceColor',[88/255, 90/255, 92/255])
        nexttile
        g = bar(allgen,'stacked');
        for i = 1:5
            g(i).FaceColor = coptpm(i,:);
        end
        nexttile
        bar(allsmax,'FaceColor',coptpm(6,:))


        clear tmpcost tmpgen allgen allsmax
    end
end

%% Visualization

%plot settings
colb = brewermap(9,'YlGnBu');
colm = brewermap(9,'BuPu');
coln = brewermap(9,'RdPu');
colw = brewermap(9,'YlGn');
colo = brewermap(9,'Oranges');

col = [colb(7,:); colm(7,:); coln(6,:); colw(7,:); colo(7,:)];

figure
set(gcf,'Units','Inches','Position',[0.5,0.5,10,5])
tf = tiledlayout(2,5);
tf.Padding = 'compact';
tf.TileSpacing = 'compact';
tf.InnerPosition = [0.11,0.18,0.85,0.65];
fs = 12;
c = 1;
X = categorical(analysis.loadcases);
X = reordercats(X,analysis.loadcases);
for l = 1:5
    ax(c) = nexttile;
    h = bar(X,analysis.mincostvar(l,:));
    h.FaceColor = col(l,:);
    %xlabel('Load Case','Interpreter','Latex')
    if c == 1
        ylh = ylabel('$\frac{\sigma}{\mu}$','Interpreter','Latex');
        set(ylh,'Rotation',0,'Units','Normalized','Position',[-.6 .5 -1], ...
            'VerticalAlignment','middle', ...
            'HorizontalAlignment','center','FontSize',16)
    end
    title(locdisp{l},'Interpreter','Latex','FontSize',fs,'Units','Normalized','Position',[0.5 1.2 -1])
    ax(c).TickLabelInterpreter = 'latex';
    ax(c).FontSize = fs;
    grid on
    box on

    ax(c+5) = nexttile(c+5);
    g = bar(X,analysis.hybridbenefit(l,:));
    g.FaceColor = col(l,:);
    xlabel('Load Case','Interpreter','Latex','FontSize',fs)
    if c == 1
        ylh = ylabel("$\frac{C_{H} - C_{S}}{C_S}$",'Interpreter','Latex');
        set(ylh,'Rotation',0,'Units','Normalized','Position',[-.6 .5 -1], ...
            'VerticalAlignment','middle', ...
            'HorizontalAlignment','center','FontSize',16)
    end
    ax(c+5).TickLabelInterpreter = 'latex';
    ax(c+5).FontSize = fs;
    grid on
    box on
    c = c + 1;
end
linkaxes(ax(1:5),'y')
linkaxes(ax(6:10),'y')


% if savefile
%     save(fullfile(selectedfolder,'SimplifiedHybridResults.mat'),"results")
% end