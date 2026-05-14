% Find optimal FFA settings

%load all data
folder = 'C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\FFA6DSensitivity_SW0';
allContents = dir(folder);
numfold = max(size(allContents));
c = 1;
for n = 1:numfold
    if ~strcmp(allContents(n).name,'.') && ~strcmp(allContents(n).name,'..')
        filefold = strcat(folder,'\',allContents(n).name);
        allsens(c) = load(fullfile(filefold,'SimplifiedResults.mat'));

        tmpnms = split(allContents(n).name,'_');
        name{c} = tmpnms{2};
        c = c + 1;
    end
end
%% choose an approx location of optimal settings
ind = 65; %This looks like a good spot for all cases

% find the FFA settings for one file
ffaopt = [allsens(1).ffa(ind).max,allsens(1).ffa(ind).pop,allsens(1).ffa(ind).gamma,allsens(1).ffa(ind).beta0,allsens(1).ffa(ind).alpha,allsens(1).ffa(ind).adamp];

%% find index of all files corresponding to that FFA setting
indopt = nan(1,length(allsens));
for i = 1:length(allsens)
    tmpffa = allsens(i).ffa;
    clear ffamat allcosts
    for j = 1:length(tmpffa)
        ffamat(j,:) = [tmpffa(j).sens{1}, tmpffa(j).sens{2}, tmpffa(j).sens{3}, tmpffa(j).sens{4}, tmpffa(j).sens{5}, tmpffa(j).sens{6}];
        allcosts(j) = allsens(i).cost(j);
    end
    [~,indopt(i)] = ismember(ffaopt, ffamat, 'row');
    costopt(i) = allsens(i).cost(indopt(i));
    errorcost(i) = (costopt(i) - min(allcosts))/min(allcosts);
end

%plot the difference between this point's cost and min cost for all cases

X = categorical({name{1},name{2},name{3},name{4},name{5},name{6}});
X = reordercats(X,{name{1},name{2},name{3},name{4},name{5},name{6}});

figure
bar(X,errorcost.*100)
ylabel('Cost Error [/%]')
title({'Cost Error vs min cost with FFA:', ...
    strcat("I: ",string(ffaopt(1))," P: ",string(ffaopt(2))," G: ",string(ffaopt(3))," B: ",string(ffaopt(4)), ...
    " A: ",string(ffaopt(5))," Ad: ",string(ffaopt(6)))})