%visualize the output from the hyperparameter sensitivity study on the 6D
%firefly algorithm

%set the folder with the sensitivity results
selectedfolder = uigetdir();

fileList = dir(fullfile(selectedfolder,'*.mat'));
numfiles = max(size(fileList));

% textList = dir(fullfile(selectedfolder,'*.txt'));
% 
% for i = 1:numfiles
%     tmp = split(fileList(i).name,'_');
%     am = tmp(5);
%     arraymat(i) = str2double(am{1});
% end
% 
% for j = 1:max(size(textList))
%     tmp = split(textList(j).name,'_');
%     ttmp = split(tmp(7),'.');
%     at = ttmp(1);
%     at = at{1};
%     arraytext(j) = str2double(at);
% end
% 
% [C,ia] = setdiff(arraytext,arraymat);

%either load the simplified mat file or load all files
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
            surv(i) = tmp.(nm{1}).output.min.surv;
            cost(i) = tmp.(nm{1}).output.min.cost;
            popbestcost{i} = tmp.(nm{1}).output.popBestCost;
            
            gen(i,:) = [tmp.(nm{1}).output.min.kWwi{1}, tmp.(nm{1}).output.min.kWi{1}, tmp.(nm{1}).output.min.kWwa{1}, tmp.(nm{1}).output.min.kWd{1}, tmp.(nm{1}).output.min.kWc{1}];
            smax(i) = tmp.(nm{1}).output.min.Smax{1};
            ffa(i) = tmp.(nm{1}).opt.ffa;
            failsurv(i) = tmp.(nm{1}).opt.failsurv;
        
            clear tmp
            savefile = 1;
            i = i + 1;
        end
    end
end

%% extract ffa parameters
numfiles = max(length(cost)); %updated numfiles in case some where unreadable

for i = 1:numfiles
    bb(i) = ffa(i).beta0;
    popn(i) = ffa(i).pop;
    gamn(i) = ffa(i).gamma;
    betn(i) = ffa(i).beta0;
    alphan(i) = ffa(i).alpha;
    adn(i) = ffa(i).adamp;
end
indgood = find(bb ~= 1);
%indgood = 1:1:numfiles;
[pk,ind] = findpeaks(cost(indgood),'MinPeakProminence',5000);
for i = 1:length(ind)
    it = string(ffa(indgood(ind(i))).max);
    pop = string(ffa(indgood(ind(i))).pop);
    gamma = string(ffa(indgood(ind(i))).gamma);
    beta = string(ffa(indgood(ind(i))).beta0);
    alpha = string(ffa(indgood(ind(i))).alpha);
    adamp = string(ffa(indgood(ind(i))).adamp);

    dispnm(i) = strcat("I: ",it, " P: ",pop," G: ",gamma," B: ",beta," A: ",alpha," Ad: ",adamp);
end

%% plots
% var = popn;
% opts = unique(var);
% for p = 1:length(opts)
%     indv{p} = find(var == opts(p));
% end
% cMapb = cmasherImport('arctic',300);
% cMapb = cMapb(50:end,:);
% cMapo = cmasherImport('sunburst',300);
% cMapo = cMapo(50:end,:);
% cMapg = cmasherImport('jungle',300);
% cMapg = cMapg(50:end,:);
% 
% figure
% hold on
% c = 1; d = 1; e = 1;
% for i = 1:numfiles
%     if sum(i == indv{1})
%         plot(popbestcost{i},'color',cMapb(c,:))
%         c = c + 1;
%     elseif sum(i == indv{2})
%         plot(popbestcost{i},'color',cMapo(d,:))
%         d = d + 1;
%     else
%         plot(popbestcost{i},'color',cMapg(d,:))
%         e = e + 1;
%     end
% end
% ylabel('Population Best Cost')
% xlabel('Iteration')

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

fs = 14;

figure
set(gcf,'Units','Inches','Position',[0.5,0.5,10,5])
tf = tiledlayout(2,2);
tf.Padding = 'loose';
tf.TileSpacing = 'loose';

ax(1) = nexttile;
plot(cost(indgood)/min(cost),'k','linewidth',1.2)
ylh = ylabel("$\frac{\mathrm{cost}}{\mathrm{min(cost)}}$",'Interpreter','Latex');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.23 .5 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center','FontSize',fs)
ax(1).TickLabelInterpreter = 'latex';
ax(1).FontSize = fs;
ylim([1,1.2])
grid on
box on

ax(2) = nexttile(3);
plot(surv(indgood),'k','linewidth',1.2)
xlabel('Sensitivity Iterations','Interpreter','Latex','FontSize',fs)
ax(2).TickLabelInterpreter = 'latex';
ax(2).FontSize = fs;
ylim([0.9,1.1])
ylh = ylabel('$a_{sim}$','Interpreter','Latex');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.23 .5 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center','FontSize',fs)
grid on
box on

ax(3) = nexttile(2);
hold on
for i = 1:5
    plot(gen(indgood,i),'linewidth',1.2,'Color',coptpm(i,:))
end
legend('Wind','Solar','Wave','Dies','Curr','NumColumns',2,'Interpreter','Latex','FontSize',fs)
ylh = ylabel({'Generation', '[kW]'},'Interpreter','Latex');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.23 .5 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center','FontSize',fs)
ax(3).TickLabelInterpreter = 'latex';
ax(3).FontSize = fs;
grid on
box on

ax(4) = nexttile(4);
plot(smax(indgood),'k','linewidth',1.2,'Color',coptpm(6,:))
hold on
ylh = ylabel({'Battery', '[kWh]'},'Interpreter','Latex');
set(ylh,'Rotation',0,'Units','Normalized','Position',[-.23 .5 -1], ...
    'VerticalAlignment','middle', ...
    'HorizontalAlignment','center','FontSize',fs)
xlabel('Sensitivity Iterations','Interpreter','Latex','FontSize',fs)
ax(4).TickLabelInterpreter = 'latex';
ax(4).FontSize = fs;
grid on
box on
linkaxes(ax,'x')

%% 
figure
tf = tiledlayout(5,1);
tf.Padding = 'tight';
tf.TileSpacing = 'tight';

ax(1) = nexttile;
plot(cost(indgood),'linewidth',1.2)
hold on
% for i = 1:length(x10)
%     xline(x10(i),'r','linewidth',1.2)
% end
ylabel('min Cost')
grid on

ax(2) = nexttile;
plot(cost(indgood)/min(cost),'linewidth',1.2)
ylabel({'Fraction of', 'Min Cost'})
grid on

ax(3) = nexttile;
plot(surv(indgood),'linewidth',1.2)
ylabel('Surv')
grid on

ax(4) = nexttile;
hold on
for i = 1:5
    plot(gen(indgood,i),'linewidth',1.2)
end
legend('Wind','Solar','Wave','Dies','Curr')
% for i = 1:length(x10)
%     xline(x10(i),'r','linewidth',1.2)
% end
ylabel({'Generation', '[kW]'})
grid on

ax(5) = nexttile;
plot(smax(indgood),'linewidth',1.2)
hold on
% for i = 1:length(x10)
%     xline(x10(i),'r','linewidth',1.2)
% end
ylabel({'Battery', '[kWh]'})
xlabel('Run Count')
grid on

linkaxes(ax,'x')

if savefile
    save(fullfile(selectedfolder,'SimplifiedResults.mat'),"cost","surv","gen","smax","ffa")
end