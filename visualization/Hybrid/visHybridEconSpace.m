function [] = visHybridEconSpace(optStruct,pm)

opt = optStruct.opt;

%adjust cost to thousands
cost = optStruct.cost/1000;
mass = optStruct.mass;
%create grid
disc = 400;
Smax = linspace(opt.Smax_1,opt.Smax_n,500);

if opt.pm == 1
    kW = linspace(opt.wind.kW_1,opt.wind.kW_m,disc); 
elseif opt.pm == 2
    kW = linspace(opt.inso.kW_1,opt.inso.kW_m,disc);
elseif opt.pm == 3
    kW = linspace(opt.wave.kW_1,opt.wave.kW_m,disc); 
elseif opt.pm == 4
    kW = linspace(opt.dies.kW_1,opt.dies.kW_m,disc);   
else
    kW = linspace(opt.curr.kW_1,opt.curr.kW_m,disc); 
end

[kWgrid,Smaxgrid] = ndgrid(kW,Smax);
costgrid = reshape(cost,[disc,500]);
massgrid = reshape(mass.total,[disc,500]);
buoygrid = reshape(optStruct.indPM,[disc,500]);

testS = linspace(1,500,500);
med = 11000.*ones(size(testS));
sm = 1200.*ones(size(testS));
lg = 16000.*ones(size(testS));

massinterp = scatteredInterpolant(mass.total,reshape(Smaxgrid,[disc*500,1]),reshape(kWgrid,[disc*500,1]));


kWsm = massinterp(sm,testS);
kWmed = massinterp(med,testS);
kWlg = massinterp(lg,testS);

kWsm(kWsm<0) = nan;
kWsm(kWsm>8) = nan;
kWmed(kWmed<0) = nan;
kWmed(kWmed>8) = nan;
kWlg(kWlg<0) = nan;
kWlg(kWlg>8) = nan;

figure
s = surf(kWgrid,Smaxgrid,costgrid);
s.EdgeColor = 'none';
s.FaceColor = 'flat';

view(0,90)
ylabel('Storage Capacity [kWh]')
xlabel('Rated Power [kW]')
title('Total Cost')
c = colorbar;
c.Label.String = '[$] in thousands';

hold on
grid on

figure
s = surf(kWgrid,Smaxgrid,buoygrid);
s.EdgeColor = 'none';
s.FaceColor = 'flat';

view(0,90)
ylabel('Storage Capacity [kWh]')
xlabel('Rated Power [kW]')
title('buoy')
c = colorbar;

hold on
grid on


figure
s = surf(kWgrid,Smaxgrid,massgrid);
s.EdgeColor = 'none';
s.FaceColor = 'flat';
view(0,90)
c = colorbar;
c.Label.String = '[kg]';
hold on
col = brewermap(10,'greys'); %colors
zloc = 1.2*max(massgrid,[],'all').*ones(size(testS));
pl(1) = plot3(kWsm,testS,zloc,'-','Color',col(4,:),'LineWidth',2,'DisplayName','Small Buoy Limit');
pl(2) = plot3(kWmed,testS,zloc,'-','Color',col(8,:),'LineWidth',2,'DisplayName','Medium Buoy Limit');
pl(3) = plot3(kWlg,testS,zloc,'r-','LineWidth',2,'DisplayName','Large Buoy Limit');

ylabel('Storage Capacity [kWh]')
xlabel('Rated Power [kW]')
title('Total Payload Mass')
legend(pl,'Location','northoutside','Orientation','horizontal')


grid on



[~,id1] = min(abs(0.3-kW));
kWr(1) = kW(id1);
Smr(1) = 15;

[~,id2] = min(abs(3.5-kW));
kWr(2) = kW(id2);
Smr(2) = 200;

[~,id3] = min(abs(5-kW));
kWr(3) = kW(id3);
Smr(3) = 350;

costcomp = optStruct.costcomp;
fn = fieldnames(costcomp);
for f = 1:length(fn)
    if length(costcomp.(fn{f})) > 1
        costcompgrid.(fn{f}) = reshape(costcomp.(fn{f}),[disc,500]);
        costcompgrid.(fn{f}) = costcompgrid.(fn{f})./1000;
    end
end

figure
set(gcf,'Units','inches')
set(gcf,'Position', [1, 1, 8.5, 4])
tiledlayout(1,3)
for j = 1:3

    ind = find(kWgrid == kWr(j) & Smaxgrid == Smr(j));
    c = 1; %counter 
    for f = 1:length(fn)
        if isfield(costcompgrid,fn{f})
            costpoint(j,f) = costcompgrid.(fn{f})(ind);
            if costpoint(j,f) > 0
                plotcost{j}(c) = costpoint(j,f);
                plotX{j}{c} = fn{f};
                c = c + 1;
            end
        else
            costpoint(j,f) = 0;
        end
    end
    ax(j) = nexttile;
    X = categorical(plotX{j});
    bh = bar(X,plotcost{j});
    maxC(j) = max(plotcost{j});
    ylabel('1000 $')
    title(strcat("kW: ",string(kWr(j))," S: ",string(Smr(j))))
end
ylim([0,max(maxC)])
linkaxes(ax,'y')

% %visualize details for the 3 regions
% 
% [~,id1] = min(abs(1.5-kW));
% kWr(1) = kW(id1);
% Smr(1) = 150;
% 
% [~,id2] = min(abs(3.5-kW));
% kWr(2) = kW(id2);
% Smr(2) = 200;
% 
% [~,id3] = min(abs(5-kW));
% kWr(3) = kW(id3);
% Smr(3) = 200;
% 
% costcomp = optStruct.costcomp;
% fn = fieldnames(costcomp);
% for f = 1:length(fn)
%     if length(costcomp.(fn{f})) > 1
%         costcompgrid.(fn{f}) = reshape(costcomp.(fn{f}),[disc,disc]);
%         costcompgrid.(fn{f}) = costcompgrid.(fn{f})./1000;
%     end
% end
% indMoor = reshape(optStruct.indMoor,[disc,disc]);
% 
% figure
% colors = colormap('winter');
% for j = 1:3
% 
%     ind = find(kWgrid == kWr(j) & Smaxgrid == Smr(j));
%     indMoorPoint(j) = indMoor(ind);
%     for f = 1:length(fn)
%         if isfield(costcompgrid,fn{f})
%             costpoint(j,f) = costcompgrid.(fn{f})(ind);
%         else
%             costpoint(j,f) = 0;
%         end
%     end
% end
% colors = colormap(brewermap(length(fn),'Spectral')); %colors
% X = categorical({'low gen','med gen','high gen'});
% X = reordercats(X,{'low gen','med gen','high gen'});
% bh = bar(X,costpoint,'stacked');
% ylabel('Cost [1000$]')
% for i = 1:length(fn)
%     bh(i).FaceColor = 'flat';
%     bh(i).CData = colors(i,:);
% end
% disp('done')
end

