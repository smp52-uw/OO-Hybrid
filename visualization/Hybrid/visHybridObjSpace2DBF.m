function [] = visHybridObjSpace2DBF(optStruct)

opt = optStruct.opt;

%adjust cost to thousands
cost = optStruct.output.cost{1}/1000;
surv = optStruct.output.surv{1};

costB = cost;
costW = cost;

costB(surv<0.99) = 2*max(cost);
costW(surv<0.99) = nan;
%create grid
disc = 500;
Smax = linspace(opt.Smax_1,opt.Smax_n,500).';

if opt.pm == 1
    kW = linspace(opt.wind.kW_1,opt.wind.kW_m,disc); 
    kWmin = optStruct.output.min.kWwi{1};
    kWrun = optStruct.output.Kwi_run{1};
elseif opt.pm == 2
    kW = linspace(opt.inso.kW_1,opt.inso.kW_m,disc);
    kWmin = optStruct.output.min.kWi{1};
    kWrun = optStruct.output.Ki_run{1};
elseif opt.pm == 3
    kW = linspace(opt.wave.kW_1,opt.wave.kW_m,disc);
    kWmin = optStruct.output.min.kWwa{1};
    kWrun = optStruct.output.Kwa_run{1};
elseif opt.pm == 4
    kW = linspace(opt.dies.kW_1,opt.dies.kW_m,disc);   
    kWmin = optStruct.output.min.kWd{1};
    kWrun = optStruct.output.Kwd_run{1};
else
    kW = linspace(opt.curr.kW_1,opt.curr.kW_m,disc).'; 
    kWmin = optStruct.output.min.kWc{1};
    kWrun = optStruct.output.Kc_run{1};
end

[kWgrid,Smaxgrid] = ndgrid(kW,Smax);
costgridW = reshape(costW,[disc,500]);
costgridB = reshape(costB,[disc,500]);
survgrid = reshape(surv,[disc,500]);

%check min point
indM = find(kWgrid == kWmin & Smaxgrid == optStruct.output.min.Smax{1});

load('appleCMap.mat');
cMap = flipud(cMap); %switching it to gr -> rd
cMap = cMap(100:end,:); %cutting off the mostly white colors

figure %cost figure (white)
hold on
Srun = optStruct.output.S_run{1};
scatter3(Srun,kWrun,costW,[],costW,'filled');
colormap(gca, cMap)
pl(1) = plot3(optStruct.output.min.Smax{1},kWmin,5*max(cost),'co','LineWidth',2,'DisplayName','Optimal Point');

view(0,90)
xlabel('Storage Capacity [kWh]')
ylabel('Rated Power [kW]')
title(strcat("Total Cost: ", optStruct.loc,", LC = ",string(optStruct.uc.loadcase), ", PM = ",string(opt.pm)))
c = colorbar;
c.Label.String = '[$] in thousands';

grid on

figure %cost figure (black)
hold on
costtr = costgridB.';
s = surf(Smaxgrid.',kWgrid.',costtr);
colormap(cMap)
s.EdgeColor = 'none';
s.FaceColor = 'flat';
pl(1) = plot3(optStruct.output.min.Smax{1},kWmin,5*max(cost),'co','LineWidth',2,'DisplayName','Optimal Point');

view(0,90)
xlabel('Storage Capacity [kWh]')
ylabel('Rated Power [kW]')
title(strcat("Total Cost: ", optStruct.loc,", LC = ",string(optStruct.uc.loadcase), ", PM = ",string(opt.pm)))
c = colorbar;
c.Label.String = '[$] in thousands';

grid on


% figure %surv figure
% hold on
% survtr = survgrid.';
% s = surf(Smaxgrid.',kWgrid.',survtr);
% colormap(cMap)
% s.EdgeColor = 'none';
% s.FaceColor = 'flat';
% pl(1) = plot3(optStruct.output.min.Smax{1},kWmin,max(cost),'co','LineWidth',2,'DisplayName','Optimal Point');
% 
% view(0,90)
% xlabel('Storage Capacity [kWh]')
% ylabel('Rated Power [kW]')
% title(strcat("Total Cost: ", optStruct.loc,", LC = ",string(optStruct.uc.loadcase), ", PM = ",string(opt.pm)))
% c = colorbar;
% c.Label.String = '[$] in thousands';
% 
% grid on



% fig = figure;
% set(gcf,'Units','inches')
% set(gcf,'Position', [1, 1, 8.5, 4])
% 
% if opt.pm == 1
%     costcat = {'kWcost_wind','Icost_wind','Scost','Pmtrl','Pinst','Pmooring','vesselcost','turbrepair','battencl'};
%     costname = {'$kWcost_{wind}$','$Icost_{wind}$','Scost','Pmtrl','Pinst','Pmooring','vesselcost','turbrepair','battencl'};
% elseif opt.pm == 2
%     costcat = {'Mcost_inso','Ecost_inso','Icost_inso','Strcost_inso','Scost','Pmtrl','Pinst','Pmooring','vesselcost','solarrepair','battencl'};
%     costname = {'$Mcost_{inso}$','$Ecost_{inso}$','$Icost_{inso}$','$Strcost_{inso}$','Scost','Pmtrl','Pinst','Pmooring','vesselcost','solarrepair','battencl'};
% end
% allfields = fieldnames(optStruct.output.min);
% for i = 1:length(costcat)
%      tempind = find(strcmp(costcat{i},allfields));
%     if isempty(tempind)
%         indcost(i) = nan;
%         Y(i) = nan;
%     else
%         indcost(i) = tempind;
%         Y(i) = optStruct.output.min.(allfields{indcost(i)});
%     end
% end
% 
% X = categorical(costname);
% 
% bh = bar(X,Y);
% ylabel('$')
% axesH = findall(fig, "Type", "axes");
% set(axesH, "TickLabelInterpreter", 'latex')
% title(strcat("Optimal Point (Total Cost = $",string(round(optStruct.output.min.cost)),")"))
% 

end

