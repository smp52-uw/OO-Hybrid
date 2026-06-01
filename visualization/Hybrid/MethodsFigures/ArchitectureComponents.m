%Make pieces of archeticture figure

%load 2D space for firefly tiles
tmp = load('C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\HybridPaper\AllLUP\LC3_NoBatteryRule\AllLUP_PD2_FFASENSE_7_05202026.mat');
optStruct = tmp.AllLUP_PD2_FFASENSE_7_05202026;
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
    kWrun = optStruct.output.Kd_run{1};
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

cMap = cmasherImport('ember',1000);
cMap = flipud(cMap); %switching it to gr -> rd

figure %cost figure (black)
set(gcf,'Unit','Inches','Position',[0.7,0.7,3,2.7])
hold on
costtr = costgridB.';
s = surf(Smaxgrid.',kWgrid.',costtr);
colormap(cMap)
s.EdgeColor = 'none';
s.FaceColor = 'flat';
s.FaceAlpha = 1;

view(0,90)
xlabel('Storage Capacity [kWh]','Interpreter','latex','FontSize',11)
ylabel('Rated Power [kW]','Interpreter','latex','FontSize',11)
set(gca,"TickLabelInterpreter",'latex','FontSize',11)
xlim([0,500])
ylim([0,8])
% c = colorbar;
% c.Label.String = '[$] in thousands';
clim([min(costtr,[],'all') 1.05*max(costgridW,[],'all')])
grid on


%% load point 1kW 10 kWh 
%load("PointOptimization.mat");

colors = brewermap(11, 'Set3');
colors2 = brewermap(11, 'Set2');
colors3 = brewermap(11, 'PRGn');
colors4 = brewermap(11, 'RdPu');
colors5 = brewermap(11, 'OrRd');
colpm(1,:) = colors(1,:);
colpm(2,:) = colors(6,:);
colpm(3,:) = colors(5,:);
colpm(4,:) = colors2(8,:);
colpm(5,:) = colors3(4,:);
colb(1,:) = colors2(4,:);
colb(2,:) = colors4(4,:);
colb(3,:) = colors5(8,:);
ts = 10228;
te = 10329;
%% plot time series
figure
set(gcf,'Unit','Inches','Position',[0.7,0.7,2,1])
plot(output.min.Pinso,'color',colpm(2,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on

figure
set(gcf,'Unit','Inches','Position',[0.7,0.7,2,1])
plot(output.min.Pwind,'color',colpm(1,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on

figure
set(gcf,'Unit','Inches','Position',[0.7,0.7,2,1])
plot(output.min.Pwave,'color',colpm(3,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on

figure
set(gcf,'Unit','Inches','Position',[0.7,0.7,2,1])
plot(output.min.Pcurr,'color',colpm(5,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on

figure
set(gcf,'Unit','Inches','Position',[0.7,0.7,2,1])
plot(output.min.Pdies,'color',colpm(4,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on

%%
figure
set(gcf,'Unit','Inches','Position',[0.7,0.7,2,1])
plot(output.min.S1,'color',colb(1,:),'linewidth',1.2)
xlim([ts,te])
%ylim([0,11*1000])
xticklabels([])
yticklabels([])
grid on 
box on

figure
set(gcf,'Unit','Inches','Position',[2.7,0.7,2,1])
plot(output.min.L,'color',colb(2,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on

figure
set(gcf,'Unit','Inches','Position',[4.7,0.7,2,1])
plot(output.min.D,'color',colb(3,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on

%% Self discharge
sd1 = output.min.S1*(batt.sdr/100)*(1/(30*24)); %[Wh] self discharge
figure
set(gcf,'Unit','Inches','Position',[4.7,0.7,2,1])
plot(sd1,'color',colb(3,:),'linewidth',1.2)
xlim([ts,te])
xticklabels([])
yticklabels([])
grid on 
box on