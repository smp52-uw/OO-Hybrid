function [] = visHybridSim(optStruct)

data = optStruct.data;
opt = optStruct.opt;
wave = optStruct.wave;
atmo = optStruct.atmo;
inso = optStruct.inso;
uc = optStruct.uc;
[data, opt] = prepHybrid(data,opt,uc,wave,atmo,inso);
output = optStruct.output;
length(data.met.time)
%extend time values
if ~isequal(length(output.min.Pdies),length(data.met.time))
    orig_l = length(data.met.time);
    vecMid = datevec(data.met.time(end));
    data.met.time = [data.met.time ; zeros(length(output.min.Pdies) ...
        - length(data.met.time),1)];
    for t = orig_l+1:length(data.met.time)
        vec = vecMid;
        vec(4) = vecMid(4) + t - orig_l;
        data.met.time(t) = datenum(vec);
    end
end

%I think this is redundant
% [data.met.wind_spd,data.met.time] = ...
%     extendToLifetime(data.met.wind_spd,data.met.time, ...
%     optStruct.uc.lifetime);

figure
set(gcf, 'Position', [20, 20, 1300, 700])
%STORAGE TIME SERIES
ax(1) = subplot(8,1,1);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.S1(1:end-1)/1000,'Color',[255,69,0]/256, ... 
    'DisplayName','Battery Storage 1','LineWidth',2)
hold on
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.S2(1:end-1)/1000,'Color',[245,238,54]/256, ... 
    'DisplayName','Battery Storage 2','LineWidth',2)
legend('show')
ylabel('[kWh]')
ylim([0 inf])
yticks([0 max(output.min.S1(1:end-1)/1000)])
xticks([])
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on

%POWER TIME SERIES
ax(2) = subplot(8,1,2);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.Pdies/1000,'Color',[161,65,225]/256, ... 
    'DisplayName','Diesel Power','LineWidth',2)
legend('show')
ylabel('[kW]')
xticks([])
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on

%POWER TIME SERIES
ax(3) = subplot(8,1,3);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.Pinso/1000,'Color',[225,177,65]/256, ... 
    'DisplayName','Solar Power','LineWidth',2)
legend('show')
ylabel('[kW]')
xticks([])
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on

%POWER TIME SERIES
ax(4) = subplot(8,1,4);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.Pwave/1000,'Color',[65,105,225]/256, ... 
    'DisplayName','Wave Power','LineWidth',2)
legend('show')
ylabel('[kW]')
xticks([])
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on

%POWER TIME SERIES
ax(5) = subplot(8,1,5);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.Pwind/1000,'Color',[6,139,33]/256, ... 
    'DisplayName','Wind Power','LineWidth',2)
legend('show')
ylabel('[kW]')
xticks([])
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on

%DUMPED POWER TIME SERIES
ax(6) = subplot(8,1,6);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.D/1000,'Color',[199,4,128]/256, ... 
    'DisplayName','Power Dumped','LineWidth',2)
legend('show')
ylabel('[kW]')
xticks([])
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on

%Load Time Series
ax(7) = subplot(8,1,7);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.L/1000,'Color',[212,92,176]/256, ... 
    'DisplayName','Load','LineWidth',2)
legend('show')
ylabel('[kW]')
%ylim([0 inf])
yticks([0 1])
xticks([])
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on


%Failure Time Series
ax(8) = subplot(8,1,8);
plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
    output.min.F,'.','Color',[0,0,0], ... 
    'DisplayName','Failures','LineWidth',2)
legend('show')
ylabel('[logical]')
ylim([0,2])
yticks([0 1])
xlabel('Time')
set(gca,'FontSize',10)
set(gca,'LineWidth',2)
grid on
% %EFFICIENCY TIME SERIES
% ax(4) = subplot(4,1,4);
% plot(datetime(data.met.time,'ConvertFrom','datenum'), ...
%     output.min.eff_t,'Color',[255,20,147]/256, ... 
%     'DisplayName','Efficiency','LineWidth',2)
% legend('show')
% ylabel('[~]')
% xlabel('Time')
% ylim([0 1.25*max(output.min.eff_t)])
% set(gca,'FontSize',16)
% set(gca,'LineWidth',2)
% grid on
% 
% set(gcf, 'Position', [100, 100, 1400, 650])

linkaxes(ax,'x')
linkaxes(ax(2:6),'y')


end

