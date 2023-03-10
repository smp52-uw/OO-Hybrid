function [] = visSurvSpace(optStruct)

opt = optStruct.opt;
output = optStruct.output;

%plot settings
Sm_max = 500; %[kWh]
Gr_max = 8; %[kW]

%create grid
[Smaxgrid,kWgrid] = meshgrid(opt.Smax,opt.kW);
lw = 1.1;
fs = 12;

a_sat = output.surv; %availability satisfied
a_sat(output.surv == 0) = nan;

figure
%s = surf(Smaxgrid,kWgrid,output.surv);
s = surf(Smaxgrid,kWgrid,a_sat);
s.EdgeColor = 'none';
s.FaceColor = 'flat';
hold on
view(0,90)
xlabel('Storage Capacity [kWh]')
ylabel('Rated Power [kW]')
ylim([-inf Gr_max])
xlim([-inf Sm_max])
c = colorbar;
c.Label.String = 'Persistence';
caxis([min(a_sat(:)) max(output.surv(:))]) %to produce cartoon
%lb = (output.min.surv)/(max(output.surv(:))/2);
% AdvancedColormap('bg l w r',8000,[1*lb,lb+.05*(1-lb),lb+0.1*(1-lb),1])
% AdvancedColormap('kkgw glw lww r',1000, ...
%        [lb,lb+.05*(1-lb),lb+0.15*(1-lb),1]);
hold on
set(gca,'FontSize',fs)
set(gca,'LineWidth',lw)
grid on

end
