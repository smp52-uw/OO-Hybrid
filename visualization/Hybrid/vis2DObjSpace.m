function [] = vis2DObjSpace(optStruct,pm1,Gr1_max)

%diesel=1, inso =2, wind = 3, wave =4
set(0,'DefaultTextFontname', 'cmr10')
set(0,'DefaultAxesFontName', 'cmr10')
set(0,'defaulttextinterpreter','latex')
j=500;
n=500;

opt = optStruct.opt;
output = optStruct.output;
Smax_full = opt.Smax;
if pm1 == 4
    kW1_full = opt.dies.kW;
elseif pm1 == 2
    kW1_full = opt.inso.kW;
elseif pm1 == 1
    kW1_full = opt.wind.kW;
elseif pm1 == 3
    kW1_full = opt.wave.kW;
elseif pm1 == 5
    kW1_full = opt.curr.kW;
end

%plot settings
Sm_max = 500; %[kWh]
%Gr1_max = 16; %[kW]

%adjust cost to thousands
size_it = size(output.cost);
%figure(1)
colordata = permute(repmat([255 255 245]'./256,[1,500,500]),[3 2 1]);
for i=1: size_it(2)
    cost = output.cost{i};
    size(cost)
    cost = reshape(cost,[j,n]);
    
    %create grid
    %[Kd,Ki,Kwi,Kwa,S]
    kW1 = kW1_full{i};
    Smax = Smax_full{i};
    [kWgrid,Smaxgrid]= ndgrid(kW1,Smax);
    lw = 1.1;
    fs = 10;
    kWlist = reshape(kWgrid,[j*n 1]);
    Slist = reshape(Smaxgrid,[j*n 1]);
    %remove failure configurations
    a_sat = cost; %availability satisfied
    
    %remove all output data other than the dimensions we want
    surv_p = output.surv{i};
    
    %reshape into 3x3x3 grid
    surv_p = reshape(surv_p,[j,n]);
    kW_2d = reshape(kWlist,[j,n]);
    S_2d = reshape(Slist,[j,n]);
    costdata = reshape(a_sat,[j,n]);
    
    costdata(surv_p < 0.99) = nan;
    %find minimum
    [m,m_ind] = min(costdata(:));
    
%     s = surf(Smaxgrid,kWgrid,cost,1.*ones(length(Smaxgrid), ... 
%         length(kWgrid),3)); %white
%     s = surf(Smaxgrid,kWgrid,cost,colordata); %off white
    s.EdgeColor = 'none';
    s.FaceColor = 'flat';
    hold on
    s = surf(Smaxgrid,kWgrid,costdata);
    s.EdgeColor = 'none';
    s.FaceColor = 'flat';
    hold on
    grid on
%     c = colorbar;
%     c.Label.String = '[kg]';
end
output.min.cost = output.min.cost;
scatter3(Smaxgrid(m_ind),kWgrid(m_ind),m*2,50,'ko', ...
    'LineWidth',1.5,'MarkerEdgeColor','k')
view(0,90)
xlabel('$S_{max}$ [kWh]')


ylim([-inf Gr1_max])
xlim([-inf Sm_max])

c = colorbar;
c.Label.String = '[kg]';
%caxis([0 max(costdata(:))/1.6]) %to produce cartoon
%lb = (output.min.cost)/(max(costdata(:))/1.6);
lb=0
%lb = (output.min.cost)/max(costdata(:))
% if pm1 == 5
%     lb = lb/2;
% end
% AdvancedColormap('bg l w r',8000,[1*lb,lb+.05*(1-lb),lb+0.1*(1-lb),1])
% AdvancedColormap('kkgw glw lww r',1000, ...
%         [lb,lb+.05*(1-lb),lb+0.15*(1-lb),1]);
p_max = max(max(s.ZData(s.YData(:,1) < Gr1_max, ... 
    s.XData(1,:) < Sm_max))); %plot max
a_min = min(costdata(:));
p_max_mult = p_max/a_min;

AdvancedColormap('kkgw glww lww r rk k',1000, ...
    [lb,lb+.025*(1-lb),lb+0.1*(1-lb),.35,.7,1]);
caxis([a_min(i) ceil(max(p_max_mult)*a_min(i))])
% AdvancedColormap('kkgw glw lww r',1000, ...
%         [lb,lb+.05*(1-lb),lb+0.15*(1-lb),1]);
hold on
set(gca,'FontSize',fs)
set(gca,'LineWidth',lw)

% c = colorbar;
% c.Label.String = '[kg]';
% caxis([0 max(a_sat(:))/2]) %to produce cartoon
% lb = (output.min.cost)/(max(a_sat(:))/2);
% AdvancedColormap('bg l w r',8000,[1*lb,lb+.05*(1-lb),lb+0.1*(1-lb),1])
% AdvancedColormap('kkgw glw lww r',1000, ...
%         [lb,lb+.05*(1-lb),lb+0.15*(1-lb),1]);

end


