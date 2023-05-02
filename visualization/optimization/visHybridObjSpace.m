function [] = visHybridObjSpace(optStruct,d,pm1,pm2)

%Set up a 2D or 3D ObjSpace
%If 2D then pm2 = ~
%diesel=1, inso =2, wind = 3, wave =4
j=9;
k=9;
l=9;
m=9;
n=9;

opt = optStruct.opt;
output = optStruct.output;
Smax = opt.Smax{end};
if pm1 == 1
    kW1 = opt.dies.kW{end};
elseif pm1 == 2
    kW1 = opt.inso.kW{end};
elseif pm1 == 3
    kW1 = opt.wind.kW{end};
elseif pm1 == 4
    kW1 = opt.wave.kW{end};
end
if pm2 == 1
    kW2 = opt.dies.kW{end};
elseif pm2 == 2
    kW2 = opt.inso.kW{end};
elseif pm2 == 3
    kW2 = opt.wind.kW{end};
elseif pm2 == 4
    kW2 = opt.wave.kW{end};
end
%plot settings
Sm_max = 500; %[kWh]
Gr1_max = 16; %[kW]
Gr2_max = 16; %[kW]

%adjust cost to thousands
output.cost = output.cost{end};
output.min.cost = output.min.cost;

%create grid

kW3 = opt.dies.kW{end};
kW4 = opt.inso.kW{end};
%[Kd,Ki,Kwi,Kwa,S]
[Kdgrid,Kigrid,Kwigrid,Kwagrid,Sgrid]= ndgrid(kW3, kW4,kW1,kW2,Smax);
lw = 1.1;
fs = 18;
Kdlist = reshape(Kdgrid,[j*k*l*m*n 1]);
Kilist = reshape(Kigrid,[j*k*l*m*n 1]);
Kwilist = reshape(Kwigrid,[j*k*l*m*n 1]);
Kwalist = reshape(Kwagrid,[j*k*l*m*n 1]);
Slist = reshape(Sgrid,[j*k*l*m*n 1]);
%remove failure configurations
a_sat = output.cost; %availability satisfied

%remove all output data other than the dimensions we want
surv_p = output.surv{end};
surv_p(Kdlist~=0 | Kilist~=0) = [];
Kwilist(Kdlist~=0 | Kilist~=0) = [];
Kwalist(Kdlist~=0 | Kilist~=0) = [];
Slist(Kdlist~=0 | Kilist~=0) = [];
a_sat(Kdlist~=0 | Kilist~=0) = [];
%find minimum
[m,m_ind] = min(a_sat(:));

%reshape into 3x3x3 grid
surv_p = reshape(surv_p,[j,k,l]);
Kwi_3d = reshape(Kwilist,[j,k,l]);
Kwa_3d = reshape(Kwalist,[j,k,l]);
S_3d = reshape(Slist,[j,k,l]);
cdata = reshape(a_sat,[j,k,l]);

a_sat(surv_p < 0.99) = nan;
%find minimum
[m,m_ind] = min(a_sat(:));
figure
%colordata = permute(repmat([255 255 245]'./256,[1,500,500]),[3 2 1]);
% if d == 3
%     s = surf(Smaxgrid,kW1grid,kW2grid,output.cost,1.*ones(length(Smaxgrid), ... 
%         length(kW1grid),3)); %white
%     s = surf(Smaxgrid,kW1grid,kW2grid,output.cost,colordata); %off white
%     s.EdgeColor = 'none';
%     s.FaceColor = 'flat';
    %hold on
[faces,verts,colors] = isosurface(S_3d,Kwi_3d,Kwa_3d,surv_p,0.99,cdata);
%p = patch(s);
%isocolors(S_3d,Kwi_3d,Kwa_3d,cdata,p)
p = patch('Vertices', verts, 'Faces', faces, ...
      'FaceVertexCData', colors, ...
      'FaceColor','interp', ...
      'edgecolor', 'interp')
isonormals(S_3d,Kwi_3d,Kwa_3d,surv_p,p)
view(3);
hold on
scatter3(Slist(m_ind),Kwilist(m_ind),Kwalist(m_ind),m*2,130,'k.', ...
   'LineWidth',1.5,'MarkerEdgeColor','k')

xlabel('Storage Capacity [kWh]')
ylabel('Rated Wind Power [kW]')
zlabel('Rated Wave Power [kW]')
%zlim([-inf Gr2_max])
% ylim([-inf Gr1_max])
% xlim([-inf Sm_max])
c = colorbar;
c.Label.String = '[kg]';
% caxis([0 max(a_sat(:))/2]) %to produce cartoon
% lb = (output.min.cost)/(max(a_sat(:))/2);
% AdvancedColormap('bg l w r',8000,[1*lb,lb+.05*(1-lb),lb+0.1*(1-lb),1])
% AdvancedColormap('kkgw glw lww r',1000, ...
%         [lb,lb+.05*(1-lb),lb+0.15*(1-lb),1]);
% %hold on
% set(gca,'FontSize',fs)
% set(gca,'LineWidth',lw)
grid on

end

