function [] = visPERObjSpace(optStruct)

%Set up a 2D or 3D ObjSpace
%If 2D then pm2 = ~
%diesel=1, inso =2, wind = 3, wave =4

j=129;
k=129;
l=129;
m =129;
n=129;
o=129;
%plot settings
Sm_max = 500; %[kWh]
Gr1_max = 16; %[kW]
Gr2_max = 16; %[kW]
lw = 1.1;
fs = 8;

opt = optStruct.opt;
output = optStruct.output;
%adjust cost to thousands
%output.cost = output.cost{end};
%output.min.cost = output.min.cost;

kW1 = output.Ki_run;
kW2 = output.Kwi_run;
kW4 = output.Kd_run;
kW3 = output.Kwa_run;
kW5 = output.Kc_run;
Smax = output.S_run;

figure


a_sat = output.cost{end};
surv_p = output.surv{end};

kW1(kW4>0 | kW3>0 | kW5>0) = [];
kW2(kW4>0 | kW3>0 | kW5>0) = [];
Smax(kW4>0 | kW3>0 | kW5>0) = [];
a_sat(kW4>0 | kW3>0 | kW5>0) = [];
surv_p(kW4>0 | kW3>0 | kW5>0) = [];

a_sat(surv_p < 0.99) = nan;

kW1(~isfinite(a_sat)) = [];
kW2(~isfinite(a_sat)) = [];
Smax(~isfinite(a_sat))= [];
a_sat(~isfinite(a_sat))= [];
surv_p(~isfinite(a_sat)) = [];
%Trying to get to a mesh that doesn't anger isonormal
%[Kigrid,Kwigrid,Sgrid]= ndgrid(kW1,kW2,Smax);
% [X,Y,Z]= meshgrid(opt.inso.kW{end},opt.wind.kW{end},opt.Smax{end});
% cost_plot = zeros(size(X));
% surv_plot = zeros(size(X));
% for a = 1:j
%     for b = 1:k
%         for c = 1:n
%             
%             x = X(a,b,c);
%             y = Y(a,b,c);
%             z = Z(a,b,c);
%             s_ind = find(kW1 == x & kW2 == y & Smax == z & kW3 == 0 & kW4 == 0 & kW5 == 0);
%             if isempty(s_ind)
%                 surv_plot(a,b,c) = 0;
%                 cost_plot(a,b,c) = inf;
%             else
%                 surv_plot(a,b,c) = surv_p(s_ind);
%                 cost_plot(a,b,c) = a_sat(s_ind);
%             end
%         end
%     end
% end

%find minimum
[min_m,min_ind] = min(a_sat(:));

%colordata = permute(repmat([255 255 245]'./256,[1,500,500]),[3 2 1]);
% if d == 3
%     s = surf(Smaxgrid,kW1grid,kW2grid,output.cost,1.*ones(length(Smaxgrid), ... 
%         length(kW1grid),3)); %white
%     s = surf(Smaxgrid,kW1grid,kW2grid,output.cost,colordata); %off white
%     s.EdgeColor = 'none';
%     s.FaceColor = 'flat';
    %hold on
%disp(it)

%[faces,verts,colors] = isosurface(X,Y,Z,surv_p,0.99,cdata);
%[faces,verts,colors] = isosurface(X,Y,Z,surv_plot,0.99,cost_plot);
%[faces,verts,colors] = isosurface(kW1,kW2,Smax,surv_p,0.99,a_sat);
%isosurface(Ki_3d,Kwi_3d,S_3d,surv_p,0.99);
% p = patch('Vertices', verts, 'Faces', faces, ...
%       'FaceVertexCData', colors, ...
%       'FaceColor','interp', ...
%       'edgecolor', 'interp')
% isonormals(X,Y,Z,surv_plot,p)
colormap(parula(100));
cmap = colormap;
max_c = max(a_sat); %maximum
percent_c = round(100.*a_sat./max_c); %cost as percent of maximum rounded to an integer to be used as index of cmap

scatter3(kW1,kW2,Smax,35,cmap(percent_c,:),'.')
view(3);
hold on
scatter3(kW1(min_ind),kW2(min_ind),Smax(min_ind), 45,130,'ro', ...
   'LineWidth',1.5,'MarkerEdgeColor','r')

%zlim([-inf Gr2_max])
% ylim([-inf Gr1_max])
% xlim([-inf Sm_max])
c = colorbar;
c.Label.String = '[Mass/Maximum Mass]';
% caxis([0 max(a_sat(:))/2]) %to produce cartoon
% lb = (output.min.cost)/(max(a_sat(:))/2);
% AdvancedColormap('bg l w r',8000,[1*lb,lb+.05*(1-lb),lb+0.1*(1-lb),1])
% AdvancedColormap('kkgw glw lww r',1000, ...
%         [lb,lb+.05*(1-lb),lb+0.15*(1-lb),1]);
% %hold on
set(gca,'FontSize',fs)
set(gca,'LineWidth',lw)
grid on
xlabel('Rated Solar Power [kW]')
ylabel('Rated Wind Power [kW]')
zlabel('Storage Capacity [kWh]')

end

