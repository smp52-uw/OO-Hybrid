function [] = visHybridObjSpace(optStruct,d,pm1,pm2)

%Set up a 2D or 3D ObjSpace
%If 2D then pm2 = ~
%diesel=1, inso =2, wind = 3, wave =4
j=9;
k=9;
l=9;
m =9;
n=9;

%plot settings
Sm_max = 500; %[kWh]
Gr1_max = 16; %[kW]
Gr2_max = 16; %[kW]
lw = 1.1;
fs = 8;

opt = optStruct.opt;
output = optStruct.output;
Smax = opt.Smax;
%adjust cost to thousands
%output.cost = output.cost{end};
%output.min.cost = output.min.cost;

if pm1 == 1
    kW1 = opt.dies.kW;
elseif pm1 == 2
    kW1 = opt.inso.kW;
elseif pm1 == 3
    kW1 = opt.wind.kW;
elseif pm1 == 4
    kW1 = opt.wave.kW;
end
if pm2 == 1
    kW2 = opt.dies.kW;
elseif pm2 == 2
    kW2 = opt.inso.kW;
elseif pm2 == 3
    kW2 = opt.wind.kW;
elseif pm2 == 4
    kW2 = opt.wave.kW;
end
kW3 = opt.dies.kW;
kW4 = opt.wave.kW;

it_max = size(output.cost);
figure
tiledlayout(3,2)
for it = 1:it_max(2)
    
    %create grid
    %[Kd,Ki,Kwi,Kwa,S]
    [Kdgrid,Kigrid,Kwigrid,Kwagrid,Sgrid]= ndgrid(kW3{it},kW1{it},kW2{it},kW4{it},Smax{it});
    Kdlist = reshape(Kdgrid,[j*k*l*m*n 1]);
    Kilist = reshape(Kigrid,[j*k*l*m*n 1]);
    Kwilist = reshape(Kwigrid,[j*k*l*m*n 1]);
    Kwalist = reshape(Kwagrid,[j*k*l*m*n 1]);
    Slist = reshape(Sgrid,[j*k*l*m*n 1]);
    a_sat = output.cost{it};
    surv_p = output.surv{it};
    %remove all output data other than the dimensions we want
% %     surv_p = output.surv{it};
% %     surv_p(Kdlist~=0 | Kwalist~=0) = [];
% %     Kwilist(Kdlist~=0 | Kwalist~=0) = [];
% %     Kilist(Kdlist~=0 | Kwalist~=0) = [];
% %     Slist(Kdlist~=0 | Kwalist~=0) = [];
% %     a_sat(Kdlist~=0 | Kwalist~=0) = [];
% % %     %find minimum
% % %     [min_m,min_ind] = min(a_sat(:));
% %     
% %     %reshape into 3x3x3 grid
% %     surv_p = reshape(surv_p,[j,k,l]);
% %     Kwi_3d = reshape(Kwilist,[j,k,l]);
% %     Ki_3d = reshape(Kilist,[j,k,l]);
% %     S_3d = reshape(Slist,[j,k,l]);
% %     cdata = reshape(a_sat,[j,k,l]);
    
    %Trying to get to a mesh that doesn't anger isonormal
    %[Kigrid,Kwigrid,Sgrid]= ndgrid(kW1,kW2,Smax);
    [X,Y,Z]= meshgrid(kW1{it},kW2{it},Smax{it});
    cost_plot = zeros(size(X));
    surv_plot = zeros(size(X));
    for a = 1:j
        for b = 1:k
            for c = 1:n
                
                x = X(a,b,c);
                y = Y(a,b,c);
                z = Z(a,b,c);
                s_ind = find(Kilist == x & Kwilist == y & Slist == z & Kwalist == 0 & Kdlist == 0);
                surv_plot(a,b,c) = surv_p(s_ind);
                cost_plot(a,b,c) = a_sat(s_ind);
            end
        end
    end
    
    a_sat(surv_p < 0.99) = nan;
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
    nexttile
    %[faces,verts,colors] = isosurface(X,Y,Z,surv_p,0.99,cdata);
    [faces,verts,colors] = isosurface(X,Y,Z,surv_plot,0.99,cost_plot);
    %isosurface(Ki_3d,Kwi_3d,S_3d,surv_p,0.99);
    p = patch('Vertices', verts, 'Faces', faces, ...
          'FaceVertexCData', colors, ...
          'FaceColor','interp', ...
          'edgecolor', 'interp')
    isonormals(X,Y,Z,surv_plot,p)
    view(3);
    hold on
    scatter3(Kilist(min_ind),Kwilist(min_ind),Slist(min_ind), m*5,130,'ro', ...
       'LineWidth',1.5,'MarkerEdgeColor','r')
    
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
    set(gca,'FontSize',fs)
    set(gca,'LineWidth',lw)
    grid on
    xlabel('Rated Solar Power [kW]')
    ylabel('Rated Wind Power [kW]')
    zlabel('Storage Capacity [kWh]')
end

end

