function [] = visP2TObjSpace(optStruct)

%Set up a 2D or 3D ObjSpace
%If 2D then pm2 = ~
%diesel=1, inso =2, wind = 3, wave =4

% j=129;
% k=129;
% l=129;
% m =129;
% n=129;
% o=129;
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

% kW1 = output.Ki_run;
% kW2 = output.Kwi_run;
% kW4 = output.Kd_run;
% kW3 = output.Kwa_run;
% kW5 = output.Kc_run;
% Smax = output.S_run;

%append all points to long arrays
% kW1 = [output.Ki_run{1,1}];
% kW2 = [output.Kwi_run{1,1}];
% kW3 = [output.Kd_run{1,1}];
% kW4 = [output.Kwa_run{1,1}];
% kW5 = [output.Kc_run{1,1}];
% Smax = [output.S_run{1,1}];
% a_sat = [output.cost{1,1}];
% surv_p = [ output.surv{1,1}];
% for i=2:5
%     disp(i)
%     kW1 = [kW1; output.Ki_run{1,i}'];
%     kW2 = [kW2; output.Kwi_run{1,i}'];
%     kW3 = [kW3; output.Kd_run{1,i}'];
%     kW4 = [kW4; output.Kwa_run{1,i}'];
%     kW5 = [kW5; output.Kc_run{1,i}'];
%     Smax = [Smax; output.S_run{1,i}'];
%     a_sat = [a_sat; output.cost{1,i}];
%     surv_p = [surv_p; output.surv{1,i}];
% end
figure(1)
tiledlayout(3,2)
for i=1:5
    nexttile
    
    kW1 = output.Ki_run{1,i};
    kW2 = output.Kwi_run{1,i};
    kW4 = output.Kd_run{1,i};
    kW3 = output.Kwa_run{1,i};
    kW5 = output.Kc_run{1,i};
    Smax = output.S_run{1,i};
    a_sat = output.cost{1,i};
    surv_p = output.surv{1,i};
    
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
    
    [min_m,min_ind] = min(a_sat(:));
    
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


figure(2)
tiledlayout(3,2)
for i=1:5
    figure(1+i)
    tiledlayout(3,2)
    
    kW1 = output.Ki_run{1,i};
    kW2 = output.Kwi_run{1,i};
    kW4 = output.Kd_run{1,i};
    kW3 = output.Kwa_run{1,i};
    kW5 = output.Kc_run{1,i};
    Smax = output.S_run{1,i};
    a_sat = output.cost{1,i};
    surv_p = output.surv{1,i};
    
    kW1(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    kW2(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    kW3(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    kW4(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    kW5(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    Smax(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    a_sat(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    surv_p(output.Kd_run{1,i}>0 | output.Kwa_run{1,i}>0 | output.Kc_run{1,i}>0) = [];
    
    a_sat(surv_p < 0.99) = nan;
    
    kW1(~isfinite(a_sat)) = [];
    kW2(~isfinite(a_sat)) = [];
    kW3(~isfinite(a_sat)) = [];
    kW4(~isfinite(a_sat)) = [];
    kW5(~isfinite(a_sat)) = [];
    Smax(~isfinite(a_sat))= [];
    a_sat(~isfinite(a_sat))= [];
    surv_p(~isfinite(a_sat)) = [];
    
    [min_m,min_ind] = min(a_sat(:));
    
    colormap(parula(100));
    cmap = colormap;
    max_c = max(a_sat); %maximum
    percent_c = round(100.*a_sat./max_c); %cost as percent of maximum rounded to an integer to be used as index of cmap
    kW_y = {kW2, kW3, kW4, kW5}
    for j=1:4
        nexttile
        scatter3(kW1,kW_y{j},Smax,35,cmap(percent_c,:),'.')
        view(3);
        hold on
        scatter3(kW1(min_ind),kW_y{j}(min_ind),Smax(min_ind), 45,130,'ro', ...
           'LineWidth',1.5,'MarkerEdgeColor','r')
        
        c = colorbar;
        c.Label.String = '[Mass/Maximum Mass]';
        % %hold on
        set(gca,'FontSize',fs)
        set(gca,'LineWidth',lw)
        grid on
        xlabel('Rated Solar Power [kW]')
        Gr_t = {'Wind','Wave','Diesel','Current'}
        ylabel(['Rated ',Gr_t{j},'  Power [kW]'])
        zlabel('Storage Capacity [kWh]')
    end
end
end
