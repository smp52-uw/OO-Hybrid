function visMooringCostMatrix
%visualize mooring cost matrix

close all
saveL = {'SFOMF';'PortHueneme';...
    'PacWave';'MidAtlSB';'BerSea';'altWETS';'altPISCES'};
%'PISCES' 'WETS' are obselete locations - ERA5 had bad data here
figure(1)
tf = tiledlayout(4,2,"TileSpacing","compact");
tf.Padding = 'compact';
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 6, 10])
set(gcf,'Color','w')
le = [0.01, 0.255, 0.5, 0.745, 0.01, 0.255, 0.5];

mincost = nan;
maxcost = nan;
for ll = 1:7
    %load mooring matrix
    clear MoorMat
    filestr = strcat(saveL{ll},'_Mooring.mat');
    load(filestr);

    mincost = min([mincost,min(MoorMat.WorstCase.cost,[],'all')]); %update min cost for scaling the color bar
    maxcost = max([maxcost,max(MoorMat.WorstCase.cost,[],'all')]); %update max cost for scaling the color bar
end
for ll = 1:7
    ax(ll) = nexttile;
    set(gca,'FontSize',8)

    %load mooring matrix
    clear MoorMat
    filestr = strcat(saveL{ll},'_Mooring.mat');
    load(filestr);
    
    cmap = brewermap(9,'Greens'); %colors
    cmap2 = brewermap(6,'OrRd'); %colors

    %make a grid of buoy size and payload mass
    sizemat = size(MoorMat.WorstCase.mass);

    pylmass = reshape(MoorMat.WorstCase.PLmass,[1,sizemat(1)*sizemat(2)]); %vector of payload mass
    pylmass = sort(unique(pylmass));
    buoysz = reshape(MoorMat.WorstCase.dia,[1,sizemat(1)*sizemat(2)]); %vector of dia
    buoysz = sort(unique(buoysz));
    
    [bbss, ttmm] = meshgrid(buoysz,pylmass); %grid of mass vs dia
    cost = nan(size(bbss)); %initialize cost grid

    for b = 1:length(buoysz)
        cost(:,b) = interp1(MoorMat.WorstCase.PLmass(b,:),MoorMat.WorstCase.cost(b,:),pylmass);
    end

    h = heatmap(buoysz,pylmass,cost,'Colormap',cmap,'ColorLimits',[mincost maxcost],'FontSize',6);
    h.MissingDataColor = cmap2(5,:);
    h.MissingDataLabel = 'nan';
    title(saveL(ll))
    xlabel('Buoy Size [m]')
    ylabel('Payload Mass [kg]')
end

%% Line Plot version of mooring space
clear MoorMat
figure(2)
set(gcf,'Units','inches')
%set(gcf,'Position', [1, 1, 6.5, 3.25])
set(gcf,'Position', [1, 1, 11.5, 5.5])
set(gcf,'Color','w')
le = [0.01, 0.255, 0.5, 0.745, 0.01, 0.255, 0.5];
for ll = 1:7
    if ll == 7
        subplot(2,4,7:8);
    else
        subplot(2,4,ll);
    end

    opos = get(gca,'OuterPosition');
    ipos = get(gca,'Position');
    opos(1) = le(ll);
    if ll < 7
       opos(3) = 0.22;
       innerwidth = ipos(3);
       innerLE(ll) = ipos(1);
    end
    set(gca,'OuterPosition',opos)
    %load mooring matrix
    filestr = strcat(saveL{ll},'_Mooring.mat');
    load(filestr);
    
    cmap4 = brewermap(9,'Blues'); %colors
    cmap2 = brewermap(9,'OrRd'); %colors
    %generate the 'clean' mass-cost matrices (removing the submerged
    %points)
    clearance = MoorMat.WorstCase.Sub;
    ind_fail = find(clearance < -0.1);
    ind_surv = find(clearance >= -0.1);
    plm_surv = MoorMat.WorstCase.PLmass;
    cost_surv = MoorMat.WorstCase.cost;
    plm_fail = MoorMat.WorstCase.PLmass;
    cost_fail = MoorMat.WorstCase.cost;

    plm_surv(ind_fail) = nan;
    plm_fail(ind_surv) = nan;
    
    cost_surv(ind_fail) = nan;
    cost_fail(ind_surv) = nan;

    hold on
    b1 = plot(MoorMat.WorstCase.PLmass(1,:),MoorMat.WorstCase.cost(1,:),'.-','Color',cmap4(4,:),'MarkerSize',15,'DisplayName','Small');
    b2 = plot(MoorMat.WorstCase.PLmass(2,:),MoorMat.WorstCase.cost(2,:),'.-','Color',cmap4(6,:),'MarkerSize',15,'DisplayName','Medium');
    b3 = plot(MoorMat.WorstCase.PLmass(3,:),MoorMat.WorstCase.cost(3,:),'.-','Color',cmap4(8,:),'MarkerSize',15,'DisplayName','Large');

    % plot(plm_fail(1,:),cost_fail(1,:),'x','Color',cmap2(4,:),'LineWidth',2,'MarkerSize',10,'DisplayName','Small - Fail')
    % plot(plm_fail(2,:),cost_fail(2,:),'x','Color',cmap2(6,:),'LineWidth',2,'MarkerSize',10,'DisplayName','Medium - Fail')
    % plot(plm_fail(3,:),cost_fail(3,:),'x','Color',cmap2(8,:),'LineWidth',2,'MarkerSize',10,'DisplayName','Large - Fail')
    % 
    % plot(plm_surv(1,:),cost_surv(1,:),'-','Color',cmap4(4,:),'LineWidth',2,'DisplayName','Small - Surv')
    % plot(plm_surv(2,:),cost_surv(2,:),'-','Color',cmap4(6,:),'LineWidth',2,'DisplayName','Medium - Surv')
    % plot(plm_surv(3,:),cost_surv(3,:),'-','Color',cmap4(8,:),'LineWidth',2,'DisplayName','Large - Surv')

    if ll == 7
        hL = legend('show','location','eastoutside');
        Lpos = get(hL,'Position');
        Lpos(1) = innerLE(4);
        set(hL,'Position',Lpos)
        %opos(3) = 0.22 + Lpos(3);
        %set(gca,'OuterPosition',opos)
        ipos(3) = innerwidth + 0.017;
        set(gca,'Position',ipos)
    end
    xlabel('Payload Mass [kg]')
    ylabel('Cost [$]')
    title(saveL(ll))
end
end