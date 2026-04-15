%you'll need the cmasherImport function added to path or the color map
%won't generate

selectedfolder = uigetdir();

fileList = dir(fullfile(selectedfolder,'*.mat'));
numfiles = max(size(fileList));

%separate by location, load case, power module
for i = 1:numfiles
    clear optStruct opt 
    tmp = load(fullfile(selectedfolder,fileList(i).name));
    nm = split(fileList(i).name,'.');
    nm = nm(1);
    optStruct = tmp.(nm{1});
    
    opt = optStruct.opt;
    
    %adjust cost to thousands
    cost = optStruct.output.cost{1}/1000;
    surv = optStruct.output.surv{1};
    dsurv = surv - 0.99; %this will be negative for bad points and positive for overdesigned points
    costW = cost;
   
    costW(surv<0.99) = nan;
    %create grid
    disc = 500;
    Smax = linspace(opt.Smax_1,opt.Smax_n,2*disc).';
    Srun = optStruct.output.S_run{1};
    Smaxmin = optStruct.output.min.Smax{1};
    costmin = optStruct.output.min.cost;
    if opt.pm == 1
        kW = linspace(opt.wind.kW_1,opt.wind.kW_m,2*disc); 
        kWmin = optStruct.output.min.kWwi{1};
        kWrun = optStruct.output.Kwi_run{1};
    elseif opt.pm == 2
        kW = linspace(opt.inso.kW_1,opt.inso.kW_m,2*disc);
        kWmin = optStruct.output.min.kWi{1};
        kWrun = optStruct.output.Ki_run{1};
    elseif opt.pm == 3
        kW = linspace(opt.wave.kW_1,opt.wave.kW_m,2*disc);
        kWmin = optStruct.output.min.kWwa{1};
        kWrun = optStruct.output.Kwa_run{1};
    elseif opt.pm == 4
        kW = linspace(opt.dies.kW_1,opt.dies.kW_m,2*disc);   
        kWmin = optStruct.output.min.kWd{1};
        kWrun = optStruct.output.Kd_run{1};
    else
        kW = linspace(opt.curr.kW_1,opt.curr.kW_m,2*disc).'; 
        kWmin = optStruct.output.min.kWc{1};
        kWrun = optStruct.output.Kc_run{1};
    end
    
    [kWgrid2,Smaxgrid2] = ndgrid(kW,Smax);

    %identify the + 2% cost space
    cost2 = costW;
    cost2(cost2 > min(cost2)*1.02) = nan;

    kWgrid = reshape(kWrun,[disc,500]);
    Smaxgrid = reshape(Srun,[disc,500]);
    costgrid2 = reshape(cost2,[disc,500]);
    survgrid = reshape(dsurv,[disc,500]);

    Vq = interpn(kWgrid,Smaxgrid,costgrid2,kWgrid2,Smaxgrid2,'nearest');
    %package output
    plotdata{i}.costW = costW;
    plotdata{i}.cost2 = cost2;
    plotdata{i}.costmin = costmin;
    plotdata{i}.surv = surv;
    plotdata{i}.kWmin = kWmin;
    plotdata{i}.kWrun = kWrun;
    plotdata{i}.Srun = Srun;
    plotdata{i}.Smaxmin = Smaxmin;
    plotdata{i}.kWgrid = kWgrid;
    plotdata{i}.Smaxgrid = Smaxgrid;
    plotdata{i}.costgrid2 = costgrid2;
    plotdata{i}.survgrid = survgrid;

    plotdata{i}.kWgrid2 = kWgrid2;
    plotdata{i}.Smaxgrid2 = Smaxgrid2;
    plotdata{i}.costgridint = Vq;

    pm(i) = optStruct.opt.pm;
    loc{i} = optStruct.loc;
    lc(i) = optStruct.uc.loadcase;
end

locoptions = {'PacWave','MidAtlSB','BerSea','altWETS','altPISCES'};
locdisp =  {'Newport','Mid-Atlantic Shelf Break','Bering Sea','O''ahu','Coastal Washington'};
loadoptions = [1 3 5];


cMap = cmasherImport('prinsenvlag',1000);

%% plot data
for i = 1:5
    figure
    tf = tiledlayout(length(locoptions),length(loadoptions));
    tf.Padding = 'compact';
    tf.TileSpacing = 'compact';
    set(gcf,'Units','inches')
    set(gcf, 'Position', [1, 1, 7, 10])
    title(tf,strcat("Power Module",string(i)))
    
    for ll = 1:length(locoptions)
        for uu = 1:length(loadoptions)
            indloc = strcmp(loc, locoptions{ll});
            indload = lc == loadoptions(uu);
            indpm = pm == i;
            inds = find(indloc & indload & indpm);
    
            titlestr = strcat(locdisp{ll},", LC = ",string(loadoptions(uu)));
            ax = nexttile;
            j = inds;

            s = surf(plotdata{j}.Smaxgrid,plotdata{j}.kWgrid,plotdata{j}.survgrid);
            colormap(gca, cMap)
            s.EdgeColor = "none";
            hold on
            plot3(plotdata{j}.Smaxmin,plotdata{j}.kWmin,5*max(plotdata{j}.cost2),'o','LineWidth',1.2,'MarkerSize',5,'Color','k');
            collabel = '$a_{sim}$ - 0.99';
            view(0,90)
            xlabel('Storage Capacity [kWh]')
            ylabel('Rated Power [kW]')
            c = colorbar;
            c.Label.String = collabel;
            ylim([0,8])

            title(titlestr)
            grid on
        end
    end
end

