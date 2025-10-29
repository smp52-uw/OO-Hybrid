%Plot the optimal 2D solution and the hybrid solution for each location for
%a particular load case

%Tile 1: PISCES LC 3
%Tile 2: BerSea LC 3
%Tile 3: PacWave LC 5
clear
close all

set(0,'defaulttextinterpreter','none')
%set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'calibri')
set(0,'DefaultAxesFontName', 'calibri')

%% merge structures 
folder = "C:\Users\smpal\Documents\OO-Hybrid-Research\NAVFAC_ReportResults\";
task2loc = {'WETS','PISCES','PacWave','MidAtlSB','BerSea'};
pmtext = {"Wind","Solar","Wave","Current"};

temp2D = load('PISCES_PD2PM2_AllLC_FFAOpt_HotelPerFull_08172025');
optStruct2D{1} = temp2D.PISCES_PD2PM2_AllLC_FFAOpt_HotelPerFull_08172025(2);
titlestr{1} = [{"PISCES: "},{"Next-gen OO and UUV"}];

temp2D = load('BerSea_PD2PM1_AllLC_FFAOpt_HotelPerFull_08162025');
optStruct2D{2} = temp2D.BerSea_PD2PM1_AllLC_FFAOpt_HotelPerFull_08162025(2);
titlestr{2} = [{"Bering Sea: "},{"Next-gen OO and UUV"}];

temp2D = load('PacWave_PD2PM1_AllLC_FFAOpt_HotelPerFull_08162025');
optStruct2D{3} = temp2D.PacWave_PD2PM1_AllLC_FFAOpt_HotelPerFull_08162025(3);
titlestr{3} = [{"PacWave: "},{"Constant 200 W"}];

temp6D = load('PISCES_PD6_AllLC_FFAOpt_09262025');
optStruct6D{1} = temp6D.PISCES_PD6_AllLC_FFAOpt_09262025(2);

temp6D = load('BerSea_PD6_LC3_FFA_PSALG_09252025');
temp6D = temp6D.BerSea_PD6_LC3_FFA_PSALG_09252025(2,1,1);
optStruct6D{2} = temp6D;

temp6D = load('PacWave_PD6_AllLC_FFAOpt_09262025');
optStruct6D{3} = temp6D.PacWave_PD6_AllLC_FFAOpt_09262025(3);


%% plotting setup
results = figure;
set(gcf,'Units','inches')
set(gcf, 'Position', [1, 1, 7, 5])

tfig = tiledlayout(1,3);
% tfig.Padding = 'loose';
% tfig.TileSpacing = 'loose';
%colors
addpath(genpath("C:\Users\smpal\Documents\cmasher"))
%wind = green
%solar = yellow
%wave = blue
%current = purple
%hybrid = pink
col{1} = cmasherImport('jungle', 50);
col{2} = cmasherImport('amber', 70);
col{3} = cmasherImport('arctic', 50);
col{4} = cmasherImport('amethyst', 50);
col{5} = cmasherImport('flamingo', 50);

col2D{1} = [col{1}(25,:) ; col{1}(30,:); col{1}(35,:); col{1}(40,:)];
col2D{2} = [col{2}(55,:) ; col{2}(60,:); col{2}(65,:); col{2}(70,:)];
col2D{3} = [col{3}(25,:) ; col{3}(30,:); col{3}(35,:); col{3}(40,:)];
col2D{4} = [col{4}(25,:) ; col{4}(30,:); col{4}(35,:); col{4}(40,:)];
col6D = [col{5}(30,:) ; col{5}(35,:); col{5}(40,:); col{5}(45,:)];


%bar chart settings
MaxGroupWidth = 0.75;
groupOffset = MaxGroupWidth;

for tind = 1:max(size(optStruct2D))
    %unpack cost variables into costdata
    %Cost Types:
    %Gen CapEx
    %Battery CapEx
    %Mooring Cost CapEx
    %Repair OpEx
    %Vessel Cost (for installation and repairs)
    %2D
    SG = optStruct2D{tind}.output.min;

    al_cost = 103.92; %$/sqft for 1/2" 6061 AL
    al_cost = al_cost*0.0929; %convert to $/m2
    batt_vol = polyval(optStruct2D{tind}.opt.p_dev.b_size,SG.Smax{1})*optStruct2D{tind}.econ.batt.volmult;
    batt_len = (batt_vol)^(1/3); %assume cube
    battencl = 2*6*(batt_len^2)*al_cost;
    SG.Icost_batt = 0.1*(SG.Scost+battencl); %installation cost (10% is a guess from Brian)

    costdata(1,1) = SG.CapEx - (SG.Pmooring + SG.Pinst + SG.Pmtrl + SG.battencl + SG.Scost + SG.Icost_batt);
    costdata(2,1) = SG.battencl + SG.Scost + SG.Icost_batt; %
    costdata(3,1) = SG.Pmooring + SG.Pmtrl;
    costdata(4,1) = SG.OpEx - SG.vesselcost;
    %costdata(5,1) = SG.vesselcost + SG.Pinst;
    
    %6D
    HG = optStruct6D{tind}.output.min;

    batt_vol = polyval(optStruct6D{tind}.opt.p_dev.b_size,HG.Smax{1})*optStruct6D{tind}.econ.batt.volmult;
    batt_len = (batt_vol)^(1/3); %assume cube
    battencl = 2*6*(batt_len^2)*al_cost;
    HG.Icost_batt = 0.1*(HG.Scost+battencl); %installation cost (10% is a guess from Brian)

    costdata(1,2) = HG.CapEx - (HG.Pmooring + HG.Pinst + HG.Pmtrl + HG.battencl + HG.Scost + HG.Icost_batt); %
    costdata(2,2) = HG.battencl + HG.Scost + HG.Icost_batt; %
    costdata(3,2) = HG.Pmooring + HG.Pmtrl;
    costdata(4,2) = HG.OpEx - HG.vesselcost;
    %costdata(5,2) = HG.vesselcost + HG.Pinst;
    
    %convert all costs to 1000 $
    costdata = costdata ./1000;
    %plot
    
    ax(tind) = nexttile();
    set(ax,'defaulttextinterpreter','latex','fontsize',10)
    hold on
    
    X = categorical({'Power CapEx','Storage CapEx','Mooring CapEx','Repairs OpEx'});
    X = reordercats(X,{'Power CapEx','Storage CapEx','Mooring CapEx','Repairs OpEx'});
    h = bar(X, costdata, 'FaceColor','flat');
    set(h,'BarWidth',groupOffset);

    set(gca, 'XTickLabel',X,'FontSize',10, 'TickLabelInterpreter', 'latex');
    h(1).CData = col2D{optStruct2D{tind}.opt.pm};
    h(2).CData = col6D;

        
    hold off
    grid on
    box on

    title(titlestr{tind},'interpreter','latex','fontsize',10)

    if tind == 1
        ylabel('Cost [\$ 1000]')
    end

    sgleg = pmtext{optStruct2D{tind}.opt.pm};
    lgd = legend(sgleg, "Hybrid",'location','northeast','Interpreter','Latex','fontsize',10);
    lgd.ItemTokenSize = [10, 20];
end

linkaxes(ax,'y')


