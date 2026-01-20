%Time domain comparison between 2D and Hybrid
% clear
% clc
% close all

addpath(genpath("C:\Users\smpal\Documents\OO-Hybrid-Research\NAVFAC_ReportResults"))
addpath(genpath("C:\Users\smpal\MREL Dropbox\Sarah Palmer\OO-Hybrid\IceModelTesting\JulStart_neg2degMelt"))

% %2D wind
% tempNI = load('BerSea_LC3_Wi257_S83_NoIce_01122026');
% optStruct{1} = tempNI;
% titlestr{1} = "2D Wind LC3";
% 
% tempFI = load('BerSea_LC3_Wi257_S83_FastMelt_01112026');
% optStruct{2} = tempFI;
% titlestr{2} = "Fast Melt";
% 
% 
% tempSI = load('BerSea_LC3_Wi257_S83_SlowMelt_01112026');
% optStruct{3} = tempSI;
% titlestr{3} = "Slow Melt";

% % 2D solar
% tempNI = load('BerSea_LC3_In8_S125_NoIce_01122026');
% optStruct{1} = tempNI;
% titlestr{1} = "2D Solar LC5";
% 
% tempFI = load('BerSea_LC5_In8_S125_FastMelt_01112026');
% optStruct{2} = tempFI;
% titlestr{2} = "Fast Melt";
% 
% 
% tempSI = load('BerSea_LC5_In8_S125_SlowMelt_01112026');
% optStruct{3} = tempSI;
% titlestr{3} = "Slow Melt";

% % Hybrid standard LC 3
% tempNI = load('BerSea_LC3_In219_Wi072_S71_NoIce_01122026');
% optStruct{1} = tempNI;
% titlestr{1} = "Standard Hybrid LC3";
% 
% tempFI = load('BerSea_LC3_In219_Wi072_S71_FastMelt_01112026');
% optStruct{2} = tempFI;
% titlestr{2} = "Fast Melt";
% 
% 
% tempSI = load('BerSea_LC3_In219_Wi072_S71_SlowMelt_01112026');
% optStruct{3} = tempSI;
% titlestr{3} = "Slow Melt";

% % Hybrid standard LC 5
tempNI = load('BerSea_LC3_In138_Wi056_S69_NoIce_01122026');
optStruct{1} = tempNI;
titlestr{1} = "Standard Hybrid LC5";

tempFI = load('BerSea_LC5_In138_Wi056_S69_FastMelt_01112026');
optStruct{2} = tempFI;
titlestr{2} = "Fast Melt";

tempSI = load('BerSea_LC5_In138_Wi056_S69_SlowMelt_01112026');
optStruct{3} = tempSI;
titlestr{3} = "Slow Melt";

% % Hybrid solar/wec LC 3
% tempNI = load('BerSea_LC3_In283_Wa092_S72_NoIce_01042026');
% optStruct{1} = tempNI;
% titlestr{1} = "Solar/WEC Hybrid LC3";
% 
% tempFI = load('BerSea_LC3_In283_Wa092_S72_FastMelt_01042026');
% optStruct{2} = tempFI;
% titlestr{2} = "Fast Melt";
% 
% tempSI = load('BerSea_LC3_In283_Wa092_S72_SlowMelt_01042026');
% optStruct{3} = tempSI;
% titlestr{3} = "Slow Melt";

%% Figure
figure
set(gcf,'Units','inches')
set(gcf, 'Position', [1, 1, 7, 5.5])
fs = 10;
szst = max(size(optStruct));
tfig = tiledlayout(6,szst);
%tfig.Padding = 'compact';
tfig.TileSpacing = 'compact';

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
col{6} = cmasherImport('neutral', 50);

col2D{1} = [col{1}(25,:) ; col{1}(30,:); col{1}(35,:); col{1}(40,:)];
col2D{2} = [col{2}(55,:) ; col{2}(60,:); col{2}(65,:); col{2}(70,:)];
col2D{3} = [col{3}(25,:) ; col{3}(30,:); col{3}(35,:); col{3}(40,:)];
col2D{4} = [col{4}(20,:) ; col{4}(30,:); col{4}(35,:); col{4}(40,:)];
colS = [col{5}(30,:) ; col{5}(35,:); col{5}(40,:); col{5}(45,:)];
colLD = [col{6}(1,:) ; col{6}(20,:); col{6}(40,:); col{6}(45,:)];

cs = 0; cwi = 0; cwa = 0; cc = 0;
for s = 1:szst

    out = optStruct{s}.output.min;
    opt = optStruct{s}.opt;
    tend = length(out.L)-1;
    tstart = 1;
    %t = tstart:1:tstart+tend;
    t = datetime(optStruct{s}.opt.wave.time,'ConvertFrom','datenum');

    ax(s) = nexttile(s); %load
    plot(t,out.L(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',colLD(1,:))
    %xlim([t(tstart),t(tstart+tend)])
    grid on
    xticks([])
    if s == 1
        yl = ylabel([{'L(t)'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yl,'Units','inches');
        ypos = get(yl,'Position');
        ypos(1) = ypos(1)*3.1; %left (2.9 for fig 1, 3.1 for fig 2)
        ypos(2) = ypos(2)*0.3; %height
        set(yl,'Position',ypos)
    end
    title(titlestr{s},'interpreter','latex','fontsize',fs)
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+3) = nexttile(s+3); %ice
    if ~isfield(opt,'ice_ts') %if the sim was run before the ice model then the equivalent is all ones
        ice = ones(size(out.L));
    else
        ice = opt.ice_ts;
    end
    plot(t,ice(tstart:tstart+tend),'LineWidth',1.5,'Color',colLD(1,:))
    %xlim([t(tstart),t(tstart+tend)])
    grid on
    xticks([])
    if s == 1
        yli = ylabel([{'Ice'},{'[~]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yli,'Units','inches');
        yposi = get(yli,'Position');
        yposi(1) = ypos(1);
        yposi(2) = ypos(2);
        set(yli,'Position',yposi)
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+6) = nexttile(s+6); %solar
    %xlim([t(tstart),t(tstart+tend)])
    if sum(out.Pinso) > 0
        plot(t,out.Pinso(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',col2D{2}(1,:))
        cs = cs + 1;
        maxPV(cs) =  max(out.Pinso(tstart:tstart+tend)./1000);
        grid on
    else
        plot(t,nan(size(t)),'LineWidth',1.5,'Color',col2D{3}(1,:))
        text(.5,.5,'No Solar Panels','Units','Normalized', ...
            'VerticalAlignment','middle',"horizontalAlignment",'center','FontWeight','normal', ...
            'FontSize',fs,'Interpreter','latex');
        box on

    end
    xticks([])
    if s == 1
        yl2 = ylabel([{'$P_{\mathrm{PV}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yl2,'Units','inches');
        ypos2 = get(yl2,'Position');
        ypos2(1) = ypos(1);
        ypos2(2) = ypos(2);
        set(yl2,'Position',ypos2)
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+9) = nexttile(s+9); %wind or wave
    %xlim([t(tstart),t(tstart+tend)])
    if sum(out.Pwind) > 0
        plot(t,out.Pwind(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',col2D{1}(1,:))
        cwi = cwi + 1;
        maxW(cwi) = max(out.Pwind(tstart:tstart+tend)./1000);
        grid on
    elseif sum(out.Pwave) > 0
        plot(t,out.Pwave(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',col2D{3}(1,:))
        cwa = cwa + 1;
        maxW(cwa) = max(out.Pwave(tstart:tstart+tend)./1000);
        grid on
    else
        plot(t,nan(size(t)),'LineWidth',1.5,'Color',col2D{3}(1,:))
        text(.5,.5,'No Wind or Wave','Units','Normalized', ...
            'VerticalAlignment','middle',"horizontalAlignment",'center','FontWeight','normal', ...
            'FontSize',fs,'Interpreter','latex');
        box on
    end
    xticks([])
    if s== 1
        if sum(out.Pwind) > 0
            yl3 = ylabel([{'$P_{\mathrm{wind}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        elseif sum(out.Pwave) > 0
            yl3 = ylabel([{'$P_{\mathrm{wave}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        else
            yl3 = ylabel([{'$P_{\mathrm{w}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        end
        set(yl3,'Units','inches');
        ypos3 = get(yl3,'Position');
        ypos3(1) = ypos(1);
        ypos3(2) = ypos(2);
        set(yl3,'Position',ypos3)
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+12) = nexttile(s+12); %storage
    %xlim([t(tstart),t(tstart+tend)])
    plot(t,out.S1((tstart:tstart+tend))./1000,'LineWidth',1.5,'Color',colS(1,:))
    hold on
    plot(t,out.S2((tstart:tstart+tend))./1000,'LineWidth',1.5,'Color',colS(3,:))
    grid on
    xticks([])
    if s == 1
        yl5 = ylabel([{'S(t)'},{'[kWh]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yl5,'Units','inches');
        ypos5 = get(yl5,'Position');
        ypos5(1) = ypos(1);
        ypos5(2) = ypos(2);
        set(yl5,'Position',ypos5)
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+15) = nexttile(s+15); %dumped
    %xlim([t(tstart),t(tstart+tend)])
    plot(t,out.D(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',colLD(3,:))
    grid on
    if s == 1
        yl6 = ylabel([{'D(t)'},{'[kWh]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yl6,'Units','inches');
        ypos6 = get(yl6,'Position');
        ypos6(1) = ypos(1);
        ypos6(2) = ypos(2);
        set(yl6,'Position',ypos6)
    end
    xlabel('Time [hr]','interpreter','latex','fontsize',fs)
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)
end


linkaxes([ax(1),ax(2),ax(3)], 'xy')
linkaxes([ax(4),ax(5),ax(6)], 'xy')
linkaxes([ax(7),ax(8),ax(9)], 'xy')
linkaxes([ax(10),ax(11),ax(12)], 'xy')
linkaxes([ax(13),ax(14),ax(15)], 'xy')
linkaxes([ax(16),ax(17),ax(18)], 'xy')
linkaxes(ax,'x')

% ylim(ax(4),[0,max(maxPV)])
% ylim(ax(7),[0,max(maxW)])
% ylim(ax(1),[0,max(out.L(tstart:tstart+tend)./1000)])
