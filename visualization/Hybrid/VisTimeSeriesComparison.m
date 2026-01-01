%Time domain comparison between 2D and Hybrid

addpath(genpath("C:\Users\smpal\Documents\OO-Hybrid-Research\WECLCMResults\6DWEC_LCM"))
addpath(genpath("C:\Users\smpal\Documents\OO-Hybrid-Research\NAVFAC_ReportResults"))
addpath(genpath("C:\Users\smpal\Documents\OO-Hybrid-Research\WECLCMResults\5D_NoMC"))

%Figure 1: Bering Sea - LC=3
% temp2D = load('BerSea_PD2PM1_AllLC_FFAOpt_HotelPerFull_08162025');
% optStruct{1} = temp2D.BerSea_PD2PM1_AllLC_FFAOpt_HotelPerFull_08162025(2);
% titlestr{1} = "Optimal single generator";
% 
% temp6D = load('BerSea_PD6_LC3_FFA_PSALG_09252025');
% optStruct{2} = temp6D.BerSea_PD6_LC3_FFA_PSALG_09252025(2,1,1);
% titlestr{2} = "Base hybrid";
% 
% 
% temp6DLCM = load('BerSea_PD6_AllLC_FFAOpt_Low1010MassCostWEC_09302025');
% optStruct{3} = temp6DLCM.BerSea_PD6_AllLC_FFAOpt_Low1010MassCostWEC_09302025(2);
% titlestr{3} = "Low WEC mass and cost hybrid";

% %Figure 2: MidAtlSB - LC=5
temp2D = load('MidAtlSB_PD2PM2_AllLC_BF500_HotelPerFull_07302025');
optStruct{1} = temp2D.MidAtlSB_PD2PM2_AllLC_BF500_HotelPerFull_07302025(3);
titlestr{1} = "Optimal single generator";

temp6D = load('MidAtlSB_PD6_LC5_FFA_PSALG_09232025');
optStruct{2} = temp6D.MidAtlSB_PD6_LC5_FFA_PSALG_09232025(2,1,1,2);
titlestr{2} = "Base hybrid";


temp5D = load('MidAtlSB_PD5PM2_AllLC_FFAOpt_10032025');
optStruct{3} = temp5D.MidAtlSB_PD5PM2_AllLC_FFAOpt_10032025(3);
titlestr{3} = "Low surface expression hybrid";

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
    tend = 24*90;
    tstart = 250*24;
    t = tstart:1:tstart+tend;

    ax(s) = nexttile(s); %load
    plot(t,out.L(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',colLD(1,:))
    xlim([tstart,tstart+tend])
    grid on
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

    ax(s+3) = nexttile(s+3); %solar
    xlim([tstart,tstart+tend])
    if sum(out.Pinso) > 0
        plot(t,out.Pinso(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',col2D{2}(1,:))
        cs = cs + 1;
        maxPV(cs) =  max(out.Pinso(tstart:tstart+tend)./1000);
        grid on
    else
        text(.5,.5,'No Solar Panels','Units','Normalized', ...
            'VerticalAlignment','middle',"horizontalAlignment",'center','FontWeight','normal', ...
            'FontSize',fs,'Interpreter','latex');
        box on

    end
    if s == 1
        yl2 = ylabel([{'$P_{\mathrm{PV}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yl2,'Units','inches');
        ypos2 = get(yl2,'Position');
        ypos2(1) = ypos(1);
        ypos2(2) = ypos(2);
        set(yl2,'Position',ypos2)
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+6) = nexttile(s+6); %wind
    xlim([tstart,tstart+tend])
    if sum(out.Pwind) > 0
        plot(t,out.Pwind(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',col2D{1}(1,:))
        cwi = cwi + 1;
        maxWI(cwi) = max(out.Pwind(tstart:tstart+tend)./1000);
        grid on
    else
        text(.5,.5,'No Wind Turbine','Units','Normalized', ...
            'VerticalAlignment','middle',"horizontalAlignment",'center','FontWeight','normal', ...
            'FontSize',fs,'Interpreter','latex');
        box on
    end
    if s== 1
        yl3 = ylabel([{'$P_{\mathrm{wind}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yl3,'Units','inches');
        ypos3 = get(yl3,'Position');
        ypos3(1) = ypos(1);
        ypos3(2) = ypos(2);
        set(yl3,'Position',ypos3)
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)


    % ax(s+9) = nexttile(s+9); %wave
    % if sum(out.Pwave) > 0
    %     plot(t,out.Pwave(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',col2D{3}(1,:))
    %     cwa = cwa + 1;
    %     maxWA(cwa) = max(out.Pwave(tstart:tstart+tend)./1000);
    %     grid on
    % 
    % else
    %     text(.5,.5,'No WEC','Units','Normalized', ...
    %         'VerticalAlignment','middle',"horizontalAlignment",'center','FontWeight','normal', ...
    %         'FontSize',fs,'Interpreter','latex');
    %     box on
    % end
    % if s == 1
    %     yl4 = ylabel([{'$P_{\mathrm{wave}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
    %     set(yl4,'Units','inches');
    %     ypos4 = get(yl4,'Position');
    %     ypos4(1) = ypos(1);
    %     ypos4(2) = ypos(2);
    %     set(yl4,'Position',ypos4)
    % end
    % set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+9) = nexttile(s+9); %current
    xlim([tstart,tstart+tend])
    if out.kWc{1} > 0.005
        plot(t,out.Pcurr(tstart:tstart+tend)./1000,'LineWidth',1.5,'Color',col2D{4}(1,:))
        cc = cc + 1;
        maxCC(cc) = max(out.Pcurr(tstart:tstart+tend)./1000);
        grid on

    else
        text(.5,.5,'No Current Turbine','Units','Normalized', ...
            'VerticalAlignment','middle',"horizontalAlignment",'center','FontWeight','normal', ...
            'FontSize',fs,'Interpreter','latex');
        box on
    end
    if s == 1
        yl4 = ylabel([{'$P_{\mathrm{current}}(t)$'},{'[kW]'}],'Interpreter','latex','rotation',0,'HorizontalAlignment','center','FontSize',fs);
        set(yl4,'Units','inches');
        ypos4 = get(yl4,'Position');
        ypos4(1) = ypos(1);
        ypos4(2) = ypos(2);
        set(yl4,'Position',ypos4)
    end
    set(gca,'TickLabelInterpreter','latex','fontsize',fs)

    ax(s+12) = nexttile(s+12); %storage
    xlim([tstart,tstart+tend])
    plot(t,out.S1((tstart:tstart+tend))./1000,'LineWidth',1.5,'Color',colS(1,:))
    grid on
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
    xlim([tstart,tstart+tend])
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
linkaxes([ax(1),ax(4),ax(7),ax(10),ax(16)],'x')

ylim(ax(4),[0,max(maxPV)])
ylim(ax(7),[0,max(maxWI)])
ylim(ax(1),[0,max(out.L(tstart:tstart+tend)./1000)])
ylim(ax(10),[0,max(maxCC)])
%ylim(ax(10),[0,max(maxWA)])
