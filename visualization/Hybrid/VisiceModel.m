function VisiceModel(data,opt)
colT = colormap(brewermap(12,'Set3')); %colors
colW = colormap(brewermap(9,'YlGn')); %colors
colS = colormap(brewermap(9,'YlOrBr')); %colors

    swsoice = data.swso.*opt.ice_ts;
    Uice = data.met.wind_spd.*opt.ice_ts;
    U = data.met.wind_spd;
    time = data.met.time;
    T = data.temperature.T2mw;
    %ice models
    clear ax
    figure
    set(gcf,'Position', [100, 100, 400, 600])
    tf = tiledlayout(3,1);
    tf.Padding = 'compact';
    tf.TileSpacing = 'compact';
   
    ax(1) = nexttile(1); %slow melt
    plot(datetime(time,'convertfrom','datenum'),T,'linewidth',1.2,'Color',colT(5,:),'DisplayName','Air Temp')
    hold on
    yline(-2,'k','linewidth',1.2,'DisplayName','Ice Trigger')
    ylabel('Air Temperature [C]','Interpreter','Latex')
    legend('Location','BestOutside','Interpreter','Latex')
    linkaxes(ax,'x')
    
    grid on

    ax(2) = nexttile(2); %slow melt
    plot(datetime(time,'convertfrom','datenum'),U,'linewidth',1.2,'color',colW(8,:),'DisplayName','Without Ice Model')
    hold on
    plot(datetime(time,'convertfrom','datenum'),Uice,'linewidth',1.2,'color',colW(5,:),'DisplayName','With Ice Model')
    yline(12,'k','linewidth',1.2,'DisplayName','Ice Trigger')
    ylabel('Wind Speed [m/s]','Interpreter','Latex')
    legend('Location','BestOutside','Interpreter','Latex')
    grid on

    ax(3) = nexttile(3); %slow melt
    plot(datetime(time,'convertfrom','datenum'),data.swso,'linewidth',1.2,'color',colS(5,:),'DisplayName','Without Ice Model')
    hold on
    plot(datetime(time,'convertfrom','datenum'),swsoice,'linewidth',1.2,'color',colS(3,:),'DisplayName','With Ice Model')
    ylabel('Solar Radiation [$W/m^2$]','Interpreter','Latex')
    legend('Location','BestOutside','Interpreter','Latex')
    grid on
    xlabel('Time','Interpreter','Latex')
    linkaxes(ax,'x')
    
end