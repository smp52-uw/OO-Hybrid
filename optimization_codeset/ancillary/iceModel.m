function [Iceslow,Icefast] = iceModel(T,U,time,opt)
%compute and apply the ice model
%Inputs
% T = air temp [C]
% U = wind speed [m/s]

%melting temp
Tmelt = -2;
%find indices of temp below -2 C and wind above 12 m/s
logicaln2 = T < -2; %ones when temp is below -2
logical12 = U > 12; %ones when wind is above 12 m/s
iceind = find(logicaln2 & logical12); %indices of both occurances

Icefast = zeros(size(time)); %logical command for ice or not
Iceslow = zeros(size(time)); %logical command for ice or not

%check for no ice
if isempty(iceind)
    return
end
if iceind(1) == 1 %ice on first time step
    Icefast(1) = 1;
    Iceslow(1) = 1;
    accumulate = 1;
else
    accumulate = 0;
end
for t = 2:length(time)
    %are we adding ice this time step
    if sum(iceind == t) %if this index is an icing index
        Icefast(t) = 1;
        Iceslow(t) = 1;
        accumulate = accumulate + 1; %tracking the number of hours of accumulation
    end
    if Icefast(t-1) == 1 && T(t) <= Tmelt %are we maintaining ice this time step?
        Icefast(t) = 1;
    end
    if Iceslow(t-1) == 1 && T(t) <= Tmelt %are we maintaining ice this time step?
        Iceslow(t) = 1;
    end

    if T(t) > Tmelt && Iceslow(t-1) == 1 %melting
        accumulate = accumulate - 1; %remove 1 hour of ice
        if accumulate >= 1 %if not all ice has melted
            Iceslow(t) = 1;
        end
    end
end

%Calculate monthly availability
dtime = datetime(time,'convertfrom','datenum');
c = 0;
for y = 2015:2020 %for year indices
    c = c +1;
    for m = 1:12 %for month indices
        indm{m} = find(month(dtime) == m & year(dtime) == y);
        meanavail{c}(m,1) = sum(1-Iceslow(indm{m}))/length(indm{m});
        meanavail{c}(m,2) = sum(1-Icefast(indm{m}))/length(indm{m});

        moyr{c}{m} = char(strcat(string(y),'-',string(m)));
    end
end




%Ice figure
if opt.pltdebug
    col = colormap(brewermap(9,'Set2')); %colors
    
    %monthly availability
    figure
    tb = tiledlayout(3,2)
    set(gcf,'Position', [100, 100, 1000, 600])
    tb.Padding = 'compact';
    tb.TileSpacing = 'compact';
    for y = 1:6 %for year indices
        nexttile
        X = categorical(moyr{y});
        X = reordercats(X,moyr{y});
        h = bar(X, meanavail{y}, 'FaceColor','flat');
    end
    %ice starting
    figure
    yyaxis left
    plot(U,'linewidth',1.2,'color',col(1,:))
    ylabel('[m/s]')
    title('Ice Accumulation Indices')
    
    yyaxis right
    plot(T,'linewidth',1.2,'Color',col(6,:))
    hold on
    xline(iceind,'r')
    ylabel('[C]')
    xlabel('Time')
    
    %ice models
    clear ax
    figure
    set(gcf,'Position', [100, 100, 1000, 600])
    tf = tiledlayout(2,2);
    tf.Padding = 'compact';
    tf.TileSpacing = 'compact';
    title(tf,"Ice Modeling")
    ax(1) = nexttile; %fast melt
    plot(datetime(time,'convertfrom','datenum'),U,'linewidth',1.2,'color',col(1,:),'DisplayName','Wind Speed')
    hold on
    yline(12,'r','DisplayName','Ice Limit')
    plot(datetime(time,'convertfrom','datenum'),Icefast.*10,'k','linewidth',1.2,'DisplayName','Ice Model')
    ylabel('[m/s]')
    title('Ice Fast Melting Model')
    legend
    
    ax(3) = nexttile(3); %fast melt
    plot(datetime(time,'convertfrom','datenum'),T,'linewidth',1.2,'Color',col(6,:),'DisplayName','Air Temp')
    hold on
    yline(-2,'r','DisplayName','Ice Limit')
    plot(datetime(time,'convertfrom','datenum'),Icefast.*10,'k','linewidth',1.2,'DisplayName','Ice Model')
    xlabel('Time')
    ylabel('[C]')
    legend
    
    
    ax(2) = nexttile(2); %slow melt
    plot(datetime(time,'convertfrom','datenum'),U,'linewidth',1.2,'color',col(1,:),'DisplayName','Wind Speed')
    hold on
    yline(12,'r','DisplayName','Ice Limit')
    plot(datetime(time,'convertfrom','datenum'),Iceslow.*10,'k','linewidth',1.2,'DisplayName','Ice Model')
    ylabel('[m/s]')
    title('Ice Slow Melting Model')
    legend
    
    ax(4) = nexttile(4); %slow melt
    plot(datetime(time,'convertfrom','datenum'),T,'linewidth',1.2,'Color',col(6,:),'DisplayName','Air Temp')
    hold on
    yline(-2,'r','DisplayName','Ice Limit')
    plot(datetime(time,'convertfrom','datenum'),Iceslow.*10,'k','linewidth',1.2,'DisplayName','Ice Model')
    ylabel('[C]')
    xlabel('Time')
    legend
    linkaxes(ax,'x')
end
end