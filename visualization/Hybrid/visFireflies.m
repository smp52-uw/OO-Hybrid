function visFireflies(results,filename)
%make animation of fireflies
c1 = brewermap(9,'Set3');
coptpm(1,:) = c1(1,:);

coptpm(2,:) = c1(6,:);
coptpm(3,:) = c1(5,:);

c4 = brewermap(9,'Set2');
coptpm(4,:) = c4(8,:);

c5 = brewermap(9,'PRGn');
coptpm(5,:) = c5(4,:);
coptpm(6,:) = c4(4,:);


%make real populations (what I saved is newpop)
pop(1).cost = results.output.cost{1};
pop(1).Kwi = results.output.Kwi_run{1};
pop(1).Ki = results.output.Ki_run{1};
pop(1).Kwa = results.output.Kwa_run{1};
pop(1).Kd = results.output.Kd_run{1};
pop(1).Kc = results.output.Kc_run{1};
pop(1).S = results.output.S_run{1};
for it = 2:length(results.output.cost)
    %append new and old population
    tmpcost = [pop(it-1).cost results.output.cost{it}];

    tmpkwi = [pop(it-1).Kwi results.output.Kwi_run{it}];
    tmpKi = [pop(it-1).Ki results.output.Ki_run{it}];
    tmpKwa = [pop(it-1).Kwa results.output.Kwa_run{it}];
    tmpKd = [pop(it-1).Kd results.output.Kd_run{it}];
    tmpKc = [pop(it-1).Kc results.output.Kc_run{it}];
    tmpS = [pop(it-1).S results.output.S_run{it}];
    % Sort
    [~, SortOrder]=sort([tmpcost]);
    tmpcost = tmpcost(SortOrder);
    tmpkwi = tmpkwi(SortOrder);
    tmpKi = tmpKi(SortOrder);
    tmpKwa = tmpKwa(SortOrder);
    tmpKd = tmpKd(SortOrder);
    tmpKc = tmpKc(SortOrder);
    tmpS = tmpS(SortOrder);
    
    % Truncate
    pop(it).cost = tmpcost(1:100);
    pop(it).Kwi = tmpkwi(1:100);
    pop(it).Ki = tmpKi(1:100);
    pop(it).Kwa = tmpKwa(1:100);
    pop(it).Kd = tmpKd(1:100);
    pop(it).Kc = tmpKc(1:100);
    pop(it).S = tmpS(1:100);
end
figure
tf = tiledlayout(3,1);
title(tf,filename)

ax(1) = nexttile(1);
hold on; grid on
ylabel(ax(1),'Cost [1000$]')
xlim([1,251])
ylim([0 5E3])

ax(2) = nexttile(2);
hold on; grid on
ylabel('Storage [kWh]')
ylim([0,500])
xlim([1,251])

ax(3) = nexttile(3);
hold on; grid on
ylabel(ax(3),'Generation [kW]')
xlabel(ax(3),'Iteration')
ylim([0,8])
xlim([1,251])

for it = 1:length(results.output.cost)

    %check for existence of min cost
    tmpind = find(results.output.cost{it} == results.output.min.cost);
    if exist("indopt",'var')
        if ~isempty(tmpind)
            indopt = tmpind;
        end
    else
        indopt = tmpind;
    end
    if ~isempty(tmpind)
        itopt = it;
    end
    x = it.*ones(size(results.output.cost{it}));
    
    indinf = results.output.cost{it} > 4E6;
    indfi = ~indinf;

    % plot(ax(1),x,results.output.cost{it}./1000,'o','Color',[88/255, 90/255, 92/255],'linewidth',1.2,'markersize',4)
    % 
    % plot(ax(2),x(indfi),results.output.S_run{it}(indfi),'o','Color',coptpm(6,:),'linewidth',1.2,'markersize',4)
    % 
    % 
    % plot(ax(3),x(indfi),results.output.Kwi_run{it}(indfi),'o','Color',coptpm(1,:),'linewidth',1.2,'markersize',4)
    % plot(ax(3),x(indfi),results.output.Ki_run{it}(indfi),'o','Color',coptpm(2,:),'linewidth',1.2,'markersize',4)
    % plot(ax(3),x(indfi),results.output.Kwa_run{it}(indfi),'o','Color',coptpm(3,:),'linewidth',1.2,'markersize',4)
    % plot(ax(3),x(indfi),results.output.Kd_run{it}(indfi),'o','Color',coptpm(4,:),'linewidth',1.2,'markersize',4)
    % plot(ax(3),x(indfi),results.output.Kc_run{it}(indfi),'o','Color',coptpm(5,:),'linewidth',1.2,'markersize',4)
    plot(ax(1),x,pop(it).cost./1000,'o','Color',[88/255, 90/255, 92/255],'linewidth',1.2,'markersize',4)

    plot(ax(2),x,pop(it).S,'o','Color',coptpm(6,:),'linewidth',1.2,'markersize',4)
    plot(ax(3),x,pop(it).Kwi,'o','Color',coptpm(1,:),'linewidth',1.2,'markersize',4)
    plot(ax(3),x,pop(it).Ki,'o','Color',coptpm(2,:),'linewidth',1.2,'markersize',4)
    plot(ax(3),x,pop(it).Kwa,'o','Color',coptpm(3,:),'linewidth',1.2,'markersize',4)
    plot(ax(3),x,pop(it).Kd,'o','Color',coptpm(4,:),'linewidth',1.2,'markersize',4)
    plot(ax(3),x,pop(it).Kc,'o','Color',coptpm(5,:),'linewidth',1.2,'markersize',4)
    %drawnow;
        
    % Optional: Pause for visibility (remove or adjust as needed)
    %exportgraphics(gcf,strcat(filename,'.gif'),'Append',true);
end

plot(ax(1),itopt,results.output.cost{itopt}(indopt)./1000,'o','Color','r','linewidth',1.2,'markersize',7)
%exportgraphics(gcf,strcat(filename,'.gif'),'Append',true);