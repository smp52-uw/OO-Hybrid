function Vis2DFFASens(optStruct,tabOutput, plottype)
%visualize the output from the hyperparameter sensitivity study on the 2D
%firefly algorithm

%Input: the output structure from the sensitivity study
%Output: plots of the cost vs iteration for each run. The color coding for
%each plot will be based on a different variable to determine which
%parameters are most important

%tabOutput = readtable("MidAtlSB_PD2PM1_LC5_FFA_FSALG_06302025.csv");
if ~isempty(tabOutput)
    it1 = find(strcmp(tabOutput.Column1,'Iteration 1'));
    for it = 1:length(it1)
        if strcmp(tabOutput.Column1(it1(it)+50), 'Iteration 51') %100 iteration run
            pop{it}.Iteration = 1:1:100;
            tempcost = tabOutput.Column2(it1(it):it1(it)+99);
            splitcost = split(tempcost,' ');
            pop{it}.BestCost = str2double(splitcost(:,end));
        else
            pop{it}.Iteration = 1:1:50;
            tempcost = tabOutput.Column2(it1(it):it1(it)+49);
            splitcost = split(tempcost,' ');
            pop{it}.BestCost = str2double(splitcost(:,end));
        end
    end
end

strsz = size(optStruct);

%colors
cola = brewermap(9,'Blues'); %colors
colb = brewermap(9,'Greens'); %colors
colc = brewermap(9,'Purples'); %colors
cold = brewermap(9,'Oranges'); %colors
cole = brewermap(9,'Reds'); %colors
colf = brewermap(9,'RdPu'); %colors

col1{1} = cola([8,6,3],:);
col1{2} = colb([8,6,3],:);
col1{3} = colc([8,6,3],:);
col2{1} = cold([8,6,3],:);
col2{2} = cole([8,6,3],:);
col2{3} = colf([6,4,3],:);

% % %all results figures
count = 1; %counter
for p = 1:strsz(1)
    figure %one figure for each number of iterations
    set(gcf,'Position', [100, 100, 900, 500])
    tf = tiledlayout(strsz(2),strsz(3));
    title(tf,strcat(" Iterations: ",string(optStruct(p,1,1,1,1,1).opt.ffa.max)))
    for q = 1:strsz(2)
        for r = 1:strsz(3)
            ax(q,r) = nexttile; %tile for pop(q), gamma(r)
            hold on
            grid on
            xlabel('Iteration')
            ylabel('Cost [$]')


            for s = 1:strsz(4)
                for u = 1:strsz(5)
                    if s == 1
                        col = col1{u};
                    else
                        col = col2{u};
                    end

                    for v = 1:strsz(6) 
                        itStruct = optStruct(p,q,r,s,u,v);

                        if strcmp(plottype,'newpop')
                            for it = 1:itStruct.opt.ffa.max + 1
                                mincost(it) = min(itStruct.output.cost{it});
                            end
                            itvec = 1:1:itStruct.opt.ffa.max + 1;
                        elseif ~isempty(tabOutput)
                            mincost = pop{count}.BestCost';
                            itvec = pop{count}.Iteration;
                        else
                            mincost = itStruct.output.popBestCost;
                            itvec = 1:1:length(mincost);
                        end

                        %pack cost, runtime, and convergence into arrays
                        absmincost(count) = itStruct.output.min.cost; %min cost for this case
                        runtime(count) = str2double(itStruct.comptime);
                        [f, ~] = fit(itvec', mincost', 'exp1');
                        expdecay_a(count) = f.a;
                        expdecay_b(count) = f.b;
                        count = count + 1;

                        plot(itvec,mincost,'linewidth',1.5,'Color',col(v,:))
                        titlestr = strcat("Population: ",string(itStruct.opt.ffa.pop)," Gamma: ",string(itStruct.opt.ffa.gamma));
                        title(titlestr)

                    end
                end
            end
        end
    end
    linkaxes(ax,'y')
end

% %removing bad points - alphad colormap
% count = 1; %counter
% for p = 1:strsz(1)
%     figure %one figure for each number of iterations
%     set(gcf,'Position', [100, 100, 900, 500])
%     tf = tiledlayout(strsz(2),strsz(3));
%     title(tf,strcat(" Iterations: ",string(optStruct(p,1,1,1,1,1).opt.ffa.max)))
%     for q = 1:strsz(2)
%         for r = 1:strsz(3)
%             ax(q,r) = nexttile; %tile for pop(q), gamma(r)
%             hold on
%             grid on
%             xlabel('Iteration')
%             ylabel('Cost [$]')
% 
% 
%             for s = 1:strsz(4)
%                 for u = 1:strsz(5)
%                     if s == 1
%                         col = col1{u};
%                     else
%                         col = col2{u};
%                     end
% 
%                     for v = 1:strsz(6) 
%                         itStruct = optStruct(p,q,r,s,u,v);
% 
%                         if strcmp(plottype,'newpop')
%                             for it = 1:itStruct.opt.ffa.max + 1
%                                 mincost(it) = min(itStruct.output.cost{it});
%                             end
%                             itvec = 1:1:itStruct.opt.ffa.max + 1;
%                         elseif ~isempty(tabOutput)
%                             mincost = pop{count}.BestCost';
%                             itvec = pop{count}.Iteration;
%                         else
%                             mincost = itStruct.output.popBestCost;
%                             itvec = 1:1:length(mincost);
%                         end
% 
%                         %pack cost, runtime, and convergence into arrays
%                         absmincost(count) = itStruct.output.min.cost; %min cost for this case
%                         runtime(count) = str2double(itStruct.comptime);
%                         [f, ~] = fit(itvec', mincost', 'exp1');
%                         expdecay_a(count) = f.a;
%                         expdecay_b(count) = f.b;
%                         count = count + 1;
% 
%                         if v ~=3 %no light
%                             plot(itvec,mincost,'linewidth',1.5,'Color',col(v,:))
%                             titlestr = strcat("Population: ",string(itStruct.opt.ffa.pop)," Gamma: ",string(itStruct.opt.ffa.gamma));
%                             title(titlestr)
%                         end
% 
%                     end
%                 end
%             end
%         end
%     end
%     linkaxes(ax,'y')
% end


%Removing bad points - Gamma Color Map
cold = brewermap(9,'YlOrRd'); %colors
col2{1} = cold([4,3,2],:);
count = 1; %counter
figure %one figure for each number of iterations
set(gcf,'Position', [100, 100, 900, 500])
tf = tiledlayout(strsz(1),strsz(2));
for p = 1:strsz(1)
    for q = 1:strsz(2) %pop
        ax(p,q) = nexttile;
        hold on
        grid on
        xlabel('Iteration')
        ylabel('Cost [$]')
        for r = 1:strsz(3) %gamma
            for s = 1:strsz(4) 
                for u = 1:strsz(5)
                    if s == 1
                        col = col1{u};
                    else
                        col = col2{u};
                    end
                    for v = 1:strsz(6) 
                        itStruct = optStruct(p,q,r,s,u,v);

                        if strcmp(plottype,'newpop')
                            for it = 1:itStruct.opt.ffa.max + 1
                                mincost(it) = min(itStruct.output.cost{it});
                            end
                            itvec = 1:1:itStruct.opt.ffa.max + 1;
                        elseif ~isempty(tabOutput)
                            mincost = pop{count}.BestCost';
                            itvec = pop{count}.Iteration;
                        else
                            mincost = itStruct.output.popBestCost;
                            itvec = 1:1:length(mincost);
                        end

                        %pack cost, runtime, and convergence into arrays
                        absmincost(count) = itStruct.output.min.cost; %min cost for this case
                        runtime(count) = str2double(itStruct.comptime);
                        [f, ~] = fit(itvec', mincost', 'exp1');
                        expdecay_a(count) = f.a;
                        expdecay_b(count) = f.b;
                        count = count + 1;

                        lintyp = ['-',':'];
                        if  v ~= 3 %dark and medium

                            plot(itvec,mincost,'linewidth',1.5,'Color',col(r,:),'LineStyle',lintyp(v))
                            titlestr = strcat(" Iterations: ",string(optStruct(p,1,1,1,1,1).opt.ffa.max)," Population: ",string(itStruct.opt.ffa.pop));
                            title(titlestr)
                        end

                    end
                end
            end
        end
    end
    linkaxes(ax,'y')
end



min(absmincost)
min(runtime)/60
max(runtime)/60
min(expdecay_b)

end