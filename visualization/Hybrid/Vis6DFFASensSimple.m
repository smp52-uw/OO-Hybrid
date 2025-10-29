function Vis6DFFASensSimple(optStruct)
%Simple visualization of FFA sensitivity results
%just min cost vs point

cola = brewermap(9,'Blues'); %colors
colb = brewermap(9,'Greens'); %colors
colc = brewermap(9,'Purples'); %colors
cold = brewermap(9,'Oranges'); %colors
cole = brewermap(9,'Reds'); %colors
colf = brewermap(9,'RdPu'); %colors

strsz = size(optStruct);

figure
hold on

count = 0;
for i = 1:strsz(1)
    for j = 1:strsz(2)
        for k = 1:strsz(3)
            if length(strsz) == 3
                titlestr = strcat("PM: ",string(optStruct(1,1,1).opt.pm), " LC: ",string(optStruct(1,1,1).uc.loadcase));
                count = count + 1;
                ffavar(count,:) = [optStruct(i,j,k).opt.ffa.sens{1}(i) optStruct(i,j,k).opt.ffa.sens{2}(j) optStruct(i,j,k).opt.ffa.sens{3}(k)...
                    optStruct(i,j,k).opt.ffa.sens{4}(1) optStruct(i,j,k).opt.ffa.sens{5}(1) optStruct(i,j,k).opt.ffa.sens{6}(1)];
                scatter(count, optStruct(i,j,k).output.min.cost/1000,'o','MarkerEdgeColor',cola(4,:),'linewidth',2)
            else
                titlestr = strcat("PM: ",string(optStruct(1,1,1,1).opt.pm), " LC: ",string(optStruct(1,1,1,1).uc.loadcase));
                %ind = find(strsz > 1,1,'last');
                for l = 1:strsz(5)
                    for m = 1:strsz(6)
                        count = count + 1;
                        %scatter(count, optStruct(i,j,k,1,l,m).output.min.cost/1000,'o','MarkerEdgeColor',cola(4,:),'linewidth',2)
    
                        scatter(count, optStruct(i,j,k,1,l,m).output.min.surv,'o','MarkerEdgeColor',cola(4,:),'linewidth',2)
                    end
                end
            end
        end
    end
end

title(titlestr)
ylabel("Best Cost [$ 1000]")
xlabel("Sensitivity Run")
end