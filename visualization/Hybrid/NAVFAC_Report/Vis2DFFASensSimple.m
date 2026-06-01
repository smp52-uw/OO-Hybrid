function Vis2DFFASensSimple(optStruct)
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
            for l = 1:strsz(4)
                for m = 1:strsz(5)
                    if m == 1
                        col = cola;
                    elseif m == 2
                        col = colb;
                    else
                        col = colc;
                    end
                    for n = 1:strsz(6)
                        count = count + 1;
                        if n == 1
                            scatter(count, optStruct(i,j,k,l,m,n).output.min.cost/1000,'o','MarkerEdgeColor',col(4,:),'linewidth',2)
                        else
                            scatter(count, optStruct(i,j,k,l,m,n).output.min.cost/1000,'x','MarkerEdgeColor',col(4,:),'linewidth',2)
                        end
                    end
                end
            end
        end
    end
end
titlestr = strcat("PM: ",string(optStruct(1,1,1,1,1,1).opt.pm), " LC: ",string(optStruct(1,1,1,1,1,1).uc.loadcase))
title(titlestr)
ylabel("Best Cost [$ 1000]")
xlabel("Sensitivity Run")
end