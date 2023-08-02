%Visualize Firefly Optimization

%Choose 3D
d1 = 1;
d2 = 2;
d3 = 3;

%Convert populations into arrays
for j=1:MaxIt+1
    hold on
    for i=1:nPop
        if ~isempty(pop_save{j}(i).Position)
            p1(i) = pop_save{j}(i).Position(d1);
            p2(i) = pop_save{j}(i).Position(d2);
            p3(i) = pop_save{j}(i).Position(d3);
        end
    end
    scatter3(p1, p2, p3)
    drawnow
    pause(2)
end
