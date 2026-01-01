function [j_1,j_max] = check_edges2(j, j_min)
%2 box telescoping
    if j_min <= 2 %left side near grid minimum
        j_1 = 1;
        j_max = j_min + 2;
    elseif j_min >= j-1 %right side near grid maximum
        j_1 = j_min-2;
        j_max = j;
    else %in the middle
        j_1 = j_min -2;
        j_max = j_min + 2;
    end
end