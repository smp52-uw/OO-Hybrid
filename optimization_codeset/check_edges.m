function [j_1,j_max] = check_edges(j, j_min)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if j_min == 1
        j_1 = 1;
        j_max = j_min + 1;
    elseif j_min == j
        j_1 = j_min-1;
        j_max = j;
    else
        j_1 = j_min -1;
        j_max = j_min + 1;
    end
end