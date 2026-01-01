function [Kd, Ki, Kwi, Kwa, Kc, S,sim_run] = loopyloops_v2(S_t, S_old, Kc_t, Kc_old, Kwa_t, Kwa_old, Kwi_t, Kwi_old, Ki_t, Ki_old, Kd_t, Kd_old,opt)

%set loop constants
% i = 1; %index of old grid style array
% dx = round(mode(diff(Kd_t)),4); %for this to work the dx must be constant
% ds = round(mode(diff(S_t)),4);
dx = mode(diff(Kd_t)); %for this to work the dx must be constant
ds = mode(diff(S_t));
% Kd = [];
% Ki = [];
% Kwi = [];
% Kwa = [];
% Kc = [];
% S = [];

%matrix of old points before adjustments
%old_dim = size(Kd_old)
% if old_dim(1) > 1
%     old_points = [Kd_old Ki_old Kwi_old Kwa_old Kc_old S_old];
% else
%     old_points = [Kd_old' Ki_old' Kwi_old' Kwa_old' Kc_old' S_old'];
% end
%grid boundaries
sgrid_max = opt.bf.N;
xgrid_max = opt.bf.M;
%make permutation matrix
% syms a b c
% p_mat = [a b c];
per_mat = [1 2 3];
per_mat = permn(per_mat,6)';
p_mat = per_mat;
disp('Loopy Loops start ...')

%everything will break if D = 0;
Kd_old(Kd_old == 0) = 0.1;
Ki_old(Ki_old == 0) = 0.1;
Kwi_old(Kwi_old == 0) = 0.1;
Kwa_old(Kwa_old == 0) = 0.1;
Kc_old(Kc_old == 0) = 0.1;
S_old(S_old == 0) = 0.1;

%Big_Matrix = zeros[5,(3^5)*length(Kd_old)];
Big_Matrix = [];
parfor (i = 1 : length(Kd_old),opt.bf.maxworkers)
    %i = 1;
    d_mat = p_mat(1,:);
    d_mat(d_mat == 2) = 1 + dx/Kd_old(i);
    d_mat(d_mat == 3) = 1 - dx/Kd_old(i);
    
    i_mat = p_mat(2,:);
    i_mat(i_mat == 2) = 1+ dx/Ki_old(i);
    i_mat(i_mat == 3) = 1+ -1*dx/Ki_old(i);
    
    wi_mat = p_mat(3,:);
    wi_mat(wi_mat== 2) = 1+ dx/Kwi_old(i);
    wi_mat(wi_mat == 3) = 1+ -1*dx/Kwi_old(i);
    
    wa_mat = p_mat(4,:);
    wa_mat(wa_mat == 2) = 1+ dx/Kwa_old(i);
    wa_mat(wa_mat == 3) = 1+ -1*dx/Kwa_old(i);

    c_mat = p_mat(5,:);
    c_mat(c_mat == 2) = 1+ dx/Kc_old(i);
    c_mat(c_mat == 3) = 1+ -1*dx/Kc_old(i);
    
    s_mat = p_mat(6,:);
    s_mat(s_mat == 2) = 1+ ds/S_old(i);
    s_mat(s_mat == 3) = 1+ -1*ds/S_old(i);
    
    step_mat = [d_mat;i_mat;wi_mat;wa_mat;c_mat;s_mat];
    
    D = diag([Kd_old(i) Ki_old(i) Kwi_old(i) Kwa_old(i) Kc_old(i) S_old(i)]);
    D_2 = D*step_mat;
    %adjust values that should have been zero
    D_2(D_2 == 0.1-dx) = -1*dx;
    D_2(D_2 == 0.1+dx) = 1*dx;
    D_2(D_2 == 0.1) = 0;
    D_2(D_2 == 0.1-ds) = -1*ds;
    D_2(D_2 == 0.1+ds) = 1*ds;
    %no negative value
    %remove any cols where one points is outside the grid boundaries

    edge_i = D_2(1,:) > xgrid_max;
    edge_i = edge_i + (D_2(2,:) > xgrid_max);
    edge_i = edge_i + (D_2(3,:) > xgrid_max);
    edge_i = edge_i + (D_2(4,:) > xgrid_max);
    edge_i = edge_i + (D_2(5,:) > xgrid_max);
    edge_i = edge_i + (D_2(6,:) > sgrid_max);
    edge_i = edge_i + (D_2(6,:) < 1);
    edge_i = edge_i + (sum(D_2 < 0,1));
    rm_i = edge_i >= 1;
    D_2(:,rm_i) = [];
    Big_Matrix = [Big_Matrix D_2];
end

Big_Matrix = round(Big_Matrix,4);
new_points = unique(Big_Matrix','rows','stable');
%new_points = new_points;
%old_points = round(old_points,4);
%[extra_points,in,io] = intersect(new_points, old_points,'stable','rows');
%new_points(in,:) = [];
new_points = new_points';
Kd = new_points(1,:);
Ki = new_points(2,:);
Kwi = new_points(3,:);
Kwa = new_points(4,:);
Kc = new_points(5,:);
S = new_points(6,:);
%sim_run is 1 for all points added to the new grid like arrays
sim_run = ones(length(S),1);
%sim_run(in) = 0; %don't re-run the old points - I want to rerun them for
%the per/tel mix

disp('Loopy Loops end ...')
end