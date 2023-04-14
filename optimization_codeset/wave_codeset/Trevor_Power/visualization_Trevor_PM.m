clear all
wsr_2 = load('2m_5kW_WEC_Power_Frame.mat');
wsr_2 = wsr_2.('WEC_pframe');
%names = get(wsr_2,'columnname');
names = fieldnames(wsr_2);
names = names(3:end-3);
names = erase(names, 'I_w_a_v_e -');
wsr_2 = table2array(wsr_2);
wsr_2 = wsr_2(1:end-1,:);

% P_500 = wsr_2(:,12)*48;
% duplicates = sum(~unique(P_500));

[H,T] = meshgrid(0.25:0.5:8.75,1:22);
grid_sz = size(H);
H = reshape(H,[grid_sz(1)*grid_sz(2),1]);
T = reshape(T,[grid_sz(1)*grid_sz(2),1]);

%plot(wsr_2(:,1),wsr_2(:,2))
tbl_size = size(wsr_2);
plotsz = ceil((tbl_size(2)-2)/3);
for i = 3:tbl_size(2)
%for i = tbl_size-12:tbl_size
    I = griddata(wsr_2(:,1),wsr_2(:,2),wsr_2(:,i),H,T);
    I(~ismember(H,wsr_2(:,1))) = NaN;
    I_grid = reshape(I,[grid_sz(1),grid_sz(2)]);
    P_grid = -48.*I_grid; %power
    [H_grid,T_grid] = meshgrid(0.25:0.5:8.75,1:22);
    subplot(3,plotsz,(i-2))
    surf(T_grid,H_grid,P_grid)
    xlabel('Wave Height')
    ylabel('Period')
    title(names(i-2))
    view(2)
    colorbar
    %caxis([0,2500])
end