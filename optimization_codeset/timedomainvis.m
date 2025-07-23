S14 = load('Test_P6253_S14.mat');
S13 = load('Test_P6253_S13.mat');
S14H = load('Test_P6253_S14_HP.mat');
S13H = load('Test_P6253_S13_HP.mat');

% NLALL = load('Test_6D6253_S13_AllNLBattDeg.mat');
% CalB2 = load('Test_6D6253_S13_CalBattDeg.mat');

figure
tiledlayout(3,2)

ax(1) = nexttile;
hold on
plot(S13.S1,'.','linewidth',2,'DisplayName','13')
plot(S14.S1,'.','linewidth',2,'DisplayName','14')
ylabel('B1')
title('No Hotel Rule')
legend

ax(2) = nexttile;
hold on
plot(S13H.S1,'.','linewidth',2,'DisplayName','13')
plot(S14H.S1,'.','linewidth',2,'DisplayName','14')
title('Hotel Rule')
ylabel('B1')

ax(3) = nexttile;
hold on
plot(S13.S2,'.','linewidth',2,'DisplayName','13');
plot(S14.S2,'.','linewidth',2,'DisplayName','14');
ylabel('B2')

ax(4) = nexttile;
hold on
plot(S13H.S2,'.','linewidth',2,'DisplayName','13')
plot(S14H.S2,'.','linewidth',2,'DisplayName','14')
ylabel('B2')

%Fails
ax(5) = nexttile;
hold on
plot(S13.F,'.','linewidth',2,'DisplayName','13');
plot(S14.F,'.','linewidth',2,'DisplayName','14');
hr = 1:1:length(S14.F);
plot(hr(S14.F>S13.F), S14.F(S14.F>S13.F),'x','linewidth',2,'DisplayName','Nonlinear');
ylabel('F')

ax(6) = nexttile;
hold on
plot(S13H.F,'.','linewidth',2,'DisplayName','13')
plot(S14H.F,'.','linewidth',2,'DisplayName','14')
hr = 1:1:length(S14H.F);
plot(hr(S14H.F>S13H.F), S14.F(S14H.F>S13H.F),'x','linewidth',2,'DisplayName','Nonlinear');
ylabel('F')




linkaxes(ax,'x')