% S2 = load('Test_P6253_S14.mat');
% S1 = load('Test_P6253_S13.mat');
% S14H = load('Test_P6253_S14_HP.mat');
% S13H = load('Test_P6253_S13_HP.mat');

% T1 = load('Test_P7856_S18_HPF.mat');
% T1.Name = "P7856_S18";
% T2 = load('Test_P7856_S19_HPF.mat');
% T2.Name = "P7856_S19";

T1 = load('Test_P7856_S18_HPF_TA6.mat');
T1.Name = "P7856_S18";
T2 = load('Test_P7856_S19_HPF_TA6.mat');
T2.Name = "P7856_S19";

figure
tiledlayout(5,1)

ax(1) = nexttile;
hold on
plot(T1.Pcurr,'.','linewidth',2,'DisplayName',T1.Name)
plot(T2.Pcurr,'.','linewidth',2,'DisplayName',T2.Name)
ylabel('Pc')

ax(2) = nexttile;
hold on
plot(T1.S1,'.','linewidth',2,'DisplayName',T1.Name)
plot(T2.S1,'.','linewidth',2,'DisplayName',T2.Name)
ylabel('B1')
legend

ax(3) = nexttile;
hold on
plot(T1.S2,'.','linewidth',2,'DisplayName',T1.Name);
plot(T2.S2,'.','linewidth',2,'DisplayName',T2.Name);
ylabel('B2')

ax(4) = nexttile;
hold on
plot(T1.L,'linewidth',2,'DisplayName',T1.Name);
plot(T2.L,'linewidth',2,'DisplayName',T2.Name);
ylabel('L')

%Fails
ax(5) = nexttile;
hold on
plot(T1.F,'.','linewidth',2,'DisplayName',T1.Name);
plot(T2.F,'.','linewidth',2,'DisplayName',T2.Name);
hr = 1:1:length(T2.F);
plot(hr(T2.F>T1.F), T2.F(T2.F>T1.F),'x','linewidth',2,'DisplayName','T2 Extra Fail');
plot(hr(T2.F<T1.F), T1.F(T2.F<T1.F),'*','linewidth',2,'DisplayName','T1 Extra Fail');
ylabel('F')
legend

linkaxes(ax,'x')

T1surv = 1 - sum(T1.F)/length(T1.F)
T2surv = 1 - sum(T2.F)/length(T2.F)