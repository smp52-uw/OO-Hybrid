tf = tiledlayout(3,1);
title(tf,'PacWave: Wi = 3, D = 1, S = 70')

ax(1) = nexttile;
plot(optStruct.output.min.L)
ylabel('L(t)')

ax(2) = nexttile;
hold on
plot(optStruct.output.min.Pdies)
plot(optStruct.output.min.Pwind)
ylabel('Power')
legend('Dies','Wind')

ax(3) = nexttile;
hold on
plot(optStruct.output.min.S1)
plot(optStruct.output.min.S2)
ylabel('S(t)')

linkaxes(ax,'x')