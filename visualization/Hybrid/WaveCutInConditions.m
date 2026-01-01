%Wave parameters for cut in conditions

g = 9.81;
rho = 1025;

%Hs range
Hs = 0.02:0.01:1; %2 cm jumps to 1 m

%Tp range
Tp = 0.1:0.1:15;

[hh,tt] = meshgrid(Hs,Tp);

k = (1./(64.*pi)).*rho.*(g.^2).*(hh.^2).*tt; %Hs2*Tp;

%limit based on 1 kw
k1 = k;
i1 = find(k1 > 30);
i2 = find(k1 < 25);
k1(i1) = nan;
k1(i2) = nan;

%limit based on 2kw
k2 = k;
i1 = find(k2 > 50);
i2 = find(k2 < 42);
k2(i1) = nan;
k2(i2) = nan;

figure
pcolor(hh,tt,k1);
hold on
pcolor(hh,tt,k2);
ylabel('Period [s]')
xlabel('Wave Height [m]')
c = colorbar;
c.Label.String = 'Wave Power Flux [W/m]';
c.Label.FontSize = 11;
