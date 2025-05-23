%Fischer Panda AGT DC 4000 12 PMS (email and website)
DiesPPI = 1.2; %Adjusting cost from 2020 ->2024 via PPI
diesLib(1).kW = 3.2;
diesLib(1).cost = 15000*DiesPPI;
diesLib(1).d = 0.56; %m
diesLib(1).V = diesLib(1).d^3; %[m3]
diesLib(1).m = 135.60; %kg
diesLib(1).c = .757; %[l/h]

%Polar Power Volvo Penta D1-13 - 8080
diesLib(2).kW = 5.5;
diesLib(2).cost = 14260*DiesPPI;
diesLib(2).d = 0.53; %m
diesLib(2).V = 0.491*0.53*0.532; %[m^3]
diesLib(2).m = 111.58; %kg
diesLib(2).c = .32*diesLib(2).kW; %[l/h]

%Polar Power Volvo Penta D1-20 - 8080
diesLib(3).kW = 8;
diesLib(3).cost = 15595*DiesPPI;
diesLib(3).d = 0.60; %m
diesLib(3).V = 0.491*0.590*0.552; %[m^3]
diesLib(3).m = 131.54; %kg
diesLib(3).c = .29*diesLib(3).kW; %[l/h]

%Polar Power Volvo Penta D1-20 - 8220
diesLib(4).kW = nan;
diesLib(4).cost = nan;
diesLib(4).d = nan; %m
diesLib(4).V = nan; %m
diesLib(4).m = nan; %kg
diesLib(4).c = nan; %[l/h]

%Polar Power Volvo Penta D1-30 - 8220
diesLib(5).kW = 14;
diesLib(5).cost = 18160*DiesPPI;
diesLib(5).d = 0.66; %m
diesLib(5).V = 0.502*0.657*0.582; %[m3]
diesLib(5).m = 152.86; %kg
diesLib(5).c = .29*diesLib(5).kW; %[l/h]

%Polar Power Volvo Penta D1-40 - 8220
diesLib(6).kW = nan;
diesLib(6).cost = nan;
diesLib(6).d = nan; %m
diesLib(6).V = nan; %m
diesLib(6).m = nan; %kg
diesLib(6).c = nan; %[l/h]

%Polar Power Volvo Penta D1-40 - 8340
diesLib(7).kW = 20;
diesLib(7).cost = 22165*DiesPPI;
%diesLib(7).d = 0.53; %m (wrong?)
diesLib(7).d = 0.74; %m (PDC 8340VP-40 D2-40)
diesLib(7).V = 0.502*0.740*0.603; %[m3]
diesLib(7).m = 179.62; %kg
diesLib(7).c = .29*diesLib(7).kW; %[l/h]

%Hamilton Ferris FCD 150ME
diesLib(8).kW = 1.8;
diesLib(8).cost = 7200*DiesPPI;
diesLib(8).d = 0.76; %m
diesLib(8).V = 0.76*0.46*0.42; %[m3]
diesLib(8).m = 56.70; %kg
diesLib(8).c = 0.568; %[l/h]

%ALTEN D1
diesLib(9).kW = 1.4;
diesLib(9).cost = 8048*DiesPPI;
diesLib(9).d = 0.6096; %m
diesLib(9).V = diesLib(9).d^3; %[m3]
diesLib(9).m = 56.70; %kg
diesLib(9).c = 0.95; %[l/h]

%ALTEN D2
diesLib(10).kW = 1.9;
diesLib(10).cost = 10000*DiesPPI;
diesLib(10).d = nan; %m
diesLib(10).V = nan;
diesLib(10).m = nan; %kg
diesLib(10).c = 1.36; %[l/h]

%ALTEN D3
diesLib(11).kW = 3.6;
diesLib(11).cost = 12000*DiesPPI;
diesLib(11).d = nan; %m
diesLib(11).V = nan;
diesLib(11).m = nan; %kg
diesLib(11).c = 2.00; %[l/h]

