%Identify the size vs rated power curves for each generation type to
%identify the minimum rated power that we should use in the design space

optInputs %get all the power constants

kW = linspace(0,8,1000); %dummy rated power vector

%wind
Dwi = 2*sqrt(1000*2*kW/(turb.eta*atmo.rho_a*pi*turb.ura^3)); %[m]

%current
Dc = sqrt(1000*2*kW/(cturb.eta*atmo.rho_w*cturb.ura^3)); %[m]

%solar
A = kW/(inso.eff*inso.rated);
Di = sqrt(A);

%Diesel
opt.p_dev.d_size = calcDeviceVal('dieselsize',[],econ.diessize_n);
Dd = polyval(opt.p_dev.d_size,kW);

figure
tf = tiledlayout(4,1);
title(tf,'Generator Diameters')

nexttile
plot(kW,Dwi,'linewidth',2)
ylabel('Wind Dia')
xlim([0,1])
grid on

nexttile
plot(kW,Dc,'linewidth',2)
ylabel('Current Dia')
xlim([0,1])
grid on

nexttile
plot(kW,Di,'linewidth',2)
ylabel('Inso Dia')
xlim([0,1])
grid on

nexttile
plot(kW,Dd,'linewidth',2)
ylabel('Diesel Dia')
xlim([0,1])
grid on
xlabel('Rated Power [kW]')