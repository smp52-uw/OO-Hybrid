%Effect of diesel error

optInputs

%user inputs
kW_dies = linspace(0,8,100);
runtime_tot = 250*ones(size(kW_dies));
opt.p_dev.d_burn = calcDeviceVal('dieselburn',[],econ.diesburn_n);
lph = polyval(opt.p_dev.d_burn,kW_dies); %[l/h], burn rate

for i = 1:length(kW_dies)
    if lph(i)*runtime_tot(i) > 800
        runtime_tot(i) = 800/lph(i);
    end
    
    %fuel mass
    dies_dens = 0.85; %[g/cm3] = [kg/L] from a chevron report
    mass_fuel(i,1) = runtime_tot(i)*lph(i)*dies_dens; %[kg]
    mass_fuel(i,2) = 3*runtime_tot(i)*lph(i)*dies_dens; %[kg]
    %Fuel cost
    fuel(i,1) = runtime_tot(i)*lph(i)*econ.dies.fcost; %cost of consumed fuel
    fuel(i,2) = 3*runtime_tot(i)*lph(i)*econ.dies.fcost; %cost of consumed fuel
end

%norm by avg cost or mass of system
folder = "C:\Users\smpal\Documents\OO-Hybrid-Research\NAVFAC_ReportResults\";
%pull hybrid data
filename = '*_PD6_*.mat';
filename = dir(strcat(folder, filename));
for j = 1:max(size(filename))
    xn = split(filename(j).name,"_");
    temp = load(fullfile(folder,filename(j).name));
    fn = fieldnames(temp);

    if sum(strcmp(xn,"AllLC")) == 1
        optStruct6D{j} = temp.(fn{1});

    elseif length(size(temp.(fn{1}))) == 3
        tempopt = temp.(fn{1})(2,1,1);
    else
        tempopt = temp.(fn{1})(2,1,1,2);
    end

    if sum(strcmp(xn,"LC1")) == 1
        optStruct6D{j} = tempopt;
    elseif sum(strcmp(xn,"LC3")) == 1
        optStruct6D{j} = tempopt;
    elseif sum(strcmp(xn,"LC5")) == 1
        optStruct6D{j} = tempopt;
    end

    %pack up cost and mass
    if max(size(optStruct6D{j})) >1
        cost(j,1) = optStruct6D{j}(1,1).output.min.cost;
        cost(j,2) = optStruct6D{j}(1,2).output.min.cost;
        cost(j,3) = optStruct6D{j}(1,3).output.min.cost;
    
        mass(j,1) = optStruct6D{j}(1,1).output.min.Pmtrl./2./econ.platform.steel;
        mass(j,2) = optStruct6D{j}(1,2).output.min.Pmtrl./2./econ.platform.steel;
        mass(j,3) = optStruct6D{j}(1,2).output.min.Pmtrl./2./econ.platform.steel;
    else
        cost(j,1) = optStruct6D{j}(1,1).output.min.cost;
        cost(j,2) = nan;
        cost(j,3) = nan;
    
        mass(j,1) = optStruct6D{j}(1,1).output.min.Pmtrl./2./econ.platform.steel;
        mass(j,2) = nan;
        mass(j,3) = nan;
    end
end

opt.p_dev.d_mass = calcDeviceVal('dieselmass',[],econ.diesmass_n);
mass_dies = polyval(opt.p_dev.d_mass,kW_dies); %mass of 1 generator [kg]

mcost = mean(cost,'all', 'omitnan');
mmass = mean(mass,'all', 'omitnan');


%%
figure
tiledlayout(2,1)

nexttile
plot(kW_dies,mass_fuel(:,2)./mass_dies','LineWidth',2,'DisplayName','Wrong')
hold on
plot(kW_dies,mass_fuel(:,1)./mass_dies','LineWidth',2,'DisplayName','Correct')
title('Fuel Mass / Generator Mass')
grid on

nexttile
plot(kW_dies,mass_fuel(:,2)./mmass,'LineWidth',2,'DisplayName','Wrong')
hold on
plot(kW_dies,mass_fuel(:,1)./mmass,'LineWidth',2,'DisplayName','Correct')
title('Fuel Mass / Avg Platform Mass')
grid on
xlabel("Diesel Generator Size [kW]")
legend
