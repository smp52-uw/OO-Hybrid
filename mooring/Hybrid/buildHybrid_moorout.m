%convert mooring output into a mat file with useful information
%Sarah May Palmer
%10-14-2024
clear
clc
%load data
moordat = readtable("Orcaflex mooring costs preliminary - Sheet1.csv"); %Initial mooring data from Curty (only for medium buoy size)

%locations
locs = {'WETS';'SFOMF';'Hueneme';'PISCES';...
    'PacWave';'Mid Atlantic';'Bering'};
saveL = {'WETS';'SFOMF';'PortHueneme';'PISCES';...
    'PacWave';'MidAtlSB';'BerSea'};
%buoysz = {'Small';'Medium';'Large'};
buoysz = {'Medium'}; %for now we only have medium data

figure
tiledlayout(2,1)
for l = 1:max(size(locs))
    for b = 1:max(size(buoysz))
        isL = strcmp(moordat.Site,locs{l}); %find rows of this location
        %indL = find(isL);
        isB = strcmp(moordat.BuoySize,buoysz{b}); %find rows of this size
        %indB = find(isB);
        isW = strcmp(moordat.Waves,'Max'); %rows of max wave conditions
        %indW = find(isW);
        isC = strcmp(moordat.Current,'Max'); %rows of max current conditions
        %indC = find(isC);

        indLBW = find(isL & isB & isW); %both location and buoy size and max waves
        indLBC = find(isL & isB & isC); %both location and buoy size and max current

        for m = 1:length(indLBW)
            MoorMat.MaxW.mass(b,m) = moordat.TotalMass(indLBW(m));
            MoorMat.MaxW.cost(b,m) = moordat.totalMooringCost(indLBW(m));
            MoorMat.MaxW.Size(b,m) = moordat.BuoySize(indLBW(m));
        end
        for m = 1:length(indLBC)
            MoorMat.MaxC.mass(b,m) = moordat.TotalMass(indLBC(m));
            MoorMat.MaxC.cost(b,m) = moordat.totalMooringCost(indLBC(m));
            MoorMat.MaxC.Size(b,m) = moordat.BuoySize(indLBC(m));
        end

        nexttile(1)
        hold on
        plot(MoorMat.MaxC.mass,MoorMat.MaxC.cost,'DisplayName',locs{l},'LineWidth',2)
        title('Max Current/Average Wave')
        ylabel('Cost $')
        xlabel('Mass [kg]')
        nexttile(2)
        hold on
        plot(MoorMat.MaxW.mass,MoorMat.MaxW.cost,'DisplayName',locs{l},'LineWidth',2)
        title('Max Wave/Average Current')
        xlabel('Mass [kg]')
        ylabel('Cost $')


    end
    filenm = strcat(saveL{l},'_Mooring','.mat');
    save(filenm, 'MoorMat')
end

legend('location','eastoutside')

