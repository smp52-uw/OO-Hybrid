%convert mooring output into a mat file with useful information
%Sarah May Palmer
%10-14-2024

clear
clc

%% load data
%moordat = readtable("Orcaflex mooring costs preliminary - Sheet1.csv"); %Initial mooring data from Curty (only for medium buoy size)
moordat = readtable("OrcaFlex Mooring model results - Regular Wave Jan 3rd 2025 results.csv"); %final mooring data from Curty (for SM/MD/LG buoys)

%locations
locs = {'WETS';'SFOMF';'Hueneme';'PISCES';...
    'PacWave';'Mid Atlantic';'Bering'};
saveL = {'WETS';'SFOMF';'PortHueneme';'PISCES';...
    'PacWave';'MidAtlSB';'BerSea'};
buoysz = {'Small';'Medium';'Large'};
buoydia = [2,4,6];

%% process and save data
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
            MoorMat.MaxW.dia(b,m) = buoydia(b);
            MoorMat.MaxW.mass(b,m) = moordat.TotalMass(indLBW(m));
            MoorMat.MaxW.PLmass(b,m) = moordat.Payload(indLBW(m));

            costtemp = moordat.totalMooringCost(indLBW(m));
            MoorMat.MaxW.cost(b,m) = str2double(costtemp{1}(2:end));
            MoorMat.MaxW.Size(b,m) = moordat.BuoySize(indLBW(m));
            MoorMat.MaxW.Sub(b,m) = moordat.MinClearance(indLBW(m));
        end
        for m = 1:length(indLBC)
            MoorMat.MaxC.dia(b,m) = buoydia(b);
            MoorMat.MaxC.mass(b,m) = moordat.TotalMass(indLBC(m));
            MoorMat.MaxC.PLmass(b,m) = moordat.Payload(indLBC(m));
            costtemp = moordat.totalMooringCost(indLBC(m));
            MoorMat.MaxC.cost(b,m) = str2double(costtemp{1}(2:end));
            MoorMat.MaxC.Size(b,m) = moordat.BuoySize(indLBC(m));
            MoorMat.MaxC.Sub(b,m) = moordat.MinClearance(indLBC(m));
        end

        %find the worst case mooring for a specific location and buoy
        if length(indLBW) == length(indLBC)
            for ii = 1:length(indLBC)
                jj = find(MoorMat.MaxW.mass(b,ii) == MoorMat.MaxC.mass(b,:)); %find current case with same buoy mass
                
                %%alternate method
                if MoorMat.MaxW.Sub(b,ii) < -0.1 || MoorMat.MaxC.Sub(b,jj) < -0.1 %fails the submergence requirement
                    MoorMat.WorstCase.dia(b,ii) = buoydia(b);
                    MoorMat.WorstCase.mass(b,ii) = MoorMat.MaxW.mass(b,ii);
                    MoorMat.WorstCase.PLmass(b,ii) = MoorMat.MaxW.PLmass(b,ii);

                    MoorMat.WorstCase.cost(b,ii) = nan; %cost is undefined for the failed moorings
                    MoorMat.WorstCase.Size(b,ii) = MoorMat.MaxW.Size(b,ii);
                    MoorMat.WorstCase.Sub(b,ii) = min([MoorMat.MaxW.Sub(b,ii),MoorMat.MaxC.Sub(b,jj)]);

                elseif strcmp(locs{l},'Hueneme') && buoydia(b) == 4 && MoorMat.MaxC.mass(b,jj) == 11000 %check for the self propelling buoy case
                    MoorMat.WorstCase.dia(b,ii) = buoydia(b);
                    MoorMat.WorstCase.mass(b,ii) = MoorMat.MaxW.mass(b,ii);
                    MoorMat.WorstCase.PLmass(b,ii) = MoorMat.MaxW.PLmass(b,ii);

                    MoorMat.WorstCase.cost(b,ii) = MoorMat.MaxW.cost(b,ii);
                    MoorMat.WorstCase.Size(b,ii) = MoorMat.MaxW.Size(b,ii);
                    MoorMat.WorstCase.Sub(b,ii) = MoorMat.MaxW.Sub(b,ii);

                elseif MoorMat.MaxW.cost(b,ii) > MoorMat.MaxC.cost(b,jj)
                    MoorMat.WorstCase.dia(b,ii) = buoydia(b);
                    MoorMat.WorstCase.mass(b,ii) = MoorMat.MaxW.mass(b,ii);
                    MoorMat.WorstCase.PLmass(b,ii) = MoorMat.MaxW.PLmass(b,ii);

                    MoorMat.WorstCase.cost(b,ii) = MoorMat.MaxW.cost(b,ii);
                    MoorMat.WorstCase.Size(b,ii) = MoorMat.MaxW.Size(b,ii);
                    MoorMat.WorstCase.Sub(b,ii) = MoorMat.MaxW.Sub(b,ii);
                else
                    MoorMat.WorstCase.dia(b,ii) = buoydia(b);
                    MoorMat.WorstCase.mass(b,ii) = MoorMat.MaxC.mass(b,jj);
                    MoorMat.WorstCase.PLmass(b,ii) = MoorMat.MaxC.PLmass(b,jj);

                    MoorMat.WorstCase.cost(b,ii) = MoorMat.MaxC.cost(b,jj);
                    MoorMat.WorstCase.Size(b,ii) = MoorMat.MaxC.Size(b,jj);
                    MoorMat.WorstCase.Sub(b,ii) = MoorMat.MaxC.Sub(b,jj);
                end
            end
        end
        
    end

    %Find lowest cost mooring cost for a given payload mass
    [mincost,minind] = min(MoorMat.WorstCase.cost, [], 1); %lowest cost payload (will ignore nan)
    loopsz = size(MoorMat.WorstCase.cost);
    for ii = 1:loopsz(2)
        MoorMat.SimMoor.dia(ii) = MoorMat.WorstCase.dia(minind(ii),ii);
        MoorMat.SimMoor.mass(ii) = MoorMat.WorstCase.mass(minind(ii),ii);
        MoorMat.SimMoor.PLmass(ii) = MoorMat.WorstCase.PLmass(minind(ii),ii);
        MoorMat.SimMoor.cost(ii) = MoorMat.WorstCase.cost(minind(ii),ii);
        MoorMat.SimMoor.Size(ii) = MoorMat.WorstCase.Size(minind(ii),ii);
        MoorMat.SimMoor.Sub(ii) = MoorMat.WorstCase.Sub(minind(ii),ii);
    end

    filenm = strcat(saveL{l},'_Mooring','.mat');
    save(filenm, 'MoorMat')
end
