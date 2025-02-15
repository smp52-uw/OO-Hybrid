%convert mooring output into a mat file with useful information
%Sarah May Palmer
%10-14-2024
clear
clc
%load data
%moordat = readtable("Orcaflex mooring costs preliminary - Sheet1.csv"); %Initial mooring data from Curty (only for medium buoy size)
moordat = readtable("OrcaFlex Mooring model results - Regular Wave Jan 3rd 2025 results.csv"); %Initial mooring data from Curty (only for medium buoy size)

%locations
locs = {'WETS';'SFOMF';'Hueneme';'PISCES';...
    'PacWave';'Mid Atlantic';'Bering'};
saveL = {'WETS';'SFOMF';'PortHueneme';'PISCES';...
    'PacWave';'MidAtlSB';'BerSea'};
buoysz = {'Small';'Medium';'Large'};
buoydia = [2,4,6];
%buoysz = {'Medium'}; %for now we only have medium data

% figure
% tiledlayout(3,1)
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

                % %%worst performing mooring method
                % if MoorMat.MaxW.Sub(b,ii) < MoorMat.MaxC.Sub(b,jj)
                %     MoorMat.WorstCase.dia(b,ii) = buoydia(b);
                %     MoorMat.WorstCase.mass(b,ii) = MoorMat.MaxW.mass(b,ii);
                %     MoorMat.WorstCase.cost(b,ii) = MoorMat.MaxW.cost(b,ii);
                %     MoorMat.WorstCase.Size(b,ii) = MoorMat.MaxW.Size(b,ii);
                %     MoorMat.WorstCase.Sub(b,ii) = MoorMat.MaxW.Sub(b,ii);
                % else
                %     MoorMat.WorstCase.dia(b,ii) = buoydia(b);
                %     MoorMat.WorstCase.mass(b,ii) = MoorMat.MaxC.mass(b,jj);
                %     MoorMat.WorstCase.cost(b,ii) = MoorMat.MaxC.cost(b,jj);
                %     MoorMat.WorstCase.Size(b,ii) = MoorMat.MaxC.Size(b,jj);
                %     MoorMat.WorstCase.Sub(b,ii) = MoorMat.MaxC.Sub(b,jj);
                % end

                
                %%alternate method
                if MoorMat.MaxW.Sub(b,ii) < -0.1 || MoorMat.MaxC.Sub(b,jj) < -0.1
                    MoorMat.WorstCase.dia(b,ii) = buoydia(b);
                    MoorMat.WorstCase.mass(b,ii) = MoorMat.MaxW.mass(b,ii);
                    MoorMat.WorstCase.PLmass(b,ii) = MoorMat.MaxW.PLmass(b,ii);

                    MoorMat.WorstCase.cost(b,ii) = MoorMat.MaxW.cost(b,ii);
                    MoorMat.WorstCase.Size(b,ii) = MoorMat.MaxW.Size(b,ii);
                    MoorMat.WorstCase.Sub(b,ii) = min([MoorMat.MaxW.Sub(b,ii),MoorMat.MaxC.Sub(b,jj)]);

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


        % nexttile(1)
        % hold on
        % plot(MoorMat.MaxC.mass,MoorMat.MaxC.cost,'DisplayName',locs{l},'LineWidth',2)
        % title('Max Current/Average Wave')
        % ylabel('Cost $')
        % xlabel('Mass [kg]')
        % nexttile(2)
        % hold on
        % plot(MoorMat.MaxW.mass,MoorMat.MaxW.cost,'DisplayName',locs{l},'LineWidth',2)
        % title('Max Wave/Average Current')
        % xlabel('Mass [kg]')
        % ylabel('Cost $')
        % 
        % nexttile(3)
        % hold on
        % plot(MoorMat.WorstCase.mass,MoorMat.WorstCase.cost,'DisplayName',locs{l},'LineWidth',2)
        % title('Max Wave/Average Current')
        % xlabel('Mass [kg]')
        % ylabel('Cost $')


    end
    filenm = strcat(saveL{l},'_Mooring','.mat');
    save(filenm, 'MoorMat')
end

%legend('location','eastoutside')

