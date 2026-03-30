%create data files for Hybrid Optimization

%update to use the current structure that LoadSites creates rather than
%duplicating code (12/31/25)
%More updates to pull in the new summer start data which was processed
%differently (1/8/25)

%More updates to pull in new WETS and PISCES

%location
clear locout
site_case = 3;

switch site_case
    case 1
        loc = 'PacWave';
        cd PacWave\
    case 2 %Old WETS - do not use, the surface roughness values are wrong
        % loc = 'WETS';
        % cd WETS\
    case 3
        loc = 'Mid-Atlantic';
        locout = 'MidAtlSB';
        cd MidAtlSB\
    case 4 %old PISCES - do not use, the surface roughness values are wrong
        % loc = 'PISCES';
        % cd PISCES\
    case 5
        loc = 'Port_Hueneme';
        cd Port_Hueneme\
    case 6
        loc = "BerSea";
        cd BerSea\
    case 7
        loc = 'SFOMF';
        cd \SFOMF
    case 8
        loc = 'altPISCES';
        cd altPISCES;
    case 9
        loc = 'altWETS';
        cd altWETS;
end

if ~exist('locout','var')
    locout = loc; %makes the out file name the same as the input file except for mid atlantic
end

matfiles = dir('*.mat') ; 
N = length(matfiles) ; 
for i = 1:N
    if ~strcmp(matfiles(i).name, strcat(loc,'.mat')) %If I've already run createAllDataHybrid I don't want to load in that file
        load(matfiles(i).name)
    end
end

data.title = loc;
data.depth = site_depth;
data.dist = 'nan';
data.wave = wave;
data.wind = wind;
data.solar = solar;

data.lat = lat;
data.lon = lon;
data.curr.vmag = current.U_z;
%data.curr.v = v;
data.curr.depth = current.z;
data.curr.time = current.time;

data.temperature = temperature;

save(locout,"data")