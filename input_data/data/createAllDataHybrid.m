%create data files for Hybrid Optimization

%location
loc = "BerSea";
cd BerSea\
matfiles = dir('*.mat') ; 
N = length(matfiles) ; 
for i = 1:N
    load(matfiles(i).name)
end

data.title = loc;
data.depth = site_depth;
data.dist = 'nan';
data.wave = wave;
data.wind = wind;
data.solar = solar;

data.lat = lat;
data.lon = lon;
data.curr.u = u;
data.curr.v = v;
data.curr.depth = depth;
data.curr.time = time;

save(loc,"data")