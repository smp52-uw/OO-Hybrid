function [K_ts] = getPowerDensity(dataStruct,type)

%%%%%%%%%%% SET VALUES %%%%%%%%%%%%%%
test1 = isfield(dataStruct,'met'); %Trent's data
test2 = isfield(dataStruct,'wind'); %Task 2 data
if test1
    met_pts = 1:length(dataStruct.met.time);
    wave_pts = 1:length(dataStruct.wave.time);
    met_time = dataStruct.met.time(met_pts);
    wind_time = met_time;
    solar_time = met_time;
    wave_time = dataStruct.wave.time(wave_pts);
    wind = dataStruct.met.wind_spd(met_pts);
    inso = dataStruct.met.shortwave_irradiance(met_pts);
    hs = dataStruct.wave.significant_wave_height(wave_pts);
    tp = dataStruct.wave.peak_wave_period(wave_pts);
    curr_pts = 1:length(dataStruct.curr.time);
    curr_time = dataStruct.curr.time(curr_pts)';
    c_vel = sqrt(dataStruct.curr.u(curr_pts).^2 + dataStruct.curr.v(curr_pts).^2);
    c_vel = c_vel';
elseif test2
    wind_pts = 1:length(dataStruct.wind.time);
    wave_pts = 1:length(dataStruct.wave.time);
    solar_pts = 1:length(dataStruct.solar.time);
    curr_pts = 1:length(dataStruct.curr.time);
    wind_time = datenum(dataStruct.wind.time(wind_pts));
    wave_time = datenum(dataStruct.wave.time(wave_pts));
    solar_time = datenum(dataStruct.solar.time(solar_pts));
    curr_time = dataStruct.curr.time(curr_pts);
    wind = dataStruct.wind.U(wind_pts);
    inso = dataStruct.solar.swso(solar_pts);
    hs = dataStruct.wave.Hs(wave_pts);
    tp = dataStruct.wave.Tp(wave_pts);
    c_vel = sqrt(dataStruct.curr.u(curr_pts).^2 + dataStruct.curr.v(curr_pts).^2);
    c_vel = c_vel';
    c_vel(c_vel>4) = nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(type,'wave')
    K_ts(:,1) = wave_time;

    rho = 1020; %[kg/m^3]
    g = 9.81; %[m/s^2]
    K_ts(:,2) = (1/(16*4*pi))*rho*g^2.*hs(:).^2.*tp(:);
end

if isequal(type,'wind')

    K_ts(:,1) = wind_time;
    rho = 1.225; %[kg/m^3]
    K_ts(:,2) = (1/2)*rho.*wind.^3;
end

if isequal(type,'inso')
    K_ts(:,1) = solar_time;
    K_ts(:,2) = inso;
end

if isequal(type,'curr')
    K_ts(:,1) = curr_time;
    c_sz = size(c_vel);
    rho = 1025; %[kg/m^3]
    K_ts(:,2) = (1/2)*rho.*c_vel.^3;

end