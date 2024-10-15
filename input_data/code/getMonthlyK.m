function [K_avg] = getMonthlyK(dataStruct,type)

%%%%%%%%%%% SET VALUES %%%%%%%%%%%%%%
test1 = isfield(dataStruct,'met');
test2 = isfield(dataStruct,'swso');
test3 = isfield(dataStruct,'wind');
if test1
    met_pts = 1:length(dataStruct.met.time);
    met_time = dataStruct.met.time(met_pts);
    wind_time = met_time;
    solar_time = met_time;
    wind = dataStruct.met.wind_spd(met_pts);
    if test2
        inso = dataStruct.swso(met_pts);
    else
        inso = dataStruct.met.shortwave_irradiance(met_pts);
    end
elseif test3 %task2 data
    wind_pts = 1:length(dataStruct.wind.time);
    solar_pts = 1:length(dataStruct.solar.time);
    wind_time = datenum(dataStruct.wind.time(wind_pts));
    solar_time = datenum(dataStruct.solar.time(solar_pts));
    wind = dataStruct.wind.U(wind_pts);
    inso = dataStruct.solar.swso(solar_pts);
end


test = isfield(dataStruct.wave,'significant_wave_height');
if test
    wave_pts = 1:length(dataStruct.wave.time);
    wave_time = dataStruct.wave.time(wave_pts);
    hs = dataStruct.wave.significant_wave_height(wave_pts);
    tp = dataStruct.wave.peak_wave_period(wave_pts);
else
    wave_pts = 1:length(dataStruct.wave.time);
    wave_time = datenum(dataStruct.wave.time(wave_pts));
    hs = dataStruct.wave.Hs(wave_pts);
    tp = dataStruct.wave.Tp(wave_pts);
end
test = isfield(dataStruct.curr,'time');
test1 = isfield(dataStruct.curr,'u');
test2 = isfield(dataStruct.curr,'speed6a');
if test
    curr_pts = 1:length(dataStruct.curr.time);
    curr_time = dataStruct.curr.time(curr_pts)';
end
if test2
    c_vel = dataStruct.curr.speed6a;
elseif test1
    c_vel = (dataStruct.curr.u.^2 + dataStruct.curr.v.^2).^0.5; %vel mag
    c_vel = c_vel';
    c_vel(c_vel>4) = nan;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isequal(type,'wave')

    dvall = datevec(wave_time);
    dvymu = unique(dvall(:,1:2),'rows');
    K_avg = zeros(length(dvymu),2);
    K_avg(:,1) = datenum(num2str(dvymu(:,1:2)));
    rho = 1025; %[kg/m^3]
    g = 9.81; %[m/s^2]
    K = (1/(16*4*pi))*rho*g^2.*hs(:).^2 ...
        .*tp(:);
    for i = 1:length(dvymu)
        %disp(dvall(i,1))
        pts = find(dvall(:,1) == dvymu(i,1) & ...
            dvall(:,2) == dvymu(i,2));
        K_avg(i,2) = mean(K(pts),'omitnan');
    end
end

if isequal(type,'wind')
    dvall = datevec(wind_time);
    dvymu = unique(dvall(:,1:2),'rows');
    K_avg = zeros(length(dvymu),2);
    K_avg(:,1) = datenum(num2str(dvymu(:,1:2)));
    rho = 1.2; %[kg/m^3]
    K = (1/2)*rho.*wind.^3;
    for i = 1:length(dvymu)
        pts = find(dvall(:,1) == dvymu(i,1) & ...
            dvall(:,2) == dvymu(i,2));
        K_avg(i,2) = mean(K(pts),'omitnan');
    end
end

if isequal(type,'inso')
    dvall = datevec(solar_time);
    dvymu = unique(dvall(:,1:2),'rows');
    K_avg = zeros(length(dvymu),2);
    K_avg(:,1) = datenum(num2str(dvymu(:,1:2)));
    K = inso;
    for i = 1:length(dvymu)
        pts = find(dvall(:,1) == dvymu(i,1) & ...
            dvall(:,2) == dvymu(i,2));
        K_avg(i,2) = mean(K(pts),'omitnan');
    end
end

if isequal(type,'curr')
    dvall = datevec(curr_time);
    dvymu = unique(dvall(:,1:2),'rows');
    c_sz = size(c_vel);
    K_avg = zeros(length(dvymu),c_sz(2)+1);
    K_avg(:,1) = datenum(num2str(dvymu(:,1:2)));
    rho = 1025; %[kg/m^3]
    for i=1:c_sz(2)
        K = (1/2)*rho.*c_vel(:,i).^3;
        for j = 1:length(dvymu)
            pts = find(dvall(:,1) == dvymu(j,1) & ...
                dvall(:,2) == dvymu(j,2));
            %disp('test')
            %disp(length(pts))
            K_avg(j,i+1) = mean(K(pts),'omitnan');
        end
    end
end