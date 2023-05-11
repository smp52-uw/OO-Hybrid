function [K_avg] = getMonthlyK(dataStruct,type)

%%%%%%%%%%% SET VALUES %%%%%%%%%%%%%%

met_pts = 1:length(dataStruct.met.time);
wave_pts = 1:length(dataStruct.wave.time);
curr_pts = 1:length(dataStruct.curr.time);
met_time = dataStruct.met.time(met_pts);
wave_time = dataStruct.wave.time(wave_pts);
wind = dataStruct.met.wind_spd(met_pts);
inso = dataStruct.met.shortwave_irradiance(met_pts);
hs = dataStruct.wave.significant_wave_height(wave_pts);
tp = dataStruct.wave.peak_wave_period(wave_pts);
c_vel = (dataStruct.curr.u.^2 + dataStruct.curr.v.^2).^0.5; %vel mag
curr_time = dataStruct.curr.time(curr_pts)';
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
        pts = find(dvall(:,1) == dvymu(i,1) & ...
            dvall(:,2) == dvymu(i,2));
        K_avg(i,2) = mean(K(pts),'omitnan');
    end
end

if isequal(type,'wind')
    dvall = datevec(met_time);
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
    dvall = datevec(dataStruct.met.time);
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
    c_vel = c_vel';
    dvall = datevec(curr_time);
    dvymu = unique(dvall(:,1:2),'rows');
    c_sz = size(c_vel)
    K_avg = zeros(length(dvymu),c_sz(2)+1);
    K_avg(:,1) = datenum(num2str(dvymu(:,1:2)));
    rho = 1025; %[kg/m^3]
    for i=1:c_sz(2)
        K = (1/2)*rho.*c_vel(:,i).^3;
        for j = 1:length(dvymu)
            pts = find(dvall(:,1) == dvymu(j,1) & ...
                dvall(:,2) == dvymu(j,2));
            K_avg(j,i+1) = mean(K(pts),'omitnan');
        end
    end
end