function [monthlyArray, K_avg] = getMonthlyK(dataStruct,type)
%Run after prepHybrid
%all time vectors will be the same and in datenum format

time = datetime(dataStruct.data.met.time,'ConvertFrom','datenum');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startDate = time(1);
endDate = time(end);
monthlyArray = startDate : calmonths(1) : endDate;

if isequal(type,'wave')
    K_avg = zeros(length(monthlyArray),1);
    K = dataStruct.opt.wave.wavepower_ts; %[kW/m]
    for i = 1:length(K_avg)
        pts = find(month(time) == month(monthlyArray(i)) & year(time) == year(monthlyArray(i)));

        K_avg(i) = mean(K(pts),'omitnan');
    end
end

if isequal(type,'wind')
    K_avg = zeros(length(monthlyArray),1);
    rho = dataStruct.data.wind.rho; %[kg/m^3]
    K = (1/2)*rho.*dataStruct.data.wind.U.^3; %This is the one without ice
    for i = 1:length(K_avg)
        pts = find(month(time) == month(monthlyArray(i)) & year(time) == year(monthlyArray(i)));
        K_avg(i) = mean(K(pts),'omitnan');
    end
end

if isequal(type,'inso')
    K_avg = zeros(length(monthlyArray),1);
    K = dataStruct.data.swso_noice;
    for i = 1:length(K_avg)
        pts = find(month(time) == month(monthlyArray(i)) & year(time) == year(monthlyArray(i)));
        K_avg(i) = mean(K(pts),'omitnan');
    end
end

if isequal(type,'curr')
    K_avg = zeros(length(monthlyArray),1);
    rho = 1025; %[kg/m^3]
    c_vel = dataStruct.data.curr.speed6a;
    K = (1/2)*rho.*c_vel(:,1).^3;
    for i = 1:length(K_avg)
        pts = find(month(time) == month(monthlyArray(i)) & year(time) == year(monthlyArray(i)));
        K_avg(i) = mean(K(pts),'omitnan');
    end
end