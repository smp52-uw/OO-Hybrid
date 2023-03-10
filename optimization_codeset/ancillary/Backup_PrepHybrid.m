%% Backup

%PrepHybrid is a combination of PrepDies/PrepInso/PrepWind/PrepWave
%Written by Sarah Palmer (based on Trent's code)

%Nothing from prepDies is active
%Prep Inso
%make time series data adequately long
[data.met.time, ind_wi] = unique(data.met.time); %need a time series without duplicates
[data.swso,data.met.insotime] = ...
    extendToLifetime(data.met.shortwave_irradiance(ind_wi),data.met.time,uc.lifetime);
%PrepWind
%make time series data adequately long
[data.met.wind_spd,data.met.windtime] = ...
    extendToLifetime(data.met.wind_spd(ind_wi),data.met.time,uc.lifetime);
%PrepWave
opt.wave.wavepower_ra = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    (wave.Hs_ra)^2*(wave.Tp_ra); %[W], wave power at rated
%extract data
Hs = data.wave.significant_wave_height; %[m]
Tp = data.wave.peak_wave_period; %[s]
opt.wave.wavepower_ts = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    Hs.^2.*Tp./1000; %[kW/m] %timeseries of wavepower
opt.wave.L = atmo.g.*Tp.^2/(2*pi); %wavelength timeseri

%extend wavepower, time, hs, tp and wavelength timeseries
[data.wave.time, ind_wa] = unique(data.wave.time);
[opt.wave.wavepower_ts,opt.wave.time] =  ...
    extendToLifetime(opt.wave.wavepower_ts(ind_wa),data.wave.time,uc.lifetime);
[opt.wave.Hs] = extendToLifetime(Hs(ind_wa),data.wave.time,uc.lifetime);
[opt.wave.Tp] = extendToLifetime(Tp(ind_wa),data.wave.time,uc.lifetime);
[opt.wave.L] = extendToLifetime(opt.wave.L(ind_wa),data.wave.time,uc.lifetime);

%% Make Clean time series
%regularly spaced time series
if opt.wave.time(1) == data.met.windtime(1) && opt.wave.time(1) == data.met.insotime(1)
    converted_time = datetime(opt.wave.time,'ConvertFrom','datenum');
    reg_time = [converted_time(1):mode(diff(converted_time)):converted_time(end)]';   
    %interpolate to same time series
    data.met.time = reg_time; %might have to convert back to datenum
    data.met.windtime = datetime(data.met.windtime,'ConvertFrom','datenum');
    data.met.insotime = datetime(data.met.insotime,'ConvertFrom','datenum');
    opt.wave.time = datetime(opt.wave.time,'ConvertFrom','datenum');
    opt.wave.wavepower_ts = interp1(opt.wave.time,opt.wave.wavepower_ts,reg_time,'linear');
    opt.wave.Hs  = interp1(opt.wave.time,opt.wave.Hs,reg_time,'linear');
    opt.wave.Tp  = interp1(opt.wave.time,opt.wave.Tp,reg_time,'linear');
    opt.wave.L  = interp1(opt.wave.time,opt.wave.L,reg_time,'linear');
    data.met.wind_spd  = interp1(data.met.windtime,data.met.wind_spd,reg_time,'linear');
    data.swso  = interp1(data.met.insotime,data.swso,reg_time,'linear');
else %If time series don't start at the same point then need to add data to the late one
    disp('Test: start not matched')
    converted_time_wave = datetime(opt.wave.time,'ConvertFrom','datenum');  
    converted_time_wind = datetime(data.met.windtime,'ConvertFrom','datenum');
    converted_time_inso = datetime(data.met.insotime,'ConvertFrom','datenum');
    %temporary clean hour spaced time series for each resource - need to be
    %exactly as many years as the lifetime of the system
    regwave_time = [converted_time_wave(1):hours(1):converted_time_wave(end)]'; 
    regwind_time = [converted_time_wind(1):hours(1):converted_time_wind(end)]'; 
    reginso_time = [converted_time_inso(1):hours(1):converted_time_inso(end)]'; 
    %interpolate data to clean time series
    opt.wave.wavepower_ts = interp1(converted_time_wave,opt.wave.wavepower_ts,regwave_time,'linear');
    opt.wave.Hs  = interp1(converted_time_wave,opt.wave.Hs,regwave_time,'linear');
    opt.wave.Tp  = interp1(converted_time_wave,opt.wave.Tp,regwave_time,'linear');
    opt.wave.L  = interp1(converted_time_wave,opt.wave.L,regwave_time,'linear');
    data.met.wind_spd  = interp1(converted_time_wind,data.met.wind_spd,regwind_time,'linear');
    data.swso  = interp1(converted_time_inso,data.swso,reginso_time,'linear');

    %temporary clean times- need all timestamps to be on the hour so
    %finding indecies works between timeseries
    wave_time = dateshift(regwave_time,'start','hour');
    wind_time = dateshift(regwind_time,'start','hour');
    inso_time = dateshift(reginso_time,'start','hour');
    %actual start time
    time_start = min([wave_time(1),wind_time(1),inso_time(1)]); %earliest time start
    if time_start == wave_time(1)
        display('Test: wave first')
        length(data.swso)
        data.met.time = wave_time; %might have to convert back to datenum
        data.met.wind_spd = align_timeseries(data.met.time,wind_time,data.met.wind_spd);
        data.swso = align_timeseries(data.met.time,inso_time,data.swso);
        %[data_out] = align_timeseries(data.met.time,wind_time,data.met.wind_spd)
        %[data_out] = align_timeseries(reg_time,t_in,a_in)

    elseif time_start == wind_time(1)
        display('Test: wind first')
        data.met.time = wind_time; %might have to convert back to datenum
        opt.wave.wavepower_ts = align_timeseries(data.met.time,wave_time,opt.wave.wavepower_ts);
        opt.wave.Hs = align_timeseries(data.met.time,wave_time,opt.wave.Hs);
        opt.wave.Tp = align_timeseries(data.met.time,wave_time,opt.wave.Tp);
        opt.wave.L = align_timeseries(data.met.time,wave_time,opt.wave.L);
        data.swso = align_timeseries(data.met.time,inso_time,data.swso);
    elseif time_start == inso_time(1)
        display('Test: inso first')
        data.met.time = inso_time; %might have to convert back to datenum
        data.met.wind_spd = align_timeseries(data.met.time,wind_time,data.met.wind_spd);
        opt.wave.wavepower_ts = align_timeseries(data.met.time,wave_time,opt.wave.wavepower_ts);
        opt.wave.Hs = align_timeseries(data.met.time,wave_time,opt.wave.Hs);
        opt.wave.Tp = align_timeseries(data.met.time,wave_time,opt.wave.Tp);
        opt.wave.L = align_timeseries(data.met.time,wave_time,opt.wave.L);
    end
end