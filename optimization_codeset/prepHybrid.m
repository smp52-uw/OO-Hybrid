function [data, opt] = prepHybrid(data,opt,uc,wave,atmo,inso)
%PrepHybrid is a combination of PrepDies/PrepInso/PrepWind/PrepWave
%Written by Sarah Palmer (based on Trent's code)

%Nothing from prepDies is active

%Get wave data
opt.wave.wavepower_ra = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    (wave.Hs_ra)^2*(wave.Tp_ra); %[W], wave power at rated
%extract data
Hs = data.wave.significant_wave_height; %[m]
Tp = data.wave.peak_wave_period; %[s]
opt.wave.wavepower_ts = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    Hs.^2.*Tp./1000; %[kW/m] %timeseries of wavepower
opt.wave.L = atmo.g.*Tp.^2/(2*pi); %wavelength timeseries

%Clean Time series
[data.met.time, ind_wi] = unique(data.met.time); %need a time series without duplicates
[data.wave.time, ind_wa] = unique(data.wave.time);
converted_time_wave = datetime(data.wave.time,'ConvertFrom','datenum');  
converted_time_wind = datetime(data.met.time,'ConvertFrom','datenum');
converted_time_inso = datetime(data.met.time,'ConvertFrom','datenum');
%temporary clean hour spaced time series for each resource - need to be
%exactly as many years as the lifetime of the system
regwave_time = [dateshift(converted_time_wave(1),'end','hour'):hours(1):dateshift(converted_time_wave(end),'end','hour')]'; 
regwind_time = [dateshift(converted_time_wind(1),'end','hour'):hours(1):dateshift(converted_time_wind(end),'end','hour')]'; 
reginso_time = [dateshift(converted_time_inso(1),'end','hour'):hours(1):dateshift(converted_time_inso(end),'end','hour')]'; 

%interpolate data to clean time series
%The clean time series will have a time slightly before the first data
%point and the interpolation makes that value go to zero
data.met.wind_spd  = interp1(converted_time_wind,data.met.wind_spd(ind_wi),regwind_time,'linear');
data.met.wind_spd = fillmissing(data.met.wind_spd,'nearest');
data.swso = fillmissing(data.met.shortwave_irradiance,'linear'); %[W/m^2]
data.swso  = interp1(converted_time_inso,data.swso(ind_wi),reginso_time,'linear');
data.swso = fillmissing(data.swso,'nearest');

%Prep Inso
%make time series data adequately long
[data.swso,data.met.insotime] = ...
    extendToLifetime(data.swso,datenum(reginso_time),uc.lifetime);
%PrepWind
%make time series data adequately long
[data.met.wind_spd,data.met.windtime] = ...
    extendToLifetime(data.met.wind_spd,datenum(regwind_time),uc.lifetime);
%PrepWave
%Interpolate wave data to clean time series
opt.wave.wavepower_ts = interp1(converted_time_wave,opt.wave.wavepower_ts(ind_wa),regwave_time,'linear');
opt.wave.Hs  = interp1(converted_time_wave,Hs(ind_wa),regwave_time,'linear');
opt.wave.Tp  = interp1(converted_time_wave,Tp(ind_wa),regwave_time,'linear');
opt.wave.L  = interp1(converted_time_wave,opt.wave.L(ind_wa),regwave_time,'linear');
opt.wave.wavepower_ts = fillmissing(opt.wave.wavepower_ts,'nearest');
opt.wave.Hs = fillmissing(opt.wave.Hs,'nearest');
opt.wave.Tp = fillmissing(opt.wave.Tp,'nearest');
opt.wave.L = fillmissing(opt.wave.L,'nearest');

%extend wavepower, time, hs, tp and wavelength timeseries
[opt.wave.wavepower_ts,opt.wave.time] =  ...
    extendToLifetime(opt.wave.wavepower_ts,datenum(regwave_time),uc.lifetime);
[opt.wave.Hs] = extendToLifetime(opt.wave.Hs,datenum(regwave_time),uc.lifetime);
[opt.wave.Tp] = extendToLifetime(opt.wave.Tp,datenum(regwave_time),uc.lifetime);
[opt.wave.L] = extendToLifetime(opt.wave.L,datenum(regwave_time),uc.lifetime);

%% Make Clean time series
%regularly spaced time series
if opt.wave.time(1) == data.met.windtime(1) && opt.wave.time(1) == data.met.insotime(1)
    disp('Time series are aligned')
%     converted_time = datetime(opt.wave.time,'ConvertFrom','datenum');
%     reg_time = [converted_time(1):mode(diff(converted_time)):converted_time(end)]';   
%     %interpolate to same time series
%     data.met.time = reg_time; %might have to convert back to datenum
%     data.met.windtime = datetime(data.met.windtime,'ConvertFrom','datenum');
%     data.met.insotime = datetime(data.met.insotime,'ConvertFrom','datenum');
%     opt.wave.time = datetime(opt.wave.time,'ConvertFrom','datenum');
%     opt.wave.wavepower_ts = interp1(opt.wave.time,opt.wave.wavepower_ts,reg_time,'linear');
%     opt.wave.Hs  = interp1(opt.wave.time,opt.wave.Hs,reg_time,'linear');
%     opt.wave.Tp  = interp1(opt.wave.time,opt.wave.Tp,reg_time,'linear');
%     opt.wave.L  = interp1(opt.wave.time,opt.wave.L,reg_time,'linear');
%     data.met.wind_spd  = interp1(data.met.windtime,data.met.wind_spd,reg_time,'linear');
%     data.swso  = interp1(data.met.insotime,data.swso,reg_time,'linear');
else %If time series don't start at the same point then need to add data to the late one
    disp('Test: start not matched')

    %temporary clean times- need all timestamps to be on the hour so
    %finding indecies works between timeseries
%     wave_time = dateshift(regwave_time,'start','hour');
%     wind_time = dateshift(regwind_time,'start','hour');
%     inso_time = dateshift(reginso_time,'start','hour');
    %actual start time
    time_start = min([opt.wave.time(1),data.met.insotime(1), data.met.windtime(1)]); %earliest time start
    if time_start == opt.wave.time(1)
        display('Test: wave first')
        data.met.wind_spd = align_timeseries(opt.wave.time,data.met.windtime,data.met.wind_spd);
        data.swso = align_timeseries(opt.wave.time,data.met.insotime,data.swso);
        data.met.time = datenum(opt.wave.time); %might have to convert back to datenum
        %[data_out] = align_timeseries(data.met.time,wind_time,data.met.wind_spd)
        %[data_out] = align_timeseries(reg_time,t_in,a_in)

    elseif time_start == data.met.windtime(1)
        display('Test: wind first')
        data.met.time = data.met.windtime;
        opt.wave.wavepower_ts = align_timeseries(data.met.time,opt.wave.time,opt.wave.wavepower_ts);
        opt.wave.Hs = align_timeseries(data.met.time,opt.wave.time,opt.wave.Hs);
        opt.wave.Tp = align_timeseries(data.met.time,opt.wave.time,opt.wave.Tp);
        opt.wave.L = align_timeseries(data.met.time,opt.wave.time,opt.wave.L);
        data.met.time = datenum(data.met.time);
        opt.wave.time = datenum(data.met.time); %might have to convert back to datenum
        %data.swso = align_timeseries(data.met.time,inso_time,data.swso);
%     elseif time_start == inso_time(1)
%         display('Test: inso first')
%         data.met.time = inso_time; %might have to convert back to datenum
%         data.met.wind_spd = align_timeseries(data.met.time,wind_time,data.met.wind_spd);
%         opt.wave.wavepower_ts = align_timeseries(data.met.time,wave_time,opt.wave.wavepower_ts);
%         opt.wave.Hs = align_timeseries(data.met.time,wave_time,opt.wave.Hs);
%         opt.wave.Tp = align_timeseries(data.met.time,wave_time,opt.wave.Tp);
%         opt.wave.L = align_timeseries(data.met.time,wave_time,opt.wave.L);
%     end
end

%winter cleaning
if inso.cleanstrat == 3 || inso.cleanstrat == 4 %winter cleaning
    %winter cleaning (if applicable)
    if data.lat < 0 %southern hemisphere
        wint_clean_mo = 5; %may
    else
        wint_clean_mo = 11; %november
    end
    dv = datevec(data.met.time);
    data.wint_clean_ind = find(dv(:,2) == wint_clean_mo & dv(:,3) == 1 ...
        & dv(:,4) == 0);
    if inso.cleanstrat == 4 %every other winter
        data.wint_clean_ind = data.wint_clean_ind(2:2:end);
    end
end

% %set mooring system
% if econ.platform.inso.boundary == 1
%     econ.platform.inso.depth(:,4) = econ.platform.inso.depth(:,1);
%     econ.platform.inso.diameter(:,4) =  ...
%         econ.platform.inso.boundary_di.*ones(5,1);
%     econ.platform.inso.cost(:,4) = ...
%         econ.platform.inso.cost(:,3).*econ.platform.inso.boundary_mf;
% end

%PrepWave
if wave.method == 1 %divide by B methodology
    disp('Error - Divide by B method chosen')
    
elseif wave.method == 2 %3d interpolation methodology
    
    %load wec sim results into structure
    wsr_1 = load('struct1m_opt');
    wsr_1 = wsr_1.('struct1m_opt');
    wsr_2 = load('struct2m_opt');
    wsr_2 = wsr_2.('struct2m_opt');
    wsr_3 = load('struct3m_opt');
    wsr_3 = wsr_3.('struct3m_opt');
    wsr_4 = load('struct4m_opt');
    wsr_4 = wsr_4.('struct4m_opt');
    wsr_5 = load('struct5m_opt');
    wsr_5 = wsr_5.('struct5m_opt');
    wsr_6 = load('struct6m_opt');
    wsr_6 = wsr_6.('struct6m_opt');
    s(6) = struct();
    s(1).wsr = wsr_1;
    s(2).wsr = wsr_2;
    s(3).wsr = wsr_3;
    s(4).wsr = wsr_4;
    s(5).wsr = wsr_5;
    s(6).wsr = wsr_6;
    %preallocate scatter arrays
    H_scat = [];
    T_scat = [];
    B_scat = [];
    CWR_scat = [];
    for b = 1:length(s)
        n = length(s(b).wsr.H);
        if ~isequal(n,length(s(b).wsr.T))
            error('Tp and Hs vectors are not equal in length.')
        end
        H = s(b).wsr.H;
        T = s(b).wsr.T;
        B = b*ones(n,1);       
        J = (1/(64*pi))*atmo.rho_w*atmo.g^2.*H.^2.*T; %find wave power
        P = reshape(s(b).wsr.mat',n,1); %find wec power (use mat not P)
        CWR = P./(J.*B); %find cwr
        %populate scatter arrays
        T_scat = [T_scat ; T];
        H_scat = [H_scat ; H];    
        B_scat = [B_scat ; B];
        CWR_scat = [CWR_scat ; CWR];
    end
    %create scattered interpolant
    opt.wave.F =  scatteredInterpolant(T_scat,H_scat,B_scat,CWR_scat);
    %create width-rated power function
    opt.wave.B_func = zeros(2,wave.B_func_n); %preallocation function
    opt.wave.B_func(1,:) = linspace(min(B_scat),max(B_scat), ...
        wave.B_func_n); %B [m]
    for i = 1:wave.B_func_n
        %find Gr for each B value
        opt.wave.B_func(2,i) = (opt.wave.B_func(1,i)*wave.eta_ct* ...
            opt.wave.F(wave.Tp_ra,wave.Hs_ra,opt.wave.B_func(1,i))* ...
            opt.wave.wavepower_ra)/((1+wave.house)*(1000));
    end
    
end
end