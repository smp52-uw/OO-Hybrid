function [data, opt] = prepHybrid(data,opt,uc,wave,atmo,inso,cturb)
%PrepHybrid is a combination of PrepDies/PrepInso/PrepWind/PrepWave
%Written by Sarah Palmer (based on Trent's code)
%Updates for any location: 6-26-23
%More updates for task 2 locations: 10-8-2024

%Figure out if you're working with Trent's locations or task data
input1 = isfield(data,'met'); %Trent's data
input2 = isfield(data,'wind'); %Task 2 data

%get current data
data.curr.vmag = (data.curr.u.^2 + data.curr.v.^2).^0.5; %velocity magnitude
data.curr.vmag = data.curr.vmag';
data.curr.vmag(data.curr.vmag>4) = nan; %Remove any non-physical current speeds

if opt.pltdebug && input2 %diagnostic plot of input data for Task 2 locations (raw data)
    col = colormap(brewermap(9,'Set2')); %colors
    figure(1)
    tf = tiledlayout(5,1);
    title(tf,strcat(data.title," - Unprocessed Data"))
    
    ax(1) = nexttile; %solar
    plot(data.solar.time,data.solar.swso,'linewidth',1.5,'Color',col(6,:))
    ylabel('[W/m^2]')
    title('Inso')

    ax(2) = nexttile; %wind
    plot(data.wind.time,data.wind.U,'linewidth',1.5,'color',col(1,:))
    ylabel('[m/s]')
    title('Wind Speed')

    ax(3) = nexttile; %current
    temptime = datetime(data.curr.time,'ConvertFrom','datenum');
    temptime.TimeZone = "UTC";
    plot(temptime,data.curr.vmag(:,1),'linewidth',1.5,'color',col(3,:))
    ylabel('[m/s]')
    title('Surface Current Speed')


    ax(4) = nexttile; %wave height
    plot(data.wave.time,data.wave.Hs,'linewidth',1.5,'color',col(4,:))
    ylabel('[m]')
    title('Wave Height')

    ax(5) = nexttile; %wave period
    plot(data.wave.time,data.wave.Tp,'linewidth',1.5,'color',col(4,:))
    ylabel('[s]')
    title('Wave Period')
    xlabel('Time')

    linkaxes(ax,'x')
end

%Get wave data
opt.wave.wavepower_ra = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    (wave.Hs_ra)^2*(wave.Tp_ra); %[W], wave power at rated
%extract data
if input1
    Hs = data.wave.significant_wave_height; %[m]
    Tp = data.wave.peak_wave_period; %[s]
elseif input2
    Hs = data.wave.Hs; %[m]
    Tp = data.wave.Tp; %[s]
end
opt.wave.wavepower_ts = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    Hs.^2.*Tp./1000; %[kW/m] %timeseries of wavepower
opt.wave.L = atmo.g.*Tp.^2/(2*pi); %wavelength timeseries

%Clean Time series
if input1
    [data.met.time, ind_wi] = unique(data.met.time); %need a time series without duplicates
    [data.wave.time, ind_wa] = unique(data.wave.time);
    converted_time_wind = datetime(data.met.time,'ConvertFrom','datenum');
    converted_time_inso = datetime(data.met.time,'ConvertFrom','datenum');
    converted_time_wave = datetime(data.wave.time,'ConvertFrom','datenum'); 
elseif input2
    [data.wind.time, ind_wi] = unique(data.wind.time); %need a time series without duplicates
    [data.solar.time, ind_in] = unique(data.solar.time);
    [data.wave.time, ind_wa] = unique(data.wave.time);
    converted_time_wind = data.wind.time;
    converted_time_inso = data.solar.time;
    converted_time_wave = data.wave.time;
end
[data.curr.time, ind_c] = unique(data.curr.time);
curr_vec = datevec(data.curr.time);
%curr_vec(:,1) = year(converted_time_wave(1)); % %shift the time to 2015 for currents-just works for AB
converted_time_curr = datetime(curr_vec);
if input2
    converted_time_curr.TimeZone = 'UTC'; %verified the time zone on the HYCOM user group
end

%temporary clean hour spaced time series for each resource
regwave_time = [dateshift(converted_time_wave(1),'end','hour'):hours(1):dateshift(converted_time_wave(end),'end','hour')]'; 
regwind_time = [dateshift(converted_time_wind(1),'end','hour'):hours(1):dateshift(converted_time_wind(end),'end','hour')]'; 
reginso_time = [dateshift(converted_time_inso(1),'end','hour'):hours(1):dateshift(converted_time_inso(end),'end','hour')]'; 
regcurr_time = [dateshift(converted_time_curr(1),'end','hour'):hours(1):dateshift(converted_time_curr(end),'end','hour')]'; 

time_start = min([regwave_time(1),regwind_time(1), reginso_time(1), regcurr_time(1)]); %earliest time start
year1 = year(time_start(1)); %starting year for cleaned/aligned time series
time_final = datetime([year1 01 01]):hours(1):datetime([year1+uc.lifetime 01 01]); %cleaned/aligned time series
time_final = time_final(1:end-1)'; %need one less than the entire time series so it ends at midnight 12/31


%% Fill big gaps
num_d = size(data.curr.vmag); 
if input2
    bad_wind = sum(isnan(data.wind.U));
    data.wind.U  = interp1(converted_time_wind,data.wind.U(ind_wi),regwind_time,'nearest');
    if bad_wind > 24
        data.wind.U = fillmiss_phaseavg(data.wind.U, regwin_time);
    end
    data.wind.U = fillmissing(data.wind.U,'nearest');
    disp('extend wind')
    [data.wind.U,data.wind.time] = extendToLifetime(data.wind.U,datenum(regwind_time),uc.lifetime); %make time series data adequately long

    bad_solar = sum(isnan(data.solar.swso));
    data.solar.swso  = interp1(converted_time_inso,data.solar.swso(ind_in),reginso_time,'nearest');
    if bad_solar > 24
        data.solar.swso = fillmiss_phaseavg(data.solar.swso, reginso_time);
    end
    data.solar.swso = fillmissing(data.solar.swso,'nearest');
    swso_neg_count = sum(data.solar.swso < 0);
    if swso_neg_count > 0
        disp('removing negative solar data')
        data.solar.swso(data.solar.swso < 0) = 0; %remove any negative swso points
    end
    disp('extend solar')
    [data.solar.swso,data.solar.time] = extendToLifetime(data.solar.swso,datenum(reginso_time),uc.lifetime); %make time series data adequately long
     
    for i = 1:min(num_d)
        bad_curr = sum(isnan(data.curr.vmag(:,i)));
        data.curr.speed(:,i) = interp1(converted_time_curr,data.curr.vmag(ind_c,i),regcurr_time,'nearest','extrap');
        if bad_curr > 24
            data.curr.speed(:,i) = fillmiss_phaseavg(data.curr.speed(:,i), regcurr_time);
        end
        data.curr.speed(:,i) = fillmissing(data.curr.speed(:,i),'nearest');
    end
    
    bad_wave = sum(isnan(Hs));
    opt.wave.wavepower_ts = interp1(converted_time_wave,opt.wave.wavepower_ts(ind_wa),regwave_time,'nearest');
    opt.wave.Hs  = interp1(converted_time_wave,Hs(ind_wa),regwave_time,'nearest');
    opt.wave.Tp  = interp1(converted_time_wave,Tp(ind_wa),regwave_time,'nearest');
    opt.wave.L  = interp1(converted_time_wave,opt.wave.L(ind_wa),regwave_time,'nearest');
    if bad_wave > 24
        opt.wave.wavepower_ts = fillmiss_phaseavg(opt.wave.wavepower_ts, regwave_time);
        nanwave = isnan(opt.wave.wavepower_ts);
        opt.wave.Hs = fillmiss_phaseavg(opt.wave.Hs, regwave_time);
        opt.wave.Tp  = fillmiss_phaseavg(opt.wave.Tp, regwave_time);
        opt.wave.L = fillmiss_phaseavg(opt.wave.L, regwave_time);
    else
        nanwave = [];
    end
    opt.wave.wavepower_ts = fillmissing(opt.wave.wavepower_ts,'nearest');
    opt.wave.Hs = fillmissing(opt.wave.Hs,'nearest');
    opt.wave.Tp = fillmissing(opt.wave.Tp,'nearest');
    opt.wave.L = fillmissing(opt.wave.L,'nearest');
end
%interpolate data to clean time series
%The clean time series will have a time slightly before the first data
%point and the interpolation makes that value go to zero - fixed with the
%fill missing
if input1
    data.met.wind_spd  = interp1(converted_time_wind,data.met.wind_spd(ind_wi),regwind_time,'linear');
    data.met.wind_spd = fillmissing(data.met.wind_spd,'nearest');
    data.swso = fillmissing(data.met.shortwave_irradiance,'linear'); %[W/m^2]
    data.swso  = interp1(converted_time_inso,data.swso(ind_wi),reginso_time,'linear');
    data.swso = fillmissing(data.swso,'nearest');
    swso_neg_count = sum(data.swso < 0);
    if swso_neg_count > 0
        disp('removing negative solar data')
        data.swso(data.swso < 0) = 0; %remove any negative swso points
    end
    disp('extend solar')
    [data.swso,data.met.insotime] = extendToLifetime(data.swso,datenum(reginso_time),uc.lifetime); %make time series data adequately long
    disp('extend wind')
    [data.met.wind_spd,data.met.windtime] = extendToLifetime(data.met.wind_spd,datenum(regwind_time),uc.lifetime); %make time series data adequately long

    for i=1:num_d(2)
        data.curr.speed(:,i) = interp1(converted_time_curr,data.curr.vmag(ind_c,i),regcurr_time,'nearest','extrap');
        data.curr.speed(:,i) = fillmissing(data.curr.speed(:,i),'nearest');
    end

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
end

%extend wavepower, time, hs, tp and wavelength timeseries
disp('extend wave')
[opt.wave.wavepower_ts,opt.wave.time] =  ...
    extendToLifetime(opt.wave.wavepower_ts,datenum(regwave_time),uc.lifetime);
[opt.wave.Hs] = extendToLifetime(opt.wave.Hs,datenum(regwave_time),uc.lifetime);
[opt.wave.Tp] = extendToLifetime(opt.wave.Tp,datenum(regwave_time),uc.lifetime);
[opt.wave.L] = extendToLifetime(opt.wave.L,datenum(regwave_time),uc.lifetime);

%extend current speed, and time to life
disp('extend current')
[data.curr.speed6(:,1),data.curr.time] = extendToLifetime(data.curr.speed(:,1),datenum(regcurr_time),uc.lifetime);
for i=2:num_d(2)
    [data.curr.speed6(:,i)] = extendToLifetime(data.curr.speed(:,i),datenum(regcurr_time),uc.lifetime);
end
%% Make Clean time series
if input1
    if opt.wave.time(1) == data.met.windtime(1) && opt.wave.time(1) == data.met.insotime(1) && opt.wave.time(1) == data.curr.time(1)
        disp('Time series are aligned')
    else %If time series don't start at the same point then need to add data to the late one
        disp('Test: start not matched')
        data.met.time = datenum(time_final);
        
        opt.wave.wavepower_ts = align_timeseries(data.met.time,opt.wave.time,opt.wave.wavepower_ts,uc);
        opt.wave.Hs = align_timeseries(data.met.time,opt.wave.time,opt.wave.Hs,uc);
        opt.wave.Tp = align_timeseries(data.met.time,opt.wave.time,opt.wave.Tp,uc);
        opt.wave.L = align_timeseries(data.met.time,opt.wave.time,opt.wave.L,uc);
        for i=1:num_d(2)
            data.curr.speed6a(:,i) = align_timeseries(data.met.time,data.curr.time,data.curr.speed6(:,i),uc);
        end
    
        data.met.wind_spd = align_timeseries(data.met.time,data.met.windtime,data.met.wind_spd,uc);
        data.swso = align_timeseries(data.met.time,data.met.insotime,data.swso,uc);
    
        data.curr.time = datenum(data.met.time);
        opt.wave.time = datenum(data.met.time); %might have to convert back to datenum
        data.met.time = datenum(data.met.time);
    end
elseif input2
    if opt.wave.time(1) == data.wind.time(1) && opt.wave.time(1) == data.solar.time(1)
        if opt.timeadj > 1
            opt.wave.time = opt.wave.time(opt.timeadj:end); %jump forward 2 months
            data.wind.time = data.wind.time(opt.timeadj:end);
            data.solar.time = data.wind.time(opt.timeadj:end);

            [data.wind.U,data.wind.time] = extendToLifetime(data.wind.U(opt.timeadj:end),datenum(data.wind.time),uc.lifetime); %make time series data adequately long
            [data.solar.swso,data.solar.time] = extendToLifetime(data.solar.swso(opt.timeadj:end),datenum(data.solar.time),uc.lifetime); %make time series data adequately long


            [opt.wave.wavepower_ts,opt.wave.time] = extendToLifetime(opt.wave.wavepower_ts(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
            [opt.wave.Hs] = extendToLifetime(opt.wave.Hs(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
            [opt.wave.Tp] = extendToLifetime(opt.wave.Tp(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
            [opt.wave.L] = extendToLifetime(opt.wave.L(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
        end
        for i=1:num_d(2) %current data isn't aligned
            data.curr.speed6a(:,i) = align_timeseries(data.wind.time,data.curr.time,data.curr.speed6(:,i),uc);
        end
        data.curr.time = datenum(data.wind.time);
        data.met.time = data.wind.time;
        data.swso = data.solar.swso;
        data.met.wind_spd = data.wind.U;
        data.met.wind_ht = 10; %all ERA5 data is at 10 m
        disp('Time series are aligned')
    else
        disp('ERROR - Task2 data should already be aligned...')
    end
end

%%Plot processed data
if opt.pltdebug && input2 %diagnostic plot of input data for Task 2 locations (aligned data)
    clear ax
    figure(2)
    tf = tiledlayout(4,1);
    title(tf,strcat(data.title," - Processed Data"))

    ax(1) = nexttile; %solar
    plot(datetime(data.met.time,'convertfrom','datenum'), data.swso,'linewidth',1.5,'Color',col(6,:))
    ylabel('[W/m^2]')
    title('Inso')

    ax(2) = nexttile; %wind
    plot(datetime(data.met.time,'convertfrom','datenum'),data.met.wind_spd,'linewidth',1.5,'color',col(1,:))
    ylabel('[m/s]')
    title('Wind Speed')

    ax(3) = nexttile; %current
    plot(datetime(data.curr.time,'convertfrom','datenum'),data.curr.speed6a(:,1),'linewidth',1.5,'color',col(3,:))
    hold on
    yline(cturb.uci,'linewidth',2)
    ylabel('[m/s]')
    title('Surface Current Speed')

    ax(4) = nexttile; %wave height
    plot(datetime(data.met.time,'convertfrom','datenum'),opt.wave.wavepower_ts,'linewidth',1.5,'color',col(4,:))
    hold on
    plot(datetime(data.met.time(nanwave),'convertfrom','datenum'),opt.wave.wavepower_ts(nanwave),'kx')
    ylabel('[W/m ??]')
    title('Wave Power')
    xlabel('Time')

    %linkaxes(ax,'x')
    %xlim([min(data.met.time),max(data.met.time)])
end
%% winter cleaning
if inso.cleanstrat == 3 || inso.cleanstrat == 4 %winter cleaning
    %winter cleaning (if applicable)
    if data.lat < 0 %southern hemisphere
        wint_clean_mo = 5; %may
    else
        wint_clean_mo = 11; %november
    end
    if input1
        dv = datevec(data.met.time);
    elseif input2
        dv = datevec(data.solar.time);
    end
    data.wint_clean_ind = find(dv(:,2) == wint_clean_mo & dv(:,3) == 1 & dv(:,4) == 0);
    if inso.cleanstrat == 4 %every other winter
        data.wint_clean_ind = data.wint_clean_ind(2:2:end);
    end
end

%% Get current power
[opt.curr.F_curr] = Current_Power(cturb,atmo);
%current_Pmatrix = load('current_Pmatrix.mat');

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

elseif wave.method == 3 %Trevor Power Matrices
    
    %load wec sim results into structure
%     wsr_1 = load('struct1m_opt');
%     wsr_1 = wsr_1.('struct1m_opt');
    wsr_2 = load('2m_5kW_WEC_Power_Frame.mat');
    wsr_2.pframe = wsr_2.('WEC_pframe');
    wsr_2.H = wsr_2.pframe(:,1); %Wave Height
    wsr_2.T = wsr_2.pframe(:,2); %Wave Period
    wsr_2.I = wsr_2.pframe(:,3:end); %current
    wsr_2.P = 48.*wsr_2.I.*-1; %Power = V*I = 48 bus voltage *I
    wsr_2.b = 2; %width
    wsr_3 = load('3m_10kW_WEC_Power_Frame.mat');
    wsr_3 = wsr_3.('WEC_pframe');
    wsr_4 = load('4m_20kW_WEC_Power_Frame.mat');
    wsr_4 = wsr_4.('WEC_pframe');
    wsr_5 = load('5m_40kW_WEC_Power_Frame.mat');
    wsr_5 = wsr_5.('WEC_pframe');
%     wsr_6 = load('struct6m_opt');
%     wsr_6 = wsr_6.('struct6m_opt');
    s(4) = struct();
    %s(1).wsr = wsr_1;
    s(1).wsr = wsr_2;
    s(2).wsr = wsr_3;
    s(3).wsr = wsr_4;
    s(4).wsr = wsr_5;
    %s(6).wsr = wsr_6;
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
        B = wsr.b*ones(n,1);       
        J = (1/(64*pi))*atmo.rho_w*atmo.g^2.*H.^2.*T; %find wave power
        %P = reshape(s(b).wsr.mat',n,1); %find wec power (use mat not P)
        P = s(b).wsr.P;
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