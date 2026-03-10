function [data, opt] = prepHybrid(data,opt,uc,wave,atmo,inso,cturb)
%PrepHybrid is a combination of PrepDies/PrepInso/PrepWind/PrepWave
%Written by Sarah Palmer (based on Trent's code)
%Updates for any location: 6-26-23
%More updates for task 2 locations: 10-8-2024

%Figure out if you're working with Trent's locations or task data
input1 = isfield(data,'met'); %Trent's data
if input1
    disp('ERROR - Using Trent data but task 2 prepHybrid ran')
    return
end

%get current data
%data.curr.vmag = (data.curr.u.^2 + data.curr.v.^2).^0.5; %velocity magnitude
data.curr.vmag = data.curr.vmag';
data.curr.vmag(data.curr.vmag>4) = nan; %Remove any non-physical current speeds

if opt.pltdebug %diagnostic plot of input data for Task 2 locations (raw data)
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
    plot(data.curr.time,data.curr.vmag(:,1),'linewidth',1.5,'color',col(3,:))
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

%length of the time series originally
data.wind.lt = length(data.wind.time);
data.solar.lt = length(data.solar.time);
data.wave.lt = length(data.wave.time);
data.curr.lt = length(data.curr.time);
data.temperature.lt = length(data.temperature.time);

%Clean Time series
[data.wind.time, data.wind.indu] = unique(data.wind.time); %need a time series without duplicates
[data.solar.time, data.solar.indu] = unique(data.solar.time);
[data.wave.time, data.wave.indu] = unique(data.wave.time);
[data.curr.time, data.curr.indu] = unique(data.curr.time);
[data.temperature.time, data.temperature.indu] = unique(data.temperature.time);

%cut data to match start time
strnm = fieldnames(data); %all fields in data
datafields = {'wind','solar','curr','wave','temperature'};
for fn = 1:length(strnm)
    if sum(strcmp(strnm{fn},datafields)) %make sure this is one of the time series structures
        if month(data.(strnm{fn}).time(1))~= opt.monthstart
            %find first index where the date matches the beginning of the month
            indstart = find(month(data.(strnm{fn}).time) == opt.monthstart & day(data.(strnm{fn}).time) == 1 & hour(data.(strnm{fn}).time) == 0,1,'first');
            indend = find(month(data.(strnm{fn}).time) == opt.monthstart & day(data.(strnm{fn}).time) == 1 & hour(data.(strnm{fn}).time) == 0,1,'last') - 1;
        else %should only be active for current data
            indstart = 1;
            indend = find(month(data.(strnm{fn}).time) == month(data.(strnm{fn}).time(1)) & day(data.(strnm{fn}).time) == day(data.(strnm{fn}).time(1)) & hour(data.(strnm{fn}).time) == hour(data.(strnm{fn}).time(1)),1,'last');
        end
    
        varnm = fieldnames(data.(strnm{fn}));
        timeN = find(strcmp(varnm,'time'));
        data.(strnm{fn}).(varnm{timeN}) = data.(strnm{fn}).(varnm{timeN})(indstart:indend);
        for vn = 1:length(varnm)
            if ~strcmp(varnm{vn},'time') && ~strcmp(varnm{vn},'indu') && length(data.(strnm{fn}).(varnm{vn})) == data.(strnm{fn}).lt
                if min(size(data.(strnm{fn}).(varnm{vn}))) == 1 %vector
                    data.(strnm{fn}).(varnm{vn}) = data.(strnm{fn}).(varnm{vn})(data.(strnm{fn}).indu); %remove non unique values
                    data.(strnm{fn}).(varnm{vn}) = data.(strnm{fn}).(varnm{vn})(indstart:indend); %cut to correct start
                else %matrix of current speeds
                    data.(strnm{fn}).(varnm{vn}) = data.(strnm{fn}).(varnm{vn})(data.(strnm{fn}).indu,:); %remove non unique values
                    data.(strnm{fn}).(varnm{vn}) = data.(strnm{fn}).(varnm{vn})(indstart:indend,:); %cut to correct start
                end
            end
        end
    end
end
    
%temporary clean hour spaced time series for each resource
regwave_time = [dateshift(data.wave.time(1),'end','hour'):hours(1):dateshift(data.wave.time(end),'end','hour')]'; 
regwind_time = [dateshift(data.wind.time(1),'end','hour'):hours(1):dateshift(data.wind.time(end),'end','hour')]'; 
reginso_time = [dateshift(data.solar.time(1),'end','hour'):hours(1):dateshift(data.solar.time(end),'end','hour')]'; 
regcurr_time = [dateshift(data.curr.time(1),'end','hour'):hours(1):dateshift(data.curr.time(end),'end','hour')]'; 
regtemp_time = [dateshift(data.temperature.time(1),'end','hour'):hours(1):dateshift(data.temperature.time(end),'end','hour')]'; 

%Get wave data
opt.wave.wavepower_ra = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    (wave.Hs_ra)^2*(wave.Tp_ra); %[W], wave power at rated
%extract data
Hs = data.wave.Hs; %[m]
Tp = data.wave.Tp; %[s]
%opt.wave.L = atmo.g.*Tp.^2/(2*pi); %wavelength ts - this isn't used anymore
opt.wave.wavepower_ts = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    Hs.^2.*Tp./1000; %[kW/m] %timeseries of wavepower
opt.wave.P = data.wave.P./1000; %[kW/m] timeseries of intermediate depth wavepower
opt.wave.Pra = data.wave.Pra; %[kw/m]

%% Fill big gaps
num_d = size(data.curr.vmag); 

disp('Starting Wind')
bad_wind = sum(isnan(data.wind.U));
data.wind.U  = interp1(data.wind.time,data.wind.U,regwind_time,'nearest'); %shift to clean time series
data.wind.z0  = interp1(data.wind.time,data.wind.z0,regwind_time,'nearest'); %shift to clean time series
data.wind.rho = interp1(data.wind.time,data.wind.rho,regwind_time,'nearest');
if bad_wind > 24
    data.wind.U = fillmiss_phaseavg(data.wind.U, regwind_time);
    data.wind.z0 = fillmiss_phaseavg(data.wind.z0,regwind_time);
    data.wind.rho = fillmiss_phaseavg(data.wind.rho,regwind_time);
end

data.wind.U = fillmissing(data.wind.U,'nearest'); %fill short outages
data.wind.z0 = fillmissing(data.wind.z0,'nearest');
data.wind.rho = fillmissing(data.wind.rho,'nearest');

[data.wind.U,data.wind.time] = extendToLifetime(data.wind.U,datenum(regwind_time),uc.lifetime); %make time series data adequately long
[data.wind.z0,data.wind.time] = extendToLifetime(data.wind.z0,datenum(regwind_time),uc.lifetime);
[data.wind.rho,data.wind.time] = extendToLifetime(data.wind.rho,datenum(regwind_time),uc.lifetime);

disp('Starting Solar')
bad_solar = sum(isnan(data.solar.swso));
data.solar.swso  = interp1(data.solar.time,data.solar.swso,reginso_time,'nearest');
if bad_solar > 24
    data.solar.swso = fillmiss_phaseavg(data.solar.swso, reginso_time);
end
data.solar.swso = fillmissing(data.solar.swso,'nearest');
swso_neg_count = sum(data.solar.swso < 0);
if swso_neg_count > 0
    disp('removing negative solar data')
    data.solar.swso(data.solar.swso < 0) = 0; %remove any negative swso points
end
[data.solar.swso,data.solar.time] = extendToLifetime(data.solar.swso,datenum(reginso_time),uc.lifetime); %make time series data adequately long
 
disp('Starting Current')
for i = 1:min(num_d)
    bad_curr = sum(isnan(data.curr.vmag(:,i)));
    data.curr.speed(:,i) = interp1(data.curr.time,data.curr.vmag(:,i),regcurr_time,'nearest','extrap');
    if bad_curr > 24
        data.curr.speed(:,i) = fillmiss_phaseavg(data.curr.speed(:,i), regcurr_time);
    end
    data.curr.speed(:,i) = fillmissing(data.curr.speed(:,i),'nearest');
end
%extend current speed, and time to life
[data.curr.speed6(:,1),data.curr.time] = extendToLifetime(data.curr.speed(:,1),datenum(regcurr_time),uc.lifetime);
for i=2:num_d(2)
    [data.curr.speed6(:,i)] = extendToLifetime(data.curr.speed(:,i),datenum(regcurr_time),uc.lifetime);
end

disp('Starting Wave')
bad_wave = sum(isnan(Hs));
opt.wave.wavepower_ts = interp1(data.wave.time,opt.wave.wavepower_ts,regwave_time,'nearest');
opt.wave.Hs  = interp1(data.wave.time,Hs,regwave_time,'nearest');
opt.wave.Tp  = interp1(data.wave.time,Tp,regwave_time,'nearest');
%opt.wave.L  = interp1(data.wave.time,opt.wave.L,regwave_time,'nearest');
opt.wave.P = interp1(data.wave.time,opt.wave.P,regwave_time,'nearest');
if bad_wave > 24
    opt.wave.wavepower_ts = fillmiss_phaseavg(opt.wave.wavepower_ts, regwave_time);
    nanwave = isnan(opt.wave.wavepower_ts);
    opt.wave.Hs = fillmiss_phaseavg(opt.wave.Hs, regwave_time);
    opt.wave.Tp  = fillmiss_phaseavg(opt.wave.Tp, regwave_time);
    %opt.wave.L = fillmiss_phaseavg(opt.wave.L, regwave_time);
    opt.wave.P = fillmiss_phaseavg(opt.wave.P,regwave_time);
else
    nanwave = [];
end
opt.wave.wavepower_ts = fillmissing(opt.wave.wavepower_ts,'nearest');
opt.wave.Hs = fillmissing(opt.wave.Hs,'nearest');
opt.wave.Tp = fillmissing(opt.wave.Tp,'nearest');
%opt.wave.L = fillmissing(opt.wave.L,'nearest');
opt.wave.P = fillmissing(opt.wave.P,'nearest');

%extend wavepower, time, hs, tp and wavelength timeseries
[opt.wave.wavepower_ts,opt.wave.time] =  ...
    extendToLifetime(opt.wave.wavepower_ts,datenum(regwave_time),uc.lifetime);
[opt.wave.Hs] = extendToLifetime(opt.wave.Hs,datenum(regwave_time),uc.lifetime);
[opt.wave.Tp] = extendToLifetime(opt.wave.Tp,datenum(regwave_time),uc.lifetime);
%[opt.wave.L] = extendToLifetime(opt.wave.L,datenum(regwave_time),uc.lifetime);
[opt.wave.P] = extendToLifetime(opt.wave.P,datenum(regwave_time),uc.lifetime);

disp('Starting Temperature')
% temperature
bad_temp = sum(isnan(data.temperature.T2mw));
data.temperature.T2mw  = interp1(data.temperature.time,data.temperature.T2mw,regtemp_time,'nearest');
data.temperature.Tss  = interp1(data.temperature.time,data.temperature.Tss,regtemp_time,'nearest');
if bad_temp > 24
    data.temperature.T2mw = fillmiss_phaseavg(data.temperature.T2mw, regtemp_time);
    data.temperature.Tss = fillmiss_phaseavg(data.temperature.Tss, regtemp_time);
end
data.temperature.T2mw = fillmissing(data.temperature.T2mw,'nearest');
data.temperature.Tss = fillmissing(data.temperature.Tss,'nearest');
[data.temperature.T2mw,data.temperature.time] = extendToLifetime(data.temperature.T2mw,datenum(regtemp_time),uc.lifetime); %make time series data adequately long
[data.temperature.Tss,data.temperature.time] = extendToLifetime(data.temperature.Tss,datenum(regtemp_time),uc.lifetime); %make time series data adequately long
%% Make Clean time series
if opt.wave.time(1) == data.wind.time(1) && opt.wave.time(1) == data.solar.time(1) && opt.wave.time(1) == data.temperature.time(1)
    if opt.timeadj > 1
        opt.wave.time = opt.wave.time(opt.timeadj:end); %jump forward by the time adjustement
        data.wind.time = data.wind.time(opt.timeadj:end);
        data.solar.time = data.solar.time(opt.timeadj:end);
        data.temperature.time = data.temperature.time(opt.timeadj:end);

        [data.wind.U,data.wind.time] = extendToLifetime(data.wind.U(opt.timeadj:end),datenum(data.wind.time),uc.lifetime); %make time series data adequately long
        [data.wind.z0] = extendToLifetime(data.wind.z0(opt.timeadj:end),datenum(data.wind.time),uc.lifetime); %make time series data adequately long
        [data.wind.rho] = extendToLifetime(data.wind.rho(opt.timeadj:end),datenum(data.wind.time),uc.lifetime); %make time series data adequately long

        [data.solar.swso,data.solar.time] = extendToLifetime(data.solar.swso(opt.timeadj:end),datenum(data.solar.time),uc.lifetime); %make time series data adequately long
        [data.temperature.T2mw,data.temperature.time] = extendToLifetime(data.temperature.T2mw(opt.timeadj:end),datenum(data.temperature.time),uc.lifetime); %make time series data adequately long
        [data.temperature.Tss,data.temperature.time] = extendToLifetime(data.temperature.Tss(opt.timeadj:end),datenum(data.temperature.time),uc.lifetime); %make time series data adequately long

        [opt.wave.wavepower_ts,opt.wave.time] = extendToLifetime(opt.wave.wavepower_ts(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
        [opt.wave.Hs] = extendToLifetime(opt.wave.Hs(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
        [opt.wave.Tp] = extendToLifetime(opt.wave.Tp(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
       % [opt.wave.L] = extendToLifetime(opt.wave.L(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
        [opt.wave.P] = extendToLifetime(opt.wave.P(opt.timeadj:end),datenum(opt.wave.time),uc.lifetime);
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


%%Plot processed data
if opt.pltdebug  %diagnostic plot of input data for Task 2 locations (aligned data)
    clear ax
    figure(2)
    tf = tiledlayout(5,1);
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
    plot(datetime(data.met.time,'convertfrom','datenum'),opt.wave.wavepower_ts,'linewidth',1.5,'color',col(4,:),'DisplayName','Deep')
    hold on
    plot(datetime(data.met.time(nanwave),'convertfrom','datenum'),opt.wave.wavepower_ts(nanwave),'kx')
    plot(datetime(data.met.time,'convertfrom','datenum'),opt.wave.P,':','linewidth',1.5,'color',col(1,:),'DisplayName','Intermediate')
    ylabel('[kW/m]')
    title('Wave Power')
    xlabel('Time')
    legend

    ax(5) = nexttile; %temps
    plot(datetime(data.met.time,'convertfrom','datenum'),data.temperature.T2mw,'linewidth',1.5,'color',col(6,:),'DisplayName','T2M')
    hold on
    plot(datetime(data.met.time,'convertfrom','datenum'),data.temperature.Tss,'linewidth',1.5,'color',col(7,:),'DisplayName','SST')
    legend
    ylabel('[C]')
    title('Temperature')
    xlabel('Time')

    %linkaxes(ax,'x')
    %xlim([min(data.met.time),max(data.met.time)])
end

if opt.pltdebug  %diagnostic plot of input data for Task 2 locations (aligned data)
    figure
    plot(data.wind.rho)
    hold on
    yline(atmo.rho_a_c)
    ylabel('Density [kg/m3]')
    title(data.title)
    
    figure
    plot(data.wind.z0.*1000)
    hold on
    yline(atmo.zo_c)
    ylabel('Surface Roughness [mm]')
    title(data.title)
end
%% apply ice model
[Iceslow,Icefast] = iceModel(data.temperature.T2mw,data.met.wind_spd,data.met.time);
if strcmp(opt.ice,'fast')
    opt.ice_ts = 1 - Icefast; %need to flip the ice model so it is zero when ice occurs
elseif strcmp(opt.ice,'slow')
    opt.ice_ts = 1 - Iceslow; %need to flip the ice model so it is zero when ice occurs
else
    opt.ice_ts = ones(size(Iceslow)); %all ones makes no ice adjustment
end

data.swso = data.swso.*opt.ice_ts;
data.met.wind_spd = data.met.wind_spd.*opt.ice_ts;

%% winter cleaning - obselete
% if inso.cleanstrat == 3 || inso.cleanstrat == 4 %winter cleaning
%     %winter cleaning (if applicable)
%     if data.lat < 0 %southern hemisphere
%         wint_clean_mo = 5; %may
%     else
%         wint_clean_mo = 11; %november
%     end
%     dv = datevec(data.solar.time);
%     data.wint_clean_ind = find(dv(:,2) == wint_clean_mo & dv(:,3) == 1 & dv(:,4) == 0);
%     if inso.cleanstrat == 4 %every other winter
%         data.wint_clean_ind = data.wint_clean_ind(2:2:end);
%     end
% end

%% Get current power
%[opt.curr.F_curr] = Current_Power(cturb,atmo); %old method of calcuating
%current power - new version is in SimHybrid

%PrepWave
if ~wave.deep
    opt.wave.wavepower_ts = opt.wave.P; %use intermediate depth calculation
    opt.wave.wavepower_ra = opt.wave.Pra; %use intermediate depth calculation
    disp('Using intermediate water wave power')
else
    disp('Using deepwater wave power')
end
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