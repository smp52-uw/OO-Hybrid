%Troubleshoot Wave Data
optInputs
%optInputs %load inputs
data = load(loc,loc);
data = data.(loc);
%load current data
curr = load(cloc);
data.curr = curr; %add current data to data structure

%Get wave data
opt.wave.wavepower_ra = (1/(16*4*pi))*atmo.rho_w*atmo.g^2* ...
    (wave.Hs_ra)^2*(wave.Tp_ra); %[W], wave power at rated
%extract data
Hs = data.wave.significant_wave_height; %[m]
input_date = data.wave.time;
[data.wave.time, ind_wa] = unique(data.wave.time);
converted_time_wave = datetime(data.wave.time,'ConvertFrom','datenum');  
%converted_time_wave.TimeZone = "UTC"
regwave_time = [dateshift(converted_time_wave(1),'end','hour'):hours(1):dateshift(converted_time_wave(end),'end','hour')]'; 
%regwave.TimeZone = "UTC"
%PrepWave
%Interpolate wave data to clean time series
opt.wave.Hs  = interp1(converted_time_wave,Hs(ind_wa),regwave_time,'linear');
%opt.wave.Hs  = interp1(converted_time_wave,Hs,regwave_time,'linear');

plot(datetime(input_date,'ConvertFrom','datenum'),Hs,'.')
hold on
plot(converted_time_wave, Hs(ind_wa),'.')
plot(regwave_time, opt.wave.Hs,'.')
legend('Input','Data Ordered','Interp','Location','bestoutside')