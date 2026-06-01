%plot wave power matrix with and without the wave hotel load
close all;
clear;

optInputs
data = load(loc,'data');
data = data.('data');

[data, opt] = prepHybrid(data,opt,uc(c),wave,atmo,inso,cturb);
house = [0.10 0.05 0]; %percent of wave power (standard is 0.1)
kW_wave = linspace(0,8,4);

figure;
t = tiledlayout(3,1);
time = datenum(data.met.time);
for h = 1:3
    wave.house = house(h);
    ax(h) = nexttile;
    hold on
    for k = 1:length(kW_wave)
        wavepower = opt.wave.wavepower_ts; %wavepower timeseries
        if wave.method == 1 %divide by B methodology - OUTDATED    
            disp('ERROR- method set to old divide by B method')
        elseif wave.method == 2 %3d interpolation methodology
            %extract data
            Hs = opt.wave.Hs; %Hs timeseries
            Tp = opt.wave.Tp; %Tp timeseries
            %find width through rated power conditions
            width = interp1(opt.wave.B_func(2,:),opt.wave.B_func(1,:),kW_wave(k)); %[m], B
            if kW_wave(k) == 0 %Set width = 0 for no wave gen (in case the curve cross (0,0))
                width = 0;
            end
            cw = width.*opt.wave.F(Tp,Hs,width*ones(length(Tp),1)); %[m] cw ts
        end
        
        %compute wave power timeseries
        Pwave{k} = wave.eta_ct*cw.*wavepower - kW_wave(k)*wave.house; %[kW] 
        Pwave{k}(Pwave{k}<0) = 0; %no negative power
        Pwave{k}(Pwave{k}>kW_wave(k)) = kW_wave(k); %no larger than rated power
        Pwave{k} = Pwave{k}*1000; %convert to watts
        if kW_wave(k) < 0.2144
            Pwave{k} = zeros(1,length(time));
            %disp('Zero wave power due to WECSIM min')
        end
        
        plot(Pwave{k})
    end
    grid on
    xlabel('Index')
    ylabel('Power [W]')
    title(strcat("House Percent: ",string(house(h))))
end

linkaxes(ax,'x')