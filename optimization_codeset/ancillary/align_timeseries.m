function [data_out] = align_timeseries(reg_time,t_in,a_in)
    
    reg_time = datetime(reg_time,'ConvertFrom','datenum');
    t_in = datetime(t_in,'ConvertFrom','datenum');
    %a_in = data.met.wind_spd
    reg_vec = datevec(reg_time(1)); %vector of the first time in the clean time series
    ind = find(reg_time == t_in(1),1,'first'); %only want the first index
    if ind == 1 %if the first data point is the first time step in the alighned TS then no adjustment to data needed
        data_out = a_in;
    else
        time_vec = datevec(t_in);
        start_i = find(time_vec(:,2) == reg_vec(2) & time_vec(:,3) == reg_vec(3) & time_vec(:,4) == reg_vec(4)); %index to start replicating
        %start_i = find(time_vec(2) == reg_vec(2) && time_vec(3) == reg_vec(3) && time_vec(4) == reg_vec(4)); %index to start replicating
        %take the first instance of the same day/hour and duplicate from there
        data_out = [a_in(start_i(1):start_i(1)+ind); a_in];   
        data_out = data_out(1:length(reg_time));
    end
end