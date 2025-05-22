function a_out = fillmiss_phaseavg(a_in, t_in)

%find the beginning of each year
t_start = find(hour(t_in) == hour(t_in(1)) & day(t_in) == day(t_in(1)) & month(t_in) == month(t_in(1))); 

if sum(isnan(a_in)) == length(a_in)
    disp('All nan')
    a_out = a_in;
elseif length(t_start) == 1
    disp('Outage occurs during a leap year')
    a_out = a_in;
else
    ly = find(day(t_in) == 29 & month(t_in) == 2); %leap day
    t_nly = t_in;
    t_nly(ly) = [];
    a_nly = a_in;
    a_nly(ly) = [];

    t_start = find(hour(t_nly) == hour(t_nly(1)) & day(t_nly) == day(t_nly(1)) & month(t_nly) == month(t_nly(1))); 
    
    %pack into a matrix for averaging
    a_mat = [];
    for i = 1:length(t_start)-1
        a_mat(:,i) = a_nly(t_start(i):t_start(i+1)-1);
    end
    
    a_avg = mean(a_mat,2,'omitnan'); %phase avg values
    t_avg = t_nly(1:t_start(2)-1); %1 year of time
    
    ind_bad = find(isnan(a_in));
    
    a_out = a_in;
    for j = 1:length(ind_bad)
        timeind = find(hour(t_avg) == hour(t_in(ind_bad(j))) & day(t_avg) == day(t_in(ind_bad(j))) & month(t_avg) == month(t_in(ind_bad(j)))); %index corresponding to the missing value
        if ~isempty(timeind) 
            a_out(ind_bad(j)) = a_avg(timeind); %replace nan with phase avg
        else
            a_out(ind_bad(j)) = nan;
        end
    end
    disp(strcat("Number of nans left ", string(sum(isnan(a_out)))))
end
end