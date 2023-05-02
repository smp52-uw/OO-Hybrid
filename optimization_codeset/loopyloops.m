function [Kd, Ki, Kwi, Kwa, S] = loopyloops(S_t, S_old, Kwa_t, Kwa_old, Kwi_t, Kwi_old, Ki_t, Ki_old, Kd_t, Kd_old,opt)

%set loop constants
i = 1; %index of old grid style array
dx = round(mode(diff(Kd_t)),4); %for this to work the dx must be constant
ds = round(mode(diff(S_t)),4);
Kd = [];
Ki = [];
Kwi = [];
Kwa = [];
S = [];

disp('Starting Massive Loop...')
%Kd,Ki,Kwi,Kwa,S
parfor (a=1:length(S_t),opt.bf.maxworkers)
%for a = 1:length(S_t)
    %disp('loop 1')
    for b = 1:length(Kwa_t)
        %disp('loop 2')
        for c = 1:length(Kwi_t)
            %disp('loop 3')
            for d = 1:length(Ki_t)
                %disp('loop 4')
                for e = 1: length(Kd_t)
                    %disp('loop 5')
                    for i = 1 : length(Kd_old)
                    %RUN_LOOP = 1;
                    %parfor (i = 1 : length(Kd_old),opt.bf.maxworkers)
                        %disp('loop 6')
                        d_use = 0; i_use = 0; wi_use = 0; wa_use = 0; s_use = 0;
                        if round(abs(Kd_t(e) - Kd_old(i)),4) == dx || round(abs(Kd_t(e) - Kd_old(i)),4) == 0
                            d_use = 1;
                        end
                        if round(abs(Ki_t(d) - Ki_old(i)),4) == dx || round(abs(Ki_t(d) - Ki_old(i)),4) == 0
                            i_use = 1;
                        end
                        if round(abs(Kwi_t(c) - Kwi_old(i)),4) == dx || round(abs(Kwi_t(c) - Kwi_old(i)),4) == 0
                            wi_use = 1;
                        end
                        if round(abs(Kwa_t(b) - Kwa_old(i)),4) == dx || round(abs(Kwa_t(b) - Kwa_old(i)),4) == 0
                            wa_use = 1;
                        end
                        if round(abs(S_t(a) - S_old(i)),4) == ds || round(abs(S_t(a) - S_old(i)),4) == 0
                            s_use = 1;
                        end
                        if d_use == 1 && i_use == 1 && wi_use == 1 && wa_use == 1 && s_use == 1
                            %append point to new grid style array
                            Kd = [Kd, Kd_t(e)];
                            Ki = [Ki, Ki_t(d)];
                            Kwi = [Kwi, Kwi_t(c)];
                            Kwa = [Kwa, Kwa_t(b)];
                            S = [S, S_t(a)];
                            %RUN_LOOP = 0;
                            break %skip out of 6th loop
                        end
                    end
                end
            end
        end
    end
end


disp('End Massive Loop...')
end