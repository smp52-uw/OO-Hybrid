%Brian Polagye
%July 20, 2022
%Edits: Sarah Palmer 2/6/2023

%Description: Generate time-resolved load cases

function [load_case,load_series] = GenerateLoadCases_v3(data)

%% Setup

% %output file
% out_file = 'load_series.mat';

%simulation parameters
dt = minutes(60);   %load case time resolution
data_time = datetime(data.met.time,'ConvertFrom','datenum');
T = length(data_time);
%T = years(5);

%case 1: high-capacity UUV
load_case(1).name = 'High-capacity UUV';
load_case(1).L_constant = 5;  %constant load [W]
load_case(1).L_intermittent = (384+90)/0.75; %intermittent load [W]
%384 W average onboard charger load + 90 W hotel load / 75% wireless transfer efficiency
load_case(1).interval = days(7);   %intermittent load repetition rate
load_case(1).duration = hours(16.67);  %intermittent load duration

%case 2:high-frequency recharge UUV
load_case(2).name = 'High-frequency Recharge UUV';
load_case(2).L_constant = 5;  %constant load [W]
load_case(2).L_intermittent = (192+90)/0.75; %intermittent load [W]
%192 W average onboard charger load + 90 W hotel load / 75% wireless transfer efficiency
load_case(2).interval = days(1);   %intermittent load repetition rate
load_case(2).duration = hours(5.21);  %intermittent load duration

%case 3: AMP + high-capacity UUV
load_case(3).name = 'Next-gen Oceanographic + UUV';
load_case(3).L_constant = 322;  %constant load [W]  %actual load will be controlled using Duration algorithm with recharge cycle
%Maximum power status - controllable based on system command
load_case(3).L_intermittent = load_case(1).L_intermittent; %intermittent load [W]   %add in high-capacity UUV
load_case(3).interval = load_case(1).interval;      %intermittent load repetition rate
load_case(3).duration = load_case(1).duration;      %intermittent load duration

%case 4: high-frequency radar
load_case(4).name = 'High-frequency Radar';
load_case(4).L_constant = 155.1;  %constant load [W]
load_case(4).L_intermittent = 0; %intermittent load [W]

%% Initialization
t = [0:dt:dt*(T-1)]';  %initialize time series

L = zeros(length(load_case),length(t)); %initialize load cases
L_status = ones(length(load_case),length(t));  %initialize load status

%% Generate load cases

%loop through all load cases
for i = 1:length(load_case)
    
    L_intermittent_timer = load_case(i).interval;       %start intermittent time
    L_charge_timer = 0;                                 %initialize charge timer

    %loop through all time
    for j = 1:length(t)

        %always incur hotel load
        L(i,j) = load_case(i).L_constant;

        %start intermittent load
        if L_intermittent_timer < 0
            L_charge_timer = load_case(i).duration;   %set charge duration
            L_status(i,j) = 2;   %note intermittent charging status
            L_intermittent_timer = load_case(i).interval;    %reset intermittent load interval
        
        %intermittent load active
        elseif L_charge_timer > 0
            L_status(i,j) = 2;

        %intermittent load time charged
        elseif L_charge_timer <= 0
            L_status(i,j) = 1;   %note return to hotel load
        end

        %incur intermittent load, if active
        if L_status(i,j) == 2
            L(i,j) = L(i,j) + load_case(i).L_intermittent;
            L_charge_timer = L_charge_timer - dt;
            
        %otherwise count down to next intermittent load
        else
            L_intermittent_timer = L_intermittent_timer - dt;
        end
        
    end

end
load_series.L = L;
load_series.L_status = L_status;
load_series.t = t;
end
% 
% %% Visualize load cases
% 
% figure(1)
% tiledlayout(length(load_case),1)
% 
% for i = 1:length(load_case)
%     nexttile
%     plot(days(t),L(i,:),'-')
%     hold on
%     yline(mean(L(i,:)),'--','linewidth',2)
% 
%     title(load_case(i).name,'fontweight','b')
%     ylabel('L [W]','fontweight','bold')
%     grid on
%     if i == length(load_case)
%         xlabel('Time [days]','fontweight','bold')
%     end
% 
%     if i == 1
%         legend('L(t)','L_{avg}','location','northwest')
%     end
% 
% end
% 
% %% Save load cases
% 
% save(out_file,'load_case','L','L_status','t')