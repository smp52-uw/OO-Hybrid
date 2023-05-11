function [F_curr, P,kW_curr,c_vel] = Current_Power(cturb,atmo)
%Generating a power matrix for power from a current turbine

a = 200;
b = 200;

kW_curr = linspace(0,8,a); %range of rated power
c_vel = linspace(0,40,b); %range of current speeds
P = zeros(length(c_vel),length(kW_curr));
for i=1:length(c_vel)
    if c_vel(i) < cturb.uci %below cut out
        P(i,:) = 0; %[W]
    elseif cturb.uci < c_vel(i) && c_vel(i) <= cturb.ura %below rated
        P(i,:) = kW_curr.*1000.*c_vel(i).^3./cturb.ura.^3; %[W]
    elseif cturb.ura < c_vel(i) && c_vel(i) <= cturb.uco %above rated
        P(i,:) = kW_curr.*1000; %[W]
    else %above cut out
        P(i,:) = 0; %[W]
    end
end

%form a mesh
[kW_curr c_vel] = meshgrid(kW_curr, c_vel); %range of current speeds
kW_curr = reshape(kW_curr, [a*b,1]);
c_vel = reshape(c_vel, [a*b,1]);
P = reshape(P,[a*b,1]);
F_curr = scatteredInterpolant(kW_curr,c_vel,P);
%test = F_curr(kW_curr(100),c_vel(100))
Pcurr.kW_curr = kW_curr;
Pcurr.c_vel = c_vel;
Pcurr.P = P;
Pcurr.F_curr = F_curr;
% out_file = 'current_Pmatrix.mat'
% save(out_file,'Pcurr')
end