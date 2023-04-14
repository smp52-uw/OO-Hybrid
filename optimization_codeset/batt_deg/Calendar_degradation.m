function [L, d_cal]= Calendar_degradation(SoC, T, t)
%% calender aging
%SM Palmer: pulled the constants and equations from the linear_degredation
%function that is a part of the rainflow battery degredation model to be
%used for battery only experiencing lifetime based degredation

k_soc = 1.039; % SoC stress model coefficient, NMC and LMO
k_t = 0.0693; % Tfact from Millner article
T_ref = 25; % reference temperature in Celsius
k_soc_cal = k_soc;
k_cal = 4.1375e-10; % time aging stress model coefficient, NMC and LMO
SoC_ref = .6; % reference cycle average SoC, NMC and LMO

SoC_avg = mean(SoC);
T_avg = mean(T);

d_cal = k_cal.*exp(k_soc_cal*(SoC_avg-SoC_ref)) .* t ...
    .* exp(k_t * (T_avg-T_ref) .* (273 + T_ref)./(273 + T_avg));


L = Nonlinear_degradation( d_cal ); %I think this is needed to get the fast degredation at the beginning of life
end