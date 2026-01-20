function [L, wave_regime] = calculate_wave_props(d,T,g)

%From Brian: pulled out of LoadSites 12/31/2025

%calculate wavelength using dispersion relation
%note: this is using a regular wave equation and is inexact for irregular seas

%initialize return values
L = zeros(size(T));
wave_regime = zeros(size(T));

%loop through each sea state
for i = 1:length(T)

    %if wave data available
    if ~isnan(T(i))
        %L_approx = (g*T(i).^2/(2*pi)).*(tanh(4*pi^2./T(i).^2.*d/g)).^0.5; %analytical approximation for intermediate depth - Brian's
        L_approx = (g*T(i).^2/(2*pi)); %deep water approx
        problem.objective = @(x) x-g*T(i).^2./(2*pi).*tanh(2*pi*d/x);   %wave dispersion relation as problem objective
        problem.x0 = L_approx;     %approximate relation (starting guess)
        problem.solver = 'fzero';
        problem.options = optimset(@fzero);
        L(i) = fzero(problem);

        if d/L(i) > 0.5
            wave_regime(i) = 1;    %deep water wave
        elseif d/L(i) < 1/20
            wave_regime(i) = 3;    %shallow water wave
        else
            wave_regime(i) = 2;    %intermediate water depth
        end
    else
        L(i) = NaN;
        wave_regime(i) = NaN;
    end

end

end