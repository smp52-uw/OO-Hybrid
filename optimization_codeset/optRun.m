function [output,opt] = optRun(opt,data,atmo,batt,econ,uc,bc,inso, ...
    turb,cturb, wave,dies)

%curve-fit device scatters, find polyvals
opt.p_dev.t = calcDeviceVal('turbine',[],econ.wind_n);
opt.p_dev.d_cost = calcDeviceVal('dieselcost',[],econ.diescost_n);
opt.p_dev.d_mass = calcDeviceVal('dieselmass',[],econ.diesmass_n);
opt.p_dev.d_size = calcDeviceVal('dieselsize',[],econ.diessize_n);
opt.p_dev.d_burn = calcDeviceVal('dieselburn',[],econ.diesburn_n);
opt.p_dev.d_vol = calcDeviceVal('dieselvol',[],econ.diesvol_n);
opt.p_dev.b_size = calcDeviceVal('lfp_vol',[],econ.battsize_n);
[opt.p_dev.b,~,opt.p_dev.kWhmax] = calcDeviceVal('agm',[],econ.batt_n);

%HYBRID Prep Function Calls
[data, opt] = prepHybrid(data,opt,uc,wave,atmo,inso,cturb);
%Hybrid load case call
[uc.loaddata, loadseries] = GenerateLoadCases_v4(data); %updated to even hour loads
uc.draw = loadseries.L(uc.loadcase,:);
if opt.pltdebug
    figure(3)
    hold on
    plot(uc.draw(1:8760),'linewidth',1.5)
    yline(mean(uc.draw),'--','linewidth',1.5)
    ylabel('[W]')
    xlabel('Index')
    title(strcat("1 Year, Load Case: ",string(uc.loadcase)))
end
%HYBRID Opt Function call
[output,opt] = optHybrid(opt,data,atmo,batt,econ,uc,bc,dies,inso,turb,cturb, wave);

% Single gen prep and opt function calls
% [data,econ] = prepInso(data,inso,econ,uc);
% data = prepWind(data,uc);
% data = prepDies(data,econ,uc);
% opt = prepWave(data,opt,wave,atmo,uc);

% if pm == 1 %WIND
%     data = prepWind(data,uc);
%     [output,opt] = optWind(opt,data,atmo,batt,econ,uc,bc,turb);
% elseif pm == 2 %SOLAR
%     [data,econ] = prepInso(data,inso,econ,uc);
%     [output,opt] = optInso(opt,data,atmo,batt,econ,uc,bc,inso);
% elseif pm == 3 %WAVE
%     opt = prepWave(data,opt,wave,atmo,uc);
%     if opt.V == 1 %optimization version 1: nelder-mead
%         if opt.nm.many
%             opt.C = length(opt.nm.bgd_array);
%             compare(opt.C) = struct();
%             costcompare = zeros(1,opt.C);
%             for i = 1:opt.C
%                 opt.c = i;
%                 opt.nm.battgriddur = opt.nm.bgd_array(i);
%                 [compare(i).output,compare(i).opt] = ...
%                     optWave_nm(opt,data,atmo,batt,econ,uc,bc,wave);
%                 costcompare(i) = compare(i).output.min.cost;
%             end
%             [~,min_ind] = min(costcompare(:));
%             output = compare(min_ind).output;
%             opt = compare(min_ind).opt;
%             opt.nm.battgriddur = opt.nm.bgd_array(min_ind);
%         else
%             [output,opt] = optWave_nm(opt,data,atmo,batt,econ,uc,bc,wave);
%         end
%     elseif opt.V == 2 %optimization version 2: brute force
%         [output,opt] = optWave(opt,data,atmo,batt,econ,uc,bc,wave);
%     end
%     results.width = output.min.width;
%     results.cw_avg = output.min.cw_avg;
%     results.cwr_avg = output.min.cwr_avg;
% elseif pm == 4 %DIESEL
%     data = prepDies(data,econ,uc);
%     [output,opt] = optDies(opt,data,atmo,batt,econ,uc,bc,dies);
% end

%print min values
results.width = output.min.width;
results.cw_avg = output.min.cw_avg;
results.cwr_avg = output.min.cwr_avg;
results.kWd = output.min.kWd;
results.kWi = output.min.kWi;
results.kWwi = output.min.kWwi;
results.kWwa = output.min.kWwa;
results.kWc = output.min.kWc;
results.Smax = output.min.Smax;
results.cost = output.min.cost;
results.CapEx = output.min.CapEx;
results.OpEx = output.min.OpEx;
results.nvi = output.min.nvi;
results.CFd = output.min.CFd;
results.CFi = output.min.CFi;
results.CFwi = output.min.CFwi;
results.CFwa = output.min.CFwa;
results.CFc = output.min.CFc;
%results.batt_L_max = max(output.min.batt_L);
results.batt_L1_max = max(output.min.batt_L1);
results.batt_L2_max = max(output.min.batt_L2);
results.nfr = output.min.nfr;
results.noc = output.min.noc;
results

end

