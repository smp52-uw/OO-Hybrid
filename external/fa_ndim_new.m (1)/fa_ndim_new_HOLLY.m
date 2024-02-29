% ---------------------------------------------------------------------- %
% The Firefly Algorithm (FA) for unconstrained function optimization     %
% by Xin-She Yang (Cambridge University) @2008-2009                      %
% Programming dates: 2008-2009, then revised and updated in Oct 2010     %
% ---------------------------------------------------------------------- %

% References -- citation details: -------------------------------------- %
% (1) Xin-She Yang, Nature-Inspired Metaheuristic Algorithms,            %
%     Luniver Press, First Edition, (2008).                              %
% (2) Xin-She Yang, Firefly Algorithm, Stochastic Test Functions and     %
%     Design Optimisation, Int. Journal of Bio-Inspired Computation,     %
%     vol. 2, no. 2, 78-84 (2010).                                       %
% ---------------------------------------------------------------------- %

%Changes made: Sarah May Palmer 9/11/2023



% -------- Start the Firefly Algorithm (FA) main loop ------------------ % 
function [pop,objfunval] = fa_ndim_new_HOLLY 

%%Call Hybrid Opt functions needed
%tTot = tic;
optInputs

data = load(loc,loc);
data = data.(loc);
%load current data
curr = load(cloc);
data.curr = curr; %add current data to data structure

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
[data, opt] = prepHybrid(data,opt,uc(c),wave,atmo,inso,cturb);
%Hybrid load case call
[uc(c).loaddata, loadseries] = GenerateLoadCases_v4(data); %updated to even hour loads
uc(c).draw = loadseries.L(uc(c).loadcase,:);
uc = uc(c);
opt.fmin = false;

%n=opt.ffa.pop;          % Population size (number of fireflies)
n = 25;
%alpha=1.0;              % Randomness strength 0--1 (highly random)
alpha = 0.2;
%beta0=1.0;              % Attractiveness constant
beta0 = 2;
% gamma=0.01;             % Absorption coefficient
gamma = 1;
theta=0.97;             % Randomness reduction factor theta=10^(-5/tMax) 
d=opt.pd;               % Number of dimensions
%tMax=opt.ffa.max;       % Maximum number of iterations
tMax=10;       % Maximum number of iterations
%Lb=-10*ones(1,d);       % Lower bounds/limits
Lb = [0 0 0 0 0 1];
%Ub=10*ones(1,d);        % Upper bounds/limits
Ub = [8 8 8 8 8 500];
dmax = norm(Ub-Lb);
% Generating the initial locations of n fireflies
ns = nan(n,d);
Lightn = Inf(1,n);
pop{1,tMax} = [];
objfunval{1,tMax} = [];
for i=1:n
   ns(i,:)=Lb+(Ub-Lb).*rand(1,d);         % Randomization
   %Lightn(i)=ObjFun(ns(i,:),opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);               % Evaluate objectives
end
%parfor (i = 1:n,opt.bf.maxworkers)
for i=1:n
    x = ns(i,:);
    Kd = x(1);
    Ki = x(2);
    Kwi = x(3);
    Kwa = x(4);
    Kc = x(5);
    S = x(6);
    [C_temp(i),S_temp(i)] = simHybrid(Kd, Ki, Kwi, Kwa, Kc, S,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);

    if S_temp(i)< 0.99
       C_temp(i) = inf;
    end
end
for i=1:n
   Lightn(i)= C_temp(i);
end

%%%%%%%%%%%%%%%%% Start the iterations (main loop) %%%%%%%%%%%%%%%%%%%%%%%
for k=1:tMax
    k
    alpha=alpha*theta;     % Reduce alpha by a factor theta
    %scale=abs(Ub-Lb);      % Scale of the optimization problem
    scale=0.05.*abs(Ub-Lb);      % Scale of the optimization problem
    % Two loops over all the n fireflies
    for i=1:n
        for j=1:n
          % Evaluate the objective values of current solutions
          Lightn(i)=ObjFun(ns(i,:),opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);           % Call the objective
          if ~isreal(Lightn(i))
              disp("Warning - IMAG VALUE")
          end
          % Update moves
          if Lightn(i)>=Lightn(j)           % Brighter/more attractive
          %if Lightn(i)>Lightn(j)           % I want the best point each iteration to stick around
             r=sqrt(sum((ns(i,:)-ns(j,:)).^2))/dmax;
             beta=beta0*exp(-gamma*r.^2);    % Attractiveness
             %steps=alpha.*(rand(1,d)-0.5).*scale;
             steps=alpha.*(rand(1,d)).*scale; %attempting to make things work???
          % The FA equation for updating position vectors
             ns(i,:)=ns(i,:)+beta.*rand(1,d).*(ns(j,:)-ns(i,:))+steps;

             nsol_tmp=ns(i,:);
             % Apply the lower bound
             I=nsol_tmp<Lb;  nsol_tmp(I)=Lb(I);
             % Apply the upper bounds
             J=nsol_tmp>Ub;  nsol_tmp(J)=Ub(J);
             ns(i,:)=nsol_tmp;
%              P_temp(j,:) = pop(i).Position ...
%                                 + beta*rand(VarSize).*(pop(j).Position-pop(i).Position) ...
%                                 + alpha*e;
             if ~isreal(ns(i,:))
                 disp("error - imag point")
             end
    
          end
       end % end for j
    end % end for i

    % Check if the new solutions/locations are within limits/bounds
    %[ns,Lightn]=findlimits(n,ns(i,:),Lb,Ub);
%     for i=1:n
%       nsol_tmp=ns(i,:);
%       % Apply the lower bound
%       I=nsol_tmp<Lb;  nsol_tmp(I)=Lb(I);
%       % Apply the upper bounds
%       J=nsol_tmp>Ub;  nsol_tmp(J)=Ub(J);
%       % Update this new move
%       ns(i,:)=nsol_tmp;
%       Lightn(i)=ObjFun(ns(i,:),opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb); 
%     end
    %% Rank fireflies by their light intensity/objectives
    [Lightn,Index]=sort(Lightn,"ascend");
    nsol_tmp=ns;
    for i=1:n
     ns(i,:)=nsol_tmp(Index(i),:);
    end
    %% Find the current best solution and display outputs
    fbest=Lightn(1)
    nbest=ns(1,:)
    pop{k} = ns;
    objfunval{k} = Lightn;
end % End of the main FA loop (up to tMax) 
end
% Make sure that new fireflies are within the bounds/limits
function [ns,Lightn]=findlimits(n,ns,Lb,Ub)
%     for i=1:n
%       nsol_tmp=ns(i,:);
%       % Apply the lower bound
%       I=nsol_tmp<Lb;  nsol_tmp(I)=Lb(I);
%       % Apply the upper bounds
%       J=nsol_tmp>Ub;  nsol_tmp(J)=Ub(J);
%       % Update this new move
%       ns(i,:)=nsol_tmp;
%       Lightn(i)=ObjFun(ns(i,:),opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb); 
%     end
end

%% Define the objective function or cost function
function z=ObjFun(x,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb)
    % The modified sphere function: z=sum_{i=1}^D (x_i-1)^2
    %z=sum((x-1).^2); % The global minimum fmin=0 at (1,1,...,1)
    Kd = x(1);
    Ki = x(2);
    Kwi = x(3);
    Kwa = x(4);
    Kc = x(5);
    S = x(6);
    [C_temp,S_temp] = simHybrid(Kd, Ki, Kwi, Kwa, Kc, S,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);
    
    z = C_temp;
    if S_temp < 0.99
       z = inf;
    end
end
% -----------------------------------------------------------------------%