%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YOEA112
% Project Title: Implementation of Firefly Algorithm (FA) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
% Changes: by SM Palmer 7-6-23 to fit hybrid optimization

% clc;
% clear;
% close all;

%%Call Hybrid Opt functions needed
tTot = tic;
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
opt.fmin = false;
%% Problem Definition

%CostFunction=@(x) Rosenbrock(x);        % Cost Function

%this needs to be connected to opt.pd
nVar=6;                 % Number of Decision Variables - number of dimensions

VarSize=[1 nVar];       % Decision Variables Matrix Size

%Min and Max of gird - needs to connect to opt.bf.M or .N
VarMin= [0 0 0 0 0 1];             % Decision Variables Lower Bound
VarMax= [8 8 8 8 8 500];             % Decision Variables Upper Bound

%adjustments for 5D opt
if opt.pd == 5
    if opt.pm == 5 %no current
        VarMax(5) = 0;
    elseif opt.pm == 4 %no diesel
        VarMax(1) = 0;
    end
end

%% Firefly Algorithm Parameters

%MaxIt=1000;         % Maximum Number of Iterations
MaxIt=10;         % Maximum Number of Iterations

%nPop=25;            % Number of Fireflies (Swarm Size)
nPop=25;            % Number of Fireflies (Swarm Size)

gamma=1;            % Light Absorption Coefficient

beta0=2;            % Attraction Coefficient Base Value

alpha=0.2;          % Mutation Coefficient

alpha_damp=0.98;    % Mutation Coefficient Damping Ratio

%added the .* to do element wise multiplication
delta=0.05.*(VarMax-VarMin);     % Uniform Mutation Range

m=2;

if isscalar(VarMin) && isscalar(VarMax)
    dmax = (VarMax-VarMin)*sqrt(nVar);
else
    dmax = norm(VarMax-VarMin);
end

%% Initialization

% Empty Firefly Structure
firefly.Position=[];
firefly.Cost=[];
firefly.Surv=[];

% Initialize Population Array
pop=repmat(firefly,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;
BestSol.Surv = [];
BestSol.Position = [];
b = 1; %counter for saving multiple best points

% Create Initial Fireflies
for i=1:nPop %(i = 1:sim_length,opt.bf.maxworkers)
   pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
end

parfor (i = 1:nPop,opt.bf.maxworkers)
   Kd = pop(i).Position(1);
   Ki = pop(i).Position(2);
   Kwi = pop(i).Position(3);
   Kwa = pop(i).Position(4);
   Kc = pop(i).Position(5);
   S = pop(i).Position(6);
   %pop(i).Cost=CostFunction(pop(i).Position);
   [C_temp(i),S_temp(i)] = ...
                simHybrid(Kd, Ki, Kwi, Kwa, Kc, S,opt,data,atmo,batt,econ,uc(c),bc,dies,inso,wave,turb,cturb);
end

for i=1:nPop
   pop(i).Cost = C_temp(i);
   pop(i).Surv = S_temp(i);
   if pop(i).Surv < 0.99
       pop(i).Cost = inf;
   end

   if pop(i).Cost<=BestSol.Cost
       BestSol.Cost = pop(i).Cost;
       BestSol.Surv{1,b}=pop(i).Surv;
       BestSol.Position{1,b} = pop(i).Position;
       b = b+1;
   end
end
%Save each population
pop_save{1} = pop;

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% Firefly Algorithm Main Loop

for it=1:MaxIt
    b = 1; %reset counter 
    newpop=repmat(firefly,nPop,1); 
    for i=1:nPop
        newpop(i).Cost = inf;
        P_temp = zeros(nPop,nVar);
        for j=1:nPop
            C_temp(j) = inf;
            S_temp(j) = 0;
            if pop(j).Cost < pop(i).Cost
                rij=norm(pop(i).Position-pop(j).Position)/dmax;
                beta=beta0*exp(-gamma*rij^m);
                %added .*
                e=delta.*unifrnd(-1,+1,VarSize);
                %e=delta*randn(VarSize);
                
                P_temp(j,:) = pop(i).Position ...
                                + beta*rand(VarSize).*(pop(j).Position-pop(i).Position) ...
                                + alpha*e;
                
                P_temp(j,:)=max(P_temp(j,:),VarMin);
                P_temp(j,:)=min(P_temp(j,:),VarMax);
                
                Kd = P_temp(j,1);
                Ki = P_temp(j,2);
                Kwi = P_temp(j,3);
                Kwa = P_temp(j,4);
                Kc = P_temp(j,5);
                S = P_temp(j,6);
                [C_temp(j),S_temp(j)] = ...
                             simHybrid(Kd, Ki, Kwi, Kwa, Kc, S,opt,data,atmo,batt,econ,uc(c),bc,dies,inso,wave,turb,cturb);
            end
        end
        for j=1:nPop
            newsol.Position = P_temp(j,:);
            newsol.Cost = C_temp(j);
            newsol.Surv = S_temp(j);
            if pop(j).Cost < pop(i).Cost
                if newsol.Surv < 0.99
                    newsol.Cost = inf;
                end
                
                if newsol.Cost <= newpop(i).Cost
                    newpop(i) = newsol;
                    if newpop(i).Cost<=BestSol.Cost
                        %BestSol=newpop(i);
                        BestSol.Cost = newpop(i).Cost;
                        BestSol.Surv{it,b}=newpop(i).Surv;
                        BestSol.Position{it,b} = newpop(i).Position;
                        b = b + 1;
                    end
                end
                
            end
        end
    end
    
    % Merge
    pop=[pop
         newpop];  %#ok
    %save populations
    pop_save{it+1} = newpop;
    % Sort
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    % Damp Mutation Coefficient
    alpha = alpha*alpha_damp;
    
end
disp(['Optimization complete after ' ...
    num2str(round(toc(tTot),2)) ' seconds.'])

%% Results

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
