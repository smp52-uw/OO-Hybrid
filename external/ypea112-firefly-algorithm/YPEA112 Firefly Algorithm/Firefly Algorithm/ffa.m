function [output] = ffa(opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb)


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
% Changes: by SM Palmer 7-10-23 to make it a function 


tTot = tic;

%% Problem Definition

%CostFunction=@(x) Rosenbrock(x);        % Cost Function

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
MaxIt=opt.ffa.max;   % Maximum Number of Iterations

%nPop=25;            % Number of Fireflies (Swarm Size)
nPop=opt.ffa.pop;    % Number of Fireflies (Swarm Size)

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

% Create Initial Fireflies
rng('shuffle','Twister') %maybe will make the cluster rand work better?
for i=1:nPop %(i = 1:sim_length,opt.bf.maxworkers)
   %pop(i).Position=unifrnd(VarMin,VarMax,VarSize); 
   pop(i).Position = VarMin + (VarMax - VarMin).*rand(VarSize); %maybe works for the cluster??
end
C_temp = nan(1,nPop);
S_temp = nan(1,nPop);
parfor (i = 1:nPop,opt.bf.maxworkers)
   Kd = pop(i).Position(1);
   Ki = pop(i).Position(2);
   Kwi = pop(i).Position(3);
   Kwa = pop(i).Position(4);
   Kc = pop(i).Position(5);
   S = pop(i).Position(6);
   %pop(i).Cost=CostFunction(pop(i).Position);
   [C_temp(i),S_temp(i)] = ...
                simHybrid(Kd, Ki, Kwi, Kwa, Kc, S,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);
end

for i=1:nPop
   pop(i).Cost = C_temp(i);
   pop(i).Surv = S_temp(i);
   if pop(i).Surv < 0.99
       pop(i).Cost = inf;
   end

   if pop(i).Cost<=BestSol.Cost
       BestSol.Cost = pop(i).Cost;
       BestSol.Surv{1}=pop(i).Surv;
       BestSol.Position{1} = pop(i).Position;
   end
end
%Save each population
for s = 1:nPop
    output.cost{1}(s) = pop(s).Cost;
    output.surv{1}(s) = pop(s).Surv;
    output.Kd_run{1}(s) = pop(s).Position(1);
    output.Ki_run{1}(s) = pop(s).Position(2);
    output.Kwi_run{1}(s) = pop(s).Position(3);
    output.Kwa_run{1}(s) = pop(s).Position(4);
    output.Kc_run{1}(s) = pop(s).Position(5);
    output.S_run{1}(s) = pop(s).Position(6);
end

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% Firefly Algorithm Main Loop

for it=1:MaxIt
    newpop=repmat(firefly,nPop,1); 
    for i=1:nPop
        newpop(i).Cost = inf;
        newpop(i).Surv = 0;
        newpop(i).Position = nan(VarSize(1),VarSize(2));
        P_temp = zeros(nPop,nVar);
        C_temp = nan(1,nPop);
        S_temp = nan(1,nPop);
        for j=1:nPop
            C_temp(j) = inf;
            S_temp(j) = 0;
            if pop(j).Cost < pop(i).Cost
                rij=norm(pop(i).Position-pop(j).Position)/dmax;
                beta=beta0*exp(-gamma*rij^m);
                %added .*
                %e=delta.*unifrnd(-1,+1,VarSize);
                e = delta.*(-1 + (1 + 1).*rand(VarSize)); %maybe works on the cluster?
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
                             simHybrid(Kd, Ki, Kwi, Kwa, Kc, S,opt,data,atmo,batt,econ,uc,bc,dies,inso,wave,turb,cturb);
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
                        BestSol.Surv{it}=newpop(i).Surv;
                        BestSol.Position{it} = newpop(i).Position;
                    end
                end
                
            end
        end
    end
    
    % Merge
    pop=[pop
         newpop];  %#ok
    %save populations
    for s = 1:nPop
        output.cost{it+1}(s) = newpop(s).Cost;
        output.surv{it+1}(s) = newpop(s).Surv;
        output.Kd_run{it+1}(s) = newpop(s).Position(1);
        output.Ki_run{it+1}(s) = newpop(s).Position(2);
        output.Kwi_run{it+1}(s) = newpop(s).Position(3);
        output.Kwa_run{it+1}(s) = newpop(s).Position(4);
        output.Kc_run{it+1}(s) = newpop(s).Position(5);
        output.S_run{it+1}(s) = newpop(s).Position(6);
    end

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

%save output for optHybrid
output.min.kWd{1} = BestSol.Position{end}(1);
output.min.kWi{1} = BestSol.Position{end}(2);
output.min.kWwi{1} = BestSol.Position{end}(3);
output.min.kWwa{1} = BestSol.Position{end}(4);
output.min.kWc{1} = BestSol.Position{end}(5);
output.min.Smax{1} = BestSol.Position{end}(6);


end