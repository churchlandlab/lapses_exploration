function [NLL,fit,yFit]=jointFit(data, model,theseCondPairs,optimality)

% Model can be '
tic;
if ~exist('optimality','var') && isempty((strfind(strcat(theseCondPairs{:}),'N')))
    optimality = 'matched';
elseif  ~exist('optimality','var')
    optimality = 'neutral';
end

%% Set of possible two-parameter constraints
allConditionLabels = {'Auditory','Multisensory','MultisensoryNeutral','AuditoryLowRel','MultisensoryLowRel','Visual'};
allConditions = {'A','M','N','a','m','V'};
allCondPairs = {'MA','MV','MN','Aa','Mm','AN','AV','VN'};
if ~exist('theseCondPairs','var')
    theseCondPairs = allCondPairs;
end

% Equal parameter for all condition pairs
Afixed = zeros(length(allCondPairs),length(allConditions));
for c = allCondPairs
    Afixed(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
    Afixed(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(2))) = -1;
end

% Equal parameter for conditions with equal salience
Asalience = zeros(length(allCondPairs),length(allConditions));
for c = {'MN'}
    Asalience(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
    Asalience(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(2))) = -1;
end

% Equal parameter for conditions of the same uncertainty
Auncertainty = zeros(length(allCondPairs),length(allConditions));
for c = {'AN'}
    Auncertainty(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
    Auncertainty(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(2))) = -1;
end

% Equal parameter for conditions of the same modality
Amodality = zeros(length(allCondPairs),length(allConditions));
for c = {'MN','Aa','Mm'}
    Amodality(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
    Amodality(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(2))) = -1;
end
%% Get parameter constraints for given model & dataset

theseConds = {data.condition};
nCondPairs = length(theseCondPairs);
nConds = length(theseConds);
condFlag = false(1,length(allConditions));
for c = 1:nConds
    condFlag = logical(condFlag + strcmp(theseConds{c},allConditionLabels));
end
condPairFlag = false(1,length(allCondPairs));
for c = 1:nCondPairs
    condPairFlag = logical(condPairFlag+strcmp(theseCondPairs{c},allCondPairs));
end
% Express linear constraints as Aeq*params' = Beq;
switch model
    
    case 'standard'
        paramLabels = {'pse','sigma','lapseProb','lapseBias'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 1, 1];
        paramDef = [12.5, 1, 0.1, 0];
        % No constraints
        constrainedParams = [];
        
    case 'fixedLapse'
        paramLabels = {'pse','sigma','lapseProb','lapseBias'};
        paramLB = [0, 0, 0, -1];
        paramUB = [20, 10, 1, 1];
        paramDef = [12.5, 1, 0.1, 0];
        % lapseProb & lapseBias constrained to be equal for all conditions
        constrainedParams = [2,3,4];
        constraints = {Auncertainty,Afixed,Afixed};
        
    case 'noLapse'
        paramLabels = {'pse','sigma','lapseProb','lapseBias'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 0, 0];
        paramDef = [12.5, 1, 0, 0];
        % lapseProb & lapseBias constrained to be equal & zero for all conditions
        constrainedParams = [2,3,4];
        constraints = {Auncertainty,Afixed,Afixed};
        
    case 'restrictedLapse'
        paramLabels = {'pse','sigma','lapseProb','lapseBias'};
        paramLB = [0, 0, 0, -0.1];
        paramUB = [20, 10, 0.1, 0.1];
        paramDef = [12.5, 1, 0, 0];
        % lapseProb & lapseBias constrained to be equal & <10% for all conditions
        constrainedParams = [2];
        constraints = {Auncertainty};
        
    case 'inattention'
        paramLabels = {'pse','sigma','pInattention','prior'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 1, 1];
        paramDef = [12.5, 1, 0.9, 0.5];
        % pInattention constrained to be equal for equal salience conditions
        constrainedParams = [2,3];
        constraints = {Auncertainty,Asalience};

        
    case 'exploration'
        paramLabels = {'pse','sigma','beta','rewardRatio'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 100, 100];
        paramDef = [12.5, 1, 5, 1];
        % Exploration (beta) allowed to vary across conditions, unit reward
        constrainedParams = [2];
        constraints = {Auncertainty};
        
    case 'explorationReparam'
        paramLabels = {'pse','sigma','rLeft','rRight'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 10, 10];
        paramDef = [12.5, 1, 1, 1];
        % Reward (rLeft/Right) allowed to vary across conditions, unit beta
        constrainedParams = [];
        constraints = {};
        
    case 'explorationThompson'
        paramLabels = {'pse','sigma','beta','rewardRatio','k0'};
        paramLB = [0, 0, 0, 0,-1];
        paramUB = [20, 10, 100, 100,1];
        paramDef = [12.5, 1, 5, 1,0];
        % Fixed scaling between sigma (uncertainty) and beta (exploration)
        constrainedParams = [2,5];
        constraints = {Auncertainty,Afixed};        
end

nConPars = length(constrainedParams);
nParams = length(paramLabels);
Aeq = zeros(nCondPairs*nConPars,nParams*nConds);
Beq = zeros(size(Aeq,1),1);

for p = 1:nConPars
    Aeq((1:nCondPairs)+(p-1)*nCondPairs,(0:nConds-1)*nParams+constrainedParams(p)) = constraints{p}(condPairFlag,condFlag);
end

Aeq = Aeq(sum(abs(Aeq),2)~=0,:);
Beq = zeros(size(Aeq,1),1);
lb = repmat(paramLB,1,nConds);
ub = repmat(paramUB,1,nConds);
params0 = repmat(paramDef,1,nConds);
options = optimset('Display', 'off') ;

%% Fitting
switch optimality
    case 'unc'
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,model,params),params0,[],[],Aeq,Beq,lb,ub,[],options);
    case 'matched'
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,model,params),params0,[],[],Aeq,Beq,lb,ub,@optimality_constraint,options);
    case 'neutral'
        A = [0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 -1 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0]; b = [0;0;0];
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,model,params),params0,[],[],Aeq,Beq,lb,ub,@optimality_constraint_neutral,options);
    case 'exploration_neutral'
        A = [0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0]; b = 0;
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,model,params),params0,[],[],Aeq,Beq,lb,ub,@optimality_constraint_neutral_exploration,options);
    case 'thompson'
        A = [0 0 1 0 0 0 0 -1 0 0 0 0 0 0 0]; b = 0;
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,model,params),params0,A,b,Aeq,Beq,lb,ub,@thompson_constraint,options);
    case 'thompson_neutral'
        A = [0 0 1 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0]; b = 0;
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,model,params),params0,A,b,Aeq,Beq,lb,ub,@thompson_constraint_neutral,options);
end
NLL=fval;
parsSE=sqrt(diag(inv(hessian)));
[~,thisNLL,yFit] = jointNLLsummary(data,model,pars);

for c = 1:nConds
    fit(c).condition = theseConds{c};
    fit(c).NLL = thisNLL(c);
    for p = 1:nParams
        fit(c).(paramLabels{p}) = pars((c-1)*nParams+p);
        fit(c).([paramLabels{p},'SE']) = parsSE((c-1)*nParams+p);
    end
end
toc
end



