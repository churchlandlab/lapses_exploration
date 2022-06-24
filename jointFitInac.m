function [NLL,fit,yFit]=jointFitInac(data, model, theseCondPairs, optimality)
inac = '';
% Model can be '
tic;

%% Set of possible two-parameter constraints
allConditionLabels = {'Auditory','Multisensory','Visual','AuditoryInac','MultisensoryInac','VisualInac'};
allConditions = {'A','M','V','a','m','v'};
allCondPairs = {'MA','MV','ma','mv','Aa','Mm','Vv'};
if ~exist('theseCondPairs','var')
    theseCondPairs = allCondPairs;
end

if ~exist('optimality','var')
    optimality = 'both';
end

% Equal parameter across conditions
Amodality = zeros(length(allCondPairs),length(allConditions));
for c = allCondPairs
    Amodality(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
    Amodality(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(2))) = -1;
end

% Unchanged parameter with inactivation
Afixed = zeros(length(allCondPairs),length(allConditions));
for c = {'Mm','Aa','Vv'}
    Afixed(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
    Afixed(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(2))) = -1;
end

% % Equal parameter within control/inactivation
Awithin = zeros(length(allCondPairs),length(allConditions));
for c = {'MA','MV','ma','mv'}
    Awithin(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
    Awithin(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(2))) = -1;
end

% % Default parameter for control
Adef = zeros(length(allCondPairs),length(allConditions));
for c = {'Mm','Aa','Vv'}
    Adef(strcmp(c,allCondPairs),strcmp(allConditions,c{1}(1))) = 1;
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
        
    case 'inattentionReward'
        modelType = 'inattention';
        paramLabels = {'pse','sigma','pInattention','bias'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 1, 1];
        paramDef = [12.5, 1, 0.1, 0.5];
        % Only bias & pse allowed to change
        constrainedParams = [2, 3];
        constraintA = {Afixed, Afixed};
        
    case 'inattentionValue'
        modelType = 'inattentionInac';
        paramLabels = {'pse','sigma','pInattention','bias','k'};
        paramLB = [0, 0, 0, 0, 0];
        paramUB = [20, 10, 1, 1,10];
        paramDef = [12.5, 1, 0.1, 1,1];
        % Only k allowed to change .
        constrainedParams = [1,2,3,4,5,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin,Adef};
        constraintB = [0,0,0,0,0,1];
        
    case 'idealObserverReward'
        modelType = 'noLapse';
        paramLabels = {'pse','sigma','lapseProb','lapseBias'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 0, 0];
        paramDef = [12.5, 1, 0, 0];
        % sigma constrained to be equal, lapseProb & lapseBias constrained to be equal & zero for all conditions
        constrainedParams = [2,3,4];
        constraintA = {Afixed,Afixed,Afixed};
        
    case 'idealObserverValue'
        modelType = 'standardInac';
        paramLabels = {'pse','sigma','lapseProb','lapseBias','k'};
        paramLB = [0, 0, 0, 0,0];
        paramUB = [20, 10, 0, 0,10];
        paramDef = [12.5, 1, 0, 0,1];
        % sigma constrained to be equal, lapseProb & lapseBias constrained to be equal & zero for all conditions
        constrainedParams = [1,2,3,4,5,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin,Adef};
        constraintB = [0,0,0,0,0,1];
        
    case 'fixedErrorReward'
        modelType = 'fixedLapse';
        paramLabels = {'pse','sigma','lapseProb','lapseBias'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 1, 1];
        paramDef = [12.5, 1, 0.1, 0];
        % sigma, lapseProb & lapseBias constrained to be equal for all conditions
        constrainedParams = [2,3,4];
        constraintA = {Afixed,Amodality,Amodality};
        
    case 'fixedErrorValue'
        modelType = 'standardInac';
        paramLabels = {'pse','sigma','lapseProb','lapseBias','k'};
        paramLB = [0, 0, 0, 0,0];
        paramUB = [20, 10, 1, 1,10];
        paramDef = [12.5, 1, 0, 0,1];
        % sigma constrained to be equal, lapseProb & lapseBias constrained to be equal & zero for all conditions
        constrainedParams = [1,2,3,3,4,4,5,5];
        constraintA = {Afixed,Afixed,Afixed,Amodality,Afixed,Amodality,Awithin,Adef};
        constraintB = [0,0,0,0,0,0,0,1];
        
    case 'explorationRewardR'
        modelType = 'exploration';
        paramLabels = {'pse','sigma','rLeft','rRight'};
        paramLB = [0, 0, 0, 0];
        paramUB = [20, 10, 100, 100];
        paramDef = [12, 2,3,3];
        % Only rRight allowed to change .
        constrainedParams = [1,2,3];
        constraintA = {Afixed,Afixed,Afixed};
        constraintB = [0,0,0,0];

        
    case 'explorationRewardL'
        modelType = 'exploration';
        paramLabels = {'pse','sigma','rLeft','rRight'};
        paramLB = [10, 0, 0, 0];
        paramUB = [15, 5, 10, 10];
        paramDef = [12.5, 3, 1, 1];
        %Only rLeft allowed to change.
        constrainedParams = [1,2,4];
        constraintA = {Afixed,Afixed,Afixed};
        
    case 'explorationValueR'
        modelType = 'explorationInac';
        paramLabels = {'pse','sigma','rLeft','rRight','kR'};
        paramLB = [10, 0, 0, 0, 0];
        paramUB = [15, 5, 10, 10, 100];
        paramDef = [12, 2,3,3, 1];
        % Only kR allowed to change .
        constrainedParams = [1,2,3,4,5,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin,Adef};
        constraintB = [0,0,0,0,0,1];
        inac = 'valueR';
        
    case 'explorationValueL'
        modelType = 'explorationInac';
        paramLabels = {'pse','sigma','rLeft','rRight','kL'};
        paramLB = [10, 0, 0, 0, 0];
        paramUB = [15, 5, 10, 10, 100];
        paramDef = [12, 2,1,1, 1];
        % Only kL allowed to change .
        constrainedParams = [1,2,3,4,5,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin,Adef};
        constraintB = [0,0,0,0,0,1];
        inac = 'valueL';
        
    case 'explorationEffort'
        modelType = 'explorationInac';
        paramLabels = {'pse','sigma','rLeft','rRight','k'};
        paramLB = [0, 0, 0, 0, -10];
        paramUB = [20, 10, 100, 100, 10];
        paramDef = [12, 2,3,3, 0];
        % Only k allowed to change .
        constrainedParams = [1,2,3,4,5,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin,Adef};
        constraintB = [0,0,0,0,0,0];
        inac = 'effort';
        
    case 'explorationEvidence'
        modelType = 'explorationInac';
        paramLabels = {'pse','sigma','rLeft','rRight','k'};
        paramLB = [0, 0, 0, 0, -10];
        paramUB = [20, 10, 100, 100, 10];
        paramDef = [12, 2,1,1, 0];
        % Only k allowed to change .
        constrainedParams = [1,2,3,4,5,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin,Adef};
        constraintB = [0,0,0,0,0,0];
        inac = 'evidence';

    case 'explorationPrior'
        modelType = 'explorationInac';
        paramLabels = {'pse','sigma','rLeft','rRight','k'};
        paramLB = [0, 0, 0, 0, 0];
        paramUB = [20, 10, 100, 100, 10];
        paramDef = [12, 2,1,1, 1];
        % Only k allowed to change .
        constrainedParams = [1,2,3,4,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin,Adef};
        constraintB = [0,0,0,0,0,1];
        inac = 'prior';
        
    case 'explorationNoise'
        modelType = 'explorationInac';
        paramLabels = {'pse','sigma','rLeft','rRight','k'};
        paramLB = [0, 0, 0, 0, 1];
        paramUB = [20, 10, 100, 100, 10];
        paramDef = [12, 2,1,1, 1];
        % Only k allowed to change .
        constrainedParams = [1,2,3,4,5];
        constraintA = {Afixed,Afixed,Afixed,Afixed,Awithin};
        constraintB = [0,0,0,0,0];
        inac = 'noise';

end

nConPars = length(constrainedParams);
nParams = length(paramLabels);
Aeq = zeros(nCondPairs*nConPars,nParams*nConds);
Beq = zeros(size(Aeq,1),1);

for p = 1:nConPars
    Aeq((1:nCondPairs)+(p-1)*nCondPairs,(0:nConds-1)*nParams+constrainedParams(p)) = constraintA{p}(condPairFlag,condFlag);
    Beq((1:nCondPairs)+(p-1)*nCondPairs)=constraintB(p);
end
inds = sum(abs(Aeq),2)~=0;
Aeq = Aeq(inds,:);
Beq = Beq(inds,:);
lb = repmat(paramLB,1,nConds);
ub = repmat(paramUB,1,nConds);
params0 = repmat(paramDef,1,nConds);
options = optimset('Display', 'off') ;

%% Fitting
% Optimality remains
switch optimality
    case 'unc'
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,modelType,params,inac),params0,[],[],Aeq,Beq,lb,ub,[],options);
    case 'both'
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,modelType,params,inac),params0,[],[],Aeq,Beq,lb,ub,@joint_optimality_constraint,options);
    case 'inac'
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,modelType,params,inac),params0,[],[],Aeq,Beq,lb,ub,@inactivation_optimality_constraint,options);
    case 'ctrl'
        [pars,fval,~,~,~,~,hessian]=fmincon(@(params)jointNLLsummary(data,modelType,params,inac),params0,[],[],Aeq,Beq,lb,ub,@optimality_constraint,options);
end
NLL=fval;
parsSE=sqrt(diag(inv(hessian)));
[~,thisNLL,yFit] = jointNLLsummary(data,modelType,pars,inac);

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



