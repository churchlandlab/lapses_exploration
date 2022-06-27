% Function to evaluate negative log likelihood for different models of lapses
% from Pisupati, Chartarifsky et al. eLife 2021. See Methods for corresponding equations
function [NLL,thisNLL,yFit]=jointNLLsummary(data,model,params,inac)

if ~exist('inac','var')
    inac = '';
end

% Extract conditions
nConds = length(data);
thisNLL = zeros(1,nConds);
xVal = 5:0.1:20;
yFit = zeros(nConds,length(xVal));

% Model parameters
for c = 1:nConds
    theseTrials = data(c).nTrials;
    theseStim = data(c).stimRates;
    theseResp = data(c).nHighResponses;
    isInac = 0;

    switch model
        case {'standard','fixedLapse','noLapse','restrictedLapse'}
            % Descriptive psychometric models with various lapse constraints
            % (fixedLapse is identical in form to motor error model)
            modelSpec = @standardModel;
            
        case 'inattention'
            % Inattention model, with a small prob. of not attending
            modelSpec = @inattentionModel;
            
        case {'exploration','explorationThompson'}
            % Softmax exploration model (variable temperature, unit reward)
            modelSpec = @explorationModel;
            
        case 'explorationReparam'
            % Reparametrized softmax exploration model (variable reward, unit temperature)
            modelSpec = @explorationReparamModel;
            
        case 'standardInac'
            % Descriptive psychometric model/motor error fit jointly to control & inactivation
            modelSpec = @standardInacModel;
            isInac = 1;
            
        case 'inattentionInac'
            % Inattention model fit jointly to control & inactivation data
            modelSpec = @inattentionInacModel;
            isInac = 1;
            
        case 'explorationInac'
            % Softmax exploration model fit jointly to control & inactivation data
            modelSpec = @explorationInacModel;
            isInac = 1;
    end
    
    if ~isInac
        % 4 effective params per condition
        pChooseHigh = modelSpec(theseStim,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
        ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
        thisNLL(c)=-sum(ll);
        yFit(c,:) = modelSpec(xVal,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
    else
        % 5 effective params per condition - 4 + 1 inactivation param
        pChooseHigh = modelSpec(theseStim,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5),inac);
        ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
        thisNLL(c)=-sum(ll);
        yFit(c,:) = modelSpec(xVal,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5),inac);
    end
end
NLL = sum(thisNLL);
end

%% ------------------------- Models --------------------------------------- 
function pChooseHigh = standardModel(theseStim,pse,sigma,lapseProb,lapseBias)
% Standard descriptive psychometric model & motor error model: 
% Cumulative normal scaled by lapse rates
pHigh=normcdf(theseStim,pse,sigma);
pChooseHigh=lapseProb*lapseBias+(1-lapseProb)*pHigh;
end

function pChooseHigh = inattentionModel(theseStim,pse,sigma,pInattention,bias)
% Inattention model: lapses arising from biased random decisions on
% pInattention fraction of trials
pHigh=normcdf(theseStim,pse,sigma);
pChooseHigh=(pInattention)*bias+(1-pInattention)*pHigh;
end

function pChooseHigh = explorationModel(theseStim,pse,sigma,beta,rewardBias)
% Exploration model with unit reward (with possible biases) and softmax
% exploration, with varying temperatures across conditions. Latent variable
% (noisy posterior) numerically integrated out to obtain pChooseHigh
rL = 1;
rR = rL*rewardBias;


warning('off','MATLAB:integral:NonFiniteValue');
warning('off','MATLAB:integral:MinStepSize')
% Analytical log likelihood conditional on latent variable 'p' (posterior)
% See eLife '21 manuscript for derivation
conditional_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
lastwarn('')
% Marginalizing out latent variable 'p' (posterior)
pChooseHigh = integral(@(p)conditional_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);

[~, msgid] = lastwarn;
if strcmp(msgid, 'MATLAB:integral:MinStepSize')
    pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
    pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
end
end

function pChooseHigh = explorationReparamModel(theseStim,pse,sigma,rL,rR)
% Reparametrized exploration model with unit temperature (beta) and varying
% rewards across conditions - to allow for manipulations to affect rewards
% without changing temperature.
beta = 1;

warning('off','MATLAB:integral:NonFiniteValue');
warning('off','MATLAB:integral:MinStepSize')
% Analytical log likelihood conditional on latent variable 'p' (posterior)
conditional_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
lastwarn('')
% Marginalizing out latent variable 'p' (posterior)
pChooseHigh = integral(@(p)conditional_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);

[~, msgid] = lastwarn;
if strcmp(msgid, 'MATLAB:integral:MinStepSize')
    pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
    pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
end
end

function pChooseHigh = standardInacModel(theseStim,pse,sigma,lapseProb,lapseBias,k,inac)
% Standard/motor error model of inactivation: Cumulative normal scaled by
% lapse rates, with inactivation biasing the p.s.e
pse = norminv(1/(1+k), pse,sigma);
pHigh=normcdf(theseStim,pse,sigma);
pChooseHigh=lapseProb*lapseBias+(1-lapseProb)*pHigh;
end

function pChooseHigh = inattentionInacModel(theseStim,pse,sigma,pInattention,bias,k,inac)
% Inattention model with inactivation biasing prior/value (i.e. affecting both lapse rates)
% note that biasing evidence is identical in inattention & motor error models
pse = norminv(1/(1+k), pse,sigma);
b = bias/(1-bias);
bias = b*k/(1+b*k);
pHigh=normcdf(theseStim,pse,sigma);
pChooseHigh=(pInattention)*bias+(1-pInattention)*pHigh;
end

function pChooseHigh = explorationInacModel(theseStim,pse,sigma,rL,rR,k,inac)
% Exploration model with inactivation biasing evidence, prior, noise
% unilateral value (multiplicative) or effort (additive)
beta = 1;
warning('off','MATLAB:integral:NonFiniteValue');
warning('off','MATLAB:integral:MinStepSize')
switch inac
    case 'evidence'
        %Additive change to evidence 
        pse = pse+sigma*k;
        
    case 'prior'
        %Change to category prior    
        pse = norminv(1/(1+k), pse,sigma);
        
    case 'noise'
        %Multiplicative change to sensory noise  
        sigma = sigma*k;
        
    case 'valueR'
        %Multiplicative change to right values        
        rR = rR*k;
        
    case 'valueL'
        %Multiplicative change to left values        
        rL = rL*k;
        
    case 'effort'
        %Additive change to one of the value (directly affects conditional_ll)
end
lastwarn('')
if strcmp(inac,'effort')
    % Effort inactivation
    conditional_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL+k)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
    % Marginalizing out latent variable 'p' (posterior)
    pChooseHigh = integral(@(p)conditional_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
    
    [~, msgid] = lastwarn;
    if strcmp(msgid, 'MATLAB:integral:MinStepSize')
        pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR+k))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
        pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL+k))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
    end
else  
    % Common conditional_ll for all other inactivations
    conditional_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
    % Marginalizing out latent variable 'p' (posterior)
    pChooseHigh = integral(@(p)conditional_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
    
    [~, msgid] = lastwarn;
    if strcmp(msgid, 'MATLAB:integral:MinStepSize')
        pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
        pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
    end
end
end
