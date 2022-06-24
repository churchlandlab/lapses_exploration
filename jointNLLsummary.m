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
    
    switch model
        case {'standard','fixedLapse','noLapse','restrictedLapse'}
            % Descriptive psychometric models with various lapse constraints
            pChooseHigh = standardModel(theseStim,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = standardModel(xVal,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            
        case 'standardInac'
            % Descriptive psychometric models fit jointly to control & inactivation
            pChooseHigh = standardInacModel(theseStim,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5));
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = standardInacModel(xVal,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5));
            
        case 'inattention'
            % Inattention model, with a small prob. of not attending
            pChooseHigh = inattentionModel(theseStim,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = inattentionModel(xVal,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            
        case 'inattentionInac'
            % Inattention model fit jointly to inactivation & control data
            pChooseHigh = inattentionInacModel(theseStim,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5));
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = inattentionInacModel(xVal,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5));
            
        case 'exploration'
            % Softmax exploration model (variable temperature, unit reward)
            pChooseHigh = explorationModel(theseStim,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = explorationModel(xVal,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            
        case 'explorationReparam'
            % Softmax exploration model (variable reward, unit temperature)
            pChooseHigh = explorationReparamModel(theseStim,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = explorationReparamModel(xVal,params((c-1)*4+1),params((c-1)*4+2),params((c-1)*4+3),params((c-1)*4+4));
            
        case 'explorationInac'
            % Softmax exploration model fit jointly to control & inac data
            pChooseHigh = explorationInacModel(theseStim,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5),inac);
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = explorationInacModel(xVal,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4),params((c-1)*5+5),inac);
            
        case 'explorationThompson'
            % Thompson sampling model, uncertainty-dependent exploration
            pChooseHigh = explorationModel(theseStim,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4));
            ll=log(pChooseHigh).*(theseResp) + log(1-pChooseHigh).*(theseTrials-theseResp);
            thisNLL(c)=-sum(ll);
            yFit(c,:) = explorationModel(xVal,params((c-1)*5+1),params((c-1)*5+2),params((c-1)*5+3),params((c-1)*5+4));
            
            
    end
end
NLL = sum(thisNLL);
end

function pChooseHigh = standardModel(theseStim,pse,sigma,lapseProb,lapseBias)
% Standard descriptive psychometric model: Cumulative normal scaled by
% lapse rates
pHigh=normcdf(theseStim,pse,sigma);
pChooseHigh=lapseProb*lapseBias+(1-lapseProb)*pHigh;
end

function pChooseHigh = standardInacModel(theseStim,pse,sigma,lapseProb,lapseBias,k)
% Standard descriptive psychometric model of inactivation: Cumulative normal scaled by
% lapse rates, with inactivation biasing the p.s.e
pse = norminv(1/(1+k), pse,sigma);
pHigh=normcdf(theseStim,pse,sigma);
pChooseHigh=lapseProb*lapseBias+(1-lapseProb)*pHigh;
end

function pChooseHigh = inattentionModel(theseStim,pse,sigma,pInattention,bias)
% Inattention model: lapses arising from biased random decisions on 
% pInattention fraction of trials
pHigh=normcdf(theseStim,pse,sigma);
pChooseHigh=(pInattention)*bias+(1-pInattention)*pHigh;
end

function pChooseHigh = inattentionInacModel(theseStim,pse,sigma,pInattention,bias,k)
% Inattention model with inactivation biasing evidence
pse = norminv(1/(1+k), pse,sigma);
b = bias/(1-bias);
bias = b*k/(1+b*k);
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
marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
lastwarn('')

pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);

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
marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
lastwarn('')

pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);

[~, msgid] = lastwarn;
if strcmp(msgid, 'MATLAB:integral:MinStepSize')
    pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
    pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
end
end


function pChooseHigh = explorationInacModel(theseStim,pse,sigma,rL,rR,k,inac)
% Exploration model with inactivation biasing evidence, prior, noise
% unilateral value (multiplicative) or effort (additive)
beta = 1;
switch inac
    case 'evidence'
        pse = pse+sigma*k;
        warning('off','MATLAB:integral:NonFiniteValue');
        warning('off','MATLAB:integral:MinStepSize')
        marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
        lastwarn('')
        
        pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:integral:MinStepSize')
            pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
            pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
        end
        
    case 'prior'
        pse = norminv(1/(1+k), pse,sigma);
        warning('off','MATLAB:integral:NonFiniteValue');
        warning('off','MATLAB:integral:MinStepSize')
        marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
        lastwarn('')
        
        pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:integral:MinStepSize')
            pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
            pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
        end
        
    case 'noise'
        sigma = sigma*k;
        warning('off','MATLAB:integral:NonFiniteValue');
        warning('off','MATLAB:integral:MinStepSize')
        marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
        lastwarn('')
        
        pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:integral:MinStepSize')
            pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
            pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
        end
        
    case 'valueR'
        rR = rR*k;
        warning('off','MATLAB:integral:NonFiniteValue');
        warning('off','MATLAB:integral:MinStepSize')
        marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
        lastwarn('')
        
        pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:integral:MinStepSize')
            pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
            pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
        end
        
    case 'valueL'
        rL = rL*k;
        warning('off','MATLAB:integral:NonFiniteValue');
        warning('off','MATLAB:integral:MinStepSize')
        marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
        lastwarn('')
        
        pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:integral:MinStepSize')
            pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
            pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
        end
        
    case 'effort'
        warning('off','MATLAB:integral:NonFiniteValue');
        warning('off','MATLAB:integral:MinStepSize')
        marginal_ll = @(p,s,pse,sigma,beta,rR,rL) (1./(1+exp(-beta*(p*(rR+rL)-rL+k)))).*(normpdf(norminv(1-p,0,sigma),pse-s,sigma))./(normpdf(norminv(1-p,0,sigma),0,sigma));
        lastwarn('')
        
        pChooseHigh = integral(@(p)marginal_ll(p,theseStim,pse,sigma,beta,rR,rL),0,1,'ArrayValued',1);
        
        [~, msgid] = lastwarn;
        if strcmp(msgid, 'MATLAB:integral:MinStepSize')
            pChooseHigh(theseStim>pse+4*sigma) =(1./(1+exp(-beta*(rR+k))))*normcdf(theseStim(theseStim>pse+4*sigma),pse,sigma);
            pChooseHigh(theseStim<pse-5*sigma) =(1./(1+exp(-beta*(-rL+k))))*(1-normcdf(theseStim(theseStim<pse-5*sigma),pse,sigma));
        end
        
end
end