function [params, handles, stats] = fitPalamedes(stimLevels, nHighResponses, nTrials, plotting, returnCIs)
%{
Sashank Pisupati
Churchland lab, Cold Spring Harbor Laboratory
April 3rd 2018

Function to fit & plot psychometric data using Prins & Kingdom's PALAMEDES toolbox.
(Prins, N & Kingdom, F. A. A. (2009) Palamedes:  Matlab routines for
analyzing psychophysical data. http://www.palamedestoolbox.org/ )

Requires the PALAMEDES toolbox: http://www.palamedestoolbox.org/
Requires the PALAMEDES folder to be in the MATLAB path


Takes an array of stimulus levels, number of "greater than" responses for each level &
number of trials for each level as inputs. For eg:
stimLevels = [ 1 2 3 4 5 6 7 8 9 10];
nHighResponses = [10,12,18,30,50,60,80,90,95,97];
nTrials = [100,100,100,100,100,100,100,100,100,100];


Returns parameter estimates (in the output structure "params") 
for a 4-param cumulative normal fit with lapse rates: 
Bias is the mean of the cumulative normal, 
sensitivity is the inverse standard deviation. 
Lapse rate and guess rate are distances of upper and lower asymptotes 
from 1 and 0 respectively, and are not constrained to be equal.

Optionally plots the data with wilson binomial CIs & the resulting fit,
and returns figure handles for all plot elements (in "handles", if plotting is 1).

Optionally returns nonparametric bootstrap re-estimates of the parameters &
95% CONFIDENCE intervals computed from these (in "stats", if returnCIs is 1), 
& plots 95% PREDICTION intervals of fit (if returnCIs is 1 and plotting is 1).
%}

%% Check arguments
%Required arguments
if ~exist('stimLevels','var')
    error('Please input stimulus levels!');
end
if ~exist('nHighResponses','var')
    error('Please input number of "greater than category boundary" responses!');
end
if ~exist('nTrials','var')
    error('Please input number of trials!');
end

%Optional arguments
if ~exist('plotting','var')
    plotting = 0;
    handles = [];
elseif ~plotting
    handles = [];
end
if ~exist('returnCIs','var')
    returnCIs = 0;
    stats = [];
elseif ~returnCIs
    handles = [];
end

%Check for PALAMEDES folder in path
if ~contains(path, ['Palamedes', pathsep])
    error('Please add the Palamedes folder to your MATLAB path!')
end

%% Setup fitting preferences
%Cumulative normal function
PF=@PAL_CumulativeNormal;
stats.logLikelihood = [];
paramLabels = {'bias','sensitivity', 'guessRate','lapseRate'};

%Parameter ranges for initial search
rangeStim = max(stimLevels)-min(stimLevels);
stats.paramRange.bias = [min(stimLevels),max(stimLevels)];
searchGrid.alpha=[min(stimLevels):rangeStim/50:max(stimLevels)]; %Bias
stats.paramRange.sensitivity = [(10^(-2))/rangeStim,(10^2)/rangeStim];
searchGrid.beta=(10.^[-2:0.1:2])/rangeStim;%Sensitivity
stats.paramRange.lapseRate = [0,1];
searchGrid.gamma=[0:0.1:0.5]; %Guess rate
stats.paramRange.guessRate = [0,1];
searchGrid.lambda=[0:0.1:0.5]; %Lapse rate

%Unconstrained asymptotes
lapseLimits = [0,1];
guessLimits = [0,1];

%Default fitting options
paramsFree=[1 1 1 1];
options=PAL_minimize('options');
lapseFit = 'nAple';




%% Fit data
%Fit parameters
[paramEstimates,stats.logLikelihood]= PAL_PFML_Fit(stimLevels,nHighResponses,nTrials,...
                                                    searchGrid,paramsFree,PF,...
                                                    'searchOptions',options,'lapseFit',lapseFit,...
                                                    'lapseLimits',lapseLimits,'guessLimits',guessLimits);
%Save parameter estimates with relevant labels
for p = 1:4
    params.(paramLabels{p}) = paramEstimates(p);
end
%Compute psychometric curve based on fit params
xout = min(stimLevels):0.1:max(stimLevels);
yout=PAL_CumulativeNormal(paramEstimates,xout);


%% Optionally plot data, fit
if plotting
    hold on
    %Get proportion chose "greater than" category
    pChoseHigh = nHighResponses./nTrials;
    %Get 95% Wilson binomial CIs for data
    z = 1.96;
    dataUpper = (pChoseHigh + z^2./(2*nTrials) + z .* ...
                sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
    dataLower = (pChoseHigh + z^2./(2*nTrials) - z .* ...
                sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
    
    %Plot datapoints
    handles.data = scatter(stimLevels, pChoseHigh,50,'MarkerFaceColor','k','MarkerEdgeColor','k');
    %Plot datapoint CIs
    for i = 1:length(pChoseHigh)
        handles.dataCI(i) = plot([stimLevels(i),stimLevels(i)],[dataUpper(i),dataLower(i)],'k-');
    end
    %Plot fitted curve
    handles.fit = plot(xout, yout,'color','k','linewidth',1);
    
    xlim([min(stimLevels)-1,max(stimLevels)+1]);
    xlabel('Stimulus strength')
    ylabel('Proportion chose high')
    hold off
end

%% Optionally get bootstrap confidence intervals
if returnCIs
    B = 400; %Number of bootstrap iterations

    %Nonparametric Bootstrap parameter re-estimates - resamples from data 
    [paramSEs, paramBootstraps]=PAL_PFML_BootstrapNonParametric(stimLevels, nHighResponses, nTrials,...
                                                                [], paramsFree, B, PF,...
                                                                'searchOptions',options, 'lapseLimits',lapseLimits,...
                                                                'guessLimits',guessLimits,'lapseFit',lapseFit,'searchGrid',searchGrid);
    %Compute 95% bootstrap CONFIDENCE intervals for parameters (paramCI)
    for p = 1:4
        stats.paramCIs.(paramLabels{p}) = 1.96*paramSEs(p);
        stats.paramBootstraps.(paramLabels{p}) = paramBootstraps(:,p);
    end
    
    %Compute and plot 95% bootstrap PREDICTION intervals for fit (fitPI)
    if plotting
        %Predictions from each bootstrap re-estimate
        hold on
        for b = 1:B
            yBootstrap(b,:)=PAL_CumulativeNormal(paramBootstraps(b,:),xout);
        end
        %Calculate 2.5% and 97.5% quantiles for predictions
        yOrdered = sort(yBootstrap,1);
        yLower = yOrdered(round(B*0.025),:);
        yUpper = yOrdered(round(B*0.975),:);
        handles.fitPI = patch([xout fliplr(xout)],[yLower fliplr(yUpper)],'k','EdgeColor','none');
        alpha(handles.fitPI,0.1);
        hold off
    end
    
end
end