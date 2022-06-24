function [handles, params, stats] = fitPalamedesReparam(stimLevels, nHighResponses, nTrials, plotting, returnCIs)
%{
Sashank Pisupati
Churchland lab, Cold Spring Harbor Laboratory
August 21st 2018

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


Returns REPARAMETRIZED parameter estimates for a 4-param cumulative normal fit with lapse
rates: 
mu is the mean of the cumulative normal, 
sigma is the standard deviation. 
Lapse probability is the sum of distances of upper and lower asymptotes
from 1 and 0 respectively, and lapse bias is their difference.

Optionally plots the data with wilson binomial CIs & the resulting fit,
and returns a figure handle (if plotting is 1).

Optionally returns bootstrap re-estimates of the parameters &
95% confidence intervals computed from these (if returnCI is 1).
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
end
if ~exist('returnCIs','var')
    returnCIs = 0;
end

%Check for PALAMEDES folder in path
if ~contains(path, ['Palamedes', pathsep])
    error('Please add the Palamedes folder to your MATLAB path!')
end

%% Setup fitting preferences
%Cumulative normal function
PF=@PAL_CumulativeNormal;
stats.logLikelihood = [];
paramLabels = {'mu','sigma', 'lapseProb','lapseBias'};

%Parameter ranges for initial search
rangeStim = max(stimLevels)-min(stimLevels);
stats.paramRange.bias = [min(stimLevels),max(stimLevels)];
searchGrid.alpha=[min(stimLevels):rangeStim/50:max(stimLevels)]; %Bias
stats.paramRange.sensitivity = [(10^(-2))/rangeStim,(10^2)/rangeStim];
searchGrid.beta=(10.^[-2:0.1:2])/rangeStim;%Sensitivity
stats.paramRange.lapseRate = [0,0.5];
searchGrid.gamma=[0:0.1:0.5]; %Guess rate
stats.paramRange.guessRate = [0,0.5];
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
%Fit parameters
[paramEstimates,stats.logLikelihood]= PAL_PFML_Fit(stimLevels,nHighResponses,nTrials,...
                                                    searchGrid,paramsFree,PF,...
                                                    'searchOptions',options,'lapseFit',lapseFit,...
                                                    'lapseLimits',lapseLimits,'guessLimits',guessLimits);
%Save reparametrized parameter estimates with relevant labels
params.(paramLabels{1}) = paramEstimates(1);
params.(paramLabels{2}) = 1./paramEstimates(2);
params.(paramLabels{3}) = paramEstimates(3)+paramEstimates(4);
params.(paramLabels{4}) = paramEstimates(3)-paramEstimates(4);

%Compute psychometric curve based on fit params
xout = min(stimLevels):0.1:max(stimLevels);
xout = [7:0.1:20];
yout=PAL_CumulativeNormal(paramEstimates,xout);


%% Optionally plot data, fit
if plotting
    hold on
    %Get proportion chose "greater than" category
    pChoseHigh = nHighResponses./nTrials;
    %Get 95% Wilson binomial CIs for data
    z = 1.96;
    dataUpper = (pChoseHigh + z^2./(2*nTrials) + z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
    dataLower = (pChoseHigh + z^2./(2*nTrials) - z .* sqrt(pChoseHigh.*(1-pChoseHigh)./nTrials + z^2./(4*nTrials.^2))) ./ (1 + z^2./nTrials);
    
    %Plot datapoints
    handles.data = scatter(stimLevels, pChoseHigh,50,'MarkerFaceColor','k','MarkerEdgeColor','k');
    %Plot datapoint CIs
    for i = 1:length(pChoseHigh)
        handles.dataCI(i) = plot([stimLevels(i),stimLevels(i)],[dataUpper(i),dataLower(i)],'k-');
    end
    %Plot fitted curve
    handles.fit = plot(xout, yout,'color','k');
    
    xlim([min(stimLevels)-1,max(stimLevels)+1]);
    xlabel('Stimulus strength')
    ylabel('Proportion chose high')
    hold off
else 
    handles = [];
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
    stats.paramBootstraps.(paramLabels{1}) = paramBootstraps(:,1);
    stats.paramBootstraps.(paramLabels{2}) = 1./paramBootstraps(:,2);
    stats.paramBootstraps.(paramLabels{3}) = paramBootstraps(:,3)+paramBootstraps(:,4);
    stats.paramBootstraps.(paramLabels{4}) = paramBootstraps(:,3)-paramBootstraps(:,4);
    
    for p = 1:4
        stats.paramCIs.(paramLabels{p}) = 1.96*std(stats.paramBootstraps.(paramLabels{p}));
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
    else
        handles = [];
    end
else
    for p = 1:4
        stats.paramBootstraps.(paramLabels{p})=NaN;
        stats.paramCIs.(paramLabels{p}) = NaN;
    end
end
end