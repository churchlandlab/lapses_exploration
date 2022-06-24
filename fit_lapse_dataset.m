% Script to plot data & generate model fits from all models for all conditions in the
% "Lapses" dataset (Pisupati, Chartarifsky et al. 2019)
% stored at http://repository.cshl.edu/id/eprint/38957/

% Dependencies for fitting:
%   jointFit.m
%   jointFitInac.m
%   jointNLLsummary.m
% Dependencies for plotting:
%   fitPalamedesReparam.m
%   axprefs.m
%   dataset.mat that contains fit parameters for the dataset above
% Palamedes toolbox: (not included, please download from https://www.palamedestoolbox.org/)
%   Prins, N. & Kingdom, F.A.A. (2018). Applying the Model-Comparison Approach to Test Specific Research Hypotheses in Psychophysical Research Using the Palamedes Toolbox. Frontiers in psychology, 9:1250
%   DOI=10.3389/fpsyg.2018.01250

%% Multisensory data - descriptive models
data = dataset.multisensory;
% Descriptive models:
% Fit a 4 parameter sigmoid function with either zero, fixed, restricted or
% variable asymptotes
models = {'noLapse','fixedLapse','restrictedLapse','standard'};
nParams = [6,8,12,12];
subjects = {'sp01','sp02','sp03','sp04','sp05','sp06','lc35','lc37','lc38','lc39','lc40','lc43','lc44','lc45','lc46','lc47','lc48'};
fit = 0;
compare = 1;

clear bic aic nTrials nPars nll
for pooled = [1,0]
    for k = 1:length(models)
        model = models{k};
        if pooled
            subject = 'metaRat';
            j = find(strcmp(subject,{data.ratName}));
            % Pooled model fitting
            if fit
                [n,f,y] = jointFit(data(j).controlSummaryData, model,{'AV','MA'},'unc');
                fits.(model) = f;
            end
            % Pooled model comparison
            if compare
                [nll(k),bic(k),aic(k)] = modelMetrics(data(j),models{k},nParams(k));
            end
            
        else
            for i = 1:length(subjects)
                subject = subjects{i};
                j = find(strcmp(subject,{data.ratName}));
                % Individual model fitting
                if fit
                    [n,f,y] = jointFit(data(j).controlSummaryData, model,{'AV','MA'},'unc');
                    dataset.multisensory(j).(model)=f;
                end
                % Individual model comparison
                if compare
                    [nll(k,j),bic(k,j),aic(k,j)] = modelMetrics(data(j),models{k},nParams(k));
                end
            end
        end
    end
    bic = bic-min(bic);
    aic = aic-min(aic);
    % Plot model comparison
    figure(2)
    subplot(2,3,3*(2-pooled))
    plotModelComparison(bic,aic,models,pooled)
end
% Plot pooled data
subplot(2,3,[1,2,4,5])
mods = [1 2 3];
colors = {[0 0.5 0], [0.8 0 0],[0 0 0.8]};
plotData(data,mods,colors)

%% Multisensory data - theoretical models
% Theoretical models:
% Fit the ideal observer (O), fixed motor error (F), inattention (I) or
% uncertainty-guided exploration (E) models to multisensory choice data
% Conditions - Auditory (A), Visual (V), Multisensory (M)
models = {'idealObserver','motorError','inattention','exploration','explorationThompson'};
nParams = [5,7,11,11,10];
subjects = {'sp01','sp02','sp03','sp04','sp05','sp06','lc35','lc37','lc38','lc39','lc40','lc43','lc44','lc45','lc46','lc47','lc48'};
pooled = 0;
fit = 0;
compare = ~fit;
plotSupplemental = 1;

clear bic aic nTrials nPars nll
for k = 1:length(models)
    model = models{k};
    if pooled
        subject = 'metaRat';
        j = find(strcmp(subject,{data.ratName}));
        % Pooled model fitting
        if fit
            if strcmp(model,'explorationThompson')
                [n,f,y] = jointFit(data(j).controlSummaryData(1:3), model,{'AV','MV','MA'},'thompson');
            else
                [n,f,y] = jointFit(data(j).controlSummaryData(1:3), model,{'AV','MV','MA'});
            end
        end
        % Pooled model comparison
        if compare
            [nll(k),bic(k),aic(k)] = modelMetrics(data(j),models{k},nParams(k));
        end
        
    else
        for i = 1:length(subjects)
            subject = subjects{i};
            j = find(strcmp(subject,{data.ratName}));
            % Individual model fitting
            if fit
                if strcmp(model,'explorationThompson')
                    [n,f,y] = jointFit(data(j).controlSummaryData(1:3), model,{'AV','MV','MA'},'thompson');
                else
                    [n,f,y] = jointFit(data(j).controlSummaryData(1:3), model,{'AV','MV','MA'});
                end
            end
            % Individual model comparison
            if compare
                [nll(k,j),bic(k,j),aic(k,j)] = modelMetrics(data(j),models{k},nParams(k));
            end
        end
    end
end
bic = bic-min(bic);
aic = aic-min(aic);

if plotSupplemental
    figure(201)
    subplot(2,3,1)
    imageModelComparison(bic,'BIC',models)
    subplot(2,3,2)
    imageModelComparison(aic,'AIC',models)
    subplot(2,3,3)
    plotModelComparison(bic,aic,models,pooled)
end
%% Congruence manipulation - theoretical models
data = dataset.neutral;
% Fit the ideal observer (O), fixed motor error (F), inattention (I) or
% uncertainty-guided exploration (E) models to the "neutral" manipulation
% Conditions - Auditory (A), Visual (V), Multisensory (M), Neutral (N)
models = {'noLapse','fixedLapse','inattention','exploration'};
nParams = [7,9,12,12];
subjects = {'lc40','lc39','lc44','lc59','lc60'};
fit = 0;
compare = ~fit;
plotSupplemental = 1;


clear bic aic nTrials nPars nll
for pooled = [1,0]
    for k = 1:length(models)
        model = models{k};
        if pooled
            subject = 'metaRat';
            j = find(strcmp(subject,{data.ratName}));
            % Pooled model fitting
            if fit
                data(j).controlSummaryData(2).condition = 'Multisensory';
                if strcmp(model,'explorationThompson')
                    [n,f,y] = jointFit(data(j).controlSummaryData(1:4), model,{'MA','MV','MN','AN','AV','VN'},'thompson_neutral');
                elseif strcmp(model,'exploration')
                    [n,f,y] = jointFit(data(j).controlSummaryData(1:4), model,{'MA','MV','MN','AN','AV','VN'},'exploration_neutral');
                else
                    [n,f,y] = jointFit(data(j).controlSummaryData(1:4), model,{'MA','MV','MN','AN','AV','VN'},'unc');
                end
                yfitsmeta.(model) = y;
                dataset.neutral(j).(model)=f;
            end
            % Pooled model comparison
            if compare
                [nll(k),bic(k),aic(k)] = modelMetrics(data(j),models{k},nParams(k));
            end
            
        else
            for i = 1:length(subjects)
                subject = subjects{i};
                j = find(strcmp(subject,{data.ratName}));
                % Individual model fitting
                if fit
                    data(j).controlSummaryData(2).condition = 'Multisensory';
                    if strcmp(model,'explorationThompson')
                        [n,f,y] = jointFit(data(j).controlSummaryData(1:4), model,{'MA','MV','MN','AN','AV','VN'},'thompson_neutral');
                    elseif strcmp(model,'exploration')
                        [n,f,y] = jointFit(data(j).controlSummaryData(1:4), model,{'MA','MV','MN','AN','AV','VN'},'exploration_neutral');
                    else
                        [n,f,y] = jointFit(data(j).controlSummaryData(1:4), model,{'MA','MV','MN','AN','AV','VN'},'neutral');
                    end
                    yfits.(model){i} = y;
                    dataset.neutral(j).(model)=f;
                end
                % Individual model comparison
                if compare
                    [nll(k,j),bic(k,j),aic(k,j)] = modelMetrics(data(j),models{k},nParams(k));
                end
            end
        end
    end
    bic = bic-min(bic);
    aic = aic-min(aic);
    
    % Plot pooled data
    figure(3)
    subplot(2,3,[1,2,4,5])
    mods = [2 3];
    colors = {[0.8 0 0],[0.8 0.5 0]};
    plotData(data,mods,colors)
    
    % Plot model comparison
    subplot(2,3,3*(2-pooled))
    plotModelComparison(bic,aic,models,pooled)
end

if plotSupplemental
    figure(301)
    subplot(2,3,1)
    imageModelComparison(bic,'BIC',models)
    subplot(2,3,2)
    imageModelComparison(aic,'AIC',models)
    subplot(2,3,3)
    plotModelComparison(bic,aic,models,pooled)
end
%% Reward manipulations - theoretical models
% Fit the ideal observer (O), fixed motor error (F), inattention (I) or
% uncertainty-guided exploration (E) models to reward manipulation data
ind = 0;
for manipulation = {'rewardInc','rewardDec','rewardProb'}
    ind = ind+1;
    models = {'noLapse','fixedLapse','inattention','exploration'};
    modelspec = {'idealObserverValue','fixedErrorValue','inattentionValue','explorationValueR'};
    nParams = [3,5,6,5];
    pooled = 1;
    fit = 0;
    compare = ~fit;
    modalities = {'Auditory','Multisensory','Visual'};
    modalitiesInac = {'AuditoryInac','MultisensoryInac','VisualInac'};
    modalityPairs = {'MA','MV','ma','mv','Aa','Mm','Vv'};
    
    
    switch manipulation{:}
        case 'rewardInc'
            % Increased reward magnitude manipulation
            % Conditions - Auditory (Equal rewards), AuditoryInac/HInc (rHigh>rLow)
            data = dataset.reward(1:4);
            mods = 1;
            j0 = 0;
            subjects = {'lc40','lc54','lc58'};
            if pooled
                subject = 'metaRatInc';
            end
        case 'rewardDec'
            % Decreased reward magnitude manipulation
            % Conditions - Auditory (Equal rewards), AuditoryInac/HDec (rHigh<rLow)
            data = dataset.reward(5:8);
            mods = 1;
            j0 = 4;
            subjects = {'sp11','sp12','dm01'};
            if pooled
                subject = 'metaRatDec';
            end
        case 'rewardProb'
            % Reward probability manipulation
            % Conditions - Visual (Reward probabilities of [1, 0; 0, 1]),
            % VisualInac/HL0.5 (Non-zero probability of reward for incorrect high
            % rate choices i.e. reward probabilities of [1, 0; 0.5, 1])
            mods = 1;
            data = dataset.reward(9:14);
            j0 = 8;
            subjects = {'sp10','sp09','sp12','sp13','sp14'};
            if pooled
                subject = 'metaRatProb';
            end
        case 'rewardIncExploit'
            % Increased reward magnitude manipulation
            % Conditions - Auditory (Equal rewards), AuditoryInac/HInc (rHigh>rLow)
            data = dataset.reward(15:18);
            mods = [1 2 3];
            j0 = 14;
            subjects = {'lc60','sp09','sp12'};
            if pooled
                subject = 'metaRatIncExploit';
            end
    end
    
    clear bic aic nTrials nPars nll
    for k = 1:length(models)
        if pooled
            j = find(strcmp(subject,{data.ratName}));
            fitParams = data(j).exploration;
            dataInac = [data(j).controlSummaryData(mods),data(j).skewSummaryData(mods)];
            [n,t,yFit] = jointNLLsummary(dataInac,'explorationReparam',reshape([[fitParams.pse]',[fitParams.sigma]',[fitParams.rLeft]',[fitParams.rRight]']',1,8),'valueR');
            
            % Pooled model fitting
            if fit
                for l = mods
                    m = find(strcmp({data(j).controlSummaryData.condition},modalities{l}));
                    data(j).skewSummaryData(m).condition = modalitiesInac{l};
                    dataInac = [data(j).controlSummaryData(m),data(j).skewSummaryData(m)];
                    [n,f,y] = jointFitInac(dataInac, modelspec{k},modalityPairs,'both');
                end
            end
            % Pooled model comparison
            if compare
                [nll(k),bic(k),aic(k)] = modelMetrics(data(j),models{k},nParams(k));
            end
            % Plot pooled data
            figure(4);
            subplot(2,2,ind)
            colors = {[0 0 0],[0 0 0],[0 0 0]};
            colorsSkew = {[0 0.5 0],[0.8 0 0],[0.8 0 0]};
            plotDataReward(data,yFit,subject,mods,colors{ind},colorsSkew{ind})
            title([manipulation,' (exploration model fits)'])
        else
            for i = 1:length(subjects)
                subject = subjects{i};
                j = find(strcmp(subject,{data.ratName}));
                % Individual model fitting
                if fit
                    % Across modality fit
                    for l = mods
                        m = find(strcmp({data(j).controlSummaryData.condition},modalities{l}));
                        data(j).skewSummaryData(m).condition = modalitiesInac{l};
                        modsData(l) = m;
                    end
                    dataInac = [data(j).controlSummaryData(m),data(j).skewSummaryData(m)];
                    [n,f,y] = jointFitInac(dataInac, modelspec{k},modalityPairs,'unc');
                    dataset.reward(j+j0).(models{k})=f;
                    fits(j+j0).(models{k}) = y;
                    
                end
                % Individual model comparison
                if compare
                    [nll(k,j),bic(k,j),aic(k,j)] = modelMetrics(data(j),models{k},nParams(k));
                end
            end
        end
    end
    if compare
        bic = bic-min(bic);
        aic = aic-min(aic);
        % Plot model comparison
        subplot(6,4,4*(ind+3))
        plotModelComparison(bic,aic,models,pooled)
        title(manipulation)
    end
    
end

%% Neural manipulations - theoretical models
data = dataset.inactivation;
% Fit the uncertainty-guided exploration model to pharmacological
% inactivation (muscimol) data and controls (saline), allowing a one-parameter
% change in either sensory evidence (S), value (V) or motor effort (E)

% Conditions - Visual, Auditory, Multisensory (Saline),
% VisualInac, AuditoryInac, MultisensoryInac (Muscimol)
% for side = {'high','low'}
for pooled = [1,0]
    ind = 0;
    for side = {'low','high'}
        switch side{:}
            case 'high'
                manipulations = {'M2Hi','pStrHi'};
                models = {'idealObserverValue','fixedErrorValue','inattentionValue','biasedEvidence','biasedValueL','biasedEffort'};
                modelspec = {'idealObserverValue','fixedErrorValue','inattentionValue','explorationEvidence','explorationValueL','explorationEffort'};
                
                
            case 'low'
                manipulations = {'M2Lo','pStrLo'};
                models = {'idealObserverValue','fixedErrorValue','inattentionValue','biasedEvidence','biasedValueR','biasedEffort'};
                modelspec = {'idealObserverValue','fixedErrorValue','inattentionValue','explorationEvidence','explorationValueR','explorationEffort'};
                
        end
        nParams = [12,12,12,11,11,11];
        fit = 0;
        compare = ~fit;
        modalities = {'Auditory','Multisensory','Visual'};
        modalitiesInac = {'AuditoryInac','MultisensoryInac','VisualInac'};
        modalityPairs = {'MA','MV','ma','mv','Aa','Mm','Vv'};
        
        for a = 1:length(manipulations)
            ind = ind+1;
            manipulation = manipulations{a};
            switch manipulation
                case 'M2Hi'
                    data = dataset.inactivation(1:4);
                    subjects = {'lc50','lc54','lc58'};
                    if pooled
                        subject = 'metaRat FOF - High';
                    end
                    j0 = 0;
                    mods = [1 2 3];
                case 'M2Lo'
                    data = dataset.inactivation(5:10);
                    subjects = {'lc47','lc48','lc50','lc54','lc58'};
                    if pooled
                        subject = 'metaRat FOF - Low';
                    end
                    j0 = 4;
                    mods = [1 2 3];
                case 'pStrHi'
                    data = dataset.inactivation(11:17);
                    subjects = {'lc53','lc56','lc57','lc45','lc46','lc55'};
                    if pooled
                        subject = 'metaRat pStr - High';
                    end
                    j0 = 10;
                    mods = [1 2 3];
                case 'pStrLo'
                    data = dataset.inactivation(18:23);
                    subjects = {'lc53','lc57','lc45','lc46','lc55'};
                    if pooled
                        subject = 'metaRat pStr - Low';
                    end
                    j0 = 17;
                    mods = [1 2 3];
            end
            clear bic aic nTrials nPars nll
            for k = 1:length(models)
                if pooled
                    j = find(strcmp(subject,{data.ratName}));
                    if ind == 1 || ind == 2
                        fitParams = data(j).biasedValueR;
                        dataInac = [data(j).controlSummaryData(mods),data(j).inactivationSummaryData(mods)];
                        [n,t,yFit] = jointNLLsummary(dataInac,'explorationInac',reshape([[fitParams.pse]',[fitParams.sigma]',[fitParams.rLeft]',[fitParams.rRight]',[fitParams.kR]']',1,30),'valueR');
                    elseif ind == 3 || ind == 4
                        fitParams = data(j).biasedValueL;
                        dataInac = [data(j).controlSummaryData(mods),data(j).inactivationSummaryData(mods)];
                        [n,t,yFit] = jointNLLsummary(dataInac,'explorationInac',reshape([[fitParams.pse]',[fitParams.sigma]',[fitParams.rLeft]',[fitParams.rRight]',[fitParams.kL]']',1,30),'valueL');
                        
                    end
                    % Pooled model fitting
                    if fit
                        for l = mods
                            m = find(strcmp({data(j).controlSummaryData.condition},modalities{l}));
                            data(j).inactivationSummaryData(m).condition = modalitiesInac{l};
                        end
                        dataInac = [data(j).controlSummaryData(mods),data(j).inactivationSummaryData(mods)];
                        [n,f,y] = jointFitInac(dataInac, modelspec{k},modalityPairs,'inac');
                        dataset.inactivation(j+j0).(models{k})=f;
                        fits(j+j0).(models{k}) = y;
                    end
                    % Pooled model comparison
                    %                 if compare
                    %                     [nll(k),bic(k),aic(k)] = modelMetrics(data(j),models{k},nParams(k));
                    %                 end
                    
                else
                    for i = 1:length(subjects)
                        subject = subjects{i};
                        j = find(strcmp(subject,{data.ratName}));
                        % Individual model fitting
                        if fit
                            for l = mods
                                m = find(strcmp({data(j).controlSummaryData.condition},modalities{l}));
                                data(j).inactivationSummaryData(m).condition = modalitiesInac{l};
                            end
                            dataInac = [data(j).controlSummaryData(mods),data(j).inactivationSummaryData(mods)];
                            [n,f,y] = jointFitInac(dataInac, modelspec{k},modalityPairs,'inac');
                            dataset.inactivation(j+j0).(models{k})=f;
                            fits(j+j0).(models{k}) = y;
                        end
                        % Individual model comparison
                        if compare
                            [nll(k,j),bic(k,j),aic(k,j)] = modelMetrics(data(j),models{k},nParams(k));
                            BIC(k,j+j0) = bic(k,j);
                            AIC(k,j+j0) = aic(k,j);
                        end
                    end
                end
            end
            if compare && ~pooled
                %             bicAll{a} = bic-min(bic);
                %             aicAll{a} = aic-min(aic);
                BIC = BIC-min(BIC);
                AIC = AIC-min(AIC);
            end
            figure(5);
            
            if pooled
                for m = 1:3
                    subplot(4,4,m+4*(ind-1))
                    plotDataInactivation(data(end),yFit([m,m+3],:),m)
                    if ind ~=4
                        xlabel('')
                    end
                    if ind == 1 && m == 2
                        title('Pooled data (explorationBiasedValue fits)')
                    end
                end
            end
        end
    end
end
subplot(4,4,4)
plotModelComparison(BIC(4:6,1:10),AIC(4:6,1:10),models(4:6),pooled)
title('M2 (individual fits)')

subplot(4,4,12)
plotModelComparison(BIC(4:6,11:22),AIC(4:6,11:22),models(4:6),pooled)
title('pStr (individual fits)')


%% Model comparison code

function [nll,bic,aic]=modelMetrics(data,model,nPars)
if isfield(data,'inactivationRawData')
    nTrials = sum([data.controlSummaryData.nTrials])+sum([data.inactivationSummaryData.nTrials]);
elseif isfield(data,'skewRawData')
    nTrials = sum([data.controlSummaryData.nTrials])+sum([data.skewSummaryData.nTrials]);
else
    nTrials = sum(data.controlRawData.completed);
end
nll = sum([data.(model)(:).NLL]);
bic = 2*nll+nPars*log(nTrials);
aic = 2*nll+nPars*2;
end

%% Plot model comparison
function plotModelComparison(bic,aic,models,pooled)
if pooled
    b = bar([bic',aic']);
    title('Fit to pooled data')
else
    b = bar([sum(bic,2),sum(aic,2)]);
    title('Fit to individual data')
end
b(1).FaceColor = 'flat';
b(1).CData = [0.8 0.6 0.6];
b(2).FaceColor = 'flat';
b(2).CData = [0.2 0.2 0.5];
legend({'BIC','AIC'},'box','off')
set(gca,'FontSize',8)
set(gca,'XTickLabel',models)
set(gca,'XTickLabelRotation',45)
ylabel('\Delta Info criterion')
set(gcf,'color','w')
end

function imageModelComparison(ic,type,models)
imagesc(ic)
colormap('pink')
title(type)
set(gca,'FontSize',8)
set(gca,'YTick',1:length(models))
set(gca,'YTickLabel',models)
end

%% Plot psychometric curves
% Control
function hdata=plotData(data,mods,colors)
for m = 1:length(mods)
    l = mods(m);
    ind = find(strcmp('metaRat',{data.ratName}));
    [hi,p,s] = fitPalamedesReparam(data(ind).controlSummaryData(l).stimRates, data(ind).controlSummaryData(l).nHighResponses, data(ind).controlSummaryData(l).nTrials,1);
    set(hi.fit,'color',colors{m},'LineWidth',2)
    set(hi.dataCI,'color',colors{m},'LineWidth',1)
    set(hi.data,'MarkerFaceColor',colors{m},'MarkerEdgeColor',colors{m})
    set(gca,'FontSize',15)
    hdata(m) = hi;
end
axprefs(gca,18,'pmf')
title('pooled data - descriptive fits')
end
% Reward manipulation
function plotDataReward(data,yFit,subject,mods,color,colorSkew)
for m = 1:length(mods)
    l = mods(m);
    ind = find(strcmp(subject,{data.ratName}));
    [hc,p,s] = fitPalamedesReparam(data(ind).controlSummaryData(l).stimRates, data(ind).controlSummaryData(l).nHighResponses, data(ind).controlSummaryData(l).nTrials,1);
    [hi,p,s] = fitPalamedesReparam(data(ind).skewSummaryData(l).stimRates, data(ind).skewSummaryData(l).nHighResponses, data(ind).skewSummaryData(l).nTrials,1);
    delete(hi.fit)
    delete(hc.fit)
    hold on
    xVal = 5:0.1:20;
    plot(xVal,yFit(1,:),'color',color,'Linewidth',2)
    plot(xVal,yFit(2,:),'color',colorSkew,'Linewidth',2)
    set(hc.dataCI,'color',color,'LineWidth',1)
    set(hi.dataCI,'color',colorSkew,'LineWidth',1)
    set(hc.data,'MarkerFaceColor',color,'MarkerEdgeColor',color)
    set(hi.data,'MarkerFaceColor',colorSkew,'MarkerEdgeColor',colorSkew)
    set(gca,'FontSize',8)
    hdata(m) = hi;
end
axprefs(gca,12,'pmf')
end
% Neural manipulation
function [h,p,s] = plotDataInactivation(data,yFit,m)
conditions = {'Auditory','Multisensory','Visual'};
colors = {[0 0.5 0], [0.8 0 0], [0 0 0.8]};
linestyles = {'-','--'};
linewidth = 1;
ciwidth = 0.5;
markersize = 20;
fontSize = 8;
xVal = 5:0.1:20;
ctrl = data.controlSummaryData(m);
inac = data.inactivationSummaryData(m);
[h.ctrl,p.ctrl,s.ctrl] = fitPalamedesReparam(ctrl.stimRates,ctrl.nHighResponses, ctrl.nTrials,1);
[h.inac,p.inac,s.inac] = fitPalamedesReparam(inac.stimRates,inac.nHighResponses, inac.nTrials,1);
delete(h.ctrl.fit)
delete(h.inac.fit)
hold on
plot(xVal,yFit(1,:),'color',colors{m},'Linewidth',2,'LineStyle',linestyles{1})
plot(xVal,yFit(2,:),'color',colors{m},'Linewidth',2,'LineStyle',linestyles{2})
set(h.ctrl.data,'SizeData',markersize,'MarkerFaceColor',colors{m},'MarkerEdgeColor',colors{m});
set(h.ctrl.dataCI,'color',colors{m},'LineWidth',ciwidth);
set(h.inac.data,'SizeData',markersize,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colors{m});
set(h.inac.dataCI,'color',colors{m},'LineWidth',ciwidth);
if m == 1
    title([data.ratName])
end
axprefs(gca,fontSize);
ylim([0,1]);
xlim([8,17]);
set(gca,'XTick',[9,16])
set(gca,'YTick',[0,1])
ylabel('p(Chose High)')
end
