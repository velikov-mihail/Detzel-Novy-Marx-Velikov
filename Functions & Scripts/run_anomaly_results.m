%% Create the return series for the 205 anomalies in Chen and Zimmerman

clear
clc

% Timekeeping
fprintf('Now working on generating anomaly results. Run started at %s.\n', char(datetime('now')));
tic;

% Load the variables we need
load dates
load ret
load gibbs_filled
load NYSE
load prc
load exchcd
load me

% make sure you download:
% https://drive.google.com/file/d/1-1RUq2wUADu_ncvQJCYY3wxDhqhHxBow/view
% https://docs.google.com/spreadsheets/d/18DvZPscKsD0_ZeeUMjyXhF1qn0emDVaj/edit#gid=70837236
[anoms, labels, anomaly_summary] = getChenZimmermanAnomalies();

% Store the names of the signals
nAnoms = height(anomaly_summary);
signal_names = anomaly_summary.Acronym;

% Store the rebalancing indicators
JuneIndicator = (dates-100*floor(dates/100)==6);
JuneDecemberIndicator = (dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==12);
quarterEndIndicator = (dates-100*floor(dates/100)==3 | dates-100*floor(dates/100)==6 | dates-100*floor(dates/100)==9 | dates-100*floor(dates/100)==12);

for i=1:nAnoms
    
    % Assign the current anomaly signal to a matrix called 'var'
    r = find(strcmp(labels, anomaly_summary.Acronym(i)));
    var = anoms(:,:,r);
        
    % Figure out the rebalancing period
    if anomaly_summary.PortfolioPeriod(i) == 12 || anomaly_summary.PortfolioPeriod(i) == 36
        var(~JuneIndicator, :) = nan;
    elseif anomaly_summary.PortfolioPeriod(i) == 6
        var(~JuneDecemberIndicator, :) = nan;
    elseif anomaly_summary.PortfolioPeriod(i) == 3
        var(~quarterEndIndicator, :) = nan;
    end
       
    % Check for sample filters
    if ~strcmp(anomaly_summary.Filter(i),'')
        filterString = regexprep(anomaly_summary.Filter(i),' ','');
        if contains(filterString, 'abs(prc)>1')
            var(abs(prc)<=1) = nan;  
        end
        if contains(filterString, 'abs(prc)>5')
            var(abs(prc)<=5) = nan;  
        end
        if contains(filterString, 'exchcd%in%c(1,2)')
            var(exchcd>2) = nan;
        end
        if contains(filterString, 'exchcd==1')
            var(exchcd>1) = nan;
        end
        if contains(filterString, 'me>me_nyse20')
            indME5 = makeUnivSortInd(me, 5, NYSE);
            var(indME5==1) = nan;
        end
    end
                        
    % Determine the sample period
    if isnan(anomaly_summary.StartMonth(i))
        sampleStart = anomaly_summary.SampleStartYear(i)*100 + 6;    
    else
        sampleStart = anomaly_summary.SampleStartYear(i)*100 + anomaly_summary.StartMonth(i);
    end
    sampleStart = max(sampleStart, dates(1));
    sampleEnd = dates(min(find(sum(isfinite(var),2)>0, 1, 'last') + anomaly_summary.PortfolioPeriod(i), length(dates)));
    dts = [sampleStart sampleEnd];
    
    % Determine the sorting breakpoints 
    if strcmp(anomaly_summary.Acronym(i), 'ChangeInRecommendation') || strcmp(anomaly_summary.Acronym(i), 'NumEarnIncrease')
        bpts = [];
        for j = 1:4
            bpts = [bpts j*100/5];
        end
        bpts = prctile(var, bpts, 2);        
        ind = 1*(var <= repmat(bpts(:,1), 1, size(var,2))) + 2*(var >= repmat(bpts(:,end), 1, size(var,2)));        
    elseif strcmp(anomaly_summary.Acronym(i), 'RDcap')
        ind = 1*(var==0) + 2*(var>0);
    elseif strcmp(anomaly_summary.Cat_Form(i), 'discrete')       
        uVals = unique(var);
        uVals(isnan(uVals)) = [];
        ind = 1*(var==min(uVals)) + 2*(var==max(uVals));
    else
        if isnan(anomaly_summary.LSQuantile(i))
            nptfs = 5;
        else
            nptfs = round(1/anomaly_summary.LSQuantile(i));
        end
        ind = makeUnivSortInd(var, nptfs);  
        if strcmp(anomaly_summary.QuantileFilter(i), 'NYSE')
            ind = makeUnivSortInd(var, nptfs, NYSE);  
        end
        ind = 1*(ind==1) + 2*(ind==nptfs);
    end
    
    % Use the weighting from the original publication
    if strcmp(anomaly_summary.StockWeight(i), 'VW')
        resAnoms(i,1) = runUnivSort(ret, ind, dates, me, 'tcosts', tcosts, ...
                                                         'weighting', 'v', ...
                                                         'timePeriod', dts, ...
                                                         'plotFigure', 0, ...
                                                         'printResults', 0);  
    else
        resAnoms(i,1) = runUnivSort(ret, ind, dates, me, 'tcosts', tcosts, ...
                                                         'weighting', 'e', ...
                                                         'timePeriod', dts, ...
                                                         'plotFigure', 0, ...
                                                         'printResults', 0);      
    end
    fprintf('\n\n\nDone with %s, which is %d out of %d.\n\n\n', char(labels(r)), i, size(anoms,3));
end

fprintf('Done running sorts at at %s.\n', char(datetime('now')));

% Store the publication dates
for i=1:nAnoms
    r = find(strcmp(anomaly_summary.Acronym, labels(i)));
    pubDates(i,1) = 100*anomaly_summary.Year(r) + 12;
end

save Results/anomReturnSeries resAnoms anomaly_summary pubDates

%%  Extract anomaly returns for the full available sample

clear
clc

fprintf('Now working on extracting anomaly returns results for the full sample. Run started at %s.\n', char(datetime('now')));

% Load the necessary workspace
load factors_dnmv
load anomReturnSeries

% Get the number of anomalies and initialize the gross and net return
% matrices
nAnoms = length(resAnoms);
gross_rets = nan(length(dates),length(resAnoms));
net_rets = nan(length(dates),length(resAnoms));

% Store the gross and net returns for the anomalies
for i=1:nAnoms
    s = find(dates==resAnoms(i).dates(1));
    e = find(dates==resAnoms(i).dates(end));
    gross_rets(s:e,i) = resAnoms(i).pret(:,end);
    net_rets(s:e,i) = resAnoms(i).netpret;
end

% Choose the sample period
s = find(dates==197201);
e = length(dates);

% Get the anomaly tcosts
tc = gross_rets - net_rets;

% Filter the factor model definitions structure
ind = ismember([factor_model_defs.label],{'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'});
factor_model_defs = factor_model_defs(ind);

% Get a structure with the factor model data (gross, net, tc's)
nFactorModels = length(factor_model_defs);
for i=1:nFactorModels
    res(i,1) = getFactorModelData(factor_model_defs(i).label, factor_model_defs, factor_struct, [s e]);
end

% Initialize the matrices we'll create
weightsStruct = struct;
gross_sharpe = nan(nAnoms, nFactorModels);
net_sharpe = nan(nAnoms, nFactorModels);

% Initialize the factor model sharpe vectors
factor_gross_sharpe = nan(1, nFactorModels);
factor_net_sharpe = nan(1, nFactorModels);

% Calculate the factor model Sharpes for 197201 - 202012
for j=1:nFactorModels      
    
    % Store the returns and trading costs;
    rets  = [res(j).gross_factors];
    tcost = [res(j).tc];
    
    % Get the gross annualized Sharpe ratio first    
    [~, tempGrossSharpe] = calcMve(rets(s:e,:));    
    factor_gross_sharpe(1,j) = sqrt(12) * tempGrossSharpe;

    % Get the net annualized Sharpe ratio next
    [~, tempNetSharpe] = calcNetMve(rets(s:e,:) - tcost(s:e,:), ...
                                   -rets(s:e,:) - tcost(s:e,:));      
    factor_net_sharpe(1,j) = sqrt(12) * tempNetSharpe;                                
end

% Create matrices
factor_gross_sharpe = repmat(factor_gross_sharpe, nAnoms, 1);
factor_net_sharpe = repmat(factor_net_sharpe, nAnoms, 1);

% Loop through the anomalies
for i=1:nAnoms
    
    % Check the first available gross return date
    b = find(dates>=197201 & isfinite(gross_rets(:,i)), 1, 'first');
    e = find(dates>=197201 & isfinite(gross_rets(:,i)), 1, 'last');
    
    % Loop through the factor models
    for j=1:nFactorModels
        
        % If start date is past 197201 we need to reestimate the factor
        % model Sharpe ratios.
        if b>s || e<length(dates)

            % Store the returns and trading costs;
            rets  = [res(j).gross_factors];
            tcost = [res(j).tc];

            % Get & store the gross annualized Sharpe ratio first    
            [~, tempGrossSharpe] = calcMve(rets(b:e,:));    
            factor_gross_sharpe(i,j) = sqrt(12) * tempGrossSharpe;

            % Get & store the net annualized Sharpe ratio next
            [~, tempNetSharpe] = calcNetMve(rets(b:e,:) - tcost(b:e,:), ...
                                           -rets(b:e,:) - tcost(b:e,:));      
            factor_net_sharpe(i,j) = sqrt(12) * tempNetSharpe;                                                                    
        end

        % Estimate the gross shape including the anomaly return
        % Store the returns and trading costs;
        rets  = [res(j).gross_factors gross_rets(:,i)];
        tcost = [res(j).tc tc(:,i)];                
        
        % Calculate and store the gross weights and annualized Sharpe
        [tempWeights, tempGrossSharpe] = calcMve(rets(b:e,:));    
        weightsStruct(i,j).gross_weights = tempWeights;
        gross_sharpe(i,j) = sqrt(12) * tempGrossSharpe;

        % Calculate and store the net weights and annualized Sharpe
        [tempWeights, tempNetSharpe] = calcNetMve(rets(b:e,:)-tcost(b:e,:), ...
                                                 -rets(b:e,:)-tcost(b:e,:));    
        weightsStruct(i,j).net_weights = tempWeights;
        net_sharpe(i,j) = sqrt(12) * tempNetSharpe;     
    end        
    fprintf('\n Done with %d/%d anomalies @ %s.\n',i,nAnoms,char(datetime('now')));
end

save Results/anomalyResults net_sharpe gross_sharpe weightsStruct factor_gross_sharpe factor_net_sharpe 


%% Extract anomaly returns for the post-publication sample

clear
clc

fprintf('Now working on extracting anomaly returns results for the post-pub sample. Run started at %s.\n', char(datetime('now')));

% Load the necessary workspace
load factors_dnmv
load anomReturnSeries

% Get the number of anomalies and initialize the gross and net return
% matrices
nAnoms = length(resAnoms);
gross_rets = nan(length(dates),length(resAnoms));
net_rets = nan(length(dates),length(resAnoms));

% Store the gross and net returns for the anomalies
for i=1:nAnoms
    s = find(dates==resAnoms(i).dates(1));
    e = find(dates==resAnoms(i).dates(end));
    gross_rets(s:e,i) = resAnoms(i).pret(:,end);
    net_rets(s:e,i) = resAnoms(i).netpret;
end

% Choose the sample period
s = find(dates==197201);
e = length(dates);

% Get the anomaly tcosts
tc = gross_rets - net_rets;

% Filter the factor model definitions structure
ind = ismember([factor_model_defs.label], {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'});
factor_model_defs = factor_model_defs(ind);

% Get a structure with the factor model data (gross, net, tc's)
nFactorModels = length(factor_model_defs);
for i=1:nFactorModels
    res(i,1) = getFactorModelData(factor_model_defs(i).label, factor_model_defs, factor_struct, [s e]);
end

% Initialize the matrices we'll create
weightsStruct = struct;
gross_sharpe = nan(nAnoms, nFactorModels);
net_sharpe = nan(nAnoms, nFactorModels);

% Initialize the factor model sharpe vectors
factor_gross_sharpe = nan(1, nFactorModels);
factor_net_sharpe = nan(1, nFactorModels);

% Calculate the factor model Sharpes for 197201 - 202012
for j=1:nFactorModels      
    
    % Store the returns and trading costs;
    rets  = [res(j).gross_factors];
    tcost = [res(j).tc];
    
    % Get the gross annualized Sharpe ratio first    
    [~, tempGrossSharpe] = calcMve(rets(s:e,:));    
    factor_gross_sharpe(1,j) = sqrt(12) * tempGrossSharpe;

    % Get the net annualized Sharpe ratio next
    [~, tempNetSharpe] = calcNetMve(rets(s:e,:) - tcost(s:e,:), ...
                                   -rets(s:e,:) - tcost(s:e,:));      
    factor_net_sharpe(1,j) = sqrt(12) * tempNetSharpe;                                
end

% Create matrices
factor_gross_sharpe = repmat(factor_gross_sharpe, nAnoms, 1);
factor_net_sharpe = repmat(factor_net_sharpe, nAnoms, 1);

% Loop through the anomalies
for i=1:nAnoms
    
    % Check the first available gross return date
    b = find(dates>=197201 & dates>=pubDates(i) & isfinite(gross_rets(:,i)), 1, 'first');
    e = find(dates>=197201 & dates>=pubDates(i) & isfinite(gross_rets(:,i)), 1, 'last');
    
    if ~isempty(b) && ~isempty(e)
        % Loop through the factor models
        for j=1:nFactorModels

            % Reestimate the factor model sharpe ratios with the new dates
            % Store the returns and trading costs;
            rets  = [res(j).gross_factors];
            tcost = [res(j).tc];
            
            % Get & store the gross annualized Sharpe ratio first    
            [~, tempGrossSharpe] = calcMve(rets(b:e,:));    
            factor_gross_sharpe(i,j) = sqrt(12) * tempGrossSharpe;

            % Get & store the net annualized Sharpe ratio next
            [~, tempNetSharpe] = calcNetMve(rets(b:e,:) - tcost(b:e,:), ...
                                           -rets(b:e,:) - tcost(b:e,:));      
            factor_net_sharpe(i,j) = sqrt(12) * tempNetSharpe;                                                                    

            % Estimate the gross shape including the anomaly return
            % Store the returns and trading costs;
            rets  = [res(j).gross_factors gross_rets(:,i)];
            tcost = [res(j).tc tc(:,i)];                

            % Calculate and store the gross weights and annualized Sharpe
            [tempWeights, tempGrossSharpe] = calcMve(rets(b:e,:));    
            weightsStruct(i,j).gross_weights = tempWeights;
            gross_sharpe(i,j) = sqrt(12) * tempGrossSharpe;

            % Calculate and store the net weights and annualized Sharpe
            [tempWeights, tempNetSharpe] = calcNetMve(rets(b:e,:)-tcost(b:e,:), ...
                                                     -rets(b:e,:)-tcost(b:e,:));    
            weightsStruct(i,j).net_weights = tempWeights;
            net_sharpe(i,j) = sqrt(12) * tempNetSharpe;     

        end        
    end
    fprintf('\n Done with %d/%d anomalies @ %s.\n', i, nAnoms, char(datetime('now')));
end

save Results/anomalyResultsPostPub net_sharpe gross_sharpe weightsStruct factor_gross_sharpe factor_net_sharpe 

