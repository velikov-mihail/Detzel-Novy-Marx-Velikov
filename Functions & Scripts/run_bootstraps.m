

%% Basic bootstraps

clear
clc

fprintf('Now working on basic bootstraps. Run started at %s.\n', char(datetime('now')));

load factors_dnmv

% Determine the sample
s = find(dates == 197201);
e = find(dates == 202112);

% Set a few constants
nMonths = e-s+1; % number of months (should be an even number)
nRuns = 1e5; % number of runs 
nFactorModels = length(factor_model_defs); % Number of factor models
pairStartMonths = [s:2:e]';

% Set seed to default
rng('default')

% Create a random binary matrix to pick IS/OS from the pairs
randomBinaryInd = randi([0 1], nMonths/2, nRuns);

% Create a random index to pick the pairs with replacement
randomPairsIndicator = randi([1 nMonths/2], nMonths/2, nRuns);

% Initialize the in- and out-of-sample indices
inSampleInd = nan(nMonths/2, nRuns);
outOfSampleInd = nan(nMonths/2, nRuns);

% Create the index for the samples following FF bootstraps
for i=1:nRuns 
    
    % Assign the months within each pair to IS/OS randomly
    isMonths = pairStartMonths + randomBinaryInd(:,i);
    oosMonths = pairStartMonths + 1 - randomBinaryInd(:,i);
    
    % Pick 
    inSampleInd(:,i) = isMonths(randomPairsIndicator(:,i));
    outOfSampleInd(:,i) = oosMonths(randomPairsIndicator(:,i));
end
index = [inSampleInd; outOfSampleInd]-s+1;

% Get the factor model structure (contains gross rets and tcosts for factors
% in each factor model)
for i=1:nFactorModels
    factor_model_struct(i) = getFactorModelData(factor_model_defs(i).label, factor_model_defs, factor_struct, [s e]);
end

%% Run the bootstrap
bootstrap_results = struct;

for i = 1:nFactorModels
    % Timekeeping
    tic;
    fprintf('Now working on %7s factor model, which is %d/%d.\n', char(factor_model_struct(i).label), i, nFactorModels);
 
    % Take the long and short net factors
    net_factors_long = factor_model_struct(i).gross_factors(s:e,:) - ...
                       factor_model_struct(i).tc(s:e,:);
    net_factors_short = - factor_model_struct(i).gross_factors(s:e,:) - ...
                          factor_model_struct(i).tc(s:e,:);
    
    % Run the bootstrap for this particular factor model
    resThisBootstrap = runSingleFactorModelBootstrap(net_factors_long, net_factors_short, index, nRuns); % 

    % Assign the bootstrap results
    bootstrap_results(i).weights = resThisBootstrap.weights;
    bootstrap_results(i).label = factor_model_struct(i).label;
    bootstrap_results(i).SR_is = sqrt(12) * resThisBootstrap.SR_is;
    bootstrap_results(i).SR_os = sqrt(12) * resThisBootstrap.SR_os;
    toc;
end

save Results/bootstrap_results bootstrap_results

% CAPM bootstrap
mkt = factor_struct(1).gross_factor;

capmIsSR = nan(nRuns,1);
capmOosSR = nan(nRuns,1);
for i=1:nRuns
    % store the in- and out-of-sample vectors
    mktIS = mkt(inSampleInd(:,i));
    mktOS = mkt(outOfSampleInd(:,i));
    
    % Calculate the in- and out-of-sample Sharpe ratios
    capmIsSR(i) = sqrt(12) * mean(mktIS) / std(mktIS);
    capmOosSR(i) = sqrt(12) * mean(mktOS) / std(mktOS);
end

save Results/bootstrap_capm capmIsSR capmOosSR

fprintf('Basic bootstraps done at %s.\n',char(datetime('now')));

%% Bootstraps with netting - new

clear
clc

fprintf('Now working on bootstraps with netting. Run started at %s.\n',char(datetime('now')));

load factors_dnmv
load gibbs_filled

% Determine the sample
s = find(dates==197201);
e = find(dates==202112);

% Set a few constants
nRuns = 2e3; % Number of runs. This takes a while, recommended to start with 1000 to evaluate how long it takes.
nMonths = e-s+1; % number of months (should be an even number)
nModels = length(factor_model_defs); % Number of factor models
nStocks = size(tcosts,2);
pairStartMonths = [s:2:e]';
nFmax = 6; % Max number of factors to determine W 3-d array size

% Set seed to default
rng('default')

% Create a random binary matrix to pick IS/OS from the pairs
randomBinaryInd = randi([0 1], nMonths/2, nRuns);

% Create a random index to pick the pairs with replacement
randomPairsIndicator = randi([1 nMonths/2], nMonths/2, nRuns);

% Create the index for the samples following FF bootstraps
for i=1:nRuns 
    % Assign the months within each pair to IS/OS randomly
    isMonths = pairStartMonths + randomBinaryInd(:,i);
    oosMonths = pairStartMonths + 1 - randomBinaryInd(:,i);

    % Pick 
    inSampleInd(:,i) = isMonths(randomPairsIndicator(:,i));
    outOfSampleInd(:,i) = oosMonths(randomPairsIndicator(:,i));
end

% Get the full sample factor and factor model results 
for i=1:nModels
   thisModel = factor_model_defs(i).label;
   res(i,1) = getFactorModelData(thisModel, factor_model_defs, factor_struct, [s e]);       
end

% Initialize a few variables
W         = nan(nRuns, nFmax, nModels);
ISSharpes = nan(nRuns, nModels);
OSSharpes = nan(nRuns, nModels);
errorFlag = zeros(nRuns, nModels);

clear factor_model_defs
clear factor_struct

for m=1:nModels
    
    % Store this model results structure
    thisModelRes = res(m);
    
    fprintf('Now working on %s. Run started at %s.\n',char(thisModelRes.label),char(datetime('now')));
    
    % Number of factors
    nF = length(thisModelRes.factor_labels);
    
    % Check if one of FF6 models that nests five-factor models
    if ismember(thisModelRes.label,{'FF6','FF6_c','FF6_sS','FF_c_sS'})
        
        % Find the same five-factor model index
        rFF5        = find( strcmp([res.label], regexprep(thisModelRes.label,'6','5')) );
        % Store the results for the five-factor model
        ff5ModelRes = res(rFF5);
        % Get the number of factors for the five-factor model
        nF_ff5      = length(ff5ModelRes.factor_labels);
        % Get the five-factor model optimal weights
        w0_ff5      = W(:, 1:nF_ff5, rFF5)';
        
        % Initite the six-factor model starting weights
        w0  = nan(nF, nRuns);
        % Find the index for the five factors in the six-factor model
        indFF5 = ismember(thisModelRes.factor_labels, ff5ModelRes.factor_labels);        
        
        % Assign the five-factor model weights to the six-factor starting
        % weights
        w0(indFF5, :)  = w0_ff5;
        % Assign zero for momentum
        w0(~indFF5, :) = 0;
    else
        % If it's not one of the six-factor models, assign the full-sample
        % weights
        w0 = repmat(thisModelRes.net_weights', 1, nRuns);
    end
    
    % Store the gross factor returns
    gross_factors = thisModelRes.gross_factors;
    % Store the change in weight matrix for the non-market factors. We
    % assume that the market is free to trade (e.g., a market ETF).    
    dW = thisModelRes.dW(:, :, 2:end);
    % Reshape it to make it a 2-d array
    dW = dW(:, :);
    
    % Initialize the weights matrix for this model
    thisModelW = nan(nRuns, nF);    
        
    % Parallel loop through the number of rooms 
    parfor i=1:nRuns
        
       % Split the tcosts
       tcostsIS = tcosts(inSampleInd(:,i), :);
       tcostsOS = tcosts(outOfSampleInd(:,i), :);
        
       % Optimization function inputs (subset here to speed up optimization)
       GrossIS = gross_factors(inSampleInd(:,i), :);
       dWFacsIS = dW(inSampleInd(:,i), :);

       % Get the optimal Sharpe and weights with trading diversification    
       % Get the default optimization options for fmincon
       opts = optimoptions('fmincon', 'Display', 'off');
       % SharpeTD calculates a Sharpe ratio with trading diversification. Fun
       % will only take the weights and pass the other inputs to SharpeTD
       fun = @(w) -SharpeTD(w, GrossIS, dWFacsIS, tcostsIS);
       
       % Allow for errors if something funky happens in the optimization
       try 
           % Find the optimal weights and Sharpe ratio           
           [wStar, SRStar] = fmincon(fun,w0(:,i),[],[],ones(1,nF),1,[],[],[],opts);
           
           % Store the optimal weights
           thisModelW(i,:) = wStar;
           
           % Store the Sharpe ratios
           ISSharpes(i,m) = -(SRStar);
           OSSharpes(i,m) = SharpeTD(wStar, gross_factors(outOfSampleInd(:,i),:), dW(outOfSampleInd(:,i),:), tcostsOS);
       catch
           errorFlag(i,m) = 1;
       end

    end
    W(:, 1:nF, m) = thisModelW;
    save results/bootstrap_netting ISSharpes OSSharpes W errorFlag 
end


% Annualize the Sharpes
ISSharpes = sqrt(12)*ISSharpes;
OSSharpes = sqrt(12)*OSSharpes;

save results/bootstrap_netting ISSharpes OSSharpes W errorFlag 

fprintf('Netting bootstraps done at %s.\n',char(datetime('now')));

