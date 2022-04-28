clear
clc

load dates
load gibbs_filled
load factors_dnmv

% Determine the sample
s = find(dates==197201);
e = find(dates==202112);

% Store a few constants
nModels = length(factor_model_defs);
nStocks = size(tcosts,2);

% Initialize a cell for the netting results
netting_sharpe = cell(nModels,1);
netting_weights = cell(nModels,1);

% Subset the tcosts 
tcosts = tcosts(s:e, :);

for i=1:nModels
    
    % Determine the model and number of factors in this iteration
    thisModel = factor_model_defs(i).label;
    nFactors = length(factor_model_defs(i).factors);    
    
    % Timekeeping
    fprintf('Now working on %s.\n',char(thisModel));
    tic;
    
    % Get the basic results for this factor model (weights, sharpe ratios,
    % etc.) that do not account for trading diversification
    res(i,1) = getFactorModelData(thisModel, factor_model_defs, factor_struct, [s e]);       
    
    % Optimization function inputs (subset here to speed up optimization)
    % Store the gross factor returns
    gross_factors = res(i).gross_factors(s:e,:);
    % Store the change in weight matrix for the non-market factors. We
    % assume that the market is free to trade (e.g., a market ETF).
    dW = res(i).dW(s:e,:,2:end);
    % Reshape it to make it a 2-d array
    dW = dW(:,:);
    
    % Start with the net weights that do not account for trading
    % diversification
    w0 = res(i).net_weights';
    
    % Get the optimal Sharpe and weights with trading diversification    
    % Get the default optimization options for fmincon
    opts = optimoptions('fmincon','Display','off');    
    % SharpeTD calculates a Sharpe ratio with trading diversification. Fun
    % will only take the weights and pass the other inputs to SharpeTD
    fun = @(w) -SharpeTD(w,gross_factors,dW,tcosts);
    % Find the optimal weights and Sharpe ratio
    [wStar, SRStar] = fmincon(fun,w0,[],[],ones(1,nFactors),1,zeros(nFactors,1),[],[],opts);
    
    % Store the weights and Sharpe ratio
    netting_weights{i} = wStar;
    netting_sharpe{i} = -sqrt(12)*SRStar;
    toc;
end

% Assign the netting results to the struct
[res.netting_sharpe] = netting_sharpe{:};
[res.netting_weights] = netting_weights{:};

% Drop the dW field from res (we don't need it and it takes space)
res = rmfield(res,'dW'); 

save Results/fullSampleResults res 
