clear
clc

fprintf('Now working on out-of-sample MVE strategy results. Run started at %s.\n', char(datetime('now')));

load factors_dnmv
load gibbs_filled

% Determine the sample
s = find(dates==197201);
e = find(dates==202112);

% Set a few constants
nModels = length(factor_model_defs); % Number of factor models
nStocks = size(tcosts,2);
nMonths = length(dates);
nFactors = arrayfun(@(x) length(x.factors),factor_model_defs);
nFmax = max(nFactors); % Max number of factors to determine W 3-d array size

for i=1:nModels
   thisModel = factor_model_defs(i).label;
   res(i,1) = getFactorModelData(thisModel,factor_model_defs,factor_struct,[s e]);       
end

% Initialize a few variables
W = nan(nMonths, nFmax, nModels);
W_netting = nan(nMonths, nFmax, nModels);
OSrets = nan(nMonths, nModels);
OSrets_netting = nan(nMonths, nModels);
errorFlag = zeros(nMonths, nModels);

% Clear the structures we don't need
clear factor_model_defs
clear factor_struct

startDate=s+120;

parfor m=1:nModels

    % Store this model results structure
    thisModelRes = res(m);
    
    % Store the start time &print it out.
    startTime = datetime('now');
    fprintf('Now working on %s. Run started at %s.\n', char(thisModelRes.label), char(startTime));
    
    % Store the number of factors
    nF = length(thisModelRes.factor_labels);

    
    % Initialize the weights matrices for this model
    thisModelW_netting = nan(nMonths,nFmax);
    thisModelW = nan(nMonths,nFmax);
    
    % Store the factor returns and tcosts
    gross_factors  = thisModelRes.gross_factors;
    net_factors    = thisModelRes.net_factors;
    tcosts_factors = gross_factors-net_factors;
    
    % Store the change in weight matrix for the non-market factors. We
    % assume that the market is free to trade (e.g., a market ETF).    
    dW = thisModelRes.dW(:,:,2:end);
    % Reshape it to make it a 2-d array
    dW = dW(:,:);
    
    % Loop through the months
    for t = startDate:nMonths
        
        % Determine the in-sample months
        inSampleInd=(s:t-1)';
        
        % Separate the tcosts
        tcostsIS = tcosts(inSampleInd,:);
        tcostsOS = tcosts(t,:);
        
        % Separate the gross returns
        GrossIS = gross_factors(inSampleInd,:);
        GrossOS = gross_factors(t,:);
        
        % Separate the dWs
        dWFacsIS = dW(inSampleInd,:);
        dWFacsOS = dW(t,:);        
        
        % Net return without netting
        [wStar, ~] = calcNetMve(net_factors(inSampleInd,:), -GrossIS-tcosts_factors(inSampleInd,:));
        
        % Calculate the out-of-sample returns without netting
        OSrets(t,m)     = (sign(wStar').*GrossOS)*wStar - tcosts_factors(t,:)*abs(wStar);
        
        % Assign the optimal weights without netting
        wStar1          = nan(nFmax,1);
        wStar1(1:nF)    = wStar;      
        thisModelW(t,:) = wStar1;
        
        
        % Figure out the starting weights for the netting results
        if t == startDate
            w0 = calcMve(GrossIS);
        else
            w0 = thisModelW_netting(t-1,:)';
            w0 = w0(1:nF);
        end

        % Get the default optimization options for fmincon        
        opts = optimoptions('fmincon','Display','off');
        % SharpeTD calculates a Sharpe ratio with trading diversification. Fun
        % will only take the weights and pass the other inputs to SharpeTD       
        fun = @(w) -SharpeTD(w,GrossIS,dWFacsIS,tcostsIS);
        
        % Allow for errors if something funky happens in the optimization
        try
            % Find the optimal weights and Sharpe ratio                       
            [wStar, ~] = fmincon(fun,w0,[],[],ones(1,nF),1,[],[],[],opts);
            
            % Store the optimal weights            
            wStar1          = nan(1,nFmax);
            wStar1(1:nF)    = wStar';
            thisModelW_netting(t,:) = wStar1;    
        catch
            errorFlag(t,m) = 1;
        end

        % Calculate the out-of-sample return
        nonMktWeights       = wStar(2:end);
        portdW              = dW(t,:) * kron(nonMktWeights, speye(nStocks)); 
        PortTCsCD           = sum(abs(portdW) .* tcostsOS, 2, 'omitnan');
        OSrets_netting(t,m) = GrossOS*wStar-PortTCsCD;
        
    end
    % Store the  weights    
    W(:,:,m) = thisModelW;
    W_netting(:,:,m) = thisModelW_netting;
    
    % Print out the time for timekeeping
    endTime = datetime('now');
    fprintf('Done with on %s. Run ended at %s and took %s.\n', char(thisModelRes.label), char(endTime), char(endTime-startTime));
end

fprintf('Done with OOS strategies at %s.\n',char(datetime('now')));

save Results/OOS_results OSrets  W OSrets_netting W_netting
