function res = runSingleFactorModelBootstrap(retsLong, retsShort, index, nRuns)
% PURPOSE: runs a simple bootstrap (no netting) for a set of factor
% returns (long and short), index with in- and out-of-sample subsample, and
% number of runs
%---------------------------------------------------
% USAGE: res = runSingleFactorModelBootstrap(rets, rets_s, index, nRuns)    
%---------------------------------------------------
% Inputs
%        -retsLong - matrix with asset returns assuming a long position
%        (nTimePeriods x nAssets)
%        -retsShort - matrix with asset returns assuming a short position
%        (nTimePeriods x nAssets)
%        -index - index indicating the in- and out-of-sample subsamples
%        -nRuns - number of bootstrap runs
% Outputs
%        -res - an array of structures with bootstrap results (optimal 
%        weights, in- and out-of-sample Sharpe ratios, model label)
% This function works by searching through all combinations of factor
% positions being long, short, or not held. That, we calculate the Sharpe
% ratios for all 3^(nFactor-1) combinations of factor directional
% positions. We use the highest Sharpe ratio that has positive weights on
% the signed factors. The only exception is that we assume that we assume
% no trading costs for the market, so any position (long, short, not held)
% is fine. 

% Initialize the output structure
res = struct;

% Define a few constants
nMonths = size(retsLong, 1)/2;
nFactors = size(retsLong,2);
nCombs = 3^(nFactors-1);

% Initialize the position index. This will tell us whether a factor is
% long, short, or not held in a given combination
posIndex = nan(nCombs, nFactors);

% We assume market is costless to trade, so always enters the MVE
% calculation as long
posIndex(:,1) = [ones(nCombs,1)];

% Find the combinations (long, short, not held) for the other factors
xt = [nCombs:-1:1]'-1; 
for i=2:(nFactors)
    ttemp = mod(xt,3);
    xt = (xt-ttemp)/3;
    posIndex(:,i) = ttemp;
end

% Initialize the structure fields
res.weights = zeros(nRuns, nFactors);
res.SR_is = zeros(nRuns, 1);
res.SR_os = zeros(nRuns, 1);

for c = 1:nRuns
    
    % STore the returns for this subsample
    thisRunRetsLong = retsLong(index(:, c), :);
    thisRunRetsShort = retsShort(index(:, c), :);
    
    % Split them into in- and out-of-sample    
    thisRunRetsLongIS = thisRunRetsLong(1:nMonths, :);
    thisRunRetsShortIS = thisRunRetsShort(1:nMonths, :);
            
    thisRunRetsLongOOS = thisRunRetsLong(nMonths+1:end, :);
    thisRunRetsShortOOS = thisRunRetsShort(nMonths+1:end, :);
    
    % Loop through the variaous combinations of long/short/not held
    for i = 1:nCombs
        
        % Store the long returns in a temp ret variable
        tempRetsIS = thisRunRetsLongIS;
        tempRetsOOS = thisRunRetsLongOOS;
        
        % Find the ones we are shorting in this combination and replace
        indShort = posIndex(i,:) == 2;
        tempRetsIS(:, indShort) = thisRunRetsShortIS(:,indShort); 
        tempRetsOOS(:, indShort) = thisRunRetsShortOOS(:,indShort);

        % Find the ones we are not using in this combination and replace
        indNotUsed = posIndex(i,:) == 0;
        tempRetsIS(:, indNotUsed) = []; 
        tempRetsOOS(:, indNotUsed) = []; 
        
        tempPosIndex = posIndex(i, :);
        tempPosIndex(:, indNotUsed) = [];

        % Calculate the MVE weights & Sharpe ratio
        [weights, SR] = calcMve(tempRetsIS);
        
        % Check if all the weights on the non-market factors are positive
        nNonMktFactors = size(tempRetsIS, 2) - 1;
        if sum(weights(2:end) > 0) == nNonMktFactors
            % Compare the resulting Sharpe ratio with the one we currently
            % have
            if SR > res.SR_is(c, 1)
                
                % Store it
                res.SR_is(c, 1) = SR;
                
                % Assign the weights
                tempWeights = zeros(1,nFactors);
                tempWeights(posIndex(i, :) == 1) = weights(tempPosIndex == 1);
                tempWeights(posIndex(i, :) == 2) = -weights(tempPosIndex == 2);
                res.weights(c,:) = tempWeights;
                
                % Calculate the OS Sharpe ratio
                oosRets = tempRetsOOS * weights;
                res.SR_os(c,1) = mean(oosRets)/std(oosRets);
            end
        end
    end
end
