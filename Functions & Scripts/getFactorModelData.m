function res = getFactorModelData(factor_model, factor_model_defs, factor_struct, startend)
% PURPOSE: extracts data on a given factor model and stores it in a
% structure
%---------------------------------------------------
% USAGE: res = getFactorModelData(factor_model, factor_model_defs, factor_struct, startend)   
%---------------------------------------------------
% Inputs
%        -factor_model - cell with factor model label
%        -factor_model_defs - structure with definitions of factor models
%                    given identified by factor model labels
%        -factor_struct - structure with data on factors
%        -startend - monthly indicators of start and end dates for sample
% Outputs
%        -res - structure with data on the factor model

% Load the dates vector
load dates

% Determine the start & end date index
s = startend(1);
e = startend(2);

% Find the definition of the model
model_labels = [factor_model_defs.label]';
r = (strcmp(model_labels,factor_model));

% Pick the factors for this model
factors = factor_model_defs(r).factors;
ind = find(ismember([factor_struct.label],factors));

% Initiate a structure and store the relevant variables for the factors
res = struct;
res.label = factor_model;
res.factor_labels = [factor_struct(ind).label];
res.gross_factors = [factor_struct(ind).gross_factor];
res.net_factors = [factor_struct(ind).net_factor];
res.tc = [factor_struct(ind).tc];
res.to = [factor_struct(ind).to];
res.dW = reshape([factor_struct(ind).dW],size(factor_struct(1).dW,1),size(factor_struct(1).dW,2),length(ind));
res.startend = startend;
res.T = e-s+1;

% Calculate gross and net Sharpe ratios
[grossW, grossSR] = calcMve(res.gross_factors(s:e,:));
[netW, netSR] = calcNetMve(res.gross_factors(s:e,:) - res.tc(s:e,:), ...
                          -res.gross_factors(s:e,:) - res.tc(s:e,:));

% Assign those to the output structure
res.gross_weights = grossW';
res.net_weights = netW';
res.gross_sharpe = sqrt(12)*grossSR;
res.net_sharpe = sqrt(12)*netSR;
