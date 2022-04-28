function dSR = SqSharpeDiff(c_scalar, long1, tc1, long2, tc2, startend)
% PURPOSE: calculates squared Sharpe ratio differences accounting for a
% scalar multiple of trading costs 
%---------------------------------------------------
% USAGE: dSR = SqSharpeDiff(c_scalar, long1, tc1, long2, tc2, startend)
%---------------------------------------------------
% Inputs
%        -c_scalar - scalar multiple
%        -long1 - matrix with asset returns (nTimePeriods x nAssets)
%        assuming long positions
%        -tc1 - matrix with trading costs (nTimePeriods x nAssets)
%        -long2 - matrix with asset returns (nTimePeriods x nAssets)
%        assuming long positions
%        -tc2 - matrix with trading costs (nTimePeriods x nAssets)
% Outputs
%        -dSR - difference of squared Sharpe ratios

% Determine the sample
s=startend(1);
e=startend(2);

% Calculate the Sharpe ratios
[~, SR1] = calcNetMve(long1(s:e,:) - c_scalar * tc1(s:e,:), ...
                     -long1(s:e,:) - c_scalar * tc1(s:e,:));

[~, SR2] = calcNetMve(long2(s:e,:) - c_scalar * tc2(s:e,:), ...
                     -long2(s:e,:) - c_scalar * tc2(s:e,:));

% Take the difference
dSR = (SR1 - SR2)^2;
