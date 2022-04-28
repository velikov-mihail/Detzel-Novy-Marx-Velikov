clear
clc

% Load the variables necessary for construction
load ret
load me
load RVOL1
load NYSE
load dates
load gibbs_filled
load ff
load ireta

% Prepare the variables 
% Sort on 1-month realized volatility with NYSE breakpoints
indRVOL1 = makeUnivSortInd(RVOL1, 2, NYSE);

% Store an indicator for low-volatility universe
lowVolInd = nan(size(ret)); 
lowVolInd(indRVOL1 == 1) = 1;

% Store the industry-relative reversals
IRR = ireta - ret;

% Store low-volatility industry-relative reversals
LVIRR = IRR .* lowVolInd;

% Get the index for the makeFactorTcosts function
indMeFrk = makeBivSortInd(me, 2, LVIRR, [30 70], 'sortType', 'unconditional', ...
                                                  'breaksFilterInd', NYSE);

 % Get the number of portfolios
nPtfs = max(max(indMeFrk));
nStocks = size(ret,2);
                                             
% Initialize the weight matrix for the Freak factor
wMeFrk = nan(size(ret));
                                          
for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(indMeFrk ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wMeFrk(indMeFrk == i) = meThisPtf(indMeFrk == i) ./ ...
                                 rptdSumMeThisPtf(indMeFrk == i);        
end     

% Create the long/short portfolio index for the factor (1=short; 2=long)
indFRK = 1*ismember(indMeFrk, [1 4]) + 2*ismember(indMeFrk, [3 6]);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [frk_tc, frk_TO, ~] = makeFactorTcosts(tcosts, me, wMeFrk/2, indFRK);
[frk_ptf_tc, frk_TO, ~] = calcTcosts(tcosts, indFRK, me, 'weighting', wMeFrk/6);

% Add the long and short portfolio trading costs
frk_tc = sum(frk_ptf_tc, 2);

% Construct the factor using the adjusted weights matrix w 
resFRK = runUnivSort(ret, indFRK, dates, wMeFrk/2, 'weighting', 'v', ...
                                                   'factorModel', 1, ...
                                                   'printResults', 0, ...
                                                   'plotFigure',0);

% Store the factor in the replication vector
frk = resFRK.pret(:,end);

save data/frk_tc frk_tc frk_TO 
save data/frk frk

% Check out its Sharpe
s = find(dates==197201);
e = find(dates==202112);
fprintf('\nGross Sharpe is %2.2f.\n Net Sharpe is %2.2f.\n',sqrt(12)*mean(frk(s:e))/std(frk(s:e)),sqrt(12)*mean(frk(s:e)-frk_tc(s:e))/std(frk(s:e)-frk_tc(s:e)));
