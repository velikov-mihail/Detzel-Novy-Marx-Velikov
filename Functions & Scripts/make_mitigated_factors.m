
%% Make HXZ factor trading costs with banding

clear
clc

% Load the variables necessary for construction
load at
load me
load ibq
load beq
load nyse
load ret
load dates
load FinFirms
load ff
load gibbs_filled

% Exclude negative quarterly book equity
BEQ(BEQ<0) = nan; 

% Lag book equity 
lagBEQ = lag(BEQ, 3, nan);

% HXZ Sorting variables: market cap first
mcap = me;
mths = dates - 100*floor(dates/100);
mcap(mths ~= 6,:) = nan; % Only June market cap
mcap = FillMonths(mcap);

% Return-on-equity uses lagged book equity
roe = IBQ./lagBEQ; 

% Investment-to-assets is same as in Fama-French's asset growth 
assetGrowth = AT./lag(AT,12,nan);
ia = FillMonths(assetGrowth);

% Exclude financials
mcap(FinFirms==1) = nan;
ia(FinFirms==1) = nan;
roe(FinFirms==1) = nan;

% Sort the signals
% 40/60 buy/hold for size
% Start by a tertile sort based on [40 60] breakpoints
indME3 = makeUnivSortInd(mcap,[40 60],NYSE);
% Slow down sales (short covers) 
indME3 = make_sS(indME3,2);
% Drop the middle portfolio
indMEsS = 1*(indME3==1) + 2*(indME3==3);

% 20/40 buy/hold for ROE
% Start with a quintile sort
indROE = makeUnivSortInd(roe,5,NYSE);  
% Slow down sales (short covers)
indROE = make_sS(indROE,2);
% Combine the middle portfolios
indROEsS = (indROE>0) + (indROE>1) + (indROE>4);

% 20/40 buy/hold for IA
% Start with a quintile sort
indIA=makeUnivSortInd(ia,5,NYSE);  
% Slow down sales (short covers)
indIA=make_sS(indIA,2);
% Combine the middle portfolios
indIAsS=(indIA>0) + (indIA>1) + (indIA>4);

% Get the index for the makeFactorTcosts function
% Assign the stocks to the 18 (= 2 x 3 x 3) portfolios. 
idx_ptfs = zeros(size(ret));
for i=1:2 % Loop over mcap bins
    for j=1:3 % Loop over roe bins
        for k=1:3 % Loop over ia bins
            n = (i-1)*9 + (j-1)*3 + k;
            idx_ptfs(indMEsS==i & indROEsS==j & indIAsS==k) = n;
        end
    end
end

% Store the number of portfolios
nPtfs = max(max(idx_ptfs));
nStocks = size(ret,2);
nMonths = size(ret,1);

% Create the weights matrix
w = nan(size(ret));
for i=1:nPtfs
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(idx_ptfs ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    thisW = meThisPtf./rptdSumMeThisPtf;
    
    % Store the weights
    w(idx_ptfs == i) = thisW(idx_ptfs == i);
end


% Create the long/short portfolio indices for the three factors (1=short; 2=long)
indME = 1*ismember(idx_ptfs,[10:18]) + 2*ismember(idx_ptfs,[1:9]);
indIA = 1*ismember(idx_ptfs,[3:3:18]) + 2*ismember(idx_ptfs,[1:3:16]);
indROE = 1*ismember(idx_ptfs,[1:3 10:12]) + 2*ismember(idx_ptfs,[7:9 16:18]);

% % Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [rme_sS_tc, rme_sS_TO, rme_sS_dW] = makeFactorTcosts(tcosts, me, w/9, indME);
% [ria_sS_tc, ria_sS_TO, ria_sS_dW] = makeFactorTcosts(tcosts, me, w/6, indIA);
% [rroe_sS_tc, rroe_sS_TO, rroe_sS_dW] = makeFactorTcosts(tcosts, me, w/6, indROE);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
[rme_sS_ptf_tc, rme_sS_TO, rme_sS_dW] = calcTcosts(tcosts, indME, me, 'weighting', w/9);
[ria_sS_ptf_tc, ria_sS_TO, ria_sS_dW] = calcTcosts(tcosts, indIA, me, 'weighting', w/6);
[rroe_sS_ptf_tc, rroe_sS_TO, rroe_sS_dW] = calcTcosts(tcosts, indROE, me, 'weighting', w/6);

% Add the long and short portfolio trading costs
rme_sS_tc = sum(rme_sS_ptf_tc, 2);
ria_sS_tc = sum(ria_sS_ptf_tc, 2);
rroe_sS_tc = sum(rroe_sS_ptf_tc, 2);

% Replicate the factors using the adjusted weights matrix w 
resME = runUnivSort(ret, indME,  dates, w/9, 'factorModel', 1, ...
                                             'printResults',0, ...
                                             'plotFigure',0);
resIA = runUnivSort(ret, indIA,  dates, w/6, 'factorModel', 1, ...
                                             'printResults',0, ...
                                             'plotFigure',0);
resROE = runUnivSort(ret, indROE, dates, w/6, 'factorModel', 1, ...
                                              'printResults',0, ...
                                              'plotFigure',0);

% Store the factors in the replication vectors
rme_sS = resME.pret(:,end);
rroe_sS = resROE.pret(:,end); 
ria_sS = resIA.pret(:,end);

save data/hxz_sS rme_sS rroe_sS ria_sS
save data/hxz_sS_tc rme_sS_tc rme_sS_TO rroe_sS_tc rroe_sS_TO ria_sS_tc ria_sS_TO
save Data/hxz_sS_dW rme_sS_dW rroe_sS_dW ria_sS_dW

%% Make sS UMD factor trading costs 

clear
clc

% Load the variables necessary for construction
load ret
load me
load bm
load R
load NYSE
load dates
load gibbs_filled
load ff

% Create the buy/hold portfolio indicators
% 40/60 buy/hold for size
% Start by a tertile sort based on [40 60] breakpoints
indME3 = makeUnivSortInd(me,[40 60],NYSE);
% Slow down sales (short covers) 
indME3 = make_sS(indME3,2);
% Drop the middle portfolio
indMEsS = 1*(indME3==1) + 2*(indME3==3);

% 20/40 buy/hold for R
% Start with a quintile sort
indR5 = makeUnivSortInd(R,5,NYSE);  
% Slow down sales (short covers)
indR5 = make_sS(indR5,2);
% Combine the middle portfolios
indRsS = (indR5>0) + (indR5>1) + (indR5>4);

% Get the indices & weight matrices for the makeFactorTcosts function
% Initialize the weight matrix for the UMD factor
wMeMom = nan(size(ret));

% Calculate the index, which has the portfolio allocations only in June
indMeMom = ptfindex_from_inds(indMEsS,indRsS);

% Store the number of portfolios
nPtfs = max(max(indMeMom));
nStocks = size(ret,2);

for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(indMeMom ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wMeMom(indMeMom == i) = meThisPtf(indMeMom == i) ./ ...
                                 rptdSumMeThisPtf(indMeMom == i);        
end

% Create the long/short portfolio indices for the three factors (1=short; 2=long)
indUMD = 1*ismember(indMeMom,[1 4])    + 2*ismember(indMeMom,[3 6]);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [umd_sS_tc, umd_sS_TO, umd_sS_dW] = makeFactorTcosts(tcosts,me,wMeMom/2,indUMD);
[umd_sS_ptf_tc, umd_sS_TO, umd_sS_dW] = calcTcosts(tcosts, indUMD, me, 'weighting', wMeMom/2);

% Add the long and short portfolio trading costs
umd_sS_tc = sum(umd_sS_ptf_tc, 2);


% Replicate the factors using the adjusted weights matrix w 
resUMD = runUnivSort(ret, indUMD, dates, wMeMom/2,    'weighting','v', ...
                                                      'factorModel',1, ...
                                                      'printResults',0, ...
                                                      'plotFigure',0);

% Store the factors in the vectors
umd_sS = resUMD.pret(:,end);

save Data/umd_sS umd_sS
save Data/umd_sS_tc umd_sS_tc  umd_sS_TO 
save Data/umd_sS_dW umd_sS_dW

%% Replicate AQR HML & make its tcosts

clear
clc

% Load the variables necessary for construction
load ret
load me
load be
load dates
load NYSE
load ff
load HMLm
load gibbs_filled

% Prepare B/M
bm = FillMonths(BE) ./ me;
% AQR also kick out negative book-equity firms 
bm(bm < 0) = nan; 


% Create the buy/hold portfolio indicators
% 40/60 buy/hold for size
% Start by a tertile sort based on [40 60] breakpoints
indME3 = makeUnivSortInd(me,[40 60],NYSE);
% Slow down sales (short covers) 
indME3 = make_sS(indME3,2);
% Drop the middle portfolio
indMEsS = 1*(indME3==1) + 2*(indME3==3);

% 20/40 buy/hold for B/M
% Start with a quintile sort
indBM5 = makeUnivSortInd(bm,5,NYSE);  
% Slow down sales (short covers)
indBM5 = make_sS(indBM5,2);
% Combine the middle portfolios
indBMsS = (indBM5>0) + (indBM5>1) + (indBM5>4);


% Initialize the weight matrix for the UMD factor
wMeBm = nan(size(ret));

% Calculate the index
indMeBm = ptfindex_from_inds(indMEsS,indBMsS);

% Store the number of portfolios
nPtfs = max(max(indMeBm));
nStocks = size(ret,2);

for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(indMeBm ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wMeBm(indMeBm == i) = meThisPtf(indMeBm == i) ./ ...
                                 rptdSumMeThisPtf(indMeBm == i);        
end


% Create the long/short portfolio index for the factor (1=short; 2=long)
indBM = 1*ismember(indMeBm, [1 4]) + 2*ismember(indMeBm, [3 6]);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [hmlm_sS_tc, hmlm_sS_TO, hmlm_sS_dW] = makeFactorTcosts(tcosts, me, wMeBm/2, indBM);
[hmlm_sS_ptf_tc, hmlm_sS_TO, hmlm_sS_dW] = calcTcosts(tcosts, indBM, me, 'weighting', wMeBm/2);

% Add the long and short portfolio trading costs
hmlm_sS_tc = sum(hmlm_sS_ptf_tc, 2);

% Construct the factor using the adjusted weights matrix
resHMLM = runUnivSort(ret, indBM, dates, wMeBm/2, 'weighting','v', ...
                                                  'factorModel',1, ...
                                                  'printResults',0, ...
                                                  'plotFigure',0);

% Store the factor in the replication vector
hmlm_sS = resHMLM.pret(:,end);

save data/hmlm_sS hmlm_sS
save Data/hmlm_sS_dW hmlm_sS_dW
save data/hmlm_sS_tc hmlm_sS_tc hmlm_sS_TO 

