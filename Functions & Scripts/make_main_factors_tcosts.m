%% Read the raw HXZ factors

clear
clc

load dates

% Store the number of months
nMonths = length(dates);

% Download them directly from global-q.org website
hxzdata = webread('http://global-q.org/uploads/1/2/2/6/122679606/q5_factors_monthly_2021.csv');
writetable(hxzdata,'Data/hxz_raw_data.csv');

% Intersect the HXZ dates with our dates
hxzdates  = 100*hxzdata.year + hxzdata.month;
[~, iDates, iHXZ] = intersect(dates,hxzdates);

% Initialize the vectors for the HXZ factors
hxzrf  = nan(nMonths,1); 
hxzmkt = nan(nMonths,1); 
rme    = nan(nMonths,1); 
ria    = nan(nMonths,1); 
rroe   = nan(nMonths,1);
reg    = nan(nMonths,1);

% Store them in vectors
hxzrf(iDates)  = hxzdata.R_F(iHXZ)/100;
hxzmkt(iDates) = hxzdata.R_MKT(iHXZ)/100;
rme(iDates)    = hxzdata.R_ME(iHXZ)/100;
ria(iDates)    = hxzdata.R_IA(iHXZ)/100;
rroe(iDates)   = hxzdata.R_ROE(iHXZ)/100;
reg(iDates)    = hxzdata.R_EG(iHXZ)/100;

% Add a constant vector
const = ones(nMonths,1)/100;

% Matrix of the four HXZ factors
hxz = [const hxzmkt rme ria rroe];

% Save the factor vectors and matrix in the /Data/ folder
save Data/hxz hxzrf hxzmkt rme ria rroe const reg hxz

%% Replicate HXZ factors and make trading costs

clear
clc

% Load the variables necessary for replication
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
indME = makeUnivSortInd(mcap, 2, NYSE);
indIA = makeUnivSortInd(ia, [30 70], NYSE);
indROE = makeUnivSortInd(roe, [30 70], NYSE);

% Assign the stocks to the 18 (= 2 x 3 x 3) portfolios. 
idx_ptfs = zeros(size(ret));
for i=1:2 % Loop over mcap bins
    for j=1:3 % Loop over roe bins
        for k=1:3 % Loop over ia bins
            n = (i-1)*9 + (j-1)*3 + k;
            idx_ptfs(indME==i & indROE==j & indIA==k) = n;
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

% % Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [rme_tc, rme_TO, rme_dW] = makeFactorTcosts(tcosts, me, w/9, indME);
% [ria_tc, ria_TO, ria_dW] = makeFactorTcosts(tcosts, me, w/6, indIA);
% [rroe_tc, rroe_TO, rroe_dW] = makeFactorTcosts(tcosts, me, w/6, indROE);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
[rme_ptf_tc, rme_TO, rme_dW] = calcTcosts(tcosts, indME, me, 'weighting', w/9);
[ria_ptf_tc, ria_TO, ria_dW] = calcTcosts(tcosts, indIA, me, 'weighting', w/6);
[rroe_ptf_tc, rroe_TO, rroe_dW] = calcTcosts(tcosts, indROE, me, 'weighting', w/6);

% Add the long and short portfolio trading costs
rme_tc = sum(rme_ptf_tc, 2);
ria_tc = sum(ria_ptf_tc, 2);
rroe_tc = sum(rroe_ptf_tc, 2);

% Store them in the replication vectors
rme_rep = resME.pret(:,end);
ria_rep = resIA.pret(:,end);
rroe_rep = resROE.pret(:,end); 

% Check the replication quality (should see >95% R2)
load data/hxz
prt(nanols(rme, [const rme_rep]))
prt(nanols(ria, [const ria_rep]))
prt(nanols(rroe, [const rroe_rep]))

% Save the replicated factors, tcosts, and dWs in the /Data/ folder
save data/hxz_rep rme_rep rroe_rep ria_rep
save data/hxz_tc rme_tc rme_TO rroe_tc rroe_TO ria_tc ria_TO
save Data/hxz_dW rme_dW rroe_dW ria_dW

%% Replicate FF factors and make trading costs 

clear
clc

% Load the variables necessary for replication
load ret
load me
load bm
load R
load NYSE
load dates
load gibbs_filled
load ff
load AT
load revt
load cogs
load xint
load xsga
load be
load FinFirms

% Operating Profitability = (REVT - COGS - XSGA - XINT)/BE
OpProf = REVT + nanmatsum(-COGS, -nanmatsum(XSGA,XINT));
OpProfToBe = OpProf./BE;

% Investment is the negative of asset growth (AT/AT_{-12})
inv = -AT./lag(AT,12,nan);

% Drop the non-positive bm observations
bm(bm <= 0) = nan;
OpProfToBe(isnan(bm)) = nan;
inv(isnan(bm))=nan;

% Get the indices & weight matrices for the makeFactorTcosts function

% Start with the FF93 factors: SMB & HML

% Initialize the weight matrix for the SMB & HML factors
wFF93 = nan(size(ret));

% Calculate the index, which has the portfolio allocations only in June
indFF93 = makeBivSortInd(me, 2, bm, [30 70], 'sortType', 'unconditional', ...
                                             'breaksFilterInd', NYSE);
% Store the number of portfolios
nPtfs = max(max(indFF93));
nStocks=size(ret,2);

% Create an index that has the portfolio allocations in each month. We need
% this to calculate the weights that we'll pass on later.
monthlyIndFF93 = indFF93;
indNoStocks = sum(monthlyIndFF93,2) == 0;
monthlyIndFF93(indNoStocks,:) = nan;
monthlyIndFF93 = FillMonths(monthlyIndFF93);

for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(monthlyIndFF93 ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wFF93(monthlyIndFF93 == i) = meThisPtf(monthlyIndFF93 == i) ./ ...
                                 rptdSumMeThisPtf(monthlyIndFF93 == i);        
end

% Do the same for RMW
% Initialize the weight matrix for the RMW factors
wMeOpProf = nan(size(ret));

% Calculate the index, which has the portfolio allocations only in June
indMeOpProf = makeBivSortInd(me, 2, OpProfToBe, [30 70], 'sortType', 'unconditional', ...
                                             'breaksFilterInd', NYSE);

% Create an index that has the portfolio allocations in each month. We need
% this to calculate the weights that we'll pass on later.
monthlyIndMeOpProf = indMeOpProf;
indNoStocks = sum(monthlyIndMeOpProf,2) == 0;
monthlyIndMeOpProf(indNoStocks,:) = nan;
monthlyIndMeOpProf = FillMonths(monthlyIndMeOpProf);

for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(monthlyIndMeOpProf ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wMeOpProf(monthlyIndMeOpProf == i) = meThisPtf(monthlyIndMeOpProf == i) ./ ...
                                 rptdSumMeThisPtf(monthlyIndMeOpProf == i);        
end

% Do the same for CMA
% Initialize the weight matrix for the RMW factors
wMeInv = nan(size(ret));

% Calculate the index, which has the portfolio allocations only in June
indMeInv = makeBivSortInd(me, 2, inv, [30 70], 'sortType', 'unconditional', ...
                                             'breaksFilterInd', NYSE);

% Create an index that has the portfolio allocations in each month. We need
% this to calculate the weights that we'll pass on later.
monthlyIndMeInv = indMeInv;
indNoStocks = sum(monthlyIndMeInv,2) == 0;
monthlyIndMeInv(indNoStocks,:) = nan;
monthlyIndMeInv = FillMonths(monthlyIndMeInv);

for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(monthlyIndMeInv ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wMeInv(monthlyIndMeInv == i) = meThisPtf(monthlyIndMeInv == i) ./ ...
                                 rptdSumMeThisPtf(monthlyIndMeInv == i);        
end


% Do the same for UMD
% Initialize the weight matrix for the RMW factors
wMeMom = nan(size(ret));

% Calculate the index, which has the portfolio allocations only in June
indMeMom = makeBivSortInd(me, 2, R, [30 70], 'sortType', 'unconditional', ...
                                             'breaksFilterInd', NYSE);

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
indSMB = 1*ismember(indFF93,[4:6])     + 2*ismember(indFF93,[1:3]);
indHML = 1*ismember(indFF93,[1 4])     + 2*ismember(indFF93,[3 6]);
indRMW = 1*ismember(indMeOpProf,[1 4]) + 2*ismember(indMeOpProf,[3 6]);
indCMA = 1*ismember(indMeInv,[1 4])    + 2*ismember(indMeInv,[3 6]);
indUMD = 1*ismember(indMeMom,[1 4])    + 2*ismember(indMeMom,[3 6]);

% % Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [smb_tc, smb_TO, smb_dW] = makeFactorTcosts(tcosts,me,wFF93/3,indSMB);
% [hml_tc, hml_TO, hml_dW] = makeFactorTcosts(tcosts,me,wFF93/2,indHML);
% [rmw_tc, rmw_TO, rmw_dW] = makeFactorTcosts(tcosts,me,wMeOpProf/2,indRMW);
% [cma_tc, cma_TO, cma_dW] = makeFactorTcosts(tcosts,me,wMeInv/2,indCMA);
% [umd_tc, umd_TO, umd_dW] = makeFactorTcosts(tcosts,me,wMeMom/2,indUMD);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
[smb_ptf_tc, smb_TO, smb_dW] = calcTcosts(tcosts, indSMB, me, 'weighting', wFF93/3);
[hml_ptf_tc, hml_TO, hml_dW] = calcTcosts(tcosts, indHML, me, 'weighting', wFF93/2);
[rmw_ptf_tc, rmw_TO, rmw_dW] = calcTcosts(tcosts, indRMW, me, 'weighting', wMeOpProf/2);
[cma_ptf_tc, cma_TO, cma_dW] = calcTcosts(tcosts, indCMA, me, 'weighting', wMeInv/2);
[umd_ptf_tc, umd_TO, umd_dW] = calcTcosts(tcosts, indUMD, me, 'weighting', wMeMom/2);

% Add the long and short portfolio trading costs
smb_tc = sum(smb_ptf_tc, 2);
hml_tc = sum(hml_ptf_tc, 2);
rmw_tc = sum(rmw_ptf_tc, 2);
cma_tc = sum(cma_ptf_tc, 2);
umd_tc = sum(umd_ptf_tc, 2);


% Replicate the factors using the adjusted weights matrix w 
resSMB = runUnivSort(ret, indSMB, dates, wFF93/3,     'weighting','v', ...
                                                      'factorModel',1, ... 
                                                      'printResults',0, ...
                                                      'plotFigure',0);
resHML = runUnivSort(ret, indHML, dates, wFF93/2,     'weighting','v', ...
                                                      'factorModel',1, ...
                                                      'printResults',0, ...
                                                      'plotFigure',0);
resRMW = runUnivSort(ret, indRMW, dates, wMeOpProf/2, 'weighting','v', ...
                                                      'factorModel',1, ...
                                                      'printResults',0, ...
                                                      'plotFigure',0);
resCMA = runUnivSort(ret, indCMA, dates, wMeInv/2,    'weighting','v', ...
                                                      'factorModel',1, ...
                                                      'printResults',0, ...
                                                      'plotFigure',0);
resUMD = runUnivSort(ret, indUMD, dates, wMeMom/2,    'weighting','v', ...
                                                      'factorModel',1, ...
                                                      'printResults',0, ...
                                                      'plotFigure',0);

% Store them in the replication vectors
smb_rep = resSMB.pret(:,end);
hml_rep = resHML.pret(:,end); 
rmw_rep = resRMW.pret(:,end);
cma_rep = resCMA.pret(:,end);
umd_rep = resUMD.pret(:,end);


% Check the replication quality (should see >95% R2)
prt(nanols(smb,[const smb_rep]))
prt(nanols(hml,[const hml_rep]))
prt(nanols(rmw,[const rmw_rep]))
prt(nanols(cma,[const cma_rep]))
prt(nanols(umd,[const umd_rep]))

% Save the replicated factors, tcosts, and dWs in the /Data/ folder
save data/ff_rep mkt smb_rep hml_rep rmw_rep cma_rep umd_rep
save data/ff_tc_dnmv hml_tc smb_tc rmw_tc cma_tc umd_tc  hml_TO smb_TO rmw_TO cma_TO umd_TO 
save data/ff_dW smb_dW hml_dW rmw_dW cma_dW umd_dW

%% Read the AQR HML(m) factor

clear
clc

load dates

% Read the factor directly from AQR's website
aqrData = webread('https://images.aqr.com/-/media/AQR/Documents/Insights/Data-Sets/The-Devil-in-HMLs-Details-Factors-Monthly.xlsx');
writetable(aqrData,'Data/hmlm_raw_data.csv');

% Check that Var1 (first column) is datetime format
aqrData = aqrData(:,{'Var1','Var25'});
if ~isdatetime(aqrData.Var1)
    aqrData.Var1 = datetime(aqrData.Var1);
end

% Turn to YYYYMM
aqrdates = 100*year(aqrData.Var1)+month(aqrData.Var1);
[~,ia,ib] = intersect(dates,aqrdates);

% Asign HML(m) factor to a vector & store
hmlm = nan(size(dates)); 
hmlm(ia) = aqrData.Var25(ib);
save data/hmlm hmlm

%% Replicate AQR HML(m) & make its tcosts

clear
clc

% Load the variables necessary for replication
load ret
load me
load be
load dates
load NYSE
load ff
load HMLm
load gibbs_filled

% Prepare the variables 
bm = FillMonths(BE)./me;
% AQR also kick out negative book-equity firms 
bm(bm < 0) = nan; 

% Initialize the weight matrix for the HMLM factor
wMeBe = nan(size(ret));

% Calculate the index, which has the portfolio allocations only in June
indMeBe = makeBivSortInd(me, 2, bm, [30 70], 'sortType', 'unconditional', ...
                                             'breaksFilterInd', NYSE);

% Get the number of portfolios
nPtfs = max(max(indMeBe));
nStocks = size(ret,2);

for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(indMeBe ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wMeBe(indMeBe == i) = meThisPtf(indMeBe == i) ./ ...
                                 rptdSumMeThisPtf(indMeBe == i);        
end

% Create the long/short portfolio index (1=short; 2=long)
indHMLm = 1*ismember(indMeBe,[1 4]) + 2*ismember(indMeBe,[3 6]);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [hmlm_tc, hmlm_TO, hmlm_dW] = makeFactorTcosts(tcosts, me, wMeBe/2, indHMLm);
[hmlm_ptf_tc, hmlm_TO, hmlm_dW] = calcTcosts(tcosts, indHMLm, me, 'weighting', wMeBe/2);

% Add the long and short portfolio trading costs
hmlm_tc = sum(hmlm_ptf_tc, 2);

% Replicate HML(m)
resHMLm = runUnivSort(ret, indHMLm, dates, wMeBe/2, 'weighting', 'v', ... 
                                                'factorModel',1, ...
                                                'printResults',0, ...
                                                'plotFigure',0);

% Store HML(m) in a vector
hmlm_rep = resHMLm.pret(:,end);

% Check the replication quality
prt(nanols(hmlm, [const hmlm_rep]))

save data/hmlm_tc hmlm_TO hmlm_tc
save data/hmlm_rep hmlm_rep
save Data/hmlm_dW hmlm_dW

%% Make RMW_C and its trading costs 

clear
clc

% Load the variables necessary for replication
load ret
load me
load bm
load R
load NYSE
load dates
load gibbs_filled
load ff
load AT
load revt
load cogs
load xint
load xsga
load be
load FinFirms
load AT
load ACT
load RECT
load AP
load INVT
load XPP
load DRC
load DRLT
load XACC

% Prepare the variables 
% CP=(REVT - COGS - XSGA - XINT -dRECT - dXPP +dAP - dINVT +dR +dXACC)/BE
OP = REVT - nanmatsum(COGS, nanmatsum(XSGA, XINT));

% Deferred revenue
DR = DRC + DRLT;

% Exlude zero or negative book value
BE(BE <= 0) = nan;

% Create changes in variables that enter accruals calculation
dRECT = RECT-lag(RECT,12,nan); % Receivables
dXPP  = XPP-lag(XPP,12,nan);   % Prepaid expenses
dAP   = AP-lag(AP,12,nan);     % Account payable
dINVT = INVT-lag(INVT,12,nan); % Inventories
dDR   = DR-lag(DR,12,nan);     % Deferred revenue
dXACC = XACC-lag(XACC,12,nan); % Accued expenses

% Create accruals
Accruals = nanmatsum(zeros(size(dRECT)), ... 
                     nanmatsum(dRECT,  ...
                     nanmatsum(dXPP, ...
                     nanmatsum(-dAP, ...
                     nanmatsum(dINVT, ...
                     nanmatsum(-dDR, -dXACC))))));

% Create cash-profitability
CP = OP - Accruals;

% Scale it by book equity
CpToBe = CP ./ BE;


% Do the same for RMW
% Initialize the weight matrix for the RMW factors
wMeCpProf = nan(size(ret));

% Calculate the index, which has the portfolio allocations only in June
indMeCpProf = makeBivSortInd(me, 2, CpToBe, [30 70], 'sortType', 'unconditional', ...
                                                     'breaksFilterInd', NYSE);

% Get the number of portfolios
nPtfs = max(max(indMeCpProf));
nStocks = size(ret,2);

% Create an index that has the portfolio allocations in each month. We need
% this to calculate the weights that we'll pass on later.
monthlyIndMeCpProf = indMeCpProf;
indNoStocks = sum(monthlyIndMeCpProf,2) == 0;
monthlyIndMeCpProf(indNoStocks,:) = nan;
monthlyIndMeCpProf = FillMonths(monthlyIndMeCpProf);

for i=1:nPtfs
    
    % Store the market cap for this portfolio
    meThisPtf = me;
    meThisPtf(monthlyIndMeCpProf ~= i) = nan;
    
    % Calculate the sum of market cap for this portfolio and store in a
    % matrix
    sumMeThisPtf = sum(meThisPtf,2,'omitnan');
    rptdSumMeThisPtf = repmat(sumMeThisPtf, 1, nStocks);   
    
    % Calculate the weights for this portfolio
    wMeCpProf(monthlyIndMeCpProf == i) = meThisPtf(monthlyIndMeCpProf == i) ./ ...
                                 rptdSumMeThisPtf(monthlyIndMeCpProf == i);        
end


% Create the long/short portfolio index (1=short; 2=long)
indRMWC = 1*ismember(indMeCpProf,[1 4]) + 2*ismember(indMeCpProf,[3 6]);

% Calculate the tcosts; Use adjusted w (see Equations (5) and (6) in Appendix)
% [rmw_c_tc, rmw_c_TO, rmw_c_dW] = makeFactorTcosts(tcosts, me, wMeCpProf/2, indRMWC);
[rmw_c_ptf_tc, rmw_c_TO, rmw_c_dW] = calcTcosts(tcosts, indRMWC, me, 'weighting', wMeCpProf/2);

% Add the long and short portfolio trading costs
rmw_c_tc = sum(rmw_c_ptf_tc, 2);


% Replicate HML(m)
resRMWC = runUnivSort(ret, indRMWC, dates, wMeCpProf/2, 'weighting', 'v', ... 
                                                'factorModel',1, ...
                                                'printResults',0, ...
                                                'plotFigure',0);

% Store it in the replication vector
rmw_c = resRMWC.pret(:,end);

% Save the replicated factor, tcosts, and dWs in the /Data/ folder
save data/ff_c rmw_c rmw_c_tc rmw_c_TO
save data/ff_c_dW rmw_c_dW

