function [alphahat,alphacov,b,covb] = gmmAlphas(lhf,rhf);

% function gmmAlphas does jointly estimated factor regressions with gmm
% corrected standard errors and returns estimated vector of alphas and
% covariance matrix. This codes borrows heavily from John Cochrane's olsgmm
% function used in Cochrane Piazzesi (AER, 2005)
% 
% Inputs:
%  lhv T x N vector, left hand variable data 
%  rhv T x K matrix, right hand variable data
%  If N > 1, this runs N regressions of the left hand columns on all the (same) right hand variables. 

%  NOTE: you must make one column of rhv a vector of ones if you want a constant. 
%        should the covariance matrix estimate take out sample means?
% Output:
%  b: regression coefficients K x 1 vector of coefficients
%  seb: K x N matrix standard errors of parameters. 
%      (Note this will be negative if variance comes out negative) 
%  v: variance covariance matrix of estimated parameters. If there are many y variables, the vcv are stacked vertically

if size(rhf,1) ~= size(lhf,1);
   disp('olsgmm: left and right sides must have same number of rows. Current rows are');
   size(lhf)
   size(rhf)
end;

T = size(lhf,1);
N = size(lhf,2);
K = size(rhf,2);
Sxxinv = inv((rhf'*rhf)/T);
b = rhf\lhf;
covb=zeros(N*K);


errv = lhf-rhf*b; %TxN

%     s2 = mean(errv.^2);
%     vary = lhf - ones(T,1)*mean(lhf);
%     vary = mean(vary.^2);

%     R2v = (1-s2./vary)';
%     R2vadj= (1 - (s2./vary)*(T-1)/(T-K))';

%compute GMM standard errors; 
for indx = 1:N;
    for jindx=(1:N);
        erri=errv(:,indx);
        errj=errv(:,jindx);
	    S = (rhf.*(erri*ones(1,K)))'*(rhf.*(errj*ones(1,K)))/T;
        covb((indx-1)*K+1:indx*K,(jindx-1)*K+1:jindx*K)=1/T*Sxxinv*S*Sxxinv;
    end
end

alphaindx=mod((1:N*K)',K)==1;
alphahat=b(1,:)';
alphacov=covb(alphaindx,alphaindx);
    