function out1 = WolakAlphaTest(alphahat,omegahat,reps)
% Function to implement Wolak's (1989, JoE) test of
%
%       H0*: alpha<=0,d2>=0,...,alpha<=0
%
%  vs.  H1*: (alpha1,alpha2,...,alphaK) in R^K  (ie: general alternative)
%
% This code borrows heavily from the replication code provided by
% Andrew Patton for: 
% "Monotonicity in Asset Returns: New Tests with Applications to the Term
% Structure, the CAPM and Portfolio Sorts", (JFE, 2010)
% Paper available at:
% http://www.economics.ox.ac.uk/members/andrew.patton/research.html
%
% INPUTS:   alphahat = [Nx1] vector of estimated alphas from GMM_alphas
%           Omegahat = [NxN] covariance matrix of estimated alphas from GMM
%           reps, a scalar, the number of simulations to use to estimate the weight function in the weighted-sum of chi2 variables (default=100)
%
% OUTPUTS:  out1, a 2x1 vector, the p-values from Wolak's test 1 (null is monotonic) and test 2 (alternative is monotonic)

% [T,N] = size(data);
N=size(alphahat,1);

if nargin<3 || isempty(reps);
    reps=1000;
end

options = optimset('Display','off','TolCon',10^-8,'TolFun',10^-8,'TolX',10^-8); %Same as Patton and Timmerman (2010)
Aineq = eye(N);
Bineq = zeros(N,1);

% warning off;  % trying to speed up the code - eliminating warning messages printed to the screen
alphatilda = fmincon('constrained_mean',alphahat,Aineq,Bineq,[],[],[],[],[],options,alphahat,omegahat);
IU = (alphahat-alphatilda)'*inv(omegahat)*(alphahat-alphatilda);  % the first test stat, equation 16
EI = alphatilda'*inv(omegahat)*alphatilda;  % the second test statistic, see just after equation 18 of Wolak (1989, JoE)
warning on;

% next: use monte carlo to obtain the weights for the weighted sum of chi-squareds
weights = zeros*ones(1+size(omegahat,1),1);
% test=[]
for jj=1:reps;
    tempdata = mvnrnd(zeros(1,size(omegahat,1)),omegahat,1);  % simulating iid Normal data
%     warning off;
    mutilda1 = fmincon('constrained_mean',tempdata',Aineq,Bineq,[],[],[],[],[],options,tempdata',omegahat);
    warning on;
%     test=[test; mutilda1'];
    temp = sum(mutilda1<-10^-7);  % counting how many elements of mutilda are greater than zero;
    %10^-8 is the tolerance from fmincon, so |mutilda1|<10^-7 is "0".
    weights(1+temp) = weights(1+temp) + 1/reps;  % adding one more unit of weight to this element of the weight vector
end
out1 = 1-weighted_chi2cdf(IU,weights(end:-1:1));  % pvalue from wolak's first test

% test
% weights
end