%% Figure 1 - Freak factor

clear
clc

load factors_dnmv

% Choose the sample
s = find(dates==197201);
e = find(dates==202112);

% Add the CAPM, FF3, and FRK model definitions
n = length(factor_model_defs);
factor_model_defs(n+1).label = 'CAPM';
factor_model_defs(n+1).factors = {'mkt'};
factor_model_defs(n+2).label = 'FF3';
factor_model_defs(n+2).factors = {'mkt','smb','hml'};
factor_model_defs(n+3).label = 'FRK';
factor_model_defs(n+3).factors = {'frk'};

% Choose the three models
factor_models = {'CAPM','FF3','FRK'};
nFactorModels = length(factor_models);

% Initialize the Sharpes matrix
sharpes=nan(nFactorModels,2);

for i=1:nFactorModels
    res(i,1) = getFactorModelData(factor_models(i),factor_model_defs,factor_struct,[s e]);
    sharpes(i,1) = (res(i).gross_sharpe);
    sharpes(i,2) = (res(i).net_sharpe);
end

% Overwrite the LVIRR net Sharpe to showcase that it is negative
lvirrGrossRets = res(3).gross_factors(s:e);
lvirrTcosts = res(3).tc(s:e);
lvirrNetRets = lvirrGrossRets - lvirrTcosts;
sharpes(3,2) = sqrt(12) * mean(lvirrNetRets) / std(lvirrNetRets);

% Plot the figure
figure;
bar(sign(sharpes) .* (sharpes.^2),'grouped');

% Format the figure
legendSize = 40;
legend(' Gross SR^2',' Net SR^2','Location','northwest','box','off');
xlim([0.5 3.5]);
ylim([-0.5 5]);
xticklabels({'CAPM','FF3','LV-IRR'});
set(gca,'FontSize',legendSize);
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Position',[0 0 8.5 11])
orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  

% Export it
export_fig('Figures/figure1.pdf','-transparent');
% print -depsc2 'Figures/figure1.eps';


%% Figure 2 - grouped bar chart

clear
clc

load factors_dnmv

% Choose the sample
s = find(dates==197201);
e = find(dates==202112);

% Choose the three models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};
nFactorModels = length(factor_models);

% Initialize the Sharpes matrix
sharpes=nan(nFactorModels,2);

% Get the Sharpe ratios
for i=1:nFactorModels
    res(i,1) = getFactorModelData(factor_models(i),factor_model_defs,factor_struct,[s e]);
    sharpes(i,1) = res(i).gross_sharpe;
    sharpes(i,2) = res(i).net_sharpe;
end

% Plot the figure
figure;
bar(sharpes.^2,'grouped');

% Format the figure
legendSize = 40;
legend(' Gross SR^2',' Net SR^2','Location','northwest','box','off');
xlim([0.5 6.5]);
ylim([0 3]);
xticklabels({'FF5','FF6','HXZ4','BS6','FF5_C','FF6_C'});
set(gca,'FontSize',legendSize);
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'Position',[0 0 8.5 11])
orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  

% Export it
export_fig('Figures/figure2.pdf','-transparent');

%% Figure 3 - MVE OOS Strategies 

clear
clc

% Load data
load factors_dnmv
load OOS_results
load dates
load pdates

% Store a few constants
fontSize = 40;
lineWidth = 4.5;
grayscale = 'n';
markers = repmat({'-',':'},1,3);

% Select factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_C','FF6_C'};

% Find the model locations
modelLoc = (ismember(upper([factor_model_defs.label]), upper(factor_models)));

% Store the MVE returns
rets = OSrets(:,modelLoc);

% Determine the starting date
s = find(isnan(sum(rets,2)),1,'last')+1;

% Prep the dates
x = pdates(s:end);
x = [2*x(1)-x(2); x];

% Cumulate the returns
y = cumprod(1+[zeros(1,size(rets,2)); rets(s:end,:)]);

% Store the default color order
defColorOrder = colororder;

% Select the color scheme
if strcmp(grayscale,'y')
    colors = [repmat([0.9,0.9,0.9],2,1); repmat([0.8,0.8,0.8],2,1); repmat([0.6,0.6,0.6],2,1)];
else
    colors = [repmat(defColorOrder(1,:),2,1); repmat(defColorOrder(2,:),2,1); repmat(defColorOrder(3,:),2,1)];
end

% plot the figure
figure;
for i=1:size(y,2)
    plot(x,y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)))              
    hold on;
end


% Format the figure
h = regexprep(factor_models, '_C', 'C');
h = cellfun(@(x) strcat([' ',x]),h,'UniformOutput',0);
xlim([1980 2022])
leg = legend(h,'Location','northwest','box','off');
leg.ItemTokenSize = [80, 50];
ylabel('Performance of $1');
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
set(gca,'LooseInset',get(gca,'TightInset'))
box on;

% Export it
export_fig('Figures/figure3.pdf','-transparent');

%% Figure 4 - Anomaly Square Sharpe improvement percentiles 

clear
clc

% Load the reults
load anomalyResults
load factors_dnmv

% Select factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5C','FF6C'};

% Get the performance metrics
gross_improvement = (gross_sharpe.^2) ./ (factor_gross_sharpe.^2) - 1;
net_improvement = (net_sharpe.^2) ./ (factor_net_sharpe.^2) - 1;

% Round down to zero if <10e-4
gross_improvement(gross_improvement<10e-4) = 0;
net_improvement(net_improvement<10e-4) = 0;

% Define some figure constants
lineWidth = 4.5;
fontSize = 20;
ylimit = 250;
xlowlimit = 50;
grayscale = 'n';
markers = repmat({'-',':'}, 1, 3);

% Define the x values (percentiles) 
x = 1:100;

% Store the default color order
defColorOrder = colororder;

% Select the color scheme
if strcmp(grayscale,'y')
    colors = [repmat([0.9,0.9,0.9],2,1); repmat([0.8,0.8,0.8],2,1); repmat([0.6,0.6,0.6],2,1)];
else
    colors = [repmat(defColorOrder(1,:),2,1); repmat(defColorOrder(2,:),2,1); repmat(defColorOrder(3,:),2,1)];
end

% Plot the first panel - Gross SR improvement

% Get the y values - percentiles defined by x from the improvement
y = prctile(gross_improvement, x, 1);

% Plot the figure
figure;
for i=1:size(y,2)
    plot(x,100*y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)));
    hold on;
end

% Format the figure
h = cellfun(@(x) strcat([' ',x]),factor_models,'UniformOutput',0);
leg=legend(h, 'Location', 'northwest', 'box', 'off');
leg.ItemTokenSize = [60, 50];
ylabel('Max SR^2 improvement (%)');
xlabel('Anomaly percentile rank by impact');
xlim([xlowlimit 100])
ylim([0 ylimit])
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
pbaspect([1 1 1])
set(gca,'LooseInset',get(gca,'TightInset'))

% Export it
if strcmp(grayscale,'y')
    export_fig('Figures/figure4a_gray.pdf','-transparent');
else
    export_fig('Figures/figure4a_color.pdf','-transparent');
end


% Get the y values - percentiles defined by x from the improvement
y = prctile(net_improvement, x, 1);

% Plot the figure
figure;
for i=1:size(y,2)
    plot(x,100*y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)));
    hold on;
end

% Format the figure
leg=legend(h, 'Location', 'northwest', 'box', 'off');
leg.ItemTokenSize = [60, 50];
ylabel('Max SR^2 improvement (%)');
xlabel('Anomaly percentile rank by impact');
xlim([xlowlimit 100])
ylim([0 ylimit])
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
pbaspect([1 1 1])
set(gca,'LooseInset',get(gca,'TightInset'))

% Export it
if strcmp(grayscale,'y')
    export_fig('Figures/figure4b_gray.pdf','-transparent');
else
    export_fig('Figures/figure4b_color.pdf','-transparent');
end

%% Figure 5 - Mitigated grouped bar chart

clear
clc

% Load the data
load factors_dnmv
load fullSampleResults

% Select the model groups
baseFactorModels    = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};
bandingFactorModels = {'FF5_sS','FF6_sS','HXZ4_sS','BS6_sS','FF5_c_sS','FF6_c_sS'};

% Store the results for each
res_base    = res(ismember([res.label], baseFactorModels));
res_banding = res(ismember([res.label], bandingFactorModels));

% Store the Sharpe ratios
sharpes = [res_base.net_sharpe; 
           res_base.netting_sharpe;    
           res_banding.net_sharpe;
           res_banding.netting_sharpe]';
    
% Plot the figure
figure;
bar(sharpes.^2,'grouped');

% Format
legendSize = 40;
legend(' No mitigation',' Netting',' Banding',' Netting+Banding','Location','northwest','box','off');
title('Net SR^2');
xlim([0.5 6.5]);
ylim([0 1.75]);
lbls = upper(regexprep(baseFactorModels,'_',''));
xticklabels(lbls);
set(gca, 'FontSize', legendSize);
set(gca, 'FontName', 'Times New Roman')
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(gca, 'LooseInset', get(gca, 'TightInset'))
set(gcf, 'Position', [0 0 8.5 11])
orient(gcf,'landscape')
set(gcf,'PaperSize',[20 10]); %set the paper size to what you want  

% Export it
export_fig('Figures/figure5.pdf','-transparent');

%% Figure 6 - MVE OOS Strategies with netting and banding

clear
clc

% Load data
load factors_dnmv
load OOS_results
load dates
load pdates

% Store a few constants
fontSize = 40;
lineWidth = 4.5;
grayscale = 'n';
markers = repmat({'-',':'},1,3);

% Select factor models
factor_models = {'FF5_sS','FF6_sS','HXZ4_sS','BS6_sS','FF5_C_sS','FF6_C_sS'};

% Find the model locations
modelLoc = (ismember(upper([factor_model_defs.label]), upper(factor_models)));

% Store the MVE returns
rets = OSrets_netting(:,modelLoc);

% Determine the starting date
s = find(isnan(sum(rets,2)),1,'last')+1;

% Prep the dates
x = pdates(s:end);
x = [2*x(1)-x(2); x];

% Cumulate the returns
y = cumprod(1+[zeros(1,size(rets,2)); rets(s:end,:)]);

% Store the default color order
defColorOrder = colororder;

% Select the color scheme
if strcmp(grayscale,'y')
    colors = [repmat([0.9,0.9,0.9],2,1); repmat([0.8,0.8,0.8],2,1); repmat([0.6,0.6,0.6],2,1)];
else
    colors = [repmat(defColorOrder(1,:),2,1); repmat(defColorOrder(2,:),2,1); repmat(defColorOrder(3,:),2,1)];
end

% plot the figure
figure;
for i=1:size(y,2)
    plot(x,y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)))              
    hold on;
end


% Format the figure
h = regexprep(factor_models, '_C', 'C');
h = regexprep(h, '_sS', '');
h = cellfun(@(x) strcat([' ',x]),h,'UniformOutput',0);
xlim([1980 2022])
leg = legend(h,'Location','northwest','box','off');
leg.ItemTokenSize = [80, 50];
ylabel('Performance of $1');
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
set(gca,'LooseInset',get(gca,'TightInset'))
box on;

% Export it
export_fig('Figures/figure6.pdf','-transparent');


%% Figure B.1 - Post-publication anomaly squared Sharpe ratio improvement percentiles

clear
clc

load anomalyResultsPostPub
load factors_dnmv

% Select factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5C','FF6C'};

% Get the performance metrics
gross_improvement = (gross_sharpe.^2) ./ (factor_gross_sharpe.^2) - 1;
net_improvement = (net_sharpe.^2) ./ (factor_net_sharpe.^2) - 1;

% Round down to zero if <10e-4
gross_improvement(gross_improvement<10e-4) = 0;
net_improvement(net_improvement<10e-4) = 0;

% Define some figure constants
lineWidth = 4.5;
fontSize = 20;
ylimit = 250;
xlowlimit = 50;
grayscale = 'n';
markers = repmat({'-',':'}, 1, 3);

% Define the x values (percentiles) 
x = 1:100;

% Store the default color order
defColorOrder = colororder;

% Select the color scheme
if strcmp(grayscale,'y')
    colors = [repmat([0.9,0.9,0.9],2,1); repmat([0.8,0.8,0.8],2,1); repmat([0.6,0.6,0.6],2,1)];
else
    colors = [repmat(defColorOrder(1,:),2,1); repmat(defColorOrder(2,:),2,1); repmat(defColorOrder(3,:),2,1)];
end

% Plot the first panel - Gross SR improvement

% Get the y values - percentiles defined by x from the improvement
y = prctile(gross_improvement, x, 1);

% Plot the figure
figure;
for i=1:size(y,2)
    plot(x,100*y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)));
    hold on;
end


% Format the figure
h = cellfun(@(x) strcat([' ',x]),factor_models,'UniformOutput',0);
leg=legend(h, 'Location', 'northwest', 'box', 'off');
leg.ItemTokenSize = [60, 50];
ylabel('Max SR^2 improvement (%)');
xlabel('Anomaly percentile rank by impact');
xlim([xlowlimit 100])
ylim([0 ylimit])
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
pbaspect([1 1 1])
set(gca,'LooseInset',get(gca,'TightInset'))

% Export it
if strcmp(grayscale,'y')
    export_fig('Figures/figureB1a_gray.pdf','-transparent');
else
    export_fig('Figures/figureB1a_color.pdf','-transparent');
end


% Get the y values - percentiles defined by x from the improvement
y = prctile(net_improvement, x, 1);

% Plot the figure
figure;
for i=1:size(y,2)
    plot(x,100*y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)));
    hold on;
end

% Format the figure
leg=legend(h, 'Location', 'northwest', 'box', 'off');
leg.ItemTokenSize = [60, 50];
ylabel('Max SR^2 improvement (%)');
xlabel('Anomaly percentile rank by impact');
xlim([xlowlimit 100])
ylim([0 ylimit])
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
pbaspect([1 1 1])
set(gca,'LooseInset',get(gca,'TightInset'))

% Export it
if strcmp(grayscale,'y')
    export_fig('Figures/figureB1b_gray.pdf','-transparent');
else
    export_fig('Figures/figureB1b_color.pdf','-transparent');
end

%% Figure E.1 - MVE OOS Strategies with banding


clear
clc

% Load data
load factors_dnmv
load OOS_results
load dates
load pdates

% Store a few constants
fontSize = 40;
lineWidth = 4.5;
grayscale = 'n';
markers = repmat({'-',':'},1,3);

% Select factor models
factor_models = {'FF5_sS','FF6_sS','HXZ4_sS','BS6_sS','FF5_C_sS','FF6_C_sS'};

% Find the model locations
modelLoc = (ismember(upper([factor_model_defs.label]), upper(factor_models)));

% Store the MVE returns
rets = OSrets(:,modelLoc);

% Determine the starting date
s = find(isnan(sum(rets,2)),1,'last')+1;

% Prep the dates
x = pdates(s:end);
x = [2*x(1)-x(2); x];

% Cumulate the returns
y = cumprod(1+[zeros(1,size(rets,2)); rets(s:end,:)]);

% Store the default color order
defColorOrder = colororder;

% Select the color scheme
if strcmp(grayscale,'y')
    colors = [repmat([0.9,0.9,0.9],2,1); repmat([0.8,0.8,0.8],2,1); repmat([0.6,0.6,0.6],2,1)];
else
    colors = [repmat(defColorOrder(1,:),2,1); repmat(defColorOrder(2,:),2,1); repmat(defColorOrder(3,:),2,1)];
end

% plot the figure
figure;
for i=1:size(y,2)
    plot(x,y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)))              
    hold on;
end


% Format the figure
h = regexprep(factor_models, '_C', 'C');
h = regexprep(h, '_sS', '');
h = cellfun(@(x) strcat([' ',x]),h,'UniformOutput',0);
xlim([1980 2022])
leg = legend(h,'Location','northwest','box','off');
leg.ItemTokenSize = [80, 50];
ylabel('Performance of $1');
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
set(gca,'LooseInset',get(gca,'TightInset'))
box on;

% Export it
export_fig('Figures/figureE1.pdf','-transparent');

%% Figure E2 - MVE OOS Strategies with netting

clear
clc

% Load data
load factors_dnmv
load OOS_results
load dates
load pdates

% Store a few constants
fontSize = 40;
lineWidth = 4.5;
grayscale = 'n';
markers = repmat({'-',':'},1,3);

% Select factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_C','FF6_C'};

% Find the model locations
modelLoc = (ismember(upper([factor_model_defs.label]), upper(factor_models)));

% Store the MVE returns
rets = OSrets_netting(:,modelLoc);

% Determine the starting date
s = find(isnan(sum(rets,2)),1,'last')+1;

% Prep the dates
x = pdates(s:end);
x = [2*x(1)-x(2); x];

% Cumulate the returns
y = cumprod(1+[zeros(1,size(rets,2)); rets(s:end,:)]);

% Store the default color order
defColorOrder = colororder;

% Select the color scheme
if strcmp(grayscale,'y')
    colors = [repmat([0.9,0.9,0.9],2,1); repmat([0.8,0.8,0.8],2,1); repmat([0.6,0.6,0.6],2,1)];
else
    colors = [repmat(defColorOrder(1,:),2,1); repmat(defColorOrder(2,:),2,1); repmat(defColorOrder(3,:),2,1)];
end

% plot the figure
figure;
for i=1:size(y,2)
    plot(x,y(:,i),char(markers(i)),'LineWidth',lineWidth,'color',(colors(i,:)))              
    hold on;
end


% Format the figure
h = regexprep(factor_models, '_C', 'C');
h = regexprep(h, '_sS', '');
h = cellfun(@(x) strcat([' ',x]),h,'UniformOutput',0);
xlim([1980 2022])
leg = legend(h,'Location','northwest','box','off');
leg.ItemTokenSize = [80, 50];
ylabel('Performance of $1');
set(gca,'FontSize',fontSize)
set(gca,'FontName','Times New Roman')
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
orient(gcf,'landscape')
set(gca,'LooseInset',get(gca,'TightInset'))
box on;

% Export it
export_fig('Figures/figureE2.pdf','-transparent');

