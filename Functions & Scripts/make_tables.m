%% Start a diary (will have all the table results)
clear
clc

diary Results/DetzelNovy-MarxVelikovTablesOutput.txt

fprintf('\n\n\n\nTable printing started @ %s\n\n\n',char(datetime('now')));

%% Print Table II - descriptive stats

fprintf('\n\n\nTable 2 output:\n\n\n');

clear
clc

load factors_dnmv

% Create a constant for the regressions below
const = 0.01 * ones(size(dates));

% Drop the non-basic factors
indToDrop = (contains([factor_struct.label], {'rep','_sS','frk'}));
factor_struct(indToDrop) = [];

% Store the number of factors
nFactors=length(factor_struct);

% Choose the sample
s = find(dates==197201);
e = find(dates==202112);

% Create the table output matrix
outputMatrix = nan(nFactors,6);

for i=1:nFactors
    % Gross return 
    grossRet = factor_struct(i).gross_factor(s:e); 
    res = ols(grossRet, const(s:e));
    outputMatrix(i,1) = res.beta;
    outputMatrix(i,2) = res.tstat;
    
    % Net return
    netRet = factor_struct(i).net_factor(s:e);
    res = ols(netRet, const(s:e));
    outputMatrix(i,3) = res.beta;
    outputMatrix(i,4) = res.tstat;    
    
    % Turnover
    TO = factor_struct(i).to(s:e,:);
    res = ols(TO, const(s:e));
    outputMatrix(i,5) = res.beta;
    outputMatrix(i,6) = (outputMatrix(i,1) - outputMatrix(i,3));
end

% Store the headers for the factors
h = upper([factor_struct.label]');
h = regexprep(h,'_C','\\textsubscript{C}');

% Clean up the headers 
h(strcmp(h,{'UMD'}))={'MOM'};
h(strcmp(h,{'RME'}))={'ME'};
h(strcmp(h,{'RIA'}))={'IA'};
h(strcmp(h,{'RROE'}))={'ROE'};
h(strcmp(h,{'HMLM'}))={'HML(m)'};

% Add a few empty columns
outputMatrix = [outputMatrix(:,1:2) nan(length(factor_struct),1) outputMatrix(:,3:4) nan(length(factor_struct),1) outputMatrix(:,5:6)];

% Print the output table
colFormatSpecification = [{'%2.2f'},{'%2.2f'},{'%2.2f'},{'%2.2f'},{'%2.2f'},{'%2.2f'},{'%2.1f'},{'%2.2f'}];
mat2Tex(outputMatrix,outputMatrix,'rowHeaders',h,'colFormatSpec',colFormatSpecification);

%% Table III - MVE weights

fprintf('\n\n\nTable 3 output:\n\n\n');

clear
clc

load factors_dnmv
load fullSampleResults

% Drop the non-basic factors
indToDrop = (contains([factor_struct.label], {'rep','_sS','frk'}));
factor_struct(indToDrop) = [];

% Define the factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};

% Leave only the ones specified by factor_models above
indToKeep = ismember([res.label], factor_models);
res = res(indToKeep);

% Store the number of factors models & factors
nFactorModels = length(factor_models);
nFactors = length(factor_struct);

% Figure out the headers for each column
colheads = [factor_struct.label];

% Initialize the output matrix
outputMatrix = nan(nFactorModels, nFactors+1);

for i=1:nFactorModels
    nFactorsThisFactorModel = length(res(i).factor_labels);
    for j=1:nFactorsThisFactorModel
        % Find the position of this factor in the table
        c = strcmp(colheads, res(i).factor_labels(j));
        % Assign it
        outputMatrix(i,c) = 100 * res(i).gross_weights(j);
    end
    outputMatrix(i,end) = res(i).gross_sharpe .^ 2;    
end

% Define the headers
h = [res.label];
h = regexprep(h, '_c', '\\textsubscript{C}');

% Specify the format for each column to be passed to mat2Tex
colFormatSpecification = [repmat({'%2.0f'},1,11) {'%2.2f'}];

fprintf('Panel A:\n');
mat2Tex(outputMatrix, outputMatrix, h, 'colFormatSpec', colFormatSpecification);


% Initialize the output matrix
outputMatrix = nan(nFactorModels, nFactors+1);

for i=1:nFactorModels
    
    % Store the number of factors in this factor model
    nFactorsThisFactorModel = length(res(i).factor_labels);
    
    for j=1:nFactorsThisFactorModel
        % Find the position of this factor in the table
        c = strcmp(colheads, res(i).factor_labels(j));
        % Assign it
        outputMatrix(i,c) = 100 * res(i).net_weights(j);
    end
    outputMatrix(i,end) = res(i).net_sharpe .^ 2;
end

% Define the headers
h = [res.label];
h = regexprep(h, '_c', '\\textsubscript{C}');

% Specify the format for each column to be passed to mat2Tex
colFormatSpecification = [repmat({'%2.0f'},1,11) {'%2.2f'}];

fprintf('Panel B:\n');
mat2Tex(outputMatrix, outputMatrix, h, 'colFormatSpec', colFormatSpecification);



%% Table IV - bootrstraped results
 
fprintf('\n\n\nTable 4 output:\n\n\n');

clear
clc

% Load the stored results
load bootstrap_results
load bootstrap_capm

% Filter the bootstrap results to keep the basic models
indToKeep = ismember([bootstrap_results.label],{'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'});
res = bootstrap_results(indToKeep);

% Store the number of runs
nRuns = length(res(1).SR_is);

% Append CAPM results
res = [res(1) res];
res(1).weights = ones(nRuns,1);
res(1).SR_is = capmIsSR(1:min(size(capmIsSR,1),nRuns));
res(1).SR_os = capmOosSR(1:min(size(capmIsSR,1),nRuns));
res(1).label = {'CAPM'};

% Store the number of factor omdels
nFactorModels = length(res);

% Store the in- and out-of-sample Sharpes
all_sharpes_is = ([res.SR_is]);
all_sharpes_os = ([res.SR_os]);

% Initialize the output matrices 
panelAOutput = nan(nFactorModels, 8);
panelBOutput = nan(nFactorModels, 8);

% Store the number of factors for each factor model
for i=1:nFactorModels
    nFactors(1,i) = size(res(i).weights,2);
end

for i=1:nFactorModels
    
    % Store the average signed squared Sharpe ratios
    panelAOutput(i,1) = mean( sign(all_sharpes_is(:,i)) .* (all_sharpes_is(:,i) .^2) );
    panelBOutput(i,1) = mean( sign(all_sharpes_os(:,i)) .* (all_sharpes_os(:,i) .^2) );
    
    % Run the comparisons
    for j=2:nFactorModels
        if i~=j
                        
            % Compare the number of factors to reward parsimony & count the
            % number of times model i beats model j
            if nFactors(i) >= nFactors(j)
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) > 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) > 0;
            else
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) >= 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) >= 0;
            end
            panelAOutput(i,j) = 100 * sum(iBetterThanJ_IS) / nRuns;
            panelBOutput(i,j) = 100 * sum(iBetterThanJ_OS) / nRuns;        
        end
    end
end

% Store the max in-sample Sharpe ratios for each sample in a matrix
maxISSharpes = max(all_sharpes_is,[],2);
maxISSharpesMat = repmat(maxISSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxIS = maxISSharpesMat == all_sharpes_is;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxIS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxIS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxIS matrix
    isMaxIS(r,:) = isMaxIS(r,:) & (nFactors == n);
end


% Store the max out-of-sample Sharpe ratios for each sample in a matrix
maxOSSharpes = max(all_sharpes_os,[],2);
maxOSSharpesMat = repmat(maxOSSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxOS = maxOSSharpesMat == all_sharpes_os;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxOS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxOS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxOS matrix
    isMaxOS(r,:) = isMaxOS(r,:) & (nFactors == n);
end

% Store the % in which each model has the max Sharpe
panelAOutput(:,end) = 100 * sum(isMaxIS,1)' / nRuns;
panelBOutput(:,end) = 100 * sum(isMaxOS,1)' / nRuns;

% Add a couple of empty columns
panelAOutput = [panelAOutput(:,1) nan(nFactorModels,1) panelAOutput(:,2:end-1)  nan(nFactorModels,1) panelAOutput(:,end)];
panelBOutput = [panelBOutput(:,1) nan(nFactorModels,1) panelBOutput(:,2:end-1 ) nan(nFactorModels,1) panelBOutput(:,end)];

% Define the headers
h = [res.label];
h = regexprep(h, '_c', '\\textsubscript{C}');

% Specify the format for each column to be passed to mat2Tex
colFormatSpecification = [{'%2.2f'} repmat({'%2.1f'},1,9)];

% Print the two panels
fprintf('Panel A:\n');
mat2Tex(panelAOutput, panelAOutput, h, 'colFormatSpec',colFormatSpecification);

fprintf('Panel B:\n');
mat2Tex(panelBOutput, panelBOutput, h, 'colFormatSpec',colFormatSpecification);

        
%% Print Table V - descriptive stats for mitigated factors

fprintf('\n\n\nTable 5 output:\n\n\n');

clear
clc


load factors_dnmv

% Create a constant for the regressions below
const = 0.01*ones(size(dates));

% Leave only the sS and mkt factors
indSsMkt = contains([factor_struct.label], {'_sS', 'mkt'});
indMktRep = strcmp([factor_struct.label], 'mkt_rep');
indToDrop = find(~indSsMkt | indMktRep);
factor_struct(indToDrop) = [];

% Store the number of factors
nFactors=length(factor_struct);

% Choose the sample
s = find(dates==197201);
e = find(dates==202112);

% Create the table output matrix
outputMatrix = nan(nFactors,6);

for i=1:nFactors
    % Gross return 
    grossRet = factor_struct(i).gross_factor(s:e); 
    res = ols(grossRet, const(s:e));
    outputMatrix(i,1) = res.beta;
    outputMatrix(i,2) = res.tstat;
    
    % Net return
    netRet = factor_struct(i).net_factor(s:e);
    res = ols(netRet, const(s:e));
    outputMatrix(i,3) = res.beta;
    outputMatrix(i,4) = res.tstat;    
    
    % Turnover
    TO = factor_struct(i).to(s:e,:);
    res = ols(TO, const(s:e));
    outputMatrix(i,5) = res.beta;
    outputMatrix(i,6) = (outputMatrix(i,1) - outputMatrix(i,3));
end

% Stor ethe headers for the factors
h = upper([factor_struct.label]');
h = regexprep(h, '_C', '\\textsubscript{C}');
h = regexprep(h, '_SS', '');

% Clean up the headers 
h(strcmp(h,{'UMD'}))={'MOM'};
h(strcmp(h,{'RME'}))={'ME'};
h(strcmp(h,{'RIA'}))={'IA'};
h(strcmp(h,{'RROE'}))={'ROE'};
h(strcmp(h,{'HMLM'}))={'HML(m)'};

% Add a few empty columns
outputMatrix = [outputMatrix(:,1:2) nan(length(factor_struct),1) outputMatrix(:,3:4) nan(length(factor_struct),1) outputMatrix(:,5:6)];

% Print the output table
mat2Tex(outputMatrix, outputMatrix, h, 2);

%% Table VI - mitigated results
 
fprintf('\n\n\nTable 6 output:\n\n\n');

clear
clc

load factors_dnmv
load fullSampleResults

% Panel A, full-sample results
fprintf('Panel A:\n');


% Find the location of the base models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};
loc_base = (ismember([res.label],factor_models));
res_base = res(loc_base);

% Find the location of the models with banding
factor_models = {'FF5_sS','FF6_sS','HXZ4_sS','BS6_sS','FF5_c_sS','FF6_c_sS'};
loc_banding = (ismember([res.label], factor_models));
res_banding = res(loc_banding);

% Get the output matrix and the header for the first line of Panel A
outputMatrix = [nan(1,2) [res_banding.netting_sharpe].^2];
h = {'SR$^2$ (Netting + Banding)'};

% Print the first line (squared Sharpes with banding & netting)
mat2Tex(outputMatrix, outputMatrix, h, 2);

% Store the output matrix for the rest of the table
outputMatrix = 100*([ ([res_banding.netting_sharpe].^2) ./ ([res_base.net_sharpe].^2);
                      ([res_base.netting_sharpe].^2)    ./ ([res_base.net_sharpe].^2);
                      ([res_banding.net_sharpe].^2)     ./ ([res_base.net_sharpe].^2) ] ...
                      -1);

% Store the headers
h={'Gain from Netting + Banding','Gain from Netting','Gain from Banding'};

% Printf the output matrix with % following the numbers
for i=1:size(outputMatrix,1)
    fprintf('\\multicolumn{3}{l}{%s} &  ',char(h(i)));
    for j=1:size(outputMatrix,2)
        fprintf(' %d\\%%% ',round(outputMatrix(i,j)));
        if j<size(outputMatrix,2)
            fprintf(' & ');
        end
    end
    fprintf('\\\\[2pt]\n');
end

fprintf('\n\n\n');

% Panels B & C

% Load the bootstrap results
load bootstrap_netting
load bootstrap_capm

% Find the bootstrap results for the mitigated models
mitModelsInd = find(contains([factor_model_defs.label],'_sS'));

% Store the number of runs
nRuns = min(length(capmIsSR), length(ISSharpes));

% Add the CAPM Sharpes
all_sharpes_is = [capmIsSR(1:nRuns) ISSharpes(1:nRuns, mitModelsInd)];
all_sharpes_os = [capmOosSR(1:nRuns) OSSharpes(1:nRuns, mitModelsInd)];

% Store the number of factor models
nFactorModels = size(all_sharpes_is, 2);

% Initialize the output matrices 
panelBOutput = nan(nFactorModels, 8);
panelCOutput = nan(nFactorModels, 8);

% Store the number of factors for each factor model
nFactors(1) = 1; % CAPM
for i=2:nFactorModels
    nFactors(1,i) = length( factor_model_defs( mitModelsInd(i-1) ).factors );
end

% Loop through the factor models
for i=1:nFactorModels
    
    % Store the average signed squared Sharpe ratios
    panelBOutput(i,1) = mean( sign(all_sharpes_is(:,i)) .* (all_sharpes_is(:,i) .^2) );
    panelCOutput(i,1) = mean( sign(all_sharpes_os(:,i)) .* (all_sharpes_os(:,i) .^2) );

    % Run the comparisons
    for j=2:nFactorModels
        if i~=j
                      
            % Compare the number of factors to reward parsimony & count the
            % number of times model i beats model j
            if nFactors(i) >= nFactors(j)
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) > 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) > 0;
            else
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) >= 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) >= 0;
            end
            panelBOutput(i,j) = 100 * sum(iBetterThanJ_IS) / nRuns;
            panelCOutput(i,j) = 100 * sum(iBetterThanJ_OS) / nRuns;        
        end
    end
end


% Store the max in-sample Sharpe ratios for each sample in a matrix
maxISSharpes = max(all_sharpes_is,[],2);
maxISSharpesMat = repmat(maxISSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxIS = maxISSharpesMat == all_sharpes_is;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxIS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxIS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxIS matrix
    isMaxIS(r,:) = isMaxIS(r,:) & (nFactors == n);
end



% Store the max out-of-sample Sharpe ratios for each sample in a matrix
maxOSSharpes = max(all_sharpes_os,[],2);
maxOSSharpesMat = repmat(maxOSSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxOS = maxOSSharpesMat == all_sharpes_os;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxOS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxOS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxOS matrix
    isMaxOS(r,:) = isMaxOS(r,:) & (nFactors == n);
end

% Store the % in which each model has the max Sharpe
panelBOutput(:,end) = 100 * sum(isMaxIS,1)' / nRuns;
panelCOutput(:,end) = 100 * sum(isMaxOS,1)' / nRuns;

% Add a couple of empty columns
panelBOutput = [panelBOutput(:,1) nan(nFactorModels,1) panelBOutput(:,2:end-1)  nan(nFactorModels,1) panelBOutput(:,end)];
panelCOutput = [panelCOutput(:,1) nan(nFactorModels,1) panelCOutput(:,2:end-1 ) nan(nFactorModels,1) panelCOutput(:,end)];

% Define the headers
h = [{'CAPM'}, [res.label]];
h = regexprep(h, '_c', '\\textsubscript{C}');

% Specify the format for each column to be passed to mat2Tex
colFormatSpecification = [{'%2.2f'} repmat({'%2.1f'},1,9)];

% Print the two panels
fprintf('Panel B:\n');
mat2Tex(panelBOutput, panelBOutput, h, 'colFormatSpec',colFormatSpecification);

fprintf('Panel C:\n');
mat2Tex(panelCOutput, panelCOutput, h, 'colFormatSpec',colFormatSpecification);

%% Online Appendix results below

%% Table IA.I - Anomaly frontier expansion summary statistics

fprintf('\n\n\nTable IA.1 output:\n\n\n');

clear
clc

% Load the results
load anomalyResults
load factors_dnmv

% Choose the factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};

% Create the headers
h = upper(factor_models');
h = regexprep(h,'_C','\\textsubscript{C}');

% Table B.1, Panel A
fprintf('Panel A: \n\n');
improvement = (gross_sharpe.^2) ./ (factor_gross_sharpe.^2) - 1;
outputMatrix = 100*[mean(improvement,1)' std(improvement,[],1)' prctile(improvement,[25 50 75 90 95 99],1)'];
mat2Tex(outputMatrix, outputMatrix, h, 1);
fprintf('\n\n\n');

% Table B.1, Panel B
fprintf('Panel B: \n\n');
improvement = (net_sharpe.^2) ./ (factor_net_sharpe.^2) - 1;
outputMatrix = 100*[mean(improvement,1)' std(improvement,[],1)' prctile(improvement,[25 50 75 90 95 99],1)'];
mat2Tex(outputMatrix, outputMatrix, h, 1);

fprintf('\n\n\n');

%% Table IA.II - bootrstraped banding results
 
fprintf('\n\n\nTable IA.2 output:\n\n\n');

clear
clc

load fullSampleResults

% Find the location of the models with banding
factor_models = {'FF5_sS','FF6_sS','HXZ4_sS','BS6_sS','FF5_c_sS','FF6_c_sS'};
loc_banding = find(ismember([res.label], factor_models));
res_banding = res(loc_banding);

% Get the output matrix and the header for the first line of Panel A
outputMatrix = [nan(1,2) [res_banding.net_sharpe].^2];
h = {'SR$^2$'};

% Print the first line (squared Sharpes with banding & netting)
mat2Tex(outputMatrix, outputMatrix, h, 2);
fprintf('\n\n\n');

% Panels B & C
% Load the stored results
load bootstrap_results
load bootstrap_capm

% Filter the bootstrap results to keep the basic models
indToKeep = ismember([bootstrap_results.label] ,factor_models);
res = bootstrap_results(indToKeep);

% Store the number of runs
nRuns = length(res(1).SR_is);

% Append CAPM results
res = [res(1) res];
res(1).weights = ones(nRuns,1);
res(1).SR_is = capmIsSR(1:min(size(capmIsSR,1),nRuns));
res(1).SR_os = capmOosSR(1:min(size(capmIsSR,1),nRuns));
res(1).label = {'CAPM'};


% Store the number of factor omdels
nFactorModels = length(res);

% Store the in- and out-of-sample Sharpes
all_sharpes_is = ([res.SR_is]);
all_sharpes_os = ([res.SR_os]);

% Initialize the output matrices 
panelBOutput = nan(nFactorModels, 8);
panelCOutput = nan(nFactorModels, 8);

% Store the number of factors for each factor model
for i=1:nFactorModels
    nFactors(1,i) = size(res(i).weights,2);
end

for i=1:nFactorModels
    
    % Store the average signed squared Sharpe ratios
    panelBOutput(i,1) = mean( sign(all_sharpes_is(:,i)) .* (all_sharpes_is(:,i) .^2) );
    panelCOutput(i,1) = mean( sign(all_sharpes_os(:,i)) .* (all_sharpes_os(:,i) .^2) );
    
    % Run the comparisons
    for j=2:nFactorModels
        if i~=j
                        
            % Compare the number of factors to reward parsimony & count the
            % number of times model i beats model j
            if nFactors(i) >= nFactors(j)
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) > 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) > 0;
            else
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) >= 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) >= 0;
            end
            panelBOutput(i,j) = 100 * sum(iBetterThanJ_IS) / nRuns;
            panelCOutput(i,j) = 100 * sum(iBetterThanJ_OS) / nRuns;        
        end
    end
end

% Store the max in-sample Sharpe ratios for each sample in a matrix
maxISSharpes = max(all_sharpes_is,[],2);
maxISSharpesMat = repmat(maxISSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxIS = maxISSharpesMat == all_sharpes_is;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxIS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxIS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxIS matrix
    isMaxIS(r,:) = isMaxIS(r,:) & (nFactors == n);
end


% Store the max out-of-sample Sharpe ratios for each sample in a matrix
maxOSSharpes = max(all_sharpes_os,[],2);
maxOSSharpesMat = repmat(maxOSSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxOS = maxOSSharpesMat == all_sharpes_os;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxOS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxOS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxOS matrix
    isMaxOS(r,:) = isMaxOS(r,:) & (nFactors == n);
end

% Store the % in which each model has the max Sharpe
panelBOutput(:,end) = 100 * sum(isMaxIS,1)' / nRuns;
panelCOutput(:,end) = 100 * sum(isMaxOS,1)' / nRuns;

% Add a couple of empty columns
panelBOutput = [panelBOutput(:,1) nan(nFactorModels,1) panelBOutput(:,2:end-1)  nan(nFactorModels,1) panelBOutput(:,end)];
panelCOutput = [panelCOutput(:,1) nan(nFactorModels,1) panelCOutput(:,2:end-1 ) nan(nFactorModels,1) panelCOutput(:,end)];

% Define the headers
h = [res.label];
h = regexprep(h, '_c', '\\textsubscript{C}');
h = regexprep(h, '_sS', '');

% Specify the format for each column to be passed to mat2Tex
colFormatSpecification = [{'%2.2f'} repmat({'%2.1f'},1,9)];

% Print the two panels
fprintf('Panel A:\n');
mat2Tex(panelBOutput, panelBOutput, h, 'colFormatSpec',colFormatSpecification);

fprintf('Panel C:\n');
mat2Tex(panelCOutput, panelCOutput, h, 'colFormatSpec',colFormatSpecification);


%% Table IA.III - bootstraped netting results
 
 
fprintf('\n\n\nTable IA.3 output:\n\n\n');

clear
clc

load fullSampleResults
load factors_dnmv

% Panel A, full-sample results
fprintf('Panel A:\n');

% Find the location of the base models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};
loc_base = find(ismember([res.label],factor_models));
res_base = res(loc_base);

% Get the output matrix and the header for the first line of Panel A
outputMatrix = [nan(1,2) [res_base.netting_sharpe].^2];
h = {'SR$^2$'};

% Print the first line (squared Sharpes with banding & netting)
mat2Tex(outputMatrix, outputMatrix, h, 2);
fprintf('\n\n\n');

% Load the bootstrap results
load bootstrap_netting
load bootstrap_capm

% Find the bootstrap results for the mitigated models
modelsInd = find(ismember([factor_model_defs.label], factor_models));

% Store the number of runs
nRuns = min(length(capmIsSR), length(ISSharpes));

% Add the CAPM Sharpes
all_sharpes_is = [capmIsSR(1:nRuns) ISSharpes(1:nRuns,modelsInd)];
all_sharpes_os = [capmOosSR(1:nRuns) OSSharpes(1:nRuns,modelsInd)];

% Store the number of factor models
nFactorModels = size(all_sharpes_is, 2);

% Initialize the output matrices 
panelBOutput = nan(nFactorModels, 8);
panelCOutput = nan(nFactorModels, 8);

% Store the number of factors for each factor model
nFactors(1) = 1; % CAPM
for i=2:nFactorModels
    nFactors(1,i) = length(factor_model_defs(modelsInd(i-1)).factors);
end

% Loop through the factor models
for i=1:nFactorModels
    
    % Store the average signed squared Sharpe ratios
    panelBOutput(i,1) = mean( sign(all_sharpes_is(:,i)) .* (all_sharpes_is(:,i) .^2) );
    panelCOutput(i,1) = mean( sign(all_sharpes_os(:,i)) .* (all_sharpes_os(:,i) .^2) );

    % Run the comparisons
    for j=2:nFactorModels
        if i~=j
                      
            % Compare the number of factors to reward parsimony & count the
            % number of times model i beats model j
            if nFactors(i) >= nFactors(j)
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) > 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) > 0;
            else
                iBetterThanJ_IS = all_sharpes_is(:,i) - all_sharpes_is(:,j) >= 0;
                iBetterThanJ_OS = all_sharpes_os(:,i) - all_sharpes_os(:,j) >= 0;
            end
            panelBOutput(i,j) = 100 * sum(iBetterThanJ_IS) / nRuns;
            panelCOutput(i,j) = 100 * sum(iBetterThanJ_OS) / nRuns;        
        end
    end
end


% Store the max in-sample Sharpe ratios for each sample in a matrix
maxISSharpes = max(all_sharpes_is,[],2);
maxISSharpesMat = repmat(maxISSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxIS = maxISSharpesMat == all_sharpes_is;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxIS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxIS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxIS matrix
    isMaxIS(r,:) = isMaxIS(r,:) & (nFactors == n);
end



% Store the max out-of-sample Sharpe ratios for each sample in a matrix
maxOSSharpes = max(all_sharpes_os,[],2);
maxOSSharpesMat = repmat(maxOSSharpes, 1, nFactorModels);

% Compare the model Sharpes to the max in each sample
isMaxOS = maxOSSharpesMat == all_sharpes_os;

% Find the cases where we have more than one max
multipleMaxInd = find(sum(isMaxOS,2)>1);
nMultipleMax = length(multipleMaxInd);

for i=1:nMultipleMax
    % Find this sample with multiple max Sharpes
    r = multipleMaxInd(i);
    
    % Identify the models 
    maxSharpeModelInd = find(isMaxOS(r,:)==1);
    
    % Get the minimum number of factors across the max Sharpe models
    n = min(nFactors(maxSharpeModelInd));
    
    % Correct the isMaxOS matrix
    isMaxOS(r,:) = isMaxOS(r,:) & (nFactors == n);
end

% Store the % in which each model has the max Sharpe
panelBOutput(:,end) = 100 * sum(isMaxIS,1)' / nRuns;
panelCOutput(:,end) = 100 * sum(isMaxOS,1)' / nRuns;

% Add a couple of empty columns
panelBOutput = [panelBOutput(:,1) nan(nFactorModels,1) panelBOutput(:,2:end-1)  nan(nFactorModels,1) panelBOutput(:,end)];
panelCOutput = [panelCOutput(:,1) nan(nFactorModels,1) panelCOutput(:,2:end-1 ) nan(nFactorModels,1) panelCOutput(:,end)];

% Define the headers
h = [{'CAPM'}, [res.label]];
h = regexprep(h, '_c', '\\textsubscript{C}');

% Specify the format for each column to be passed to mat2Tex
colFormatSpecification = [{'%2.2f'} repmat({'%2.1f'},1,9)];

% Print the two panels
fprintf('Panel B:\n');
mat2Tex(panelBOutput, panelBOutput, h, 'colFormatSpec',colFormatSpecification);

fprintf('Panel C:\n');
mat2Tex(panelCOutput, panelCOutput, h, 'colFormatSpec',colFormatSpecification);



%% Table IA.IV - Break-even analysis

fprintf('\n\n\nTable IA.4 output:\n\n\n');

clear
clc

% This script finds the mulitple (c_multiple) applied to the transaction
% costs measure that would cause the Sharpe ratios of two competing models
% to be equal.
% The output is a matrix whos (i,j) entry is the c_multiple that makes
% model i's Sharpe ratio equal to model j's if Sharpe(i)<Sharpe(j) before
% discount.
% c_multiple = 0 if there is no reduction in costs that makes Sharpe(i) 
% as high as Sharpe(j). 
% c_multiple=1 if Sharpe(i) > Sharpe(j) even before discounting costs.

% Load the factors & factor model data
load factors_dnmv

% Select the sample
s = find(dates==197201);
e = find(dates==202112);

% Choose the factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};

% Store the number of factor models
nFactorModels = length(factor_models);

% Get the factor model data
for i=1:nFactorModels
    res(i,1) = getFactorModelData(factor_models(i), factor_model_defs, factor_struct, [s e]);
end

%If row Sharpe < column Sharpe, what c_mult needed to make a tie?
BreakEvenCmultipliers = nan(nFactorModels, nFactorModels);

% Loop over the factor models
for i=1:nFactorModels
    for j=1:nFactorModels
        if i~=j
            [c_scalar, success] = BreakEvenC(res(i).gross_factors, res(i).tc, res(j).gross_factors, res(j).tc, dates, [s e]);
            
            if res(i).net_sharpe < res(j).net_sharpe
                BreakEvenCmultipliers(i,j) = c_scalar;
            else
                BreakEvenCmultipliers(i,j) = 1;
            end
        end
    end
end

% Clean up the headers 
h = regexprep(factor_models,'_c','\\textsubscript{C}');

% Print the output
mat2Tex(BreakEvenCmultipliers, BreakEvenCmultipliers, h, 2);

%% Table IA.V - Ledoit and Wolf Sharpe ratio comparisons

clear
clc

fprintf('\n\n\nTable IA.5 output:\n\n\n');


% Load the OOS results
load OOS_results
load factors_dnmv
load dates

% Determine the sample
s = find(dates==197201);
e = find(dates==202112);

% Choose the factor models
factor_models = {'FF5','FF6','HXZ4','BS6','FF5_c','FF6_c'};

% Store the numer of factor models
nFactorModels = length(factor_models);

% Initialize the results matrices
diffs = nan(nFactorModels, nFactorModels); %Gross squared "Information ratio"
pvals = nan(nFactorModels, nFactorModels); %pvals for gross Spanning tests

diffs_banding = nan(nFactorModels, nFactorModels); %Gross squared "Information ratio"
pvals_banding = nan(nFactorModels, nFactorModels); %pvals for gross Spanning tests

diffs_netting = nan(nFactorModels, nFactorModels); %Gross squared "Information ratio"
pvals_netting = nan(nFactorModels, nFactorModels); %pvals for gross Spanning tests

diffs_netting_banding = nan(nFactorModels, nFactorModels); %Gross squared "Information ratio"
pvals_netting_banding = nan(nFactorModels, nFactorModels); %pvals for gross Spanning tests

% Loop over the models
for i=1:nFactorModels 
    for j=i:nFactorModels
        if i < j
            % Plain models
            rets = OSrets(s+120:e,[i j]);
            [~,pval,~,~,~,SR1hat,SR2hat]=sharpeHAC(rets);
            diffs(i,j) = 12*(SR2hat^2-SR1hat^2);
            pvals(i,j) = pval;
            
            % Banding
            rets = OSrets(s+120:e,[i j]+6);
            [~,pval,~,~,~,SR1hat,SR2hat]=sharpeHAC(rets);
            diffs_banding(i,j) = 12*(SR2hat^2-SR1hat^2);
            pvals_banding(i,j) = pval;
            
            % Netting
            rets = OSrets_netting(s+120:e,[i j]);
            [~,pval,~,~,~,SR1hat,SR2hat]=sharpeHAC(rets);
            diffs_netting(i,j) = 12*(SR2hat^2-SR1hat^2);
            pvals_netting(i,j) = pval;
            
            % Netting + Banding
            rets = OSrets_netting(s+120:e,[i j]+6);
            [~,pval,~,~,~,SR1hat,SR2hat]=sharpeHAC(rets);
            diffs_netting_banding(i,j) = 12*(SR2hat^2-SR1hat^2);
            pvals_netting_banding(i,j) = pval;
        end
    end
end
    
% Store and clean up the headers
h = factor_models';
h = regexprep(h,'_c','\\textsubscript{C}');
h = regexprep(h,'_sS1','');

% Print the output
fprintf('\n\n\nTable Ledoit Wolfe output w/o netting:\n\n\n');

a=[diffs nan(nFactorModels,1) diffs_banding];
pA=[pvals nan(nFactorModels,1) pvals_banding];

mat2Tex(a,pA,h,2,'(');

fprintf('\n\n\nTable Ledoit Wolfe output w netting:\n\n\n');

a=[diffs_netting nan(nFactorModels,1) diffs_netting_banding];
pA=[pvals_netting nan(nFactorModels,1) pvals_netting_banding];

mat2Tex(a,pA,h,2,'(');

%% Table IA.VI - replication summary

fprintf('\n\n\nTable IA.6 output:\n\n\n');

clear
clc

% Load the factor data
load factors_dnmv

% Define a constant to use in the regressions below
const = 0.01*ones(size(dates));

% Choose the factors we are comparing
indToChoose = ~contains([factor_struct.label],{'_rep','_sS','frk','_c'});
factorList = [factor_struct(indToChoose).label];

% Get thre replications
factorRepList = cellfun(@(x) strcat([x,'_rep']), factorList, 'UniformOutput',0);

% Store the number of factors
nFactorList = length(factorList);

% Find the locations of the two lists
factorListLoc = (ismember([factor_struct.label],factorList));
factorRepListLoc = (ismember([factor_struct.label],factorRepList));

% Get the gross returns
factorRets = [factor_struct(factorListLoc).gross_factor];
factorsRetRets = [factor_struct(factorRepListLoc).gross_factor];

% Initialize the output matrices
coeffs = nan(nFactorList,5);
tstats = nan(nFactorList,5);

% Determine the sample
s = find(dates==197201);
e = find(dates==202112);

for i = 1:nFactorList
    % Get the average return of the factor
    res1 = ols(factorRets(s:e,i),[const(s:e)]);
    
    % Get the average return of the factor replication
    res2 = ols(factorsRetRets(s:e,i),[const(s:e)]);
    
    % Regress the factor replication on the original factor
    res3 = ols(factorRets(s:e,i),[const(s:e) factorsRetRets(s:e,i)]);
    
    % Get the correlation between the two
    temp = corrcoef(factorRets(s:e,i),factorsRetRets(s:e,i));
    
    % Assign to the output matricess
    coeffs(i,:) = [res1.beta res2.beta res3.beta(2) res3.rbar temp(1,2)];
    tstats(i,:) = [res1.tstat res2.tstat res3.tstat(2) nan nan];
end

% Store the headers
h = upper(factorList);

% Clean up the headers 
h(strcmp(h,{'UMD'}))={'MOM'};
h(strcmp(h,{'RME'}))={'ME'};
h(strcmp(h,{'RIA'}))={'IA'};
h(strcmp(h,{'RROE'}))={'ROE'};
h(strcmp(h,{'HMLM'}))={'HML(m)'};

mat2Tex(coeffs,tstats,h,2);



%% Timekeeping


fprintf('\n\n\n\nDone with table printing @ %s\n\n\n',char(datetime('now')));
diary off


