
clear
clc

restoredefaultpath; % Start with the default path
matlabPackagePath='D:\Published Repos\Detzel-Novy-Marx-Velikov\';
paperCodePath='D:\Published Repos\Detzel-Novy-Marx-Velikov\Detzel-Novy-Marx-Velikov\';

% Add the relevant folders (with subfolders) to the path
addpath(genpath([matlabPackagePath,'Data']))
addpath(genpath([matlabPackagePath,'Functions']))
addpath(genpath([paperCodePath]))
cd(paperCodePath)

%% Add the directories

if ~exist([pwd,'Data'], 'dir')
    mkdir(['Data'])
end
if ~exist([pwd,'Results'], 'dir')
    mkdir(['Results'])
end
if ~exist([pwd,'Figures'], 'dir')
    mkdir(['Figures'])
end
addpath(genpath(pwd));

%% Start a log file

startLogFile([paperCodePath], 'dnmv')

%% These three scripts just make and save the main factors and tcosts

run('make_main_factors_tcosts.m')
run('make_mitigated_factors.m')
run('make_freak_factor.m')

%% This script organizes the factors in a coherent way

run('organize_factors.m');

%% This script runs the full sample results

run('run_full_sample_results.m');

%% This script runs the bootstrap results

run('run_bootstraps.m');

%% This script runs the out of sample MVE strategy results

run('run_oos_results.m');

%% This script runs the anomaly results

run('run_anomaly_results.m');

%%  This script prints all tables

run('make_tables.m');

%% This script makes all figures 

run('make_figures.m');

%% End the log file

diary('off');