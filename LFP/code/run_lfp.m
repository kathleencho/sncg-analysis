%% RUN_LFP
% Minimal example for the manuscript LFP analysis.
%
% The example MAT file contains bilateral mPFC LFP, sampling rate,
% initial-association start time, and rule-shift trial timing/outcome labels.
%
% This script performs only the per-session signal analysis. It does not
% run cohort statistics or generate manuscript figures.

clear; clc;

codeDir = fileparts(mfilename('fullpath'));
lfpDir = fileparts(codeDir);
matFile = fullfile(lfpDir, 'example_data', 'example_session.mat');

OUT = lfp_analysis(matFile);

disp(OUT.primarySummary);
