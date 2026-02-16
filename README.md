# Sncg+ CCK basket cells enable prefrontal gamma modulation and cognitive flexibility

This repository contains MATLAB code and example datasets used to analyze:
- Calcium imaging data
- Local Field Potential (LFP) recordings using Hilbert amplitude envelope analysis

## Requirements
- MATLAB R2023a
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Calcium imaging pipeline
The 'calcium' folder contains individual trial activity vectors for an example mouse.  
Each matlab file includes:

popvect: a 1x4 cell array containing population activity vectors for:
1. Correct trials within the first 5 IA trials
2. Error trials within the first 5 IA trials
3. Correct trials within the first 5 RS trials
4. Error trials within the first 5 RS trials

actsim: a 4x4 similarity matrix
The upper triangular portion contains average similarities between trial types. For example row 1 column 3 corresponds to similarity between IA correct and RS correct trials. The script popactvectsim.m demonstrates how activity vectors and similarity matrices were computed from binary activity rasters.
This analysis was previously described in Cho et al., 2023: https://doi.org/10.5281/zenodo.7709805

## LFP processing pipeline
The 'LFP' folder contains example session-level data and scripts used to compute task-evoked broadband gamma power. 
LFP processing steps:
1. Notch filtering at 50 Hz
2. Band-pass filtering (30–60 Hz)
3. Hilbert transform to obtain the analytic amplitude
4. Power estimation as squared amplitude
5. Baseline normalization using a 60 s pre-task baseline
6. Within-session task modulation computed as error – correct trial power

Example data are provided in example_data/example_lfp_summary.mat, which contains:
T_all: table containing band-level LFP summary metrics (including error, correct, and error-correct modulation values)
baselineTable: table containing baseline window timing information

The script run_lfp.m reproduces the core analysis pipeline using the example dataset.
