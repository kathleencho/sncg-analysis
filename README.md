# sncg-analysis

This repository contains MATLAB code and example data used to analyze:
- Calcium imaging data
- Local field potential (LFP) recordings

## Requirements
- MATLAB R2023a
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox

## Calcium imaging pipeline

The `calcium` folder contains individual trial activity vectors for an example mouse.

Each MATLAB file includes:

`popvect`: a 1×4 cell array containing population activity vectors for:
1. Correct trials within the first 5 IA trials
2. Error trials within the first 5 IA trials
3. Correct trials within the first 5 RS trials
4. Error trials within the first 5 RS trials

`actsim`: a 4×4 similarity matrix.

The upper triangular portion contains average similarities between trial types. For example, row 1 column 3 corresponds to similarity between IA correct and RS correct trials. The script `popactvectsim.m` demonstrates how activity vectors and similarity matrices were computed from binary activity rasters.

This analysis was previously described in Cho et al., 2023:
https://doi.org/10.5281/zenodo.7709805

## LFP analysis

The `LFP` folder contains MATLAB code and example data for the primary session-level outcome-related 30–60-Hz LFP measure used in the revised manuscript.

The primary analysis uses a fixed pre-initial-association reference window, decision-epoch power from dig onset to dig end, and outcome modulation across the first five rule-shift trials calculated as mean(Error) − mean(Correct) using all valid trials.

See `LFP/README.md` for the complete analysis description and instructions for running the example workflow.

### LFP analysis revision

During manuscript revision, the LFP analysis workflow was re-evaluated and the pre-revision implementation was superseded. The current repository contains the validated LFP analysis used for the revised manuscript.
