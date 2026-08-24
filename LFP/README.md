# LFP analysis

Minimal MATLAB code for the primary per-session LFP measure used in the manuscript.

The analysis starts from one self-contained `.mat` file containing bilateral mPFC LFP, sampling rate, the initial-association start time, and rule-shift trial timing/outcome labels.

Core analysis:

1. zero-phase 50-Hz notch filtering and harmonics below Nyquist;
2. zero-phase second-order Butterworth filtering at 30–60 Hz;
3. Hilbert power;
4. separate left/right dB normalization to the fixed pre-IA reference window (IA start −120 to −60 s);
5. mean power during each decision epoch (dig onset to dig end);
6. bilateral averaging after left/right normalization;
7. First5 outcome modulation = **mean(Error) − mean(Correct)** using all valid trials.

## Files

```text
LFP/
├── README.md
├── code/
│   ├── run_lfp.m
│   ├── lfp_analysis.m
│   └── lfp_validate.m
└── example_data/
    └── example_session.mat
```

Run `code/run_lfp.m` in MATLAB. `example_session.mat` contains a cropped segment of an experimental LFP recording provided to demonstrate the analysis workflow. The example is not the complete dataset and is not used independently for manuscript-level inference.

The manuscript recordings were sampled at 2048 Hz. Signal Processing Toolbox functions `butter`, `filtfilt`, `iirnotch`, and `hilbert` are required.
