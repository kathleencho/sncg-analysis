function OUT = lfp_analysis(matFile, varargin)
% LFP_ANALYSIS Analyze one bilateral mPFC LFP session.
%
% This function illustrates the revised fixed-baseline, simple-outcome-mean
% analysis. It intentionally contains no cohort-level inferential statistics.
%
% Required variables in the MAT file
% ----------------------------------
% rightLFP   : right-hemisphere LFP vector
% leftLFP    : left-hemisphere LFP vector
% Fs         : sampling frequency in Hz
% IAStartSec : initial-association start time in seconds from LFP onset
% trials     : rule-shift trial table or numeric matrix containing:
%              Trial, t_start, t_end, code
%              code 1 = Correct; code 2 = Error
%
% Signal processing
% -----------------
% 1. Zero-phase notch filtering at 50 Hz and harmonics below Nyquist.
% 2. Zero-phase second-order Butterworth band-pass filtering.
% 3. Hilbert power: abs(hilbert(filtered LFP)).^2.
% 4. Separate left/right normalization in dB to the fixed pre-IA reference
%    window from IAStartSec-120 to IAStartSec-60 s.
% 5. Mean normalized power during each decision epoch [t_start,t_end).
% 6. Bilateral value = mean of normalized left and right values.
%
% Within-session endpoint
% -----------------------
% For each band and trial window:
%   mean(all valid Error trials) - mean(all valid Correct trials)
%
% No equal-count matching, trial subsampling, resampling, or random draws
% are performed.
%
% Name-value options
% ------------------
% OutputDir     : output folder; empty disables file export
% ExpectedFs    : expected sampling rate; default 2048 Hz
% Include30To90 : include 30-90 Hz sensitivity band; default true
%
% Outputs
% -------
% OUT.trialBandPower : one row per trial x frequency band
% OUT.sessionSummary : Error-minus-Correct summaries by band and window
% OUT.primarySummary : First5 gamma_30_60 row
% OUT.referenceQC    : fixed-reference information for each band
% OUT.bandDefinitions: analyzed bands and roles
% OUT.validation     : compact input/output validation
% OUT.parameters     : resolved settings

p = inputParser;
p.addRequired('matFile', @(x) ischar(x) || isstring(x));
p.addParameter('OutputDir', "", @(x) ischar(x) || isstring(x));
p.addParameter('ExpectedFs', 2048, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.addParameter('Include30To90', true, ...
    @(x) islogical(x) || (isnumeric(x) && isscalar(x) && ismember(x,[0 1])));
p.parse(matFile, varargin{:});
P = p.Results;
P.Include30To90 = logical(P.Include30To90);

matFile = char(string(matFile));
[inputQC, D] = lfp_validate(matFile, [], 'ExpectedFs', P.ExpectedFs);

rightLFP = D.rightLFP;
leftLFP = D.leftLFP;
Fs = D.Fs;
IAStartSec = D.IAStartSec;
trials = D.trials;

nSamples = numel(rightLFP);
recordingDurationSec = nSamples / Fs;
referenceStartSec = IAStartSec - 120;
referenceEndSec = IAStartSec - 60;
[iRef1, iRef2] = local_time_indices( ...
    referenceStartSec, referenceEndSec, Fs, nSamples);

bands = local_band_definitions(P.Include30To90);

% Notch the continuous bilateral signals once before band-specific filters.
rightNotched = local_notch_chain(rightLFP, Fs);
leftNotched = local_notch_chain(leftLFP, Fs);

trialRows = cell(height(bands), 1);
referenceRows = cell(height(bands), 1);

for b = 1:height(bands)
    bandName = bands.BandName(b);
    lowHz = bands.BandLowHz(b);
    highHz = bands.BandHighHz(b);

    if highHz >= Fs/2
        error('Band %s exceeds Nyquist at Fs=%.12g Hz.', bandName, Fs);
    end

    [bf, af] = butter(2, [lowHz highHz] / (Fs/2), 'bandpass');
    rightFiltered = filtfilt(bf, af, rightNotched);
    leftFiltered = filtfilt(bf, af, leftNotched);

    rightPower = abs(hilbert(rightFiltered)).^2;
    leftPower = abs(hilbert(leftFiltered)).^2;

    rightReferencePower = mean(rightPower(iRef1:iRef2), 'omitnan');
    leftReferencePower = mean(leftPower(iRef1:iRef2), 'omitnan');

    if ~isfinite(rightReferencePower) || rightReferencePower <= 0 || ...
            ~isfinite(leftReferencePower) || leftReferencePower <= 0
        error('Reference-window power is invalid for band %s.', bandName);
    end

    nTrials = height(trials);
    right_dB = nan(nTrials,1);
    left_dB = nan(nTrials,1);
    bilateral_dB = nan(nTrials,1);
    decisionEpochSamples = nan(nTrials,1);
    validLFP = false(nTrials,1);

    for k = 1:nTrials
        [i1, i2] = local_time_indices( ...
            trials.t_start(k), trials.t_end(k), Fs, nSamples);

        pRight = mean(rightPower(i1:i2), 'omitnan');
        pLeft = mean(leftPower(i1:i2), 'omitnan');

        if isfinite(pRight) && pRight > 0 && ...
                isfinite(pLeft) && pLeft > 0
            right_dB(k) = 10*log10(pRight / rightReferencePower);
            left_dB(k) = 10*log10(pLeft / leftReferencePower);
            bilateral_dB(k) = mean([right_dB(k), left_dB(k)], 'omitnan');
            decisionEpochSamples(k) = i2 - i1 + 1;
            validLFP(k) = isfinite(bilateral_dB(k));
        end
    end

    trialRows{b} = table( ...
        trials.Trial, trials.TrialOrdinal, trials.Outcome, trials.code, ...
        trials.t_start, trials.t_end, trials.t_end-trials.t_start, ...
        repmat(bandName,nTrials,1), ...
        repmat(lowHz,nTrials,1), repmat(highHz,nTrials,1), ...
        repmat(bands.Role(b),nTrials,1), ...
        decisionEpochSamples, right_dB, left_dB, bilateral_dB, validLFP, ...
        trials.IsFirst5, trials.IsTrials6To10, trials.IsFirst10, ...
        'VariableNames', {'Trial','TrialOrdinal','Outcome','OutcomeCode', ...
        't_start','t_end','DecisionEpochDurationSec', ...
        'BandName','BandLowHz','BandHighHz','BandRole', ...
        'DecisionEpochSamples','Right_dB','Left_dB','Bilateral_dB', ...
        'ValidLFP','IsFirst5','IsTrials6To10','IsFirst10'});

    referenceRows{b} = table( ...
        bandName, lowHz, highHz, bands.Role(b), ...
        IAStartSec, referenceStartSec, referenceEndSec, iRef1, iRef2, ...
        rightReferencePower, leftReferencePower, ...
        'VariableNames', {'BandName','BandLowHz','BandHighHz','BandRole', ...
        'IAStartSec','ReferenceStartSec','ReferenceEndSec', ...
        'ReferenceStartIndex','ReferenceEndIndex', ...
        'RightReferencePower','LeftReferencePower'});
end

trialBandPower = vertcat(trialRows{:});
referenceQC = vertcat(referenceRows{:});
sessionSummary = local_session_summaries(trialBandPower, bands);

primaryMask = sessionSummary.BandName == "gamma_30_60" & ...
    sessionSummary.Window == "First5";
primarySummary = sessionSummary(primaryMask,:);
if height(primarySummary) ~= 1
    error('Expected exactly one primary First5 gamma_30_60 summary row.');
end

OUT = struct();
OUT.trialBandPower = trialBandPower;
OUT.sessionSummary = sessionSummary;
OUT.primarySummary = primarySummary;
OUT.referenceQC = referenceQC;
OUT.bandDefinitions = bands;
OUT.inputFile = string(matFile);
OUT.samplingFrequencyHz = Fs;
OUT.recordingDurationSec = recordingDurationSec;
OUT.parameters = P;
OUT.validation = struct('input', inputQC.input, 'output', table());

[outputQC, ~] = lfp_validate(matFile, OUT, 'ExpectedFs', P.ExpectedFs);
OUT.validation = outputQC;

outDir = string(P.OutputDir);
if strlength(outDir) > 0
    if ~isfolder(outDir)
        mkdir(outDir);
    end

    [~, baseName] = fileparts(matFile);
    prefix = fullfile(char(outDir), baseName);

    writetable(trialBandPower, [prefix '_trial_band_power.csv']);
    writetable(sessionSummary, [prefix '_session_error_minus_correct.csv']);
    writetable(primarySummary, [prefix '_primary_first5_gamma_30_60.csv']);
    writetable(referenceQC, [prefix '_reference_qc.csv']);
    writetable(OUT.validation.input, [prefix '_input_validation.csv']);
    writetable(OUT.validation.output, [prefix '_output_validation.csv']);
    save([prefix '_analysis.mat'], 'OUT', '-v7.3');
end
end

function bands = local_band_definitions(include30To90)
BandName = [ ...
    "delta_1_4"; ...
    "theta_4_8"; ...
    "alpha_8_12"; ...
    "beta_12_30"; ...
    "gamma_30_60"; ...
    "narrow40_38_42"];

BandLowHz = [1;4;8;12;30;38];
BandHighHz = [4;8;12;30;60;42];
Role = [ ...
    "exploratory"; ...
    "exploratory"; ...
    "exploratory"; ...
    "exploratory"; ...
    "primary"; ...
    "exploratory"];

if include30To90
    BandName(end+1,1) = "gamma_30_90";
    BandLowHz(end+1,1) = 30;
    BandHighHz(end+1,1) = 90;
    Role(end+1,1) = "sensitivity";
end

bands = table(BandName,BandLowHz,BandHighHz,Role);
end

function summary = local_session_summaries(T, bands)
windows = ["First5"; "Trials6To10"; "First10"; "AllRS"];
rows = cell(height(bands)*numel(windows),1);
r = 0;

for b = 1:height(bands)
    B = T(T.BandName == bands.BandName(b),:);

    for w = 1:numel(windows)
        r = r + 1;
        windowName = windows(w);

        switch windowName
            case "First5"
                keep = B.IsFirst5;
            case "Trials6To10"
                keep = B.IsTrials6To10;
            case "First10"
                keep = B.IsFirst10;
            case "AllRS"
                keep = true(height(B),1);
            otherwise
                error('Unknown analysis window: %s', windowName);
        end

        X = B(keep & B.ValidLFP,:);
        correctValues = X.Bilateral_dB(X.Outcome == "Correct");
        errorValues = X.Bilateral_dB(X.Outcome == "Error");

        nCorrect = numel(correctValues);
        nError = numel(errorValues);
        meanCorrect = mean(correctValues,'omitnan');
        meanError = mean(errorValues,'omitnan');
        validEstimate = nCorrect >= 1 && nError >= 1 && ...
            isfinite(meanCorrect) && isfinite(meanError);

        if validEstimate
            errorMinusCorrect = meanError - meanCorrect;
        else
            errorMinusCorrect = NaN;
        end

        rows{r} = table( ...
            bands.BandName(b), bands.BandLowHz(b), bands.BandHighHz(b), ...
            bands.Role(b), windowName, nCorrect, nError, ...
            meanCorrect, meanError, errorMinusCorrect, validEstimate, ...
            'VariableNames', {'BandName','BandLowHz','BandHighHz','BandRole', ...
            'Window','nCorrect','nError','MeanCorrect_dB','MeanError_dB', ...
            'ErrorMinusCorrect_dB','ValidEstimate'});
    end
end

summary = vertcat(rows{:});
end

function x = local_notch_chain(x, Fs)
% Apply zero-phase notches at 50 Hz and harmonics below Nyquist.
x = double(x(:));
for frequencyHz = 50:50:(Fs/2 - 1)
    w0 = frequencyHz / (Fs/2);
    qualityFactor = max(30, frequencyHz/1);
    [b,a] = iirnotch(w0, w0/qualityFactor);
    x = filtfilt(b,a,x);
end
end

function [i1, i2] = local_time_indices(t1, t2, Fs, nSamples)
% Convert half-open interval [t1,t2) to inclusive MATLAB indices.
if ~isfinite(t1) || ~isfinite(t2) || t1 < 0 || t2 <= t1
    error('Invalid time interval [%.12g, %.12g).', t1, t2);
end

i1 = floor(t1*Fs) + 1;
i2 = ceil(t2*Fs);

if i1 < 1 || i2 > nSamples || i2 < i1
    error(['Time interval [%.12g, %.12g) is outside the recording ' ...
        '(duration %.12g s).'], t1, t2, nSamples/Fs);
end
end
