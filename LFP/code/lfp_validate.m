function [Q, D] = lfp_validate(matFile, OUT, varargin)
% LFP_VALIDATE Validate public LFP example input and analysis output.
%
% Usage
% -----
% [Q,D] = lfp_validate(matFile, [], 'ExpectedFs', 2048)
% [Q,D] = lfp_validate(matFile, OUT, 'ExpectedFs', 2048)
%
% Q.input  : input-validation table
% Q.output : output-validation table (empty when OUT is empty)
%
% D contains normalized variables used by lfp_analysis:
%   rightLFP, leftLFP, Fs, IAStartSec, trials

p = inputParser;
p.addRequired('matFile', @(x) ischar(x) || isstring(x));
p.addRequired('OUT');
p.addParameter('ExpectedFs', 2048, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
p.parse(matFile, OUT, varargin{:});

expectedFs = double(p.Results.ExpectedFs);
matFile = char(string(matFile));

%% ========================================================================
% INPUT VALIDATION
% ========================================================================

if ~isfile(matFile)
    error('MAT file not found:\n%s', matFile);
end

S = load(matFile);

required = {'rightLFP','leftLFP','Fs','IAStartSec','trials'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('MAT file is missing required variable "%s".', required{k});
    end
end

rightLFP = double(S.rightLFP(:));
leftLFP  = double(S.leftLFP(:));
Fs = double(S.Fs);
IAStartSec = double(S.IAStartSec);

if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    error('Fs must be a finite positive scalar.');
end

if abs(Fs - expectedFs) > 1e-9
    error('Fs is %.12g Hz; expected %.12g Hz.', Fs, expectedFs);
end

if ~isscalar(IAStartSec) || ~isfinite(IAStartSec)
    error('IAStartSec must be a finite scalar.');
end

if isempty(rightLFP) || isempty(leftLFP)
    error('LFP vectors must be nonempty.');
end

if numel(rightLFP) ~= numel(leftLFP)
    error('rightLFP and leftLFP must have the same number of samples.');
end

if any(~isfinite(rightLFP)) || any(~isfinite(leftLFP))
    error('LFP vectors must contain only finite values.');
end

if Fs/2 <= 90
    error('Sampling frequency is too low for the requested frequency bands.');
end

trials = local_validate_trials(S.trials);

if height(trials) < 5
    error('At least five rule-shift trials are required.');
end

nSamples = numel(rightLFP);
recordingDurationSec = nSamples / Fs;

referenceStartSec = IAStartSec - 120;
referenceEndSec   = IAStartSec - 60;

if referenceStartSec < 0
    error('IAStartSec is too early for the IA-120 to IA-60 s reference window.');
end

local_time_indices(referenceStartSec, referenceEndSec, Fs, nSamples);

for k = 1:height(trials)
    local_time_indices(trials.t_start(k), trials.t_end(k), Fs, nSamples);
end

D = struct();
D.rightLFP = rightLFP;
D.leftLFP = leftLFP;
D.Fs = Fs;
D.IAStartSec = IAStartSec;
D.trials = trials;

inputTable = table( ...
    string(matFile), ...
    Fs, ...
    numel(rightLFP), ...
    numel(leftLFP), ...
    recordingDurationSec, ...
    IAStartSec, ...
    referenceStartSec, ...
    referenceEndSec, ...
    height(trials), ...
    true, ...
    'VariableNames', { ...
    'MATFile', ...
    'Fs', ...
    'RightSamples', ...
    'LeftSamples', ...
    'RecordingDurationSec', ...
    'IAStartSec', ...
    'ReferenceStartSec', ...
    'ReferenceEndSec', ...
    'TrialCount', ...
    'AllChecksPassed'});

Q = struct();
Q.input = inputTable;
Q.output = table();

%% ========================================================================
% OUTPUT VALIDATION
% ========================================================================

if isempty(OUT)
    return;
end

requiredFields = { ...
    'trialBandPower', ...
    'sessionSummary', ...
    'primarySummary', ...
    'referenceQC', ...
    'bandDefinitions', ...
    'samplingFrequencyHz', ...
    'recordingDurationSec'};

for k = 1:numel(requiredFields)
    if ~isfield(OUT, requiredFields{k})
        error('Analysis output is missing field "%s".', requiredFields{k});
    end
end

T = OUT.trialBandPower;
Ssum = OUT.sessionSummary;
Psum = OUT.primarySummary;
R = OUT.referenceQC;
B = OUT.bandDefinitions;

if ~istable(T) || ~istable(Ssum) || ~istable(Psum) || ...
        ~istable(R) || ~istable(B)
    error('Expected analysis outputs are not MATLAB tables.');
end

nBands = height(B);
expectedTrialRows = height(trials) * nBands;
trialRowCountOK = height(T) == expectedTrialRows;

validRowsFinite = true;
if ismember('ValidLFP', T.Properties.VariableNames) && ...
        ismember('Bilateral_dB', T.Properties.VariableNames)
    use = logical(T.ValidLFP);
    validRowsFinite = all(isfinite(T.Bilateral_dB(use)));
end

requiredWindows = ["First5","Trials6To10","First10","AllRS"];
expectedSummaryRows = nBands * numel(requiredWindows);
summaryRowCountOK = height(Ssum) == expectedSummaryRows;

referenceRowCountOK = height(R) == nBands;

primaryRowCountOK = height(Psum) == 1;

primaryBandOK = false;
primaryWindowOK = false;
primaryEstimateFinite = false;

if primaryRowCountOK
    if ismember('BandName', Psum.Properties.VariableNames)
        primaryBandOK = string(Psum.BandName(1)) == "gamma_30_60";
    end
    if ismember('Window', Psum.Properties.VariableNames)
        primaryWindowOK = string(Psum.Window(1)) == "First5";
    end
    if ismember('ErrorMinusCorrect_dB', Psum.Properties.VariableNames)
        primaryEstimateFinite = isfinite(Psum.ErrorMinusCorrect_dB(1));
    end
end

samplingRateOK = ...
    isscalar(OUT.samplingFrequencyHz) && ...
    isfinite(OUT.samplingFrequencyHz) && ...
    abs(double(OUT.samplingFrequencyHz)-Fs) <= 1e-9;

durationOK = ...
    isscalar(OUT.recordingDurationSec) && ...
    isfinite(OUT.recordingDurationSec) && ...
    abs(double(OUT.recordingDurationSec)-recordingDurationSec) <= 1/Fs;

allChecksPassed = ...
    trialRowCountOK && ...
    validRowsFinite && ...
    summaryRowCountOK && ...
    referenceRowCountOK && ...
    primaryRowCountOK && ...
    primaryBandOK && ...
    primaryWindowOK && ...
    primaryEstimateFinite && ...
    samplingRateOK && ...
    durationOK;

outputTable = table( ...
    trialRowCountOK, ...
    validRowsFinite, ...
    summaryRowCountOK, ...
    referenceRowCountOK, ...
    primaryRowCountOK, ...
    primaryBandOK, ...
    primaryWindowOK, ...
    primaryEstimateFinite, ...
    samplingRateOK, ...
    durationOK, ...
    allChecksPassed, ...
    'VariableNames', { ...
    'TrialRowCountOK', ...
    'ValidRowsFinite', ...
    'SummaryRowCountOK', ...
    'ReferenceRowCountOK', ...
    'PrimaryRowCountOK', ...
    'PrimaryBandOK', ...
    'PrimaryWindowOK', ...
    'PrimaryEstimateFinite', ...
    'SamplingRateOK', ...
    'RecordingDurationOK', ...
    'AllChecksPassed'});

Q.output = outputTable;

if ~allChecksPassed
    disp(outputTable);
    error('LFP output validation failed.');
end

end


%% ========================================================================
% LOCAL HELPERS
% ========================================================================

function trials = local_validate_trials(x)

if istable(x)
    raw = x;
elseif isnumeric(x)
    if size(x,2) == 3
        raw = table((1:size(x,1))',x(:,1),x(:,2),x(:,3), ...
            'VariableNames',{'Trial','t_start','t_end','code'});
    elseif size(x,2) >= 4
        raw = table(x(:,1),x(:,2),x(:,3),x(:,4), ...
            'VariableNames',{'Trial','t_start','t_end','code'});
    else
        error(['Numeric trials must have columns [t_start t_end code] ' ...
            'or [Trial t_start t_end code].']);
    end
else
    error('trials must be a table or numeric matrix.');
end

names = string(raw.Properties.VariableNames);
keys = lower(regexprep(names,'[^A-Za-z0-9]',''));

trialIdx = find(ismember(keys, ...
    ["trial","trialnumber","trialid","trialindex"]),1);
startIdx = find(ismember(keys, ...
    ["tstart","start","digstart","diggingstart"]),1);
endIdx = find(ismember(keys, ...
    ["tend","end","digend","diggingend"]),1);
codeIdx = find(ismember(keys, ...
    ["code","outcomecode","correcterror","outcome"]),1);

if isempty(startIdx) || isempty(endIdx) || isempty(codeIdx)
    error('trials must contain t_start, t_end, and code/outcome columns.');
end

n = height(raw);

if isempty(trialIdx)
    Trial = (1:n)';
else
    Trial = local_numeric(raw{:,trialIdx},'Trial');
end

t_start = local_numeric(raw{:,startIdx},'t_start');
t_end   = local_numeric(raw{:,endIdx},'t_end');

rawCode = raw{:,codeIdx};
if isnumeric(rawCode) || islogical(rawCode)
    code = double(rawCode(:));
else
    s = lower(strtrim(string(rawCode(:))));
    code = nan(n,1);
    code(ismember(s,["1","correct","c"])) = 1;
    code(ismember(s,["2","error","incorrect","e"])) = 2;
end

if any(~ismember(code,[1 2]))
    error('Trial outcomes must decode to 1=Correct or 2=Error.');
end

if any(t_start < 0) || any(t_end <= t_start)
    error('Trial timing contains invalid intervals.');
end

% Preserve supplied chronological order. TrialOrdinal is the RS position
% within the supplied rule-shift trial table.
TrialOrdinal = (1:n)';

Outcome = strings(n,1);
Outcome(code==1) = "Correct";
Outcome(code==2) = "Error";

IsFirst5 = TrialOrdinal <= 5;
IsTrials6To10 = TrialOrdinal >= 6 & TrialOrdinal <= 10;
IsFirst10 = TrialOrdinal <= 10;

trials = table( ...
    Trial, TrialOrdinal, Outcome, code, ...
    t_start, t_end, ...
    IsFirst5, IsTrials6To10, IsFirst10);

end


function x = local_numeric(v,label)

if isnumeric(v) || islogical(v)
    x = double(v(:));
else
    x = str2double(strtrim(string(v(:))));
end

if any(~isfinite(x))
    error('%s contains nonnumeric or nonfinite values.',label);
end

end


function [i1,i2] = local_time_indices(t1,t2,Fs,nSamples)
% Half-open interval [t1,t2) represented by inclusive MATLAB indices.

if ~isfinite(t1) || ~isfinite(t2) || t1 < 0 || t2 <= t1
    error('Invalid time interval [%.12g, %.12g).',t1,t2);
end

i1 = floor(t1*Fs)+1;
i2 = ceil(t2*Fs);

if i1 < 1 || i2 > nSamples || i2 < i1
    error(['Time interval [%.12g, %.12g) is outside the recording ' ...
        '(duration %.12g s).'],t1,t2,nSamples/Fs);
end

end
