%This script was originally developed for Cho et al., 2023.
%It is used without modification for the analyses reported in:
%"Sncg+ CCK basket cells enable prefrontal gamma modulation and cognitive
%flexibility"

% This code computes a single population activity vector for the period up
% to 10 seconds after the onset of digging on each trial

% It then computes the similarities between activity vectors corresponding
% to different trial types

% maximum duration of each time period to analyze
maxdur = 200;% 200 frames = 10 sec

% assume you have a matrix of times called 'times'
% each row represents one trial
% column 1 = trial START
% column 2 = dig START
% column 3 = dig done / trial END
% column 4 = iti START
% column 5 = iti END

% assume you also have a vector of 0s and 1s corresponding to whether each
% trial is correct (1) or incorrect (0) called 'outcomes'

% finally assume the first RS trial is in the variable 'firstrstrial'

N = size(times);
ntrials = N(1);

trialtype = zeros(1,ntrials);
trialtype(firstrstrial:ntrials) = 1;

% identify the different types of trials
iacorrect = find(outcomes & ~trialtype & (1:ntrials <= 5));
iaerror = find(~outcomes & ~trialtype & (1:ntrials <= 5));
rscorrect = find(outcomes & trialtype & (1:ntrials <= firstrstrial+4));
rserror = find(~outcomes & trialtype & (1:ntrials <= firstrstrial+4));

trialcat{1} = iacorrect;
trialcat{2} = iaerror;
trialcat{3} = rscorrect;
trialcat{4} = rserror;

N = size(raster);
ncells = N(1);

% for each trial compute the normalized population vector

clear popvect;

for i=1:4,
    for j=1:length(trialcat{i}),

        % for each post-dig period
        trialno = trialcat{i}(j);
        dur = min(maxdur, times(trialno,3) - times(trialno,2));
        frames = times(trialno,2):(times(trialno,2)+dur);
        totact = sum(raster(:,frames)');
        noact = find(~totact);
        tmp = mean(raster(:,frames)');
        
        % normalized population vector
        popvect{i}(j,:) = tmp / sqrt(sum(tmp.*tmp));
    end
end


% compute the similarity of population activity vectors and significant correlations between trials of various types

for i=1:4,
    for j=i:4,
        actsim(i,j) = 0;
        count = 0;
        for m=1:length(trialcat{i}),
            for n=1:length(trialcat{j}),
                if i==j && m==n,
                    continue;
                end
                count = count +1;
                actsim(i,j) = actsim(i,j) + sum(popvect{i}(m,:) .* popvect{j}(n,:));
            end
        end
        actsim(i,j) = actsim(i,j) / count;
    end
end
