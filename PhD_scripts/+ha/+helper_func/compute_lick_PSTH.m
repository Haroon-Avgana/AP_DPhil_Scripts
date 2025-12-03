function OUT = compute_lick_PSTH(event_times, align_times, edges, varargin)
% COMPUTEPSTH  Per-trial PSTH aligned to events, with average (counts & rate).
%
% OUT = computePSTH(event_times, align_times, edges, 'Name',Value,...)
%
% Inputs
%   event_times : [nEvents x 1] lick (or spike) timestamps (numeric or datetime)
%   align_times : [nTrials x 1] alignment timestamps (one per trial)
%   edges       : [1 x (B+1)] bin edges relative to 0 (e.g., -1:0.01:2)
%
% Name-Value (optional)
%   'SmoothSigma' : Gaussian sigma for smoothing the average, in the SAME
%                   units as edges (seconds if numeric, seconds if datetime).
%                   Default = 0 (no smoothing).
%   'Return'      : 'count' or 'rate' (affects OUT.avg). Default = 'rate'.
%
% Output struct OUT
%   .countsPerTrial  : [nTrials x B] counts for each trial
%   .avgCount        : [1 x B] mean counts across trials
%   .avgRate         : [1 x B] mean rate (counts / binWidth)
%   .avg             : [1 x B] alias of avgRate or avgCount per 'Return'
%   .binCenters      : [1 x B]
%   .binWidth        : scalar
%   .nTrialsUsed     : scalar (trials with non-NaN align time)
%
% Notes
% - Uses histcounts semantics: right edge of last bin is open.
% - If timestamps are datetime, differences are computed in seconds.

% -------------------- parse inputs --------------------
p = inputParser;
p.addParameter('SmoothSigma', 0, @(x)isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('Return', 'rate', @(s)ischar(s) || isstring(s));
p.parse(varargin{:});
sigma = double(p.Results.SmoothSigma);
ret   = lower(string(p.Results.Return));

% ensure columns
event_times = event_times(:);
align_times = align_times(:);

% drop NaN/NaT alignments (keep trial order for the rest)
if isdatetime(align_times)
    keep = ~isnat(align_times);
else
    keep = ~isnan(align_times);
end
align_times = align_times(keep);
nTrials = numel(align_times);

B = numel(edges)-1;
binWidth = diff(edges(1:2));
binCenters = edges(1:end-1) + binWidth/2;

countsPerTrial = zeros(nTrials, B);

% determine mode (datetime vs numeric)
use_datetime = isdatetime(event_times) && isdatetime(align_times);

% -------------------- compute per-trial counts --------------------
if nTrials>0
    % Pre-sort events for speed (not required, but helps if huge)
    event_times = sort(event_times);

    % vectorized loop via arrayfun → cell → mat
    if use_datetime
        diff_fun = @(t) seconds(event_times - t);
    else
        diff_fun = @(t) (event_times - t);
    end

    C = arrayfun(@(t) histcounts(diff_fun(t), edges), align_times, ...
                 'UniformOutput', false);
    countsPerTrial = vertcat(C{:});
end

% -------------------- averages --------------------
avgCount = mean(countsPerTrial, 1, 'omitnan');
avgRate  = avgCount / binWidth;

% optional smoothing (on both series)
if sigma > 0
    % sigma is in the same units as edges, convert to bins:
    sigma_bins = sigma / binWidth;
    % build Gaussian kernel (truncate at ±4σ)
    xk = -ceil(4*sigma_bins):ceil(4*sigma_bins);
    gk = exp(-0.5*(xk/sigma_bins).^2);
    gk = gk / sum(gk);
    avgCount = conv(avgCount, gk, 'same');
    avgRate  = conv(avgRate,  gk, 'same');
end

% choose default OUT.avg
if ret=="count"
    avgOut = avgCount;
else
    avgOut = avgRate;
end

% -------------------- pack output --------------------
OUT = struct( ...
    'countsPerTrial', countsPerTrial, ...
    'avgCount',       avgCount, ...
    'avgRate',        avgRate, ...
    'avg',            avgOut, ...
    'binCenters',     binCenters, ...
    'binWidth',       binWidth, ...
    'nTrialsUsed',    nTrials );
end
