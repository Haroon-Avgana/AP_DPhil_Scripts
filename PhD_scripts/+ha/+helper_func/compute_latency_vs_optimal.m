function S = compute_latency_vs_optimal(lick_times, stim_on_times, reward_times, reward_avaliblility_time)

% COMPUTE_LATENCY_VS_OPTIMAL
%   For each trial, compute:
%     - latency_to_first_lick  : first lick after stim onset
%     - diff_from_optimal      : |reward_time - reward_availability_time|
%   Then return per-day medians.

% Inputs:
%   lick_times                : [nLicks x 1] numeric or datetime
%   stim_on_times             : [nTrials x 1] numeric or datetime
%   reward_times              : [nTrials x 1] numeric or datetime (NaN allowed)
%   reward_avaliblility_time  : [nTrials x 1] numeric or datetime (NaN allowed)

% Outputs (struct S):
%   .latency_to_first_lick     [nTrials x 1]
%   .diff_from_optimal_reward  [nTrials x 1]
%   .median_latency            scalar (median of latency_to_first_lick, omitnan)
%   .median_diff               scalar (median of diff_from_optimal_reward, omitnan)

    % Ensure column vectors
    lick_times               = lick_times(:);
    stim_on_times            = stim_on_times(:);
    reward_times             = reward_times(:);
    reward_avaliblility_time = reward_avaliblility_time(:);

    % Make sure the trial number matches
    if length(reward_times) ~= length(stim_on_times)
        nTrials = min(numel(reward_times), numel(stim_on_times));
    else
        nTrials = numel(stim_on_times);
    end
    
    latency_to_first_lick    = nan(nTrials, 1);
    diff_from_optimal_reward = nan(nTrials, 1);

    for tr = 1:nTrials
        % First-lick latency after stim onset
        dt = lick_times - stim_on_times(tr);
        dt = dt(dt > 0);

        if ~isempty(dt) % calculate latency to lick from stim onset
            latency_to_first_lick(tr) = min(dt);
        end

        % Reward - availability time latency calculation
        if ~isnan(reward_times(tr)) && ~isnan(reward_avaliblility_time(tr))
            diff_from_optimal_reward(tr) = abs(reward_times(tr) - reward_avaliblility_time(tr));
        end
    end

    % Calculate medians (NaN-safe)
    median_lat  = median(latency_to_first_lick, 'omitnan');
    median_diff = median(diff_from_optimal_reward, 'omitnan');

    % Pack outputs
    S = struct( ...
        'latency_to_first_lick',    latency_to_first_lick, ...
        'diff_from_optimal_reward', diff_from_optimal_reward, ...
        'median_latency',           median_lat, ...
        'median_diff',              median_diff);
end