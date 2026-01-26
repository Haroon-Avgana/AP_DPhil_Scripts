function activity = compute_activity(traces, response_window)
    % traces: n_days × T
    % Returns: n_days × 1 (max activity in window per day)
    activity = max(traces(:, response_window), [], 2);
end
