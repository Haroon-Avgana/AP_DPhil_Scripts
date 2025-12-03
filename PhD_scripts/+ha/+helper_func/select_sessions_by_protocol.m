function [keptLocalIdx] = select_sessions_by_protocol(rd, opts)
% Return original indices to keep from this animal's recording_day array.
% rd(i).workflow : char/string
% rd(i).day      : datetime (or parseable to datetime)

    if isempty(rd)
        keptLocalIdx = [];
        return
    end

    % --- collect & sort by day ---
    workflows_animal = string({rd.workflow})';
    days_animal      = {rd.day}';
    dates_cat = datetime(vertcat(days_animal{:}), 'InputFormat','yyyy-MM-dd');
    [~, sort_idx] = sort(dates_cat);
    workflows_sorted = workflows_animal(sort_idx);

    % --- map names to a stable index using your master list ---
    [found, ~] = ismember(workflows_sorted, string(opts.protocol_list));
    if ~all(found)
        missing = unique(workflows_sorted(~found));
        msg = sprintf('Unknown workflows: %s', strjoin(cellstr(missing), ', '));
        if opts.error_on_unknown, error(msg); else, warning(msg); end
    end

    N = numel(workflows_sorted);
    keepMask = true(N,1);

    % --- drop all big_stim sessions ---
    bigMask = ismember(workflows_sorted, string(opts.big_names));
    keepMask(bigMask) = false;

    % --- optionally drop "following" protocols after the first big_stim ---
    if opts.drop_following
        firstBig = find(bigMask, 1, 'first');
        if ~isempty(firstBig)
            followMask = ismember(workflows_sorted, string(opts.following_names)) ...
                         & ((1:N)' > firstBig);
            keepMask(followMask) = false;
        end
    end

    % map kept sorted positions back to original indices
    keptLocalIdx = sort_idx(keepMask);
end