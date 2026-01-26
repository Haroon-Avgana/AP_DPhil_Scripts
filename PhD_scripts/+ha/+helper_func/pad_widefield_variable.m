function stacked_data = pad_widefield_variable(widefield_cat, field_name, verbose)
    % PAD_WIDEFIELD_VARIABLE Extracts and pads widefield data with NaNs for missing days
    %
    % Inputs:
    %   widefield_cat - concatenated widefield structure array
    %   field_name - string, name of field to extract (e.g., 'right_stim_aligned_V')
    %   verbose - (optional) boolean, whether to print diagnostics (default: true)
    %
    % Output:
    %   stacked_data - V×T×days array with NaN padding for missing days
    %
    % Example:
    %   widefield_cat = cat(2, passive_data.widefield);
    %   right_stim_data = pad_widefield_variable(widefield_cat, 'right_stim_aligned_V');
    %   left_stim_data = pad_widefield_variable(widefield_cat, 'left_stim_aligned_V', false);
    
    if nargin < 3
        verbose = true;
    end
    
    % Extract field from all days
    all_data = {widefield_cat.(field_name)};
    
    % Find first non-empty day for template dimensions
    non_empty_idx = find(~cellfun(@isempty, all_data), 1, 'first');
    
    if isempty(non_empty_idx)
        error('No valid data found in field: %s', field_name);
    end
    
    % Get V and T dimensions from non-empty example
    example_data = all_data{non_empty_idx};
    n_V = size(example_data, 1);
    n_T = size(example_data, 2);
    
    % Process each day: average across trials or pad with NaNs
    processed_data = cell(1, length(all_data));
    empty_count = 0;
    
    for day = 1:length(all_data)
        if isempty(all_data{day})
            % Pad empty days with NaN array
            processed_data{day} = nan(n_V, n_T);
            empty_count = empty_count + 1;
        else
            % Average across trials (3rd dimension)
            processed_data{day} = mean(all_data{day}, 3);
        end
    end
    
    % Concatenate - all days have consistent V×T dimensions
    stacked_data = cat(3, processed_data{:});
    
    % Print diagnostics if requested
    if verbose
        fprintf('\n===== %s =====\n', field_name);
        fprintf('Total days: %d\n', size(stacked_data, 3));
        fprintf('Dimensions: %d×%d×%d (V×T×days)\n', size(stacked_data));
        fprintf('Valid days: %d\n', length(all_data) - empty_count);
        fprintf('Empty days (padded): %d\n', empty_count);
        
        if empty_count > 0
            empty_days = find(all(all(isnan(stacked_data), 1), 2));
            fprintf('Empty day indices: %s\n', mat2str(empty_days));
        end
    end
end