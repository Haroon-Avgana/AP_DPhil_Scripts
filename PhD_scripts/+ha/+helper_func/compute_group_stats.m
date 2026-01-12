%% Hierarchical averaging: Sessions → Pre/Post → Animals → Group

function [group_mean, group_sem] = compute_group_stats(group_sessions, condition)
    % condition: 'plus_pre', 'plus_post', 'minus_pre', 'minus_post'
    
    n_animals = length(group_sessions);
    animal_means = [];
    animal_sems = [];
    
    for a = 1:n_animals
        sessions = group_sessions(a).(condition);
        
        if isempty(sessions)
            continue;
        end
        
        % Step 1: Average across sessions for this animal
        animal_mean = mean(sessions, 1, 'omitnan');
        
        % Step 2: SEM across sessions for this animal
        animal_sem = std(sessions, 0, 1, 'omitnan') / sqrt(size(sessions, 1));
        
        animal_means = [animal_means; animal_mean];
        animal_sems = [animal_sems; animal_sem];
    end
    
    % Step 3: Average across animals
    group_mean = mean(animal_means, 1, 'omitnan');
    
    % Step 4: Propagate SEM across animals
    % Pooled variance: mean of squared SEMs (converts SEM back to variance, averages, converts back)
    pooled_var = mean(animal_sems.^2, 1, 'omitnan');
    group_sem = sqrt(pooled_var) / sqrt(size(animal_means, 1));
end
