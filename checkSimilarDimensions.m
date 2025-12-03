function areSimilar = checkSimilarDimensions(C)
% Check if dimensions of all elements in cell array C are similar

% Get dimensions of first element
dims = size(C{1});

% Compare dimensions of all other elements
for i = 2:numel(C)
    if ~isequal(size(C{i}), dims)
        areSimilar = false;
        return;
    end
end

% If all dimensions are similar, return true
areSimilar = true;
end