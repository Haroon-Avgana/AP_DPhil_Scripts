
% a decoder that uses mPFC activity to classify whether the trial is
% classical or operant


% Download Matlab machine learning package

%%

% Variables needed:

% X - nTrials x nFeatures (mPFC ROI in a given window)
% y - nTrials x 1 : 0 classical, 1 operant
% sessionId nTrials x 1 the session number (i.e. day)
% mouseId nTrials x1 the mouse ID

% Have a mouse dummy regressor
mouseDummy = dummyvar(categorical(mouseId));   % nTrials × nMice   (full rank)
X_aug      = [X  mouseDummy];                  % append to feature set



% Build a session pair CV folds

foldStruct = {};     % each cell will hold a logical test mask
mice = unique(mouseId);

for m = mice'
    % classical & operant sessions for this mouse
    clSess = unique(sessionId(mouseId==m & y==0));
    opSess = unique(sessionId(mouseId==m & y==1));
    
    for c = clSess'
        for o = opSess'
            testMask = (sessionId==c) | (sessionId==o);
            foldStruct{end+1} = testMask;
        end
    end
end
nFolds = numel(foldStruct);

% Helper for class-balanced weights
balWeights = @(idx) ...
    (y(idx)==0) .* 0.5/sum(idx & y==0) + (y(idx)==1).*0.5/sum(idx & y==1);


% Cross validation with permutation option


doPerm   = false;          % set true for null distribution
nPerm    = 1000;           % permutations if needed
accTrue  = NaN(nFolds,1);
aucTrue  = NaN(nFolds,1);

for f = 1:nFolds
    testIdx  =  foldStruct{f};    % logical vector
    trainIdx = ~foldStruct{f};
    
    yTrain = y;
    if doPerm
        % shuffle labels *within mouse & session* in training set only
        for sid = unique(sessionId(trainIdx))'
            sub = trainIdx & sessionId==sid;
            yTrain(sub) = yTrain(sub(randperm(sum(sub))));
        end
    end
    
    w = balWeights(trainIdx);
    
    mdl = fitclinear(X_aug(trainIdx,:), yTrain(trainIdx), ...
                     'Learner','logistic', ...
                     'Regularization','ridge', ...
                     'Lambda',0, ...
                     'Weights', w);
    
    [yHat,score] = predict(mdl, X_aug(testIdx,:));
    
    % balanced accuracy
    accTrue(f) = 0.5*(sum(yHat==1 & y(testIdx)==1)/sum(y(testIdx)==1) + ...
                      sum(yHat==0 & y(testIdx)==0)/sum(y(testIdx)==0));
    
    [~,~,~,aucTrue(f)] = perfcurve(y(testIdx), score(:,2), 1);
end

fprintf('\nBalanced accuracy  : %.3f  ±  %.3f (s.e.)\n', ...
        mean(accTrue), std(accTrue)/sqrt(nFolds));
fprintf('ROC–AUC            : %.3f  ±  %.3f (s.e.)\n', ...
        mean(aucTrue), std(aucTrue)/sqrt(nFolds));