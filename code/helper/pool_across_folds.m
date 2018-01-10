function roi = pool_across_folds(roi1, roi2)

roi = roi1;
for ee = 1:length(roi1)
    % average tc properties across folds
    roi(ee).experiments = [roi1(ee).experiments; roi2(ee).experiments];
    roi(ee).trial_predsS = cellfun(@(X, Y) (X + Y) ./ 2, ...
        roi1(ee).trial_predsS, roi2(ee).trial_predsS, 'uni', false);
    roi(ee).trial_predsT = cellfun(@(X, Y) (X + Y) ./ 2, ...
        roi1(ee).trial_predsT, roi2(ee).trial_predsT, 'uni', false);
    roi(ee).trial_predsP = cellfun(@(X, Y) (X + Y) ./ 2, ...
        roi1(ee).trial_predsP, roi2(ee).trial_predsP, 'uni', false);
    switch roi(ee).model.num_channels
        case 3
            roi(ee).trial_preds = cellfun(@(X, Y, Z) X + Y + Z, ...
                roi(ee).trial_predsS, roi(ee).trial_predsT, roi(ee).trial_predsP, 'uni', false);
        case 2
            roi(ee).trial_preds = cellfun(@(X, Y) X + Y, ...
                roi(ee).trial_predsS, roi(ee).trial_predsT, 'uni', false);
        case 1
            roi(ee).trial_preds = cellfun(@(X, Y) (X + Y) ./ 2, ...
                roi1(ee).trial_preds, roi2(ee).trial_preds, 'uni', false);
    end
    roi(ee).trials = cellfun(@(X, Y) [X Y], ...
        roi1(ee).trials, roi2(ee).trials, 'uni', false);
    roi(ee).runs = vertcat(roi1(1).runs, roi2(1).runs);
    roi1(ee) = tch_noise_ceil(roi1(ee)); roi2(ee) = tch_noise_ceil(roi2(ee));
    roi(ee).noise_ceils = cellfun(@(X, Y) (X + Y) ./ 2, ...
        roi1(ee).noise_ceils, roi2(ee).noise_ceils, 'uni', false);
    % average model properties across folds
    roi(ee).model.fit_exps = [roi1(ee).model.fit_exps; roi2(ee).model.fit_exps];
    roi(ee).model.betas = cellfun(@(X, Y) (X + Y) ./ 2, ...
        roi1(ee).model.betas, roi2(ee).model.betas, 'uni', false);
    roi(ee).model.stdevs = cellfun(@(X, Y) (X + Y) ./ 2, ...
        roi1(ee).model.stdevs, roi2(ee).model.stdevs, 'uni', false);
    roi(ee).model.varexp = cellfun(@(X, Y) (X + Y) ./ 2, ...
        roi1(ee).model.varexp, roi2(ee).model.varexp, 'uni', false);
    param_names = fieldnames(roi1(ee).model.params);
    if strcmp(roi1(ee).model.type, '1ch-balloon')
        param_names = [];
    end
    for pp = 1:length(param_names)
        pn = param_names{pp};
        roi(ee).model.params.(pn) = cellfun(@(X, Y) (X + Y) ./ 2, ...
            roi1(ee).model.params.(pn), roi2(ee).model.params.(pn), 'uni', false);
    end
end
