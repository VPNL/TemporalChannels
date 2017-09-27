function [roi, model] = model_roi(name, type, fit_exps, val_exps)
% Wrapper function that fits a temporal model object (ModelTS) to a region
% time series object (ROI) and plots the fit and predictions of the 
% model. To validate the solution, include the optional fourth input
% to predict data from 'val_exps' using a model fit to 'fit_exps'. 
% Validation performance is computed separately for each row of experiment
% names in 'val_exps' such that {'Exp1' 'Exp2'; 'Exp3' 'Exp4'}
% 
% INPUTS
%   1) name: directory name of ROI to model (e.g., 'V1')
%   2) type: which type of model to use
%   3) fit_exps: set of experiments for fitting the model (e.g., {'Exp1' 'Exp2'})
%   4) val_exps: set/s of experiments for validating the fit (optional)
% 
% OUTPUTS
%   1) roi: fitted ROI object containing measured and predicted responses
%   2) model: ModelTS object used to fit and predict ROI responses
% 
% EXAMPLES
% 
% Fit a GLM to multiple experiments in V1:
% [roi, model] = model_roi('V1', 'glm', {'Exp1' 'Exp2' 'Exp3'});
% 
% Validate the fit of 2ch model across data from multiple experiments:
% [roi, model] = model_roi('V1', '2ch', {'Exp1' 'Exp2'}, {'Exp3' 'Exp4'});
%
% Validate the fit of cts-pow model on data separately for each experiment:
% [roi, model] = model_roi('V1', '2ch', {'Exp1' 'Exp2'}, {'Exp1'; 'Exp2'; 'Exp3'; 'Exp4'});
% 
% AS 2/2017


%% Setup paths and check inputs
mpath = fileparts(mfilename('fullpath'));
addpath(genpath(mpath));

% determine whether performing cross-validations
if nargin == 4
    cv_flag = 1;
elseif nargin == 3
    cv_flag = 0;
else
    error('Unexpected input arguements.');
end
smodels = {'glm' 'balloon' 'cts-pow' 'cts-div' 'dcts'}; % single-channel

%% Fit the model to fit_exps

% setup ROI object for fitting ModelTS
roi(1) = ROI(name, fit_exps);
fprintf('\nExtracting run time series for %s...\n', roi(1).nickname)
roi(1) = tc_runs(roi(1));

% setup ModelTS object to applyt to ROI
model(1) = ModelTS(type, fit_exps, roi(1).sessions);
fprintf('Coding the stimulus...\n')
model(1) = code_stim(model(1));
fprintf('Generating predictors...\n')
model(1) = pred_runs(model(1));
model(1) = pred_trials(model(1));

% fit ModelTS to ROI
fprintf('Extracting trial time series...\n')
roi(1) = tc_trials(roi(1), model(1));
fprintf('Fitting the %s model...\n', model(1).type)
[roi(1), model(1)] = tc_fit(roi(1), model(1), 1);
if model(1).num_channels > 1
    model(1) = norm_model(model(1), 1);
end
[roi(1), model(1)] = tc_fit(roi(1), model(1));
roi(1) = tc_pred(roi(1), model(1));


%% validate the model on test_exps if applicable

if cv_flag
    num_vals = size(val_exps, 1);
    for vv = 1:num_vals
        fprintf('Performing validation...\n')
        vn = 1 + vv;
        % setup ROI and ModelTS objects for validation
        roi(vn) = ROI(name, val_exps(vv, :), roi(1).sessions);
        roi(vn) = tc_runs(roi(vn));
        model(vn) = ModelTS(type, val_exps(vv, :), roi(vn).sessions);
        model(vn) = code_stim(model(vn));
        if model(vn).num_channels > 1 model(vn).normT = model(1).normT; end
        if model(vn).num_channels > 2 model(vn).normD = model(1).normD; end
        model(vn) = pred_runs(model(vn));
        model(vn) = pred_trials(model(vn));
        % setup model struct by fitting model directly to data
        roi(vn) = tc_trials(roi(vn), model(vn));
        [roi(vn), model(vn)] = tc_fit(roi(vn), model(vn), 1, fit_exps);
        roi(vn) = tc_pred(roi(vn), model(vn));
        % use model fit to data from fit_exps to predict data in val_exps
        roi(vn) = recompute(roi(vn), model(vn), roi(1).model);
    end
end


%% Plot model predictions

% plot model fit and predictions for fit_exps and val_exps
for rr = 1:length(roi)
    plot_model(roi(rr));
end

%% Save output in results
fname = [roi(1).nickname '_' roi(1).model.type '_fit' [roi(1).experiments{:}]];
if cv_flag
    for vv = 1:num_vals
        fname = [fname '_val' [roi(vv + 1).experiments{:}]];
    end
end
fname = [fname '.mat'];
fpath = fullfile(roi(1).project_dir, 'results', fname);
save(fpath, 'roi', 'model', '-v7.3');

end
