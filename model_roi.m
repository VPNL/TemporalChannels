function [roi, model] = model_roi(name, type, exps_fit, exps_val)
% Wrapper function that fits a temporal model object (ModelTS) to a region
% time series object (ROI) and plots the fit and predictions of the 
% model. To validate the solution, include the optional fourth input
% to predict data from 'exps_val' using a model fit to 'exps_fit'.
% 
% INPUTS
%   1) name: directory name of ROI to model (e.g., 'V1')
%   2) type: which type of model to use
%   3) exps_fit: list of experiments for fitting the model (e.g., {'Exp1' 'Exp2'})
%   4) exps_val: list of experiments for validating the fit (optional)
% 
% OUTPUTS
%   1) roi: fitted ROI object containing measured and predicted responses
%   2) model: ModelTS object used to fit and predict ROI responses
% 
% EXAMPLES
% 
% Fit the 'standard' model to multiple experiments in V1:
% [roi, model] = model_roi('V1', 'standard', {'Exp1' 'Exp2' 'Exp3'});
% 
% Cross-validate the fit of a 2channel model using independent data:
% [roi, model] = model_roi('V1', '2ch', {'Exp1' 'Exp2'}, {['Exp3'});
% 
% AS 2/2017


%% Setup paths and check inputs

% add paths to class objects and helper functions
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


%% Fit the model to fit_exps

% setup ROI object for fitting ModelTS
roi(1) = ROI(name, exps_fit);
fprintf('\nExtracting run time series...\n')
roi(1) = tc_runs(roi(1));

% setup ModelTS object to applyt to ROI
model(1) = ModelTS(type, exps_fit, roi(1).sessions);
fprintf('Coding the stimulus...\n')
model(1) = code_stim(model(1));
fprintf('Generating predictors...\n')
model(1) = run_preds(model(1));
model(1) = trial_preds(model(1));

% fit ModelTS to ROI
fprintf('Extracting trial time series...\n')
roi(1) = tc_trials(roi(1), model(1));
fprintf('Fitting the model...\n')
[roi(1), model(1)] = tc_fit(roi(1), model(1), 1);
roi(1) = tc_pred(roi(1), model(1));


%% validation the model on test_exps if applicable

if cv_flag
    fprintf('Performing validation...\n')
    % setup ROI and ModelTS objects for validation
    roi(2) = ROI(name, exps_val);
    roi(2) = tc_runs(roi(2));
    model(2) = ModelTS(type, exps_val, roi(2).sessions);
    model(2) = code_stim(model(2));
    model(2) = run_preds(model(2));
    model(2) = trial_preds(model(2));
    % setup model struct by fitting model directly to data
    roi(2) = tc_trials(roi(2), model(2));
    [roi(2), model(2)] = tc_fit(roi(2), model(2));
    roi(2) = tc_pred(roi(2), model(2));
    % use model fit to data from exps_fit to predict data in exps_val
    roi(2) = recompute(roi(2), model(2), roi(1).model);
end


%% Plot results

% plot model fit and predictions for fit_exps
plot_model(roi(1));

% plot cross-validated predictions for test_exps if applicable
if cv_flag; plot_model(roi(2)); end


%% Save results
fname = [roi(1).nickname '_' roi(1).model.type '_fit' [roi(1).experiments{:}]];
if cv_flag
    fname = [fname '_val' [roi(2).experiments{:}]];
end
fname = [fname '.mat'];
fpath = fullfile(roi(1).project_dir, 'results', fname);
save(fpath, 'roi', 'model', '-v7.3');

end
