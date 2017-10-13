% Stores and operates on fMRI time series of an ROI across multiple scan
% sessions. Used with tchModel object to fit and validate various temporal
% encoding models in each participant with a given ROI.
%
% CONSTRUCTOR INPUTS
%   1) name: name of ROI to model (e.g., 'V2')
%   2) exps: list of experiments to model (e.g., {'Exp1' 'Exp2'})
%   3) isessions: list of sessions to analyze (optional)
%
% METHODS
%   tch_runs -- preprocesses and stores time series of all voxels in ROI
%   tch_trials -- compiles trial-level time series sorted by experiment
%   tch_fit -- fits tchModel object to the mean time series of each ROI
%   tch_pred -- predicts trial responses using model solutoin
%   plot_runs -- plots measured vs. predicted responses for each run
%   plot_exps -- plots comparison of trial responses across experiments
%   plot_roi -- plots various components of ROI analysis
%   recompute -- validates model solution on indpendent data
%
% Example model fitting steps ("model" is a tchModel object):
%   roi = tchROI('V1', {'Exp1' 'Exp2'})
%   roi = tch_runs(roi)
%   roi = tch_trials(roi, model)
%   roi = tch_fit(roi, model)
%   roi = tch_pred(roi, model)
%   fig = plot_roi(roi, 'model')
%
% AS 2/2017

classdef tchROI
    
    properties
        name        % name of data directories for this region
        experiments % array of experiments to model
        sessions    % array of sessions that have data for this region
        model = []; % data structure of models fits for each session
        predS = {}; % predicted sustained contributions per trial type
        predT = {}; % predicted transient contributions per trial type
        predD = {}; % predicted delay channel contributions per trial type
        pred = {};  % total predicted contributions for each trial type
    end
    
    properties (Hidden)
        runs = {};        % responses for each run (TR x voxel)
        trials = {};      % responses for each trial type (TR x voxel)
        baseline = {};    % mean baseline response across all trial types
        noise_ceils = {}; % estiamte of noice ceiling for each experiment
        isessions = {};   % user-specified session list (optional)
    end
    
    properties (Constant, Hidden)
        project_dir = fileparts(fileparts(which(mfilename, 'class')));
        tr = 1; % fMRI TR (s)
    end
    
    properties (Dependent)
        run_avgs   % average timecourse across voxels for each run
        trial_avgs % average timecourse across voxels for each trial type
    end
    
    properties (Dependent, Hidden)
        all_sessions % paths to all session directories
        num_runs     % number of runs per experiment
        filenames    % paths to data files from sessions
        nickname     % ROI nickname
        session_ids  % session nicknames
        pred_sum     % sum of preds across all channels
        predS_sum    % sum of preds across S channels
        predT_sum    % sum of preds across T channels
        predD_sum    % sum of preds across D channels
    end
    
    
    methods
        
        % class constructor
        function roi = tchROI(name, exps, isessions)
            if nargin == 2
                roi.name = name;
                roi.experiments = force_cell(exps);
            elseif nargin == 3
                roi.name = name;
                roi.experiments = force_cell(exps);
                roi.isessions = force_cell(isessions);
            else
                error('Unexpected input arguments');
            end
        end
        
        % find all sessions in data directory if not specified by user
        function all_sessions = get.all_sessions(roi)
            if isempty(roi.isessions)
                all_sessions = find_sessions(roi.project_dir);
            else
                all_sessions = roi.isessions;
            end
        end
        
        % find the number of runs per experiment for each session
        function num_runs = get.num_runs(roi)
            sessions = roi.sessions; nsess = length(sessions);
            num_runs = zeros(length(roi.experiments), nsess);
            for ss = 1:nsess
                spath = fullfile(sessions{ss}, 'ROIs', roi.name);
                % find paths to data files for each experiment
                for ee = 1:length(roi.experiments)
                    d = dir(fullfile(spath, roi.experiments{ee}, 'Run*.mat'));
                    num_runs(ee, ss) = length({d.name});
                end
            end
        end
        
        % find the paths to the data files for each session
        function filenames = get.filenames(roi)
            nruns = roi.num_runs; filenames = {};
            % for each session
            for ss = 1:length(roi.sessions)
                rcnt = 0;
                spath = fullfile(roi.sessions{ss}, 'ROIs', roi.name);
                % for each experiment
                for ee = 1:length(roi.experiments)
                    % store paths to data file for each run
                    for rr = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        edir = fullfile(spath, roi.experiments{ee});
                        fname = ['Run' num2str(rr) '.mat'];
                        filenames{rcnt, ss} = fullfile(edir, fname);
                    end
                end
            end
        end
        
        % average run time series across all voxels in each ROI
        function run_avgs = get.run_avgs(roi)
            run_avgs = cellfun(@(x) mean(x, 2), roi.runs, 'uni', false);
            empty_mats = cellfun(@isempty, run_avgs);
            run_avgs(empty_mats) = {[]};
        end
        
        % average trial time series across all voxels in each ROI
        function trial_avgs = get.trial_avgs(roi)
            trial_avgs = cellfun(@(x) mean(x, 2), roi.trials, 'uni', false);
        end
        
        % replace problematic characters in ROI name
        function nickname = get.nickname(roi)
            nickname = roi.name;
            nickname = strrep(nickname, '_', '-');
            nickname = strrep(nickname, ' ', '-');
            nickname = strrep(nickname, '.', '-');
        end
        
        % label each session with an ID string
        function session_ids = get.session_ids(roi)
            sessions = roi.sessions;
            session_ids = cell(1, length(sessions));
            for ss = 1:length(sessions)
                [~, session_id] = fileparts(sessions{ss});
                session_ids{ss} = session_id;
            end
        end
        
        % sum trial predictors across all channels
        function pred_sum = get.pred_sum(roi)
            pred_sum = cellfun(@(X) sum(X, 2), roi.pred, 'uni', false);
        end
        
        % sum trial predictors across sustained channels
        function predS_sum = get.predS_sum(roi)
            predS_sum = cellfun(@(X) sum(X, 2), roi.predS, 'uni', false);
        end
        
        % sum trial predictors across transient channels
        function predT_sum = get.predT_sum(roi)
            predT_sum = cellfun(@(X) sum(X, 2), roi.predT, 'uni', false);
        end
        
        % sum trial predictors across all delay channels
        function predD_sum = get.predD_sum(roi)
            predD_sum = cellfun(@(X) sum(X, 2), roi.predD, 'uni', false);
        end
        
        % find set of all_sessions with current ROI
        function roi = select_sessions(roi)
            sessions = {}; all_sessions = roi.all_sessions; scnt = 0;
            for ss = 1:length(all_sessions)
                [~, session_id] = fileparts(all_sessions{ss});
                spath = fullfile(roi.project_dir, 'data', session_id);
                cpath = fullfile(spath, 'ROIs', roi.name);
                ecnt = 0;
                for ee = 1:length(roi.experiments)
                    if exist(fullfile(cpath, roi.experiments{ee}), 'dir') == 7
                        ecnt = ecnt + 1;
                    end
                end
                if ecnt == length(roi.experiments)
                    scnt = scnt + 1;
                    sessions{scnt} = spath;
                end
            end
            % error if no sessions with ROI are found
            if scnt == 0
                error(['No sessions found with ' roi.name '.']);
            else
                roi.sessions = sessions;
            end
        end
        
        % preprocess and store run timeseries of each voxel
        function roi = tch_runs(roi, detrend)
            if nargin == 1
                detrend = 4;
            end
            roi = select_sessions(roi); % select sessions with region
            fpaths = roi.filenames;     % find paths to data files
            raw_runs = cellfun(@(X) tch_load(X, 'tSeries'), fpaths, 'uni', false);
            roi.runs = cellfun(@(X) tch_psc(X, detrend), raw_runs, 'uni', false);
        end
        
        % check dimensionality of roi time series and model predictions
        function check_model(roi, model)
            empty_cells = cellfun(@isempty, model.run_durs);
            rds = model.run_durs; rds(empty_cells) = {0};
            comp = cellfun(@length, roi.run_avgs) ~= cell2mat(rds);
            if sum(comp(:)) > 0
                error('dimensions of data and model do not match');
            end
        end
        
        % compile time series for each trial type
        function roi = tch_trials(roi, model)
            % check model and get design parameters
            check_model(roi, model);
            sessions = roi.sessions; nsess = length(sessions);
            nconds = max(cellfun(@length, model.cond_list));
            nexps = length(roi.experiments);
            nruns_max = size(model.run_preds, 1); nruns = model.num_runs;
            % compile all trial time series for each session and experiment
            onsets = model.tonsets; offsets = model.toffsets;
            trials = cell(nconds, nsess, nexps); conds = model.tconds;
            for ss = 1:nsess
                rcnt = 0;
                for ee = 1:nexps
                    cnts = zeros(1, length(model.cond_list{ee}));
                    for rr = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        % estimate prestimulus baseline response for run
                        oframes = repmat(onsets{rcnt, ss}, model.pre_dur, 1);
                        b = repmat([1:model.pre_dur]', 1, length(oframes));
                        psf = oframes - b; bf = psf(:);
                        % calculate mean baseline response to subtract
                        bs = mean(roi.run_avgs{rcnt, ss}(bf));
                        roi.baseline{rcnt, ss} = bs;
                        % store peri-stimulus time series sorted by trial
                        for tt = 1:length(onsets{rcnt, ss})
                            % get condition number of this trial
                            cond = conds{rcnt, ss}(tt);
                            cond_idx = find(strcmp(cond, model.cond_list{ee}));
                            cnts(cond_idx) = cnts(cond_idx) + 1;
                            % get TR corresponding to onset of pre_dur
                            on_idxs = onsets{rcnt, ss}(tt) - model.pre_dur;
                            onset = on_idxs / roi.tr + 1;
                            off_idxs = offsets{rcnt, ss}(tt) - model.gap_dur / 2;
                            offset = (floor(off_idxs) + model.post_dur) / roi.tr + 1;
                            % extract the peri-stimulus time window
                            trial_avg = roi.run_avgs{rcnt, ss}(onset:offset) - bs;
                            trials{cond_idx, ss, ee}(:, cnts(cond_idx)) = trial_avg;
                        end
                    end
                end
            end
            roi.trials = trials;
        end
        
        % estimate noise ceiling for each ROI using inter-trial variability
        function roi = tch_noise_ceil(roi)
            trials = roi.trials;
            trials_avg = cellfun(@(X) mean(X, 2), trials, 'uni', false);
            trials_err = cellfun(@(X, Y) X - repmat(Y, 1, size(X, 2)), trials, trials_avg, 'uni', false);
            trials_err = cellfun(@(X) sum(sum(X .^ 2)), trials_err, 'uni', false);
            total_err = []; total_var = []; ceils = [];
            for ss = 1:size(trials, 2)
                total_err(ss) = sum([trials_err{:, ss}]);
                trial_mean = mean(mean(vertcat(trials{:, ss})));
                trials_var = vertcat(trials{:, ss}) - trial_mean;
                total_var(ss) = sum(sum(trials_var .^ 2));
                roi.noise_ceils{ss} = 1 - (total_err / total_var);
            end
        end
        
        % use GLM to fit weights for each predictor in model
        function [roi, model] = tch_fit(roi, model, optimize_flag, fit_exps)
            if nargin < 3; optimize_flag = 0; end
            if nargin < 4; fit_exps = model.experiments; end
            check_model(roi, model); sessions = roi.sessions;
            nruns_max = size(model.run_preds, 1); % max number of runs
            nruns = nruns_max - sum(cellfun(@isempty, model.run_preds));
            % concatenate data and preds across all runs in each session
            for ss = 1:length(sessions)
                run_preds = vertcat(model.run_preds{:, ss}); npreds = size(run_preds, 2);
                run_durs = model.run_durs(:, ss); num_runs = sum(cell2mat(run_durs) > 0);
                b0_cell = cellfun(@(X) zeros(X, num_runs), run_durs, 'uni', false);
                for rr = 1:num_runs; b0_cell{rr}(:, rr) = 1; end
                b0 = cell2mat(b0_cell); predictors = [run_preds b0];
                run_avgs = roi.run_avgs(:, ss); baseline = roi.baseline(:, ss);
                tc_cell = cellfun(@(X, Y) X - Y, run_avgs, baseline, 'uni', false);
                tc = vertcat(tc_cell{:}); roi.model.run_tcs{ss} = tc;
                % fit GLM and store betas, SEMs, and variance explained
                mm = tch_glm(tc, predictors);
                roi.model.run_preds{ss} = predictors * mm.betas';
                roi.model.betas{ss} = mm.betas(1:npreds);
                roi.model.stdevs{ss} = mm.stdevs(1:npreds);
                roi.model.residual{ss} = mm.residual;
                rbetas = mm.betas(npreds + 1:npreds + nruns(ss));
                % store paramters of nuisance regressors
                roi.model.rbetas{ss} = rbetas;
                rstdevs = mm.stdevs(npreds + 1:npreds + nruns(ss));
                roi.model.rstdevs{ss} = rstdevs;
                ve = 1 - (sum(mm.residual .^ 2) / sum((tc - mean(tc)) .^ 2));
                roi.model.varexp{ss} = ve;
            end
            % optimize model parameters if applicable
            omodels = {'1ch-pow' '1ch-div' '1ch-dcts' '1ch-exp' '1ch-cexp' ...
                '2ch-pow-quad' '2ch-pow-rect' '2ch-div-quad' ...
                '2ch-exp-quad' '2ch-exp-rect' ...
                '2ch-cexp-quad' '2ch-cexp-rect' ...
                '3ch-lin-quad-exp' '3ch-lin-rect-exp' ...
                '3ch-lin-quad-exp' '3ch-lin-rect-exp' ...
                '3ch-pow-quad-exp' '3ch-pow-rect-exp' ...
                '3ch-exp-quad-exp' '3ch-exp-rect-exp' ...
                '3ch-cexp-quad-exp' '3ch-cexp-rect-exp'};
            if optimize_flag && sum(strcmp(model.type, omodels))
                param_names = fieldnames(model.params);
                for ss = 1:length(sessions)
                    fname_grid = ['grid_search_results_' model.type ...
                        '_fit' [fit_exps{:}] '.mat'];
                    fpath_grid = fullfile(sessions{ss}, 'ROIs', ...
                        roi.name, fname_grid);
                    fname_grad = ['grad_desc_results_' model.type ...
                        '_fit' [fit_exps{:}] '.mat'];
                    fpath_grad = fullfile(sessions{ss}, 'ROIs', ...
                        roi.name, fname_grad);
                    % load optimization results if saved, otherwise compute
                    if exist(fpath_grad, 'file') == 2
                        fprintf('Loading gradient descent results. \n');
                        load(fpath_grad);
                        oroi = tchROI(roi.name, roi.experiments, sessions{ss});
                        oroi = tch_runs(oroi);
                        omodel = tchModel(model.type, roi.experiments, sessions{ss});
                        omodel.normT = model.normT; omodel.normD = model.normD;
                        for pp = 1:length(param_names)
                            pn = param_names{pp};
                            omodel.params.(pn){1} = params.(pn){1};
                        end
                        omodel = code_stim(omodel);
                        omodel = update_param(omodel, pn, 0);
                        omodel = pred_runs(omodel);
                        omodel = pred_trials(omodel);
                        oroi = tch_trials(oroi, omodel);
                        oroi = tch_fit(oroi, omodel);
                    elseif exist(fpath_grid, 'file') == 2
                        fprintf('Loading grid search results. \n');
                        load(fpath_grid);
                        oroi = tchROI(roi.name, roi.experiments, sessions{ss});
                        oroi = tch_runs(oroi);
                        omodel = tchModel(model.type, roi.experiments, sessions{ss});
                        omodel.normT = model.normT; omodel.normD = model.normD;
                        omodel = code_stim(omodel);
                        for mm = 1:length(params)
                            omodel(mm) = omodel(1); oroi(mm) = oroi(1);
                            for pp = 1:length(param_names)
                                pn = param_names{pp};
                                omodel(mm).params.(pn){1} = params(mm).(pn){1};
                                omodel(mm) = update_param(omodel(mm), pn, 0);                            
                            end
                            omodel(mm) = pred_runs(omodel(mm));
                            omodel(mm) = pred_trials(omodel(mm));
                            oroi(mm) = tch_trials(oroi(mm), omodel(mm));
                            oroi(mm) = tch_fit(oroi(mm), omodel(mm));
                        end
                        [oroi, omodel] = tch_optimize_fit(oroi, omodel);
                        params = omodel(1).params;
                        save(fpath_grad, 'params', '-v7.3');
                    else
                        [oroi, omodel] = tch_grid_search(roi, model, ss, 5);
                        params = omodel(1).params;
                        for mm = 2:length(omodel)
                            params(mm) = omodel(mm).params;
                        end
                        save(fpath_grid, 'params', '-v7.3');
                        [oroi, omodel] = tch_optimize_fit(oroi, omodel);
                        params = omodel(1).params;
                        save(fpath_grad, 'params', '-v7.3');
                    end
                    % copy optimized parameters for session
                    for pp = 1:length(param_names)
                        opt_params = omodel.params.(param_names{pp}){1};
                        model.params.(param_names{pp}){ss} = opt_params;
                        model = update_param(model, param_names{pp}, 0);
                    end
                end
                model = pred_runs(model);
                model = pred_trials(model);
                [roi, model] = tch_fit(roi, model, 0);
            end
            % carry over model parameters for all sessions to roi.model
            roi.model.type = model.type;
            roi.model.num_channels = model.num_channels;
            roi.model.cond_list = model.cond_list;
            roi.model.cat_list = unique([model.cats{:}]);
            roi.model.pre_dur = model.pre_dur;
            roi.model.post_dur = model.post_dur;
            roi.model.fit_exps = model.experiments;
            roi.model.params = model.params;
            roi.model.irfs = model.irfs;
        end
        
        % predict responses for each trial type using model solution
        function roi = tch_pred(roi, model)
            nexps = length(roi.experiments);
            nconds = max(cellfun(@length, model.cond_list));
            ncats = length(unique([model.cats{:}]));
            nsubs = length(roi.sessions);
            % preallocate predictor array for each trial type
            roi.pred = cell(nconds, nsubs, nexps);
            % preallocate S and T predictors if using multi-channel model
            if model.num_channels > 1
                roi.predS = cell(nconds, nsubs, nexps);
                roi.predT = cell(nconds, nsubs, nexps);
            end
            if model.num_channels > 2
                roi.predD = cell(nconds, nsubs, nexps);
            end
            % predict response for each session, experiment, and trial type
            for ss = 1:nsubs
                for ee = 1:nexps
                    for cc = 1:length(model.cond_list{ee})
                        % if using a multi-channel model
                        if model.num_channels == 1
                            amp = roi.model.betas{ss}(1:ncats);
                            % scale trial predictors by betas
                            pred = model.trial_preds.pred{cc, ss, ee};
                            pred = pred .* repmat(amp, size(pred, 1), 1);
                            % store trial predictors in roi
                            roi.pred{cc, ss, ee} = pred;
                        end
                        if model.num_channels > 1
                            ampS = roi.model.betas{ss}(1:ncats);
                            ampT = roi.model.betas{ss}(ncats + 1:2 * ncats);
                            % scale trial predictors by betas
                            predS = model.trial_preds.S{cc, ss, ee};
                            predT = model.trial_preds.T{cc, ss, ee};
                            fmriS = predS .* repmat(ampS, size(predS, 1), 1);
                            fmriT = predT .* repmat(ampT, size(predT, 1), 1);
                            % store trial predictors in roi
                            roi.predS{cc, ss, ee} = fmriS;
                            roi.predT{cc, ss, ee} = fmriT;
                            if model.num_channels == 2
                                roi.pred{cc, ss, ee} = fmriS + fmriT;
                            end
                            if model.num_channels > 2
                                ampD = roi.model.betas{ss}(2 * ncats + 1:3 * ncats);
                                predD = model.trial_preds.D{cc, ss, ee};
                                fmriD = predD .* repmat(ampD, size(predD, 1), 1);
                                roi.predD{cc, ss, ee} = fmriD;
                                roi.pred{cc, ss, ee} = fmriS + fmriT + fmriD;
                            end
                        end
                    end
                end
            end
        end
        
        % plot various components of ROI analysis
        function fig = plot_roi(roi, plot_type, save_flag)
            if nargin < 2; plot_type = 'model'; end
            if nargin < 3; save_flag = 0; end
            switch plot_type
                case 'model'
                    fig = plot_roi_model(roi, save_flag);
                case 'runs'
                    fig = plot_roi_runs(roi, save_flag);
                case 'exps'
                    fig = plot_roi_exps(roi, save_flag);
                case 'trials'
                    fig = plot_roi_trials(roi, save_flag);
            end
        end
        
        % recompute model performance given specified weights
        function roi = recompute(roi, model, fit)
            check_model(roi, model);
            [nruns, nsubs] = size(model.run_preds);
            for ss = 1:nsubs
                run_preds = vertcat(model.run_preds{:, ss}); npreds = size(run_preds, 2);
                run_durs = model.run_durs(:, ss); num_runs = sum(cell2mat(run_durs) > 0);
                b0_cell = cellfun(@(X) zeros(X, num_runs), run_durs, 'uni', false);
                for rr = 1:num_runs
                    b0_cell{rr}(:, rr) = 1;
                end
                b0 = cell2mat(b0_cell); predictors = [run_preds b0];
                run_avgs = roi.run_avgs(:, ss); baseline = roi.baseline(:, ss);
                tc_cell = cellfun(@(X, Y) X - Y, run_avgs, baseline, 'uni', false);
                tc = vertcat(tc_cell{:}); roi.model.run_tcs{ss} = tc;
                beta_vec = zeros(1, npreds + num_runs);
                beta_vec(1:npreds) = fit.betas{ss};
                beta_vec(npreds + 1:npreds + num_runs) = roi.model.rbetas{ss};
                % predict predict fMRI responses for each run
                run_pred = predictors * beta_vec';
                roi.model.run_preds{ss} = run_pred;
                res = tc - run_pred;
                % calculate model performance
                roi.model.varexp{ss} = 1 - (sum(res.^2) / sum((tc - mean(tc)).^2));
            end
            % store new fit in roi structure
            roi.model.betas = fit.betas;
            roi.model.stdevs = fit.stdevs;
            roi.model.fit_exps = fit.fit_exps;
            roi = tch_pred(roi, model);
        end
        
    end
    
end
