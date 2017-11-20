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
%   plot_roi -- plots various components of ROI analysis
%   tch_recompute -- validates model solution on indpendent data
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
        tr = 1;           % fMRI TR (s)
        project_dir = fileparts(fileparts(which(mfilename, 'class')));
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
        pred_sum     % sum of predictors across all channels
        predS_sum    % sum of predictors across S channels
        predT_sum    % sum of predictors across T channels
        predD_sum    % sum of predictors across D channels
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
            num_runs = zeros(length(roi.experiments), length(roi.sessions));
            for ss = 1:length(roi.sessions)
                spath = fullfile(roi.sessions{ss}, 'ROIs', roi.name);
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
                spath = fullfile(roi.sessions{ss}, 'ROIs', roi.name); rcnt = 0;
                % for each experiment
                for ee = 1:length(roi.experiments)
                    % store paths to data file for each run
                    for rr = 1:nruns(ee, ss)
                        edir = fullfile(spath, roi.experiments{ee});
                        fname = ['Run' num2str(rr) '.mat']; rcnt = rcnt + 1;
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
            session_ids = cell(1, length(roi.sessions));
            for ss = 1:length(roi.sessions)
                [~, session_ids{ss}] = fileparts(roi.sessions{ss});
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
            sessions = {}; scnt = 0;
            for ss = 1:length(roi.all_sessions)
                [~, session_id] = fileparts(roi.all_sessions{ss});
                spath = fullfile(roi.project_dir, 'data', session_id);
                cpath = fullfile(spath, 'ROIs', roi.name); ecnt = 0;
                for ee = 1:length(roi.experiments)
                    if exist(fullfile(cpath, roi.experiments{ee}), 'dir') == 7
                        ecnt = ecnt + 1;
                    end
                end
                if ecnt == length(roi.experiments)
                    scnt = scnt + 1; sessions{scnt} = spath;
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
        function roi = tch_runs(roi, detrend_option)
            if nargin == 1
                % use option 3 for regions with high SNR
                if sum(strcmp(roi.name, {'V1' 'V2' 'V3'}))
                    detrend_option = 3;
                else
                    detrend_option = 4;
                end
            end
            % select sessions with region and find paths to data files
            roi = select_sessions(roi); fpaths = roi.filenames;
            raw_runs = cellfun(@(X) tch_load(X, 'tSeries'), fpaths, 'uni', false);
            roi.runs = cellfun(@(X) tch_psc(X, detrend_option), raw_runs, 'uni', false);
        end
        
        % check dimensionality of roi time series and model predictions
        function check_model(roi, model)
            empty_cells = cellfun(@isempty, model.run_durs);
            rds = model.run_durs; rds(empty_cells) = {0};
            rds = cellfun(@(X) X / roi.tr, rds, 'uni', false);
            comp = cellfun(@length, roi.run_avgs) ~= cell2mat(rds);
            if sum(comp(:)) > 0
                error('Dimensions of data and model do not match');
            end
        end
        
        % compile time series for each trial type
        function roi = tch_trials(roi, model)
            % check model and get design parameters
            check_model(roi, model);
            nconds = max(cellfun(@length, model.cond_list));
            nexps = length(roi.experiments); nruns = model.num_runs;
            % compile all trial time series for each session and experiment
            onsets = model.tonsets; offsets = model.toffsets;
            trials = cell(nconds, length(roi.sessions), nexps);
            for ss = 1:length(roi.sessions)
                rcnt = 0;
                for ee = 1:nexps
                    cnts = zeros(1, length(model.cond_list{ee}));
                    for rr = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        % estimate prestimulus baseline response for run
                        npsfs = round(model.pre_dur / model.tr);
                        oframes = repmat(onsets{rcnt, ss}, npsfs, 1) / model.tr;
                        psf = oframes - repmat((1:npsfs)', 1, length(oframes));
                        % calculate mean baseline response to subtract
                        bs = mean(roi.run_avgs{rcnt, ss}(psf(:)), 1);
                        roi.baseline{rcnt, ss} = bs;
                        % store peri-stimulus time series sorted by trial
                        for tt = 1:length(onsets{rcnt, ss})
                            % get condition number of this trial
                            cond = model.tconds{rcnt, ss}(tt);
                            cond_idx = find(strcmp(cond, model.cond_list{ee}));
                            cnts(cond_idx) = cnts(cond_idx) + 1;
                            % get TR corresponding to onset of pre_dur
                            on_idx = (onsets{rcnt, ss}(tt) - model.pre_dur) / roi.tr + 1;
                            off_time = offsets{rcnt, ss}(tt) - model.gap_dur / 2;
                            off_idx = (floor(off_time) + model.post_dur + 1) / roi.tr;
                            % extract the peri-stimulus time window
                            trial_avg = roi.run_avgs{rcnt, ss}(on_idx:off_idx) - bs;
                            trials{cond_idx, ss, ee}(:, cnts(cond_idx)) = trial_avg;
                        end
                    end
                end
            end
            roi.trials = trials;
        end
        
        % estimate noise ceiling for each ROI using inter-trial variability
        function roi = tch_noise_ceil(roi)
            % find average response in for each trial type per session
            trials_avg = cellfun(@(X) mean(X, 2), roi.trials, 'uni', false);
            % calculate residual between individual trials and means
            trials_err = cellfun(@(X, Y) X - repmat(Y, 1, size(X, 2)), ...
                roi.trials, trials_avg, 'uni', false);
            trials_err = cellfun(@(X) sum(sum(X .^ 2)), trials_err, 'uni', false);
            % calculate varianced explained by mean trial responses
            for ss = 1:size(roi.trials, 2)
                total_err = sum(cell2mat(trials_err(:, ss)));
                trial_mean = mean(mean(cell2mat(roi.trials(:, ss))));
                trials_var = cell2mat(roi.trials(:, ss)) - trial_mean;
                total_var = sum(sum(trials_var .^ 2));
                roi.noise_ceils{ss} = 1 - (total_err / total_var);
            end
        end
        
        % use GLM to fit channel weights for each predictor in model
        function [roi, model] = tch_fit(roi, model, optim_proc, fit_exps)
            if nargin < 3; optim_proc = 0; end
            if nargin < 4; fit_exps = model.experiments; end
            check_model(roi, model); npreds = size(model.run_preds{1, 1}, 2);
            % subtract baseline estimates from centered time series
            tcs = cellfun(@(X, Y) X - Y, roi.run_avgs, roi.baseline, 'uni', false);
            rfs = cellfun(@(X) X / model.tr, model.run_durs, 'uni', false);
            % concatenate data and preds across all runs in each session
            for ss = 1:length(roi.sessions)
                nruns = sum(model.num_runs(:, ss));
                % construct nuisance regressors and merge with predictors
                b0 = cell2mat(cellfun(@(X, Y) code_stim_vec(zeros(X, nruns), 1:X, Y), ...
                    rfs(:, ss), num2cell(1:size(tcs, 1))', 'uni', false));
                predictors = [cell2mat(model.run_preds(:, ss)) b0];
                % fit GLM and store betas, SEMs, and variance explained
                tc = cell2mat(tcs(:, ss)); roi.model.run_tcs{ss} = tc;
                mm = tch_glm(tc, predictors);
                roi.model.run_preds{ss} = predictors * mm.betas';
                roi.model.betas{ss} = mm.betas(1:npreds);
                roi.model.stdevs{ss} = mm.stdevs(1:npreds);
                roi.model.residual{ss} = mm.residual;
                res_var = sum(mm.residual .^ 2) ./ sum((tc - mean(tc)) .^ 2);
                roi.model.varexp{ss} = 1 - res_var;
                % store paramters of nuisance regressors
                roi.model.rbetas{ss} = mm.betas(npreds + 1:npreds + nruns);
                roi.model.rstdevs{ss} = mm.stdevs(npreds + 1:npreds + nruns);
            end
            % optimize model parameters if applicable
            if model.optimize_flag == 1 && optim_proc == 1
                [roi, model] = tch_optimize_fmincon(roi, model, fit_exps);
            elseif model.optimize_flag == 1 && optim_proc == 2
                [roi, model] = tch_optimize_grid_grad(roi, model, fit_exps);
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
                        % if using a single-channel model
                        if model.num_channels == 1
                            amp = roi.model.betas{ss}(1:ncats);
                            % scale trial predictors by betas
                            pred = model.trial_preds.pred{cc, ss, ee};
                            pred = pred .* repmat(amp, size(pred, 1), 1);
                            % store trial predictors in roi
                            roi.pred{cc, ss, ee} = pred;
                        end
                        % if using a mutli-channel model
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
        function roi = tch_recompute(roi, model, fit)
            check_model(roi, model);
            % subtract baseline estimates from centered time series
            tcs = cellfun(@(X, Y) X - Y, roi.run_avgs, roi.baseline, 'uni', false);
            rfs = cellfun(@(X) X / model.tr, model.run_durs, 'uni', false);
            for ss = 1:length(roi.sessions)
                nruns = sum(model.num_runs(:, ss));
                % construct nuisance regressors and merge with predictors
                b0 = cell2mat(cellfun(@(X, Y) code_stim_vec(zeros(X, nruns), 1:X, Y), ...
                    rfs(:, ss), num2cell(1:size(tcs, 1))', 'uni', false));
                predictors = [cell2mat(model.run_preds(:, ss)) b0];
                % predict fMRI responses for each run
                tc = cell2mat(tcs(:, ss));
                beta_vec = [fit.betas{ss} roi.model.rbetas{ss}];
                run_pred = predictors * beta_vec';
                % calculate cross-validated model performance
                roi.model.run_preds{ss} = run_pred; res = tc - run_pred;
                res_var = sum(res .^ 2) / sum((tc - mean(tc)) .^ 2);
                roi.model.varexp{ss} = 1 - res_var;
            end
            % store new fit in roi structure
            roi.model.betas = fit.betas; roi.model.stdevs = fit.stdevs;
            roi.model.fit_exps = fit.fit_exps; roi = tch_pred(roi, model);
        end
        
    end
    
end
