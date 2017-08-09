% Stores and operates on fMRI time series of an ROI across multiple scan
% sessions. Used with ModelTS object to fit and validate various temporal
% encoding models in each participant with a given ROI.
%
% CONSTRUCTOR INPUTS
%   1) name: name of ROI to model (e.g., 'V2')
%   2) exps: list of experiments to model (e.g., {'Exp1' 'Exp2'})
%
% METHODS
%   tc_runs -- preprocesses and stores time series of all voxels in ROI
%   tc_trials -- compiles trial-level time series sorted by experiment
%   tc_fit -- fits ModelTS object to the mean time series of each ROI
%   tc_pred -- predicts trial responses using model solutoin
%   plot_runs -- plots measured vs. predicted responses for each run
%   plot_exps -- plots comparison of trial responses across experiments
%   plot_model -- plots measured vs. predicted reponses for each trial type
%   recompute -- validates model solution on indpendent data
%
% Example model fitting steps ("model" is a ModelTS object):
%   roi = ROI('V1', {'Exp1' 'Exp2'})
%   roi = tc_runs(roi)
%   roi = tc_trials(roi, model)
%   roi = tc_fit(roi, model)
%   roi = tc_pred(roi, model)
%   fig = plot_model(roi)
%
% AS 2/2017

classdef ROI
    
    properties
        name        % name of data directories for this ROI
        experiments % array of experiments to model
        model = []; % data structure of models fits for each session
        predS = {}; % predicted sustained contributions per trial type
        predT = {}; % predicted transient contributions per trial type
        predD = {}; % predicted delay activity contributions per trial type
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
        sessions   % paths to sessions that have current ROI
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
        function roi = ROI(name, exps, isessions)
            if nargin == 2
                roi.name = name;
                roi.experiments = force_cell(exps);
            elseif nargin == 3
                roi.name = name;
                roi.experiments = force_cell(exps);
                roi.isessions = force_cell(isessions);
            else
                error('incorrect input arguments');
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
        
        % find set of all_sessions with current ROI
        function sessions = get.sessions(roi)
            sessions = {}; scnt = 0;
            for ss = 1:length(roi.all_sessions)
                [~, session_id] = fileparts(roi.all_sessions{ss});
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
                error(['No sessions found with: ' roi.name '.']);
            end
        end
        
        % find the number of runs per experiment for each session
        function num_runs = get.num_runs(roi)
            sessions = roi.sessions; nsess = length(sessions);
            num_runs = zeros(length(roi.experiments), length(roi.sessions));
            for ss = 1:nsess
                spath = fullfile(sessions{ss}, 'ROIs', roi.name);
                % find paths to data files for each experiment
                for ee = 1:length(roi.experiments)
                    ecnt = 1;
                    while exist(fullfile(spath, roi.experiments{ee}, ['Run' num2str(ecnt) '.mat']), 'file') == 2
                        num_runs(ee, ss) = num_runs(ee, ss) + 1;
                        ecnt = ecnt + 1;
                    end
                end
            end
        end
        
        % find the paths to the data files for each session
        function filenames = get.filenames(roi)
            sessions = roi.sessions; nsess = length(sessions);
            nruns = roi.num_runs;
            filenames = {};
            % for each session
            for ss = 1:nsess
                rcnt = 0;
                spath = fullfile(sessions{ss}, 'ROIs', roi.name);
                % for each experiment
                for ee = 1:length(roi.experiments)
                    % store paths to data file for each run
                    for rr = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        filenames{rcnt, ss} = fullfile(spath, roi.experiments{ee}, ['Run' num2str(rr) '.mat']);
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
            nickname = strrep(roi.name, '_', '-');
            nickname = strrep(roi.name, ' ', '-');
            nickname = strrep(roi.name, '.', '-');
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
        
        % preprocess and store run timeseries of each voxel
        function roi = tc_runs(roi)
            fpaths = roi.filenames; % paths to data files
            raw_runs = cellfun(@(X) loadTS(X, 'tSeries'), fpaths, 'uni', false);
            roi.runs = cellfun(@(X) psc(X), raw_runs, 'uni', false);
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
        function roi = tc_trials(roi, model)
            % check model and get design parameters
            check_model(roi, model);
            sessions = roi.sessions; nsess = length(sessions);
            nconds = max(cellfun(@length, model.cond_list));
            nexps = length(roi.experiments);
            nruns_max = size(model.run_preds, 1);
            nruns = model.num_runs;
            % get onsets, offsets, and trial orders
            onsets = model.tonsets;
            offsets = model.toffsets;
            conds = model.tconds;
            % compile all trial time series for each session and experiment
            trials = cell(nconds, nsess, nexps);
            for ss = 1:nsess
                rcnt = 0;
                for ee = 1:nexps
                    cnts = zeros(1, length(model.cond_list{ee}));
                    for rr = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        % estimate prestimulus baseline response for run
                        oframes = repmat(onsets{rcnt, ss}, model.pre_dur, 1);
                        psf = oframes - repmat([1:model.pre_dur]', 1, length(oframes));
                        bf = psf(:);
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
                            onset = (onsets{rcnt, ss}(tt) - model.pre_dur) / roi.tr + 1;
                            offset = (floor(offsets{rcnt, ss}(tt) - model.gap_dur / 2) + model.post_dur) / roi.tr + 1;
                            % extract the peri-stimulus time window
                            trials{cond_idx, ss, ee}(:, cnts(cond_idx)) = roi.run_avgs{rcnt, ss}(onset:offset) - bs;
                        end
                    end
                end
            end
            roi.trials = trials;
        end
        
        % estimate noise ceiling for each ROI using inter-trial variability
        function roi = tc_noise_ceil(roi)
            trials = roi.trials;
            trials_avg = cellfun(@(X) mean(X, 2), trials, 'uni', false);
            trials_err = cellfun(@(X, Y) X - repmat(Y, 1, size(X, 2)), trials, trials_avg, 'uni', false);
            trials_err = cellfun(@(X) sum(sum(X .^ 2)), trials_err, 'uni', false);
            total_err = []; total_var = []; ceils = [];
            for ss = 1:size(trials, 2)
                total_err(ss) = sum([trials_err{:, ss}]);
                trials_var = vertcat(trials{:, ss}) - mean(mean(vertcat(trials{:, ss})));
                total_var(ss) = sum(sum(trials_var .^ 2));
                roi.noise_ceils{ss} = 1 - (total_err / total_var);
            end
        end
        
        % use GLM to fit weights for each predictor in model
        function [roi, model] = tc_fit(roi, model, optimize_flag)
            if nargin < 3
                optimize_flag = 0;
            end
            check_model(roi, model);
            nruns_max = size(model.run_preds, 1); % max number of runs
            nruns = nruns_max - sum(cellfun(@isempty, model.run_preds));
            sessions = roi.sessions;
            % concatenate data and preds across all runs in each session
            for ss = 1:length(sessions)
                predictors = []; tc = [];
                for rr = 1:nruns(ss)
                    [nframes, npreds] = size(model.run_preds{rr, ss});
                    pred = zeros(nframes, npreds + nruns(ss));
                    pred(:, 1:npreds) = model.run_preds{rr, ss};
                    pred(:, npreds + rr) = 1; % add nuisance run regressors
                    predictors = [predictors; pred];
                    tc = [tc; roi.run_avgs{rr, ss}(:) - roi.baseline{rr, ss}];
                end
                roi.model.run_tcs{ss} = tc;
                % fit GLM and store betas, SEMs, and variance explained
                mm = glmTS(tc, predictors);
                roi.model.run_preds{ss} = predictors * mm.betas';
                roi.model.betas{ss} = mm.betas(1:npreds);
                roi.model.stdevs{ss} = mm.stdevs(1:npreds);
                roi.model.rbetas{ss} = mm.betas(npreds + 1:npreds + nruns(ss)); % nuisance regressor betas
                roi.model.rstdevs{ss} = mm.stdevs(npreds + 1:npreds + nruns(ss)); % nuisance regressor stdevs
                roi.model.varexp{ss} = 1 - (sum(mm.residual.^2) / sum((tc - mean(tc)).^2));
            end
            % optimize model parameters if applicable
            omodels = {'cts' 'cts-norm' 'dcts' '2ch-cts' '2ch-dcts' '3ch'};
            if optimize_flag && sum(strcmp(model.type, omodels))
                param_names = fieldnames(model.params);
                for ss = 1:length(sessions)
                    % load grid search results if saved, otherwise compute
                    fname = ['grid_search_results_' model.type '_fit' [model.experiments{:}] '.mat'];
                    fpath = fullfile(sessions{ss}, 'ROIs', roi.nickname, fname);
                    if exist(fpath, 'file') == 2
                        fprintf('Loading grid search results. \n');
                        load(fpath);
                    else
                        [rois, models] = grid_search(roi, model, ss, 5);
                        save(fpath, 'rois', 'models', '-v7.3');
                    end
                    % load grad desc results if saved, otherwise compute
                    fname = ['grad_desc_results_' model.type '_fit' [model.experiments{:}] '.mat'];
                    fpath = fullfile(sessions{ss}, 'ROIs', roi.nickname, fname);
                    if exist(fpath, 'file') == 2
                        fprintf('Loading gradient descent results. \n');
                        load(fpath);
                    else
                        [rois, models] = optimize_fit(rois, models);
                        save(fpath, 'rois', 'models', '-v7.3');
                    end
                    for pp = 1:length(param_names)
                        model.params.(param_names{pp}){ss} = models.params.(param_names{pp}){1};
                        model = update_param(model, param_names{pp}, 0);
                    end
                    model = pred_runs(model);
                    model = pred_trials(model);
                    [roi, model] = tc_fit(roi, model, 0);
                end
            end
            % carry over model parameters to roi model struct
            roi.model.type = model.type;
            roi.model.cond_list = model.cond_list;
            roi.model.cat_list = unique([model.cats{:}]);
            roi.model.pre_dur = model.pre_dur;
            roi.model.post_dur = model.post_dur;
            roi.model.fit_exps = model.experiments;
            roi.model.params = model.params;
            roi.model.irfs = model.irfs;
        end
        
        % predict responses for each trial type using model solution
        function roi = tc_pred(roi, model)
            nexps = length(roi.experiments);
            nconds = max(cellfun(@length, model.cond_list));
            ncats = length(unique([model.cats{:}]));
            nsubs = length(roi.sessions);
            % preallocate predictor array for each trial type
            roi.pred = cell(nconds, nsubs, nexps);
            % preallocate S and T predictors if using multi-channel model
            cmodels = {'2ch' '2ch-lin' '2ch-cts' '2ch-dcts' '3ch' 'htd'};
            if sum(strcmp(model.type, cmodels))
                roi.predS = cell(nconds, nsubs, nexps);
                roi.predT = cell(nconds, nsubs, nexps);
                if strcmp(model.type, '3ch')
                    roi.predD = cell(nconds, nsubs, nexps);
                end
            end
            % predict response for each session, experiment, and trial type
            for ss = 1:nsubs
                for ee = 1:nexps
                    for cc = 1:length(model.cond_list{ee})
                        % if using multi-channel model
                        if sum(strcmp(model.type, cmodels))
                            % find beta values for each channel
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
                            % add delay channel if using a 3ch model
                            if strcmp(model.type, '3ch')
                                ampD = roi.model.betas{ss}(2 * ncats + 1:3 * ncats);
                                % scale trial predictors by betas
                                predD = model.trial_preds.D{cc, ss, ee};
                                fmriD = predD .* repmat(ampD, size(predD, 1), 1);
                                roi.predD{cc, ss, ee} = fmriD;
                                roi.pred{cc, ss, ee} = fmriS + fmriT + fmriD;
                            else
                                roi.pred{cc, ss, ee} = fmriS + fmriT;
                            end
                            % if using single-channel model
                        else
                            % find beta values
                            amp = roi.model.betas{ss}(1:ncats);
                            % scale trial predictors by betas
                            pred = model.trial_preds.pred{cc, ss, ee};
                            pred = pred .* repmat(amp, size(pred, 1), 1);
                            % store trial predictors in roi
                            roi.pred{cc, ss, ee} = pred;
                        end
                    end
                end
            end
        end
        
        % plot measurement vs. prediction for runs in each session
        function plot_runs(roi)
            % get design parameters
            sessions = roi.sessions; nsess = length(sessions);
            % setup figure
            fig_name = [roi.nickname ' run timecourses'];
            fig_pos = [.1 .3 .8 .6];
            fig = figTS(fig_name, fig_pos);
            % plot run time series and predictors for each session
            for ss = 1:nsess
                subplot(nsess, 1, ss); hold on;
                plot(roi.model.run_tcs{ss}, 'k');
                plot(roi.model.run_preds{ss}, 'r');
                if ss == 1
                    ylabel('% signal');
                end
                [~, session_id] = fileparts(sessions{ss});
                leg_str = {[session_id ': ' num2str(100 * roi.model.varexp{ss}, 2) '%'] 'pred'};
                legend(leg_str, 'Location', 'NorthEastOutside'); legend boxoff;
                axis tight; ylims = get(gca, 'ylim'); ylim([ylims(1) ceil(ylims(2))]);
                set(gca, 'XColor', 'w', 'TickDir', 'out', 'YTick', [0 ceil(ylims(2))], 'FontSize', 8);
            end
        end
        
        % plot measured responses for each trial type across experiments
        function plot_exps(roi)
            % get design parameters
            nexps = length(roi.experiments);
            cond_list = roi.model.cond_list;
            all_conds = unique([cond_list{:}], 'stable');
            nconds = length(all_conds);
            pre_dur = roi.model.pre_dur; post_dur = roi.model.post_dur;
            cols = lines(nexps);
            cond_idxs = idx_trials(roi);
            % setup figure
            fig_name = [roi.nickname ' trial responses'];
            fig_pos = [.1 .3 .8 .4];
            fig = figTS(fig_name, fig_pos);
            % plot responses to trials of the same type across experiments
            xcnt = 3; zlc = xcnt;
            y_max = 0; y_min = -1;
            for cc = 1:nconds
                % get duration of trial time window
                if nexps > 1
                    tl = length(roi.trial_avgs{cond_idxs(cc, 1), 1, 1});
                else
                    tl = length(roi.trial_avgs{cond_idxs(cc, 1), 1});
                end
                % plot custom zero line for trial
                plot([zlc - 1 zlc + tl], [0 0], 'k-');
                % plot measured response in peristimulus time window
                x = xcnt:xcnt + tl - 1;
                for ee = 1:nexps
                    if cond_idxs(cc, ee) > 0
                        y_m = [roi.trial_avgs{cond_idxs(cc, ee), :, ee}]';
                        [me(ee), cymin, cymax] = lineTS(x, y_m, 1, cols(ee, :), cols(ee, :), 'sem');
                        y_min = min([y_min cymin]); y_max = max([y_max cymax]);
                    end
                end
                % plot stimulus
                plot([xcnt + pre_dur - 1 xcnt + tl - post_dur], [-.5 -.5], 'k-', 'LineWidth', 4);
                text(xcnt + pre_dur - 1, -.8, all_conds{cc}, 'FontSize', 8);
                xcnt = xcnt + tl + 3; zlc = xcnt;
            end
            % format plot
            ylabel('fMRI (% signal)');
            ylim([floor(y_min) ceil(y_max)]);
            legend(me(:), roi.experiments, 'Location', 'NorthWest'); legend boxoff;
            title(roi.nickname, 'Interpreter', 'none');
            set(gca, 'XColor', 'w', 'TickDir', 'out', 'FontSize', 8);
        end
        
        % plot measurement vs. prediction for each trial type
        function plot_model(roi, save_flag)
            if nargin == 1
                save_flag = 0;
            elseif nargin > 2 || nargin < 1
                error('unexpected input arguements');
            end
            % get design parameters and label data
            nexps = length(roi.experiments);
            sessions = roi.sessions; nsess = length(sessions);
            npreds = length(roi.model.betas{1});
            xlabs = label_preds(roi.model);
            amps = reshape([roi.model.betas{:}], npreds, [])';
            R2 = [roi.model.varexp{:}];
            % setup figure
            fig_name = [roi.nickname ' - ' roi.model.type ' model'];
            fig_pos = [.1 .1 .8 .3 + nexps * .2];
            fig = figTS(fig_name, fig_pos);
            % plot model solution
            subplot(1 + nexps, 2, 1); hold on; barTS(amps, [0 0 0]);
            xlabel('Predictor'); ylabel('Beta (% signal)');
            tstr1 = roi.nickname;
            tstr2 = [roi.model.type ' model fit to ' strjoin(roi.model.fit_exps, '/')];
            tstr3 = ['R^{2} in ' strjoin(roi.experiments, '/') ' = ' num2str(mean(R2), 3)];
            title({tstr1; tstr2; tstr3});
            set(gca, 'XTickLabel', xlabs);
            % plot variance explained for each session
            varexp = [roi.model.varexp{:}];
            subplot(1 + nexps, 2, 2); hold on; barTS(R2, [0 0 0]); xcnt = 1;
            for ss = 1:nsess
                ypos = max([0 varexp(ss)]) + .1;
                lab = num2str(varexp(ss), 2);
                text(xcnt, ypos, lab, 'HorizontalAlignment', 'center', 'FontSize', 6);
                xcnt = xcnt + 1;
            end
            xlabel('Session'); ylabel('R^2'); ylim([0 1]);
            title('Individual Subjects');
            set(gca, 'XTickLabel', strrep(roi.session_ids, '_', '-'));
            % plot measurement vs prediction for each trial type
            pre_dur = roi.model.pre_dur;
            post_dur = roi.model.post_dur;
            y_max = 0; y_min = -2;
            for ee = 1:nexps
                ax(ee) = subplot(1 + nexps, 1, ee + 1); hold on;
                xcnt = 3; zlc = xcnt;
                for cc = 1:length(roi.trial_avgs(:, 1, ee))
                    % plot custom zero line for trial
                    tl = length(roi.trial_avgs{cc, 1, ee});
                    plot([zlc - 1 zlc + tl], [0 0], 'k-');
                    % plot measured response for peristimulus time window
                    x = xcnt:xcnt + tl - 1;
                    y_m = [roi.trial_avgs{cc, :, ee}]';
                    [me, cymin, cymax] = lineTS(x, y_m, 1, [.9 .9 .9], [.7 .7 .7], 'std');
                    y_min = min([y_min cymin]); y_max = max([y_max cymax]);
                    % plot model prediction for peristimulus time window
                    y_p = [roi.pred_sum{cc, :, ee}]';
                    [pr, cymin, cymax] = lineTS(x, y_p, 2, [0 0 0]);
                    y_min = min([y_min cymin]); y_max = max([y_max cymax]);
                    % plot separate channel contributions if applicable
                    if sum(strcmp(roi.model.type, {'2ch' '2ch-lin' '2ch-cts' '2ch-dcts' '3ch' 'htd'}))
                        y_pS = [roi.predS_sum{cc, :, ee}]';
                        [sp, cymin, cymax] = lineTS(x, y_pS, 1, [0 0 1]);
                        y_min = min([y_min cymin]); y_max = max([y_max cymax]);
                        y_pT = [roi.predT_sum{cc, :, ee}]';
                        [tp, cymin, cymax] = lineTS(x, y_pT, 1, [1 0 0]);
                        y_min = min([y_min cymin]); y_max = max([y_max cymax]);
                        if strcmp(roi.model.type, '3ch')
                            y_pD = [roi.predD_sum{cc, :, ee}]';
                            [dp, cymin, cymax] = lineTS(x, y_pD, 1, [0 1 0]);
                            y_min = min([y_min cymin]); y_max = max([y_max cymax]);
                        end
                    end
                    % plot stimulus
                    plot([xcnt + pre_dur - 1 xcnt + tl - post_dur], [-.5 -.5], 'k-', 'LineWidth', 4);
                    text(xcnt + pre_dur - 1, -1, roi.model.cond_list{ee}(cc), 'FontSize', 8);
                    xcnt = xcnt + tl + 3; zlc = xcnt;
                end
                % set legend and format plot
                if sum(strcmp(roi.model.type, {'standard' 'cts' 'cts-norm' 'dcts'}))
                    leg_str = {roi.nickname [roi.model.type ' model']};
                    legend([me pr], leg_str, 'Location', 'NorthWest');
                else
                    leg_str1 = [roi.nickname ' (N = ' num2str(nsess) ')'];
                    leg_str2 = [roi.model.type ' model'];
                    leg_str3 = {'S contribution' 'T contribution'};
                    leg_str = [leg_str1 leg_str2 leg_str3];
                    if strcmp(roi.model.type, '3ch')
                        leg_str = [leg_str 'D contribution'];
                        legend([me pr sp tp dp], leg_str, 'Location', 'NorthWest');
                    else
                        legend([me pr sp tp], leg_str, 'Location', 'NorthWest');
                    end
                end
                legend boxoff;
                title([roi.experiments{ee}], 'FontSize', 8);
                ylabel('fMRI (% signal)');
                set(gca, 'XColor', 'w', 'TickDir', 'out', 'FontSize', 8);
            end
            % norm y-axis limit across experiments
            for ee = 1:nexps
                set(ax(ee), 'YLim', [floor(y_min) ceil(y_max)], 'YTick', floor(y_min):ceil(y_max));
            end
            % save to results directory if applicable
            if save_flag
                fpath = fullfile(roi.project_dir, 'figures');
                fname = [roi.nickname '_' roi.model.type 'Model' ...
                    '_fit' [roi.model.fit_exps{:}] ...
                    '_test' [roi.experiments{:}] ...
                    '_' date '.jpg'];
                saveas(fig, fullfile(fpath, fname), 'jpg');
            end
        end
        
        % recompute model performance given indepdendently-fit weights
        function roi = recompute(roi, model, fit)
            check_model(roi, model);
            [nruns, nsubs] = size(model.run_preds);
            for ss = 1:nsubs
                run_preds = []; tc = [];
                npreds = size(model.run_preds{1, ss}, 2);
                nsruns = length(roi.model.rbetas{ss});
                for rr = 1:nsruns
                    nframes = size(model.run_preds{rr, ss}, 1);
                    pm = zeros(nframes, npreds + nsruns);
                    pm(:, 1:npreds) = model.run_preds{rr, ss};
                    pm(:, npreds + rr) = 1;
                    run_preds = [run_preds; pm];
                    tc = [tc; roi.run_avgs{rr, ss} - roi.baseline{rr, ss}];
                end
                roi.model.run_tcs{ss} = tc;
                beta_vec = zeros(1, npreds + nsruns);
                beta_vec(1:npreds) = fit.betas{ss};
                beta_vec(npreds + 1:npreds + nsruns) = roi.model.rbetas{ss};
                % predict predict fMRI responses for each run
                run_pred = run_preds * beta_vec';
                roi.model.run_preds{ss} = run_pred;
                res = tc - run_pred;
                % calculate model performance
                roi.model.varexp{ss} = 1 - (sum(res.^2) / sum((tc - mean(tc)).^2));
            end
            % store new fit in roi structure
            roi.model.betas = fit.betas;
            roi.model.stdevs = fit.stdevs;
            roi.model.fit_exps = fit.fit_exps;
            roi = tc_pred(roi, model);
        end
        
    end
    
end
