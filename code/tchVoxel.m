% Stores and operates on fMRI time series of individual voxels across
% multiple scan sessions. Used with tchModel object to fit and validate
% various temporal encoding models in each voxel. 
% 
% CONSTRUCTOR INPUTS
%   1) exps: list of experiments to model (e.g., {'Exp1' 'Exp2'})
%   2) session_list: list of sessions to model (optional)
%
% METHODS
%   tch_runs -- preprocesses and stores time series of all voxels
%   tch_trials -- compiles trial-level time series sorted by experiment
%   tch_fit -- fits tchModel object to the mean time series of each voxel
%   recompute -- validates model solution on indpendent data
% 
% Example model fitting steps ("model" is a tchModel object):
%   vox = tchVoxel({'Exp1' 'Exp2'})
%   vox = tch_runs(vox)
%   vox = tch_trials(vox, model)
%   vox = tch_fit(vox, model)
% 
% AS 7/2017

classdef tchVoxel
    
    properties
        experiments % array of experiments to model
        sessions    % array of experimental sessions
        model = []; % data structure of models fits for each session
    end
    
    properties (Hidden)
        runs = {};      % voxel responses for each run (TRs, voxels)
        trials = {};    % voxel responses for each trial type (TRs, voxels)
        baseline = {}   % mean baseline response across all trial types
        isessions = {}; % user-specified session list (optional)
    end
    
    properties (Constant, Hidden)
        project_dir = fileparts(fileparts(which(mfilename, 'class')));
        tr = 1; % fMRI TR (s)
    end
    
    properties (Dependent, Hidden)
        all_sessions % paths to all session directories
        num_runs     % number of runs per experiment
        filenames    % paths to data files from sessions
        session_ids  % session nicknames
        pred_sum     % sum of preds across all channels
        predS_sum    % sum of preds across S channels
        predT_sum    % sum of preds across T channels
        predD_sum    % sum of preds across D channels
    end
    
    
    methods
        
        % class constructor
        function vox = tchVoxel(exps, isessions)
            if nargin == 1
                vox.experiments = force_cell(exps);
            elseif nargin == 2
                vox.experiments = force_cell(exps);
                vox.isessions = force_cell(isessions);
            else
                error('incorrect input arguments');
            end
        end
        
        % find all sessions in data directory if not specified by user
        function all_sessions = get.all_sessions(vox)
            if isempty(vox.isessions)
                all_sessions = find_sessions(vox.project_dir);
            else
                all_sessions = vox.isessions;
            end
        end
                
        % find the number of runs per experiment for each session
        function num_runs = get.num_runs(vox)
            num_runs = zeros(length(vox.experiments), length(vox.sessions));
            for ss = 1:length(vox.sessions)
                spath = fullfile(vox.sessions{ss}, 'Voxels');
                % find paths to data files for each experiment
                for ee = 1:length(vox.experiments)
                    d = dir(fullfile(spath, vox.experiments{ee}, 'Run*.mat'));
                    num_runs(ee, ss) = length({d.name});
                end
            end
        end
        
        % find the paths to the data files for each session
        function filenames = get.filenames(vox)
            sessions = vox.sessions; nsess = length(sessions);
            filenames = {}; nruns = vox.num_runs;
            % for each session
            for ss = 1:nsess
                rcnt = 0;
                spath = fullfile(sessions{ss}, 'Voxels');
                % for each experiment
                for ee = 1:length(vox.experiments)
                    % store paths to data file for each run
                    for rr = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        filenames{rcnt, ss} = fullfile(spath, vox.experiments{ee}, ['Run' num2str(rr) '.mat']);
                    end
                end
            end
        end
        
        % label each session with an ID string
        function session_ids = get.session_ids(vox)
            sessions = vox.sessions;
            session_ids = cell(1, length(sessions));
            for ss = 1:length(sessions)
                [~, session_id] = fileparts(sessions{ss});
                session_ids{ss} = session_id;
            end
        end
        
        % sum trial predictors across all channels
        function pred_sum = get.pred_sum(vox)
            pred_sum = cellfun(@(X) sum(X, 2), vox.pred, 'uni', false);
        end
        
        % sum trial predictors across sustained channels
        function predS_sum = get.predS_sum(vox)
            predS_sum = cellfun(@(X) sum(X, 2), vox.predS, 'uni', false);
        end
        
        % sum trial predictors across transient channels
        function predT_sum = get.predT_sum(vox)
            predT_sum = cellfun(@(X) sum(X, 2), vox.predT, 'uni', false);
        end
        
        % sum trial predictors across all delay channels
        function predD_sum = get.predD_sum(vox)
            predD_sum = cellfun(@(X) sum(X, 2), vox.predD, 'uni', false);
        end
        
        % find set of all_sessions with experiments
        function vox = select_sessions(vox)
            sessions = {}; scnt = 0;
            for ss = 1:length(vox.all_sessions)
                [~, session_id] = fileparts(vox.all_sessions{ss});
                spath = fullfile(vox.project_dir, 'data', session_id);
                cpath = fullfile(spath, 'Voxels');
                ecnt = 0;
                for ee = 1:length(vox.experiments)
                    if exist(fullfile(cpath, vox.experiments{ee}), 'dir') == 7
                        ecnt = ecnt + 1;
                    end
                end
                if ecnt == length(vox.experiments)
                    scnt = scnt + 1;
                    sessions{scnt} = spath;
                end
            end
            % error if no sessions with specificed experiments are found
            if scnt == 0
                error('No sessions found with all experiments');
            else
                vox.sessions = sessions;
            end
        end
        
        % preprocess and store run timeseries of each voxel
        function vox = tch_runs(vox)
            vox = select_sessions(vox); % sessions with all experiments
            fpaths = vox.filenames;     % paths to data files
            raw_runs = cellfun(@(X) tch_load(X, 'tSeries'), fpaths, 'uni', false);
            vox.runs = cellfun(@(X) tch_psc(X), raw_runs, 'uni', false);
        end
        
        % check dimensionality of vox time series and model predictions
        function check_model(vox, model)
            empty_cells = cellfun(@isempty, model.run_durs);
            rds = model.run_durs; rds(empty_cells) = {0};
            comp = cellfun(@(X) size(X, 1), vox.runs) ~= cell2mat(rds);
            if sum(comp(:)) > 0
                error('dimensions of data and model do not match');
            end
        end
        
        % compile time series for each trial type
        function vox = tch_trials(vox, model)
            % check model and get design parameters
            check_model(vox, model);
            sessions = vox.sessions; nsess = length(sessions);
            nconds = max(cellfun(@length, model.cond_list));
            nexps = length(vox.experiments);
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
                        bs = mean(vox.runs{rcnt, ss}(bf, :));
                        vox.baseline{rcnt, ss} = bs;
                        % store peri-stimulus time series sorted by trial
                        for tt = 1:length(onsets{rcnt, ss})
                            % get condition number of this trial
                            cond = conds{rcnt, ss}(tt);
                            cond_idx = find(strcmp(cond, model.cond_list{ee}));
                            cnts(cond_idx) = cnts(cond_idx) + 1;
                            % get TR corresponding to onset of pre_dur
                            onset = (onsets{rcnt, ss}(tt) - model.pre_dur) / vox.tr + 1;
                            offset = (floor(offsets{rcnt, ss}(tt) - model.gap_dur / 2) + model.post_dur) / vox.tr + 1;
                            % extract the peri-stimulus time window
                            trials{cond_idx, ss, ee}(:, :, cnts(cond_idx)) = vox.runs{rcnt, ss}(onset:offset, :) - repmat(bs, length(onset:offset), 1);
                        end
                    end
                end
            end
            vox.trials = trials;
        end
        
        % use GLM to fit weights for each predictor in model
        function [vox, model] = tch_fit(vox, model, optimize_flag, fit_exps)
            if nargin < 3; optimize_flag = 0; end
            if nargin < 4; fit_exps = model.experiments; end
            check_model(vox, model); sessions = vox.sessions;
            nruns_max = size(model.run_preds, 1);
            nruns = nruns_max - sum(cellfun(@isempty, model.run_preds));
            % concatenate data and preds across all runs in each session
            for ss = 1:length(sessions)
                predictors = []; tc = [];
                for rr = 1:nruns(ss)
                    [nframes, npreds] = size(model.run_preds{rr, ss});
                    pred = zeros(nframes, npreds + nruns(ss)); 
                    pred(:, 1:npreds) = model.run_preds{rr, ss};
                    pred(:, npreds + rr) = 1; % add nuisance run regressors
                    predictors = [predictors; pred];
                    tc = [tc; vox.runs{rr, ss} - repmat(vox.baseline{rr, ss}, size(vox.runs{rr, ss}, 1), 1)];
                end
                vox.model.run_tcs{ss} = tc;
                % fit GLM and store betas, SEMs, and variance explained
                mm = tch_glm(tc, predictors);
                vox.model.run_preds{ss} = predictors * squeeze(mm.betas);
                vox.model.betas{ss} = mm.betas(:, 1:npreds, :);
                vox.model.stdevs{ss} = mm.stdevs(:, 1:npreds, :);
                vox.model.rbetas{ss} = mm.betas(:, npreds + 1:npreds + nruns(ss), :); % nuisance regressor betas
                vox.model.rstdevs{ss} = mm.stdevs(:, npreds + 1:npreds + nruns(ss), :); % nuisance regressor stdevs
                vox.model.varexp{ss} = ones(1, size(tc, 2)) - (sum(mm.residual .^ 2, 1) ./ ...
                    sum((tc - repmat(mean(tc, 1), size(tc, 1), 1)) .^ 2, 1));
            end
            % optimize model parameters if applicable
            omodels = {'cts-pow' 'cts-div' 'dcts' ...
                '2ch-pow-quad' '2ch-pow-rect' '2ch-div' '2ch-dcts' '2ch-opt' ...
                '3ch-lin-quad' '3ch-lin-rect' '3ch-pow-quad' '3ch-pow-rect' '3ch-opt'};
            if optimize_flag && sum(strcmp(model.type, omodels))
                param_names = fieldnames(model.params);
                for ss = 1:length(sessions)
                    fname_grid = ['grid_search_results_' model.type ...
                        '_fit' [fit_exps{:}] '.mat'];
                    fpath_grid = fullfile(sessions{ss}, 'Voxels', fname_grid);
                    fname_grad = ['grad_desc_results_' model.type ...
                        '_fit' [fit_exps{:}] '.mat'];
                    fpath_grad = fullfile(sessions{ss}, 'Voxels', fname_grad);
                    % load optimization results if saved, otherwise compute
                    if exist(fpath_grad, 'file') == 2
                        fprintf('Loading gradient descent results. \n');
                        load(fpath_grad);
                    elseif exist(fpath_grid, 'file') == 2
                        fprintf('Loading grid search results. \n');
                        load(fpath_grid);
                        [voxels, models] = tch_optimize_fit(voxels, models);
                        save(fpath_grad, 'voxels', 'models', '-v7.3');
                    else
                        [voxels, models] = tch_grid_search(vox, model, ss, 5);
                        save(fpath_grid, 'voxels', 'models', '-v7.3');
                        [rois, models] = tch_optimize_fit(rois, models);
                        save(fpath_grad, 'voxels', 'models', '-v7.3');
                    end
                    % copy optimized parameters for session
                    for pp = 1:length(param_names)
                        opt_params = models.params.(param_names{pp}){1};
                        model.params.(param_names{pp}){ss} = opt_params;
                        model = update_param(model, param_names{pp}, 0);
                    end
                    model = pred_runs(model);
                    model = pred_trials(model);
                    [vox, model] = tch_fit(vox, model, 0);
                end
            end
            % carry over model parameters to vox model struct
            vox.model.type = model.type;
            vox.model.cond_list = model.cond_list;
            vox.model.cat_list = unique([model.cats{:}]);
            vox.model.pre_dur = model.pre_dur;
            vox.model.post_dur = model.post_dur;
            vox.model.fit_exps = model.experiments;
            vox.model.params = model.params;
            vox.model.irfs = model.irfs;
        end
        
        % recompute model performance given indepdendently-fit weights
        function vox = recompute(vox, model, fit)
            check_model(vox, model);
            [nruns, nsubs] = size(model.run_preds);
            for ss = 1:nsubs
                run_preds = []; tc = [];
                npreds = size(model.run_preds{1, ss}, 2);
                nsruns = size(vox.model.rbetas{ss}, 2);
                for rr = 1:nsruns
                    nframes = size(model.run_preds{rr, ss}, 1);
                    pm = zeros(nframes, npreds + nsruns);
                    pm(:, 1:npreds) = model.run_preds{rr, ss};
                    pm(:, npreds + rr) = 1;
                    run_preds = [run_preds; pm];
                    tc = [tc; vox.runs{rr, ss} - repmat(vox.baseline{rr, ss}, size(vox.runs{rr, ss}, 1), 1)];
                end
                vox.model.run_tcs{ss} = tc;
                beta_vecs = zeros(1, npreds + nsruns, size(tc, 2));
                beta_vecs(1, 1:npreds, :) = fit.betas{ss};
                beta_vecs(1, npreds + 1:npreds + nsruns, :) = vox.model.rbetas{ss};
                % rectify negative weights for cross-validation
                beta_vecs(1, beta_vecs(1:npreds, :) < 0) = 0;
                % predict predict fMRI responses for each run
                run_pred = run_preds * squeeze(beta_vecs);
                vox.model.run_preds{ss} = run_pred;
                res = tc - run_pred;
                % calculate model performance
                vox.model.varexp{ss} = ones(1, size(tc, 2)) - (sum(res .^ 2, 1) ./ sum((tc - repmat(mean(tc, 1), size(tc, 1), 1)) .^ 2, 1));
            end
            % store new fit in vox structure
            vox.model.betas = fit.betas;
            vox.model.stdevs = fit.stdevs;
            vox.model.fit_exps = fit.fit_exps;
        end
        
    end
    
end
