% Stores and operates on fMRI time series of individual voxels across
% multiple scan sessions. Used with ModelTS object to fit and validate
% various temporal encoding models in each voxel. 
% 
% CONSTRUCTOR INPUTS
%   1) exps: list of experiments to model (e.g., {'Exp1' 'Exp2'})
%   2) session_list: list of sessions to model (optional)
%
% METHODS
%   tc_runs -- preprocesses and stores time series of all voxels
%   tc_trials -- compiles trial-level time series sorted by experiment
%   tc_fit -- fits ModelTS object to the mean time series of each voxel
%   recompute -- validates model solution on indpendent data
% 
% Example model fitting steps ("model" is a ModelTS object):
%   vox = Voxel({'Exp1' 'Exp2'})
%   vox = tc_runs(vox)
%   vox = tc_trials(vox, model)
%   vox = tc_fit(vox, model)
%   vox = tc_pred(vox, model)
% 
% AS 7/2017

classdef Voxel
    
    properties
        experiments % array of experiments to model
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
    
    properties (Dependent)
        sessions   % paths to sessions that have current ROI
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
        function vox = Voxel(exps, isessions)
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
        
        % find set of all_sessions with voxel data
        function sessions = get.sessions(vox)
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
            % error if no sessions with voxel data are found
            if scnt == 0
                error('No sessions found with voxel data.');
            end
        end
        
        % find the number of runs per experiment for each session
        function num_runs = get.num_runs(vox)
            sessions = vox.sessions; nsess = length(sessions);
            num_runs = zeros(length(vox.experiments), length(vox.sessions));
            for ss = 1:nsess
                spath = fullfile(sessions{ss}, 'Voxels');
                % find paths to data files for each experiment
                for ee = 1:length(vox.experiments)
                    ecnt = 1;
                    while exist(fullfile(spath, vox.experiments{ee}, ['Run' num2str(ecnt) '.mat']), 'file') == 2
                        num_runs(ee, ss) = num_runs(ee, ss) + 1;
                        ecnt = ecnt + 1;
                    end
                end
            end
        end
        
        % find the paths to the data files for each session
        function filenames = get.filenames(vox)
            sessions = vox.sessions; nsess = length(sessions);
            nruns = vox.num_runs;
            filenames = {};
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
        
        % preprocess and store run timeseries of each voxel
        function vox = tc_runs(vox)
            fpaths = vox.filenames; % paths to data files
            raw_runs = cellfun(@(X) loadTS(X, 'tSeries'), fpaths, 'uni', false);
            vox.runs = cellfun(@(X) psc(X), raw_runs, 'uni', false);
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
        function vox = tc_trials(vox, model)
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
        function [vox, model] = tc_fit(vox, model, optimize_flag)
            if nargin < 3
                optimize_flag = 0;
            end
            check_model(vox, model);
            nruns_max = size(model.run_preds, 1); % max number of runs
            nruns = nruns_max - sum(cellfun(@isempty, model.run_preds));
            % concatenate data and preds across all runs in each session
            for ss = 1:length(vox.sessions)
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
                mm = glmTS(tc, predictors);
                vox.model.run_preds{ss} = predictors * squeeze(mm.betas);
                vox.model.betas{ss} = mm.betas(:, 1:npreds, :);
                vox.model.stdevs{ss} = mm.stdevs(:, 1:npreds, :);
                vox.model.rbetas{ss} = mm.betas(:, npreds + 1:npreds + nruns(ss), :); % nuisance regressor betas
                vox.model.rstdevs{ss} = mm.stdevs(:, npreds + 1:npreds + nruns(ss), :); % nuisance regressor stdevs
                vox.model.varexp{ss} = ones(1, size(tc, 2)) - (sum(mm.residual .^ 2, 1) ./ sum((tc - repmat(mean(tc, 1), size(tc, 1), 1)) .^ 2, 1));
            end
            % optimize model parameters if applicable
            omodels = {'cts' 'cts-norm' 'dcts' '2ch-cts' '2ch-dcts' '3ch'};
            if optimize_flag && sum(strcmp(model.type, omodels))
                % load grid search results if saved, otherwise compute
                fname = ['grid_search_' vox.nickname '_' model.type '.mat'];
                fpath = fullfile(vox.project_dir, 'tmp', fname);
                if exist(fpath, 'file') == 2
                    fprintf('Loading grid search results. \n');
                    load(fpath);
                else
                    [voxs, models] = grid_search(vox, model, 3);
                    save(fpath, 'voxs', 'models', '-v7.3');
                end
                % load grad desc results if saved, otherwise compute
                fname = ['grad_desc_' vox.nickname '_' model.type '.mat'];
                fpath = fullfile(vox.project_dir, 'tmp', fname);
                if exist(fpath, 'file') == 2
                    fprintf('Loading gradient descent results. \n');
                    load(fpath);
                else
                    [voxs, models] = optimize_fit(voxs, models);
                    save(fpath, 'voxs', 'models', '-v7.3');
                end
                model.params = models.params;
                model.irfs = models.irfs;
                model = pred_runs(model);
                model = pred_trials(model);
                [vox, model] = tc_fit(vox, model, 0);
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
