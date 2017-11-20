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
%   tch_recompute -- validates model solution on indpendent data
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
        noise_ceils = {}; % estiamte of noice ceiling for each experiment
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
            filenames = {}; nruns = vox.num_runs;
            % for each session
            for ss = 1:length(vox.sessions)
                spath = fullfile(vox.sessions{ss}, 'Voxels'); rcnt = 0;
                % for each experiment
                for ee = 1:length(vox.experiments)
                    % store paths to data file for each run
                    for rr = 1:nruns(ee, ss)
                        fname = ['Run' num2str(rr) '.mat']; rcnt = rcnt + 1;
                        filenames{rcnt, ss} = fullfile(spath, vox.experiments{ee}, fname);
                    end
                end
            end
        end
        
        % label each session with an ID string
        function session_ids = get.session_ids(vox)
            session_ids = cell(1, length(vox.sessions));
            for ss = 1:length(vox.sessions)
                [~, session_ids{ss}] = fileparts(vox.sessions{ss});
            end
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
        function vox = tch_runs(vox, detrend_option)
            if nargin == 1; detrend_option = 3; end
            vox = select_sessions(vox); fpaths = vox.filenames;
            raw_runs = cellfun(@(X) tch_load(X, 'tSeries'), fpaths, 'uni', false);
            vox.runs = cellfun(@(X) tch_psc(X, detrend_option), raw_runs, 'uni', false);
        end
        
        % check dimensionality of vox time series and model predictions
        function check_model(vox, model)
            empty_cells = cellfun(@isempty, model.run_durs);
            rds = model.run_durs; rds(empty_cells) = {0};
            comp = cellfun(@(X) size(X, 1), vox.runs) ~= cell2mat(rds);
            if sum(comp(:)) > 0
                error('Dimensions of data and model do not match');
            end
        end
        
        % compile time series for each trial type
        function vox = tch_trials(vox, model)
            % check model and get design parameters
            check_model(vox, model); nruns = model.num_runs;
            nconds = max(cellfun(@length, model.cond_list));
            % get onsets, offsets, and trial orders
            onsets = model.tonsets; offsets = model.toffsets;
            % compile all trial time series for each session and experiment
            trials = cell(nconds, length(vox.sessions), length(vox.experiments));
            for ss = 1:length(vox.sessions)
                rcnt = 0;
                for ee = 1:length(vox.experiments)
                    cnts = zeros(1, length(model.cond_list{ee}));
                    for rr = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        % estimate prestimulus baseline response for run
                        oframes = repmat(onsets{rcnt, ss}, model.pre_dur, 1);
                        psf = oframes - repmat((1:model.pre_dur)', 1, length(oframes));
                        % calculate mean baseline response to subtract
                        bs = mean(vox.runs{rcnt, ss}(psf(:), :));
                        vox.baseline{rcnt, ss} = bs;
                        % store peri-stimulus time series sorted by trial
                        for tt = 1:length(onsets{rcnt, ss})
                            % get condition number of this trial
                            cond = model.tconds{rcnt, ss}(tt);
                            cond_idx = find(strcmp(cond, model.cond_list{ee}));
                            cnts(cond_idx) = cnts(cond_idx) + 1;
                            % get TR corresponding to onset of pre_dur
                            onset = (onsets{rcnt, ss}(tt) - model.pre_dur) / vox.tr + 1;
                            off_idxs = offsets{rcnt, ss}(tt) - model.gap_dur / 2;
                            offset = (floor(off_idxs) + model.post_dur) / vox.tr + 1;
                            % extract the peri-stimulus time window
                            bsm = repmat(bs, length(onset:offset), 1);
                            trial = vox.runs{rcnt, ss}(onset:offset, :) - bsm;
                            trials{cond_idx, ss, ee}(:, :, cnts(cond_idx)) = trial;
                        end
                    end
                end
            end
            vox.trials = trials;
        end
        
        % use GLM to fit weights for each predictor in model
        function [vox, model] = tch_fit(vox, model)
            check_model(vox, model); npreds = size(model.run_preds{1, 1}, 2);
            % subtract baseline estimates from centered time series
            bsm = cellfun(@(X, Y) repmat(X, size(Y, 1), 1), ...
                vox.baseline, vox.runs, 'uni', false);
            tcs = cellfun(@(X, Y) X - Y, vox.runs, bsm, 'uni', false);
            % concatenate data and preds across all runs in each session
            for ss = 1:length(vox.sessions)
                nruns = sum(model.num_runs(:, ss));
                % construct nuisance regressors and merge with predictors
                b0 = cell2mat(cellfun(@(X, Y) code_stim_vec(zeros(X, nruns), 1:X, Y), ...
                    model.run_durs(:, ss), num2cell(1:size(tcs, 1))', 'uni', false));
                predictors = [cell2mat(model.run_preds(:, ss)) b0];
                % fit GLM and store betas, SEMs, and variance explained
                tc = cell2mat(tcs(:, ss)); vox.model.run_tcs{ss} = tc;
                mm = tch_glm(tc, predictors);
                vox.model.run_preds{ss} = predictors * squeeze(mm.betas);
                vox.model.betas{ss} = mm.betas(:, 1:npreds, :);
                vox.model.stdevs{ss} = mm.stdevs(:, 1:npreds, :);
                vox.model.residual{ss} = mm.residual;
                res_var = sum(mm.residual .^ 2) ./ sum((tc - mean(tc)) .^ 2);
                vox.model.varexp{ss} = ones(1, size(tc, 2)) - res_var;
                % store parameters of nuisance regressors
                vox.model.rbetas{ss} = mm.betas(:, npreds + 1:npreds + nruns(ss), :);
                vox.model.rstdevs{ss} = mm.stdevs(:, npreds + 1:npreds + nruns(ss), :);
            end
            % carry over model parameters to vox model struct
            vox.model.type = model.type;
            vox.model.num_channels = model.num_channels;
            vox.model.cond_list = model.cond_list;
            vox.model.cat_list = unique([model.cats{:}]);
            vox.model.pre_dur = model.pre_dur;
            vox.model.post_dur = model.post_dur;
            vox.model.fit_exps = model.experiments;
            vox.model.params = model.params;
            vox.model.irfs = model.irfs;
        end
        
        % recompute model performance given indepdendently-fit weights
        function vox = tch_recompute(vox, model, fit)
            check_model(vox, model);
            bsm = cellfun(@(X, Y) repmat(X, size(Y, 1), 1), ...
                vox.baseline, vox.run_avgs, 'uni', false);
            tcs = cellfun(@(X, Y) X - Y, vox.run_avgs, bsm, 'uni', false);
            for ss = 1:length(model.sessions)
                nruns = sum(model.num_runs(:, ss));
                % construct nuisance regressors and merge with predictors
                b0 = cell2mat(cellfun(@(X, Y) code_stim_vec(zeros(X, nruns), 1:X, Y), ...
                    model.run_durs(:, ss), num2cell(1:size(tcs, 1))', 'uni', false));
                predictors = [cell2mat(model.run_preds(:, ss)) b0];
                % predict fMRI responses for each run
                tc = cell2mat(tcs(:, ss));
                beta_vecs = [fit.betas{ss} vox.model.rbetas{ss}];
                run_preds = predictors * squeeze(beta_vecs);
                % calculate cross-validated model performance
                vox.model.run_preds{ss} = run_preds; res = tc - run_preds;
                res_var = sum(res .^ 2) ./ sum((tc - mean(tc)) .^ 2);
                vox.model.varexp{ss} = ones(1, size(tc, 2)) - res_var;
            end
            % store new fit in vox structure
            vox.model.betas = fit.betas; vox.model.stdevs = fit.stdevs;
            vox.model.fit_exps = fit.fit_exps;
        end
        
    end
    
end
