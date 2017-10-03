% tchModel: Code for modeling fMRI responses to time-varying visual stimuli.
% 
% CONSTRUCTOR INPUTS
%   1) type: which model to use
%      Hemodynamic models:
%        'glm'      -- general linear model for fMRI (Boynton 1996)
%        'htd'      -- hemodynamic temporal derivative (HTD; Henson 2002)
%        'balloon'  -- nonlinear balloon model (Buxton, 1998)
%      Single-channel linear-nonlinear (LN) models:
%        'cts-pow'  -- CTS with power law (Zhou 2017)
%        'cts-div'  -- CTS with divisive normalization (Zhou 2017)
%        'dcts'     -- dynamic CTS (Zhou 2017)
%      Multi-channel models:
%        '2ch-lin-lin'  -- 2 channels with linear S and linsear T
%        '2ch-lin-quad' -- 2 channels with linear S and quadratic T
%        '2ch-lin-rect' -- 2 channels with linear S and rectified T
%        '2ch-pow-quad' -- 2 channels with CTS-pow on S and quadratic T
%        '2ch-pow-rect' -- 2 channels with CTS-pow on S and rectified T
%        '2ch-div-quad' -- 2 channels with CTS-norm on S and quadratic T
%        '2ch-opt'      -- 2 channels with optimized S and T IRFs
%        '3ch-lin-quad' -- 3 channels with linear S, quadratic T, and D
%        '3ch-lin-rect' -- 3 channels with linear S, rectified T, and D
%        '3ch-pow-quad' -- 3 channels with CTS-pow on S, quadratic T, and D
%        '3ch-pow-rect' -- 3 channels with CTS-pow on S, rectified T, and D
%        '3ch-opt'      -- 3-channels with optimized S, T, and D IRFs
%   2) exps: array of experiments for fitting model (e.g., {'Exp1' 'Exp2'})
%   3) sessions: array of paths to session data directories
%
% METHODS
%   code_stim   -- compiles onset and offset times for each stimulus
%   norm_model  -- customize predictor normalization parameters (optional)
%   pred_runs   -- generates design matrix for each run of data
%   pred_trials -- generates response predictors for each trial type
% 
% Example model generation steps:
%   exps = {'Exp1' 'Exp2' 'Exp3'};
%   model = tchModel('2ch', exps, ROI('V1', exps).sessions);
%   model = code_stim(model);
%   model = pred_runs(model);
%   model = pred_trials(model);
% 
% AS 2/2017

classdef tchModel
    
    properties
        type              % model type identifier
        experiments       % array of experiments to model
        sessions = {};    % array of paths to session directories
        params = [];      % model parameters structure
        irfs = [];        % model impulse response functions structure
        run_preds = [];   % fMRI predictors for each run
        trial_preds = []; % fMRI predictors for each trial type
    end
    
    properties (Hidden)
        num_exps        % number of experiments included in model
        tonsets = {};   % trial onset times (s)
        toffsets = {};  % trial offset times (s)
        tconds = {};    % condition number of each element in tonsets
        onsets = {};    % image onset times (s)
        offsets = {};   % image offset times (s)
        conds = {};     % condition number of each element in onsets
        cats = {};      % category of each element in onsets
        cond_list = {}; % list of condition names in each experiment
        run_durs = {};  % run durations (s)
        stim = {};      % stimulus step function
        normT = 20;     % transient channel normalization scalar
        normD = 20;     % delay channel normalization scalar
    end

    properties (Constant, Hidden)
        % path to project directory
        project_dir = fileparts(fileparts(which(mfilename, 'class')));
        % descriptors for each model implemented
        types = {'glm' 'htd' 'balloon' ... % hemodynamic models
            'cts-pow' 'cts-div' 'dcts' ... % single-channel LN models
            '2ch-lin-lin' '2ch-lin-quad' '2ch-lin-rect'  ... % analytical 2ch models
            '2ch-pow-quad' '2ch-pow-rect' '2ch-div-quad' '2ch-dcts-quad' '2ch-opt' ...  % optimized 2ch models
            '3ch-lin-quad' '3ch-lin-rect' '3ch-pow-quad' '3ch-pow-rect' '3ch-opt'}; % optimized 3ch models
        % experimental parameters
        tr = 1;         % fMRI TR (s)
        gap_dur = 1/60; % forced gap between stimuli (s)
        pre_dur = 4;    % pre-stimulus baseline period (s)
        post_dur = 12;  % post-stimulus baseline duration (s)
        fs = 1000;      % stimulus/neural sampling rate (Hz)
    end
    
    properties (Dependent, Hidden)
        num_runs     % numer of runs per experiment per session
        stimfiles    % paths to stimulus history files for each run
        num_channels % number of channels in model (1-3 per condition)
    end
    
    methods
        
        % tchModel class constructor
        function model = tchModel(type, exps, sessions)
            % standardize the input arguements
            if nargin > 1 && nargin < 4
                model.type = lower(type);
                model.experiments = force_cell(exps);
                if nargin == 2
                    model.sessions = find_sessions(model.project_dir);
                else
                    model.sessions = force_cell(sessions);
                end
            else
                error('Incorrect number of input arguments.');
            end
            % check for unexpected model type
            if sum(strcmp(model.type, model.types)) ~= 1
                error('Unexpected model type argument.');
            end
            % initialize model parameters
            nsess = length(model.sessions);
            [params, irfs] = tch_init_params(type, nsess, model.fs);
            model.params = params; model.irfs = irfs;
            model.num_exps = length(model.experiments);
        end
        
        % find the number of runs per experiment for each session
        function num_runs = get.num_runs(model)
            sessions = model.sessions; nsess = length(sessions);
            num_runs = zeros(length(model.experiments), nsess);
            for ss = 1:nsess
                spath = fullfile(sessions{ss}, 'Stimuli');
                % find paths to data files for each experiment
                for ee = 1:length(model.experiments)
                    d = dir(fullfile(spath, [model.experiments{ee} '_Run*.txt']));
                    fnames = {d.name}; num_runs(ee, ss) = length(fnames);
                end
            end
        end
        
        % Find paths to stimulus timing parameter files.
        function stimfiles = get.stimfiles(model)
            sessions = model.sessions; nsess = length(sessions);
            nruns = model.num_runs; stimfiles = cell(max(sum(nruns, 1)), nsess);
            for ss = 1:nsess
                rcnt = 0;
                for ee = 1:model.num_exps
                    ecnt = 1;
                    fstem = fullfile(sessions{ss}, 'Stimuli');
                    fname = [model.experiments{ee} '_Run'];
                    for ff = 1:nruns(ee, ss)
                        rcnt = rcnt + 1;
                        stimfiles{rcnt, ss} = fullfile(fstem, [fname num2str(ecnt) '.txt']);
                        ecnt = ecnt + 1;
                    end
                end
            end
        end
        
        % Find paths to stimulus timing parameter files.
        function num_channels = get.num_channels(model)
            switch model.type(1:3)
                case '2ch'
                    num_channels = 2;
                case 'htd'
                    num_channels = 2;
                case '3ch'
                    num_channels = 3;
                otherwise
                    num_channels = 1;
            end
        end
        
        % code onset, offset, and category of each stimulus in experiments
        function model = code_stim(model)
            % get session and run parameters
            stimfiles = model.stimfiles; [nruns_max, nsess] = size(stimfiles);
            fs = model.fs; gd = model.gap_dur;
            % get image onsets/offsets and condition orders from stimfiles
            [on, off, c, ims, ton, toff, tc, rd, cl] = cellfun(...
                @tch_stimfile, stimfiles, 'uni', false);
            % store stimulus parameters in tchModel object
            model.onsets = on;     % image onset times (s)
            model.offsets = off;   % image offset times (s)
            model.conds = c;       % condition labels per image
            model.cats = ims;      % category labels per image
            model.tonsets = ton;   % trial onset times (s)
            model.toffsets = toff; % trial offset times (s)
            model.tconds = tc;     % condition labels per trial
            model.run_durs = rd;   % duration of each run (s)
            % get list of stimulus categories sorted alphabetically
            cat_list = unique([ims{:}]);
            % get list of trial conditions from stimfile header
            rcnt = 1;
            for ee = 1:model.num_exps
                model.cond_list{ee} = cl{rcnt, 1};
                rcnt = rcnt + model.num_runs(ee, 1);
            end
            % find frame indices of image onsets and offsets for each cat
            stims = cellfun(@(X) zeros(X * fs, length(cat_list)), rd, 'uni', false);
            empty_cells = cellfun(@isempty, stimfiles); stims(empty_cells) = {[]};
            for cc = 1:length(cat_list)
                cat_idxs = cellfun(@(X) find(strcmp(cat_list(cc), X)), ims, 'uni', false);
                onset_idxs = cellfun(@(X, Y) round(fs * X(Y)) + 1, on, cat_idxs, 'uni', false);
                offset_idxs = cellfun(@(X, Y) round(fs * (X(Y) - gd / 2)), off, cat_idxs, 'uni', false);
                onset_idxs(empty_cells) = {1}; offset_idxs(empty_cells) = {1};
                % find frame indices of gap offset times
                gapoff_idxs = cellfun(@(X, Y) round(fs * (X(Y) + gd / 2)), off, cat_idxs, 'uni', false);
                gapoff_idxs(empty_cells) = {1};
                % compile all frame indices during stimuli and gaps
                stim_idxs = cellfun(@code_stim_idxs, onset_idxs, offset_idxs, 'uni', false);
                gap_idxs = cellfun(@code_stim_idxs, offset_idxs, gapoff_idxs, 'uni', false);
                % code stimulus as a step function with gaps at image offsets
                cc_c = repmat({cc}, nruns_max, nsess);
                stims = cellfun(@code_stim_vec, stims, stim_idxs, cc_c, repmat({1}, nruns_max, nsess), 'uni', false);
                stims = cellfun(@code_stim_vec, stims, gap_idxs, cc_c, repmat({0}, nruns_max, nsess), 'uni', false);
            end
            stims(empty_cells) = {[]};
            model.stim = stims;
        end
        
        % compute custom normalization parameters using run predictors
        function model = norm_model(model, custom_norm)
            if nargin == 1
                custom_norm = 0;
            end
            % if using a multi-channel model
            if custom_norm && model.num_channels > 1 
                % code run_preds for all fitting experiments
                imodel = tchModel(model.type, model.experiments, model.sessions);
                imodel = code_stim(imodel);
                imodel.normT = 1; imodel.normD = 1;
                imodel = pred_runs(imodel);
                for ss = 1:length(imodel.sessions)
                    % get predictors for all runs in this session
                    spreds = cell2mat(imodel.run_preds(:, ss));
                    npreds = size(spreds, 2);
                    % find max of predictors
                    if model.num_channels == 2
                        maxS = max(max(spreds(:, 1:npreds / 2)));
                        maxT = max(max(spreds(:, npreds / 2 + 1:npreds)));
                    elseif model.num_channels == 3
                        maxS = max(max(spreds(:, 1:npreds / 3)));
                        maxT = max(max(spreds(:, npreds / 3 + 1:2 * npreds / 3)));
                        maxD = max(max(spreds(:, 2 * npreds / 3 + 1:npreds)));
                        normDs{ss} = maxS / maxD;
                    end
                    % compute scalars to normalize max heights
                    normTs{ss} = maxS / maxT;
                end
                % set normalization constant to average across sessions
                model.normT = mean([normTs{:}]);
                if model.num_channels > 2
                    model.normD = mean([normDs{:}]);
                end
            end
        end
        
        % generate fMRI predictors for each session and run
        function model = pred_runs(model)
            switch model.type
                case 'glm'
                    model = pred_runs_glm(model);
                case 'htd'
                    model = pred_runs_htd(model);
                case 'balloon'
                    model = pred_runs_balloon(model);
                case 'cts-pow'
                    model = pred_runs_cts_pow(model);
                case 'cts-div'
                    model = pred_runs_cts_div(model);
                case 'dcts'
                    model = pred_runs_dcts(model);
                case '2ch-lin-lin'
                    model = pred_runs_2ch_lin_lin(model);
                case '2ch-lin-quad'
                    model = pred_runs_2ch_lin_quad(model);
                case '2ch-lin-rect'
                    model = pred_runs_2ch_lin_rect(model);
                case '2ch-pow-quad'
                    model = pred_runs_2ch_pow_quad(model);
                case '2ch-pow-rect'
                    model = pred_runs_2ch_pow_rect(model);
                case '2ch-div-quad'
                    model = pred_runs_2ch_div_quad(model);
                case '2ch-dcts-quad'
                    model = pred_runs_2ch_dcts_quad(model);
                case '2ch-opt'
                    model = pred_runs_2ch_opt(model);
                case '3ch-lin-quad'
                    model = pred_runs_3ch_lin_quad(model);
                case '3ch-lin-rect'
                    model = pred_runs_3ch_lin_rect(model);
                case '3ch-pow-quad'
                    model = pred_runs_3ch_pow_quad(model);
                case '3ch-pow-rect'
                    model = pred_runs_3ch_pow_rect(model);
                case '3ch-opt'
                    model = pred_runs_3ch_opt(model);
            end
        end
        
        % generate fMRI predictors for each session and trial type
        function model = pred_trials(model)
            switch model.type
                case 'glm'
                    model = pred_trials_glm(model);
                case 'htd'
                    model = pred_trials_htd(model);
                case 'balloon'
                    model = pred_trials_balloon(model);
                case 'cts-pow'
                    model = pred_trials_cts_pow(model);
                case 'cts-div'
                    model = pred_trials_cts_div(model);
                case 'dcts'
                    model = pred_trials_dcts(model);
                case '2ch-lin-lin'
                    model = pred_trials_2ch_lin_lin(model);
                case '2ch-lin-quad'
                    model = pred_trials_2ch_lin_quad(model);
                case '2ch-lin-rect'
                    model = pred_trials_2ch_lin_rect(model);
                case '2ch-pow-quad'
                    model = pred_trials_2ch_pow_quad(model);
                case '2ch-pow-rect'
                    model = pred_trials_2ch_pow_rect(model);
                case '2ch-div-quad'
                    model = pred_trials_2ch_div_quad(model);
                case '2ch-dcts-quad'
                    model = pred_trials_2ch_dcts_quad(model);
                case '2ch-opt'
                    model = pred_trials_2ch_opt(model);
                case '3ch-lin-quad'
                    model = pred_trials_3ch_lin_quad(model);
                case '3ch-lin-rect'
                    model = pred_trials_3ch_lin_rect(model);
                case '3ch-pow-quad'
                    model = pred_trials_3ch_pow_quad(model);
                case '3ch-pow-rect'
                    model = pred_trials_3ch_pow_rect(model);
                case '3ch-opt'
                    model = pred_trials_3ch_opt(model);
            end
        end
        
    end
    
end
