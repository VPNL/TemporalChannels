% tchModel: Code for modeling fMRI responses to time-varying visual stimuli.
%
% CONSTRUCTOR INPUTS
%   1) type: which model to use
%      Hemodynamic models:
%        '1ch-lin'     -- general linear model for fMRI (Boynton 1996)
%        '2ch-lin-htd' -- hemodynamic temporal derivative (HTD; Henson 2002)
%        '1ch-balloon' -- nonlinear balloon model (Buxton, 1998)
%      Single-channel models:
%        '1ch-pow'  -- CTS with power law (CTS-p; Zhou 2017)
%        '1ch-div'  -- CTS with divisive normalization (CTS-n; Zhou 2017)
%        '1ch-dcts' -- dynamic CTS (dCTS; Zhou 2017)
%        '1ch-exp'  -- stimulus-specific adaptation model (exponential)
%        '1ch-cexp' -- compressed adaptation model (CTS-p + exponential)
%      Dual-channel models:
%        '2ch-lin-lin'   -- linear susatined and linear transient
%        '2ch-lin-quad'  -- linear susatined and quadratic transient
%        '2ch-lin-rect'  -- linear susatined and rectified transient
%        '2ch-pow-quad'  -- sustained with CTS-p and quadratic transient
%        '2ch-pow-rect'  -- sustained with CTS-p and rectified transient
%        '2ch-div-quad'  -- sustained with CTS-n and quadratic transient
%        '2ch-exp-quad'  -- adapted sustained and quadratic transient
%        '2ch-exp-rect'  -- adapted sustained and rectified transient
%        '2ch-cexp-quad' -- compressed/adapted sustained and quadratic transient
%        '2ch-cexp-rect' -- compressed/adapted sustained and rectified transient
%      Multi-channel models:
%        '3ch-lin-quad-exp'  -- linear sustained, quadratic transient, and delay
%        '3ch-lin-rect-exp'  -- linear sustained, rectified transient, and delay
%        '3ch-pow-quad-exp'  -- CTS-p on sustained, quadratic transient, and delay
%        '3ch-pow-rect-exp'  -- CTS-p on sustained, rectified transient, and delay
%        '3ch-exp-quad-exp'  -- adapted sustained, quadratic transient, and delay
%        '3ch-exp-rect-exp'  -- adapted sustained, rectified transient, and delay
%        '3ch-cexp-quad-exp' -- compressed/adapted sustained, quadratic transient, and delay
%        '3ch-cexp-rect-exp' -- compressed/adapted sustained, rectified transient, and delay
%   2) exps: array of experiments for fitting model (e.g., {'Exp1' 'Exp2'})
%   3) sessions: array of names or paths to session data directories
%
% METHODS
%   code_stim   -- compiles onset and offset times for each stimulus
%   norm_model  -- customize predictor normalization parameters (optional)
%   pred_runs   -- generates design matrix for each run of data
%   pred_trials -- generates response predictors for each trial type
%
% Example model generation steps:
%   exps = {'Exp1' 'Exp2' 'Exp3'};
%   roi = tchROI('V1', exps);
%   roi = select_sessions(roi);
%   model = tchModel('2ch', exps, roi.sessions);
%   model = norm_model(model, 1);
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
        params = [];      % structure of model parameters
        irfs = [];        % structure of model impulse response functions
        run_preds = [];   % fMRI predictors for each run
        trial_preds = []; % fMRI predictors for each trial type
    end
    
    properties (Hidden)
        num_exps        % number of experiments in model
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
        delay = {};     % delay activity step function
        adapt = {};     % sustained activity without adaptation decay
        delay_act = {}; % delay activity with decay
        adapt_act = {}; % sustained activity with adaptation decay
        normT = 10;     % transient channel normalization scalar
        normD = 2;      % delay channel normalization scalar
    end
    
    properties (Constant, Hidden)
        % path to project directory
        project_dir = fileparts(fileparts(which(mfilename, 'class')));
        % descriptors for each model implemented
        types = {'1ch-lin' '2ch-lin-htd' '1ch-balloon' ...
            '1ch-pow' '1ch-div' '1ch-dcts' '1ch-exp' '1ch-cexp' ...
            '2ch-lin-lin' '2ch-lin-quad' '2ch-lin-rect' ...
            '2ch-pow-quad' '2ch-pow-rect' '2ch-div-quad' ...
            '2ch-exp-quad' '2ch-exp-rect' ...
            '2ch-cexp-quad' '2ch-cexp-rect' ...
            '3ch-lin-quad-exp' '3ch-lin-rect-exp' ...
            '3ch-pow-quad-exp' '3ch-pow-rect-exp' ...
            '3ch-exp-quad-exp' '3ch-exp-rect-exp' ...
            '3ch-cexp-quad-exp' '3ch-cexp-rect-exp' ...
            '2ch-lin-quad-opt' '3ch-lin-rect-exp-opt'};
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
        
        % Find paths to stimulus timing parameter files
        function stimfiles = get.stimfiles(model)
            nsess = length(model.sessions); nruns = model.num_runs;
            stimfiles = cell(max(sum(nruns, 1)), nsess);
            for ss = 1:nsess
                rcnt = 0;
                for ee = 1:model.num_exps
                    ecnt = 1;
                    fstem = fullfile(model.sessions{ss}, 'Stimuli');
                    fname = [model.experiments{ee} '_Run'];
                    for ff = 1:nruns(ee, ss)
                        rcnt = rcnt + 1; sn = [num2str(ecnt) '.txt'];
                        stimfiles{rcnt, ss} = fullfile(fstem, [fname sn]);
                        ecnt = ecnt + 1;
                    end
                end
            end
        end
        
        % Determine number of channels in model
        function num_channels = get.num_channels(model)
            switch model.type(1:3)
                case '1ch'
                    num_channels = 1;
                case '2ch'
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
            sfiles = model.stimfiles; [nruns_max, nsess] = size(sfiles);
            empty_cells = cellfun(@isempty, sfiles);
            fs = model.fs; gd = model.gap_dur;
            % get image onsets/offsets and condition orders from stimfiles
            [on, off, c, im, ton, toff, tc, rd, cl] = cellfun(@tch_stimfile, ...
                sfiles, 'uni', false); cat_list = unique([im{:}]);
            % store stimulus parameters in tchModel object
            model.onsets = on;     % image onset times (s)
            model.offsets = off;   % image offset times (s)
            model.conds = c;       % condition labels per image
            model.cats = im;       % category labels per image
            model.tonsets = ton;   % trial onset times (s)
            model.toffsets = toff; % trial offset times (s)
            model.tconds = tc;     % condition labels per trial
            model.run_durs = rd;   % duration of each run (s)
            % get list of trial conditions from stimfile header
            rcnt = 1;
            for ee = 1:model.num_exps
                model.cond_list{ee} = cl{rcnt, 1};
                rcnt = rcnt + model.num_runs(ee, 1);
            end
            % find frame indices of image onsets and offsets for each cat
            stims = cellfun(@(X) zeros(X * fs, length(cat_list)), ...
                rd, 'uni', false); stims(empty_cells) = {[]};
            for cc = 1:length(cat_list)
                cat_idxs = cellfun(@(X) find(strcmp(cat_list(cc), X)), ...
                    im, 'uni', false);
                on_idxs = cellfun(@(X, Y) round(fs * X(Y)) + 1, ...
                    on, cat_idxs, 'uni', false); on_idxs(empty_cells) = {1};
                off_idxs = cellfun(@(X, Y) round(fs * (X(Y) - gd / 2)), ...
                    off, cat_idxs, 'uni', false); off_idxs(empty_cells) = {1};
                % find frame indices of gap offset times
                goff_idxs = cellfun(@(X, Y) round(fs * (X(Y) + gd / 2)), ...
                    off, cat_idxs, 'uni', false); goff_idxs(empty_cells) = {1};
                % compile all frame indices during stimuli and gaps
                stim_idxs = cellfun(@code_stim_idxs, ...
                    on_idxs, off_idxs, 'uni', false);
                gap_idxs = cellfun(@code_stim_idxs, ...
                    off_idxs, goff_idxs, 'uni', false);
                % code stimulus as a step function with gaps at offsets
                cc_x = repmat({cc}, nruns_max, nsess);
                stims = cellfun(@code_stim_vec, stims, stim_idxs, ...
                    cc_x, repmat({1}, nruns_max, nsess), 'uni', false);
                stims = cellfun(@code_stim_vec, stims, gap_idxs, ...
                    cc_x, repmat({0}, nruns_max, nsess), 'uni', false);
            end
            stims(empty_cells) = {[]}; model.stim = stims;
            if ~isempty(strfind(model.type, 'ch-exp'))
                model = code_adapt_decay(model, 'exp');
            elseif ~isempty(strfind(model.type, 'ch-cexp'))
                model = code_adapt_decay(model, 'cexp');
            end
            if model.num_channels > 2 && ~isempty(strfind(model.type, '-exp'))
                model = code_delay_decay(model);
            end
        end
                
        % code decay of activity in sustained channels
        function model = code_adapt_decay(model, type)
            if nargin < 2; type = 'exp'; end
            % get session and stimulus information
            nruns_max = size(model.onsets, 1); fs = model.fs;
            empty_cells = cellfun(@isempty, model.onsets);
            % find frame indices of delay onsets and offsets for each cat
            adapt_exps = repmat(model.irfs.adapt_exp, nruns_max, 1);
            nrfS = repmat(model.irfs.nrfS, nruns_max, 1);
            adapt_exps(empty_cells) = {[]}; nrfS(empty_cells) = {[]};
            adapts = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
                model.stim, nrfS, 'uni', false);
            if strcmp(type, 'cexp')
                ce = repmat(model.params.epsilon, nruns_max, 1);
                adapts = cellfun(@(X, Y) X .^ Y, adapts, ce, 'uni', false);
            end
            adapt_acts = cellfun(@(X, Y, Z, F) code_exp_decay(X, Y, Z, F, fs), ...
                adapts, model.onsets, model.offsets, adapt_exps, 'uni', false);
            model.adapt = adapts; model.adapt_act = adapt_acts;
        end
        
        % code decay of activity in delay channels
        function model = code_delay_decay(model)
            % get session and stimulus information
            nruns_max = size(model.onsets, 1); rds = model.run_durs;
            fs = model.fs; empty_cells = cellfun(@isempty, model.onsets);
            ons = model.onsets; offs = model.offsets;
            ons = cellfun(@(X, Y) [X(2:end) Y], ons, rds, 'uni', false);
            delays = cellfun(@code_delay_act, model.stim, 'uni', false);
            delay_exps = repmat(model.irfs.delay_exp, nruns_max, 1);
            delay_exps(empty_cells) = {[]};
            delay_acts = cellfun(@(X, Y, Z, F) code_exp_decay(X, Y, Z, F, fs), ...
                delays, offs, ons, delay_exps, 'uni', false);
            model.delay = delays; model.delay_act = delay_acts;
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
                [normTs, normDs] = deal(cell(1, size(imodel.run_preds, 2)));
                cat_list = unique([imodel.cats{:}]);
                npreds = length(cat_list) * imodel.num_channels;
                for ss = 1:length(imodel.sessions)
                    % get predictors for all runs in this session
                    spreds = cell2mat(imodel.run_preds(:, ss));
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
                case '1ch-lin'
                    model = pred_runs_1ch_lin(model);
                case '2ch-lin-htd'
                    model = pred_runs_2ch_lin_htd(model);
                case '1ch-balloon'
                    model = pred_runs_1ch_balloon(model);
                case '1ch-pow'
                    model = pred_runs_1ch_pow(model);
                case '1ch-div'
                    model = pred_runs_1ch_div(model);
                case '1ch-dcts'
                    model = pred_runs_1ch_dcts(model);
                case '1ch-exp'
                    model = pred_runs_1ch_exp(model);
                case '1ch-cexp'
                    model = pred_runs_1ch_cexp(model);
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
                case '2ch-exp-quad'
                    model = pred_runs_2ch_exp_quad(model);
                case '2ch-exp-rect'
                    model = pred_runs_2ch_exp_rect(model);
                case '2ch-cexp-quad'
                    model = pred_runs_2ch_cexp_quad(model);
                case '2ch-cexp-rect'
                    model = pred_runs_2ch_cexp_rect(model);
                case '3ch-lin-quad-exp'
                    model = pred_runs_3ch_lin_quad_exp(model);
                case '3ch-lin-rect-exp'
                    model = pred_runs_3ch_lin_rect_exp(model);
                case '3ch-pow-quad-exp'
                    model = pred_runs_3ch_pow_quad_exp(model);
                case '3ch-pow-rect-exp'
                    model = pred_runs_3ch_pow_rect_exp(model);
                case '3ch-exp-quad-exp'
                    model = pred_runs_3ch_exp_quad_exp(model);
                case '3ch-exp-rect-exp'
                    model = pred_runs_3ch_exp_rect_exp(model);
                case '3ch-cexp-quad-exp'
                    model = pred_runs_3ch_cexp_quad_exp(model);
                case '3ch-cexp-rect-exp'
                    model = pred_runs_3ch_cexp_rect_exp(model);
                case '2ch-lin-quad-opt'
                    model = pred_runs_2ch_lin_quad_opt(model);
                case '3ch-lin-rect-exp-opt'
                    model = pred_runs_3ch_lin_rect_exp_opt(model);
            end
        end
        
        % generate fMRI predictors for each session and trial type
        function model = pred_trials(model)
            switch model.type
                case '1ch-lin'
                    model = pred_trials_1ch_lin(model);
                case '2ch-lin-htd'
                    model = pred_trials_2ch_lin_htd(model);
                case '1ch-balloon'
                    model = pred_trials_1ch_balloon(model);
                case '1ch-pow'
                    model = pred_trials_1ch_pow(model);
                case '1ch-div'
                    model = pred_trials_1ch_div(model);
                case '1ch-dcts'
                    model = pred_trials_1ch_dcts(model);
                case '1ch-exp'
                    model = pred_trials_1ch_exp(model);
                case '1ch-cexp'
                    model = pred_trials_1ch_cexp(model);
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
                case '2ch-exp-quad'
                    model = pred_trials_2ch_exp_quad(model);
                case '2ch-exp-rect'
                    model = pred_trials_2ch_exp_rect(model);
                case '2ch-cexp-quad'
                    model = pred_trials_2ch_cexp_quad(model);
                case '2ch-cexp-rect'
                    model = pred_trials_2ch_cexp_rect(model);
                case '3ch-lin-quad-exp'
                    model = pred_trials_3ch_lin_quad_exp(model);
                case '3ch-lin-rect-exp'
                    model = pred_trials_3ch_lin_rect_exp(model);
                case '3ch-pow-quad-exp'
                    model = pred_trials_3ch_pow_quad_exp(model);
                case '3ch-pow-rect-exp'
                    model = pred_trials_3ch_pow_rect_exp(model);
                case '3ch-exp-quad-exp'
                    model = pred_trials_3ch_exp_quad_exp(model);
                case '3ch-exp-rect-exp'
                    model = pred_trials_3ch_exp_rect_exp(model);
                case '3ch-cexp-quad-exp'
                    model = pred_trials_3ch_cexp_quad_exp(model);
                case '3ch-cexp-rect-exp'
                    model = pred_trials_3ch_cexp_rect_exp(model);
                case '2ch-lin-quad-opt'
                    model = pred_trials_2ch_lin_quad_opt(model);
                case '3ch-lin-rect-exp-opt'
                    model = pred_trials_3ch_lin_rect_exp_opt(model);
            end
        end
        
    end
    
end
