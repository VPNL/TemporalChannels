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
%      Dual-channel models:
%        '2ch-lin-lin'   -- linear susatined and linear transient
%        '2ch-lin-quad'  -- linear susatined and quadratic transient
%        '2ch-lin-rect'  -- linear susatined and rectified transient
%        '2ch-pow-quad'  -- sustained with CTS-p and quadratic transient
%        '2ch-pow-rect'  -- sustained with CTS-p and rectified transient
%        '2ch-div-quad'  -- sustained with CTS-n and quadratic transient
%        '2ch-exp-quad'  -- adapted sustained and quadratic transient
%        '2ch-exp-rect'  -- adapted sustained and rectified transient
%      Multi-channel models:
%        '3ch-lin-quad-exp'  -- linear sustained, quadratic transient, and persistent
%        '3ch-lin-rect-exp'  -- linear sustained, rectified transient, and persistent
%        '3ch-pow-quad-exp'  -- CTS-p on sustained, quadratic transient, and persistent
%        '3ch-pow-rect-exp'  -- CTS-p on sustained, rectified transient, and persistent
%        '3ch-exp-quad-exp'  -- adapted sustained, quadratic transient, and persistent
%        '3ch-exp-rect-exp'  -- adapted sustained, rectified transient, and persistent
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
%   model = tchModel('2ch-lin-quad', exps, roi.sessions);
%   model = code_stim(model); model = norm_model(model, 1);
%   model = pred_runs(model); model = pred_trials(model);
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
        types = {};     % list of all valid model types
        num_exps = [];  % number of experiments in model
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
        persist = {};   % persistent activity step function
        adapt = {};     % sustained activity without adaptation decay
        persist_act = {}; % persistent activity with decay
        adapt_act = {}; % sustained activity with adaptation decay
        normT = 20;     % transient channel normalization scalar
        normP = 2;      % persistent channel normalization scalar
        tr = 1;         % fMRI TR (s)
        gap_dur = 1/60; % forced gap between stimuli (s)
        pre_dur = 4;    % pre-stimulus baseline period (s)
        post_dur = 12;  % post-stimulus baseline duration (s)
        fs = 1000;      % stimulus/neural sampling rate (Hz)
        project_dir = fileparts(fileparts(which(mfilename, 'class')));
    end
    
    properties (Dependent, Hidden)
        num_runs      % numer of runs per experiment per session
        stimfiles     % paths to stimulus history files for each run
        num_channels  % number of channels in model (1-3 per condition)
        optimize_flag % whether model requires nonlinear optimization
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
            % check for invalid model type and initialize model parameters
            nsess = length(model.sessions); check_model_type(model.type);
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
        
        % Determine number of channels in model
        function optimize_flag = get.optimize_flag(model)
            optimize_flag = check_model_type(model.type);
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
            if model.num_channels > 2 && isfield(model.params, 'tau_pe')
                model = code_persist_decay(model);
            end
        end
                
        % code decay of activity in sustained channels
        function model = code_adapt_decay(model, type)
            if nargin < 2; type = 'exp'; end
            % get session and stimulus information
            nruns_max = size(model.onsets, 1); fs = model.fs;
            empty_cells = cellfun(@isempty, model.onsets);
            % find frame indices of persistent onsets and offsets for each cat
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
        
        % code decay of activity in persistent channels
        function model = code_persist_decay(model)
            % get session and stimulus information
            nruns_max = size(model.onsets, 1); rds = model.run_durs;
            fs = model.fs; empty_cells = cellfun(@isempty, model.onsets);
            ons = model.onsets; offs = model.offsets;
            % offs = cellfun(@(X) X + 1 / 30, offs, 'uni', false);
            ons = cellfun(@(X, Y) [X(2:end) Y], ons, rds, 'uni', false);
            persists = cellfun(@code_persist_act, model.stim, 'uni', false);
            persist_exps = repmat(model.irfs.persist_exp, nruns_max, 1);
            persist_exps(empty_cells) = {[]};
            persist_acts = cellfun(@(X, Y, Z, F) code_exp_decay(X, Y, Z, F, fs), ...
                persists, offs, ons, persist_exps, 'uni', false);
            model.persist = persists; model.persist_act = persist_acts;
        end
        
        % compute custom normalization parameters using run predictors
        function model = norm_model(model, custom_norm)
            if nargin == 1; custom_norm = 0; end
            % only normalize if using a multi-channel model
            if custom_norm == 1 && model.num_channels > 1
                % code run_preds for all fitting experiments
                imodel = tchModel(model.type, model.experiments, model.sessions);
                imodel.normT = 1; imodel.normP = 1; nch = imodel.num_channels;
                imodel = code_stim(imodel); imodel = pred_runs(imodel); 
                [normTs, normPs] = deal(zeros(1, length(imodel.sessions)));
                ncats = length(unique([imodel.cats{:}]));
                for ss = 1:length(imodel.sessions)
                    % get predictors for all runs in this session
                    spreds = cell2mat(imodel.run_preds(:, ss));
                    % find max of each set of channels
                    maxS = max(max(spreds(:, 0 * ncats + 1:1 * ncats)));
                    maxT = max(max(spreds(:, 1 * ncats + 1:2 * ncats)));
                    % compute scalars to normalize max heights to sustained
                    normTs(ss) = maxS / maxT;
                    if nch > 2
                        maxP = max(max(spreds(:, 2 * ncats + 1:3 * ncats)));
                        normPs(ss) = maxS / maxP;
                    end
                end
                % set normalization constant to average across sessions
                model.normT = mean(normTs);
                if nch > 2; model.normP = mean(normPs); end
            end
        end
        
        % generate fMRI predictors for each session and run
        function model = pred_runs(model)
            switch model.type
                case '1ch-lin'
                    model = pred_runs_1ch_lin(model);
                case '1ch-exp'
                    model = pred_runs_1ch_exp(model);
                case '1ch-rect'
                    model = pred_runs_1ch_rect(model);
                case '1ch-quad'
                    model = pred_runs_1ch_quad(model);
                case '1ch-cquad'
                    model = pred_runs_1ch_cquad(model);
                case '1ch-sig'
                    model = pred_runs_1ch_sig(model);
                case '2ch-lin-rect'
                    model = pred_runs_2ch_lin_rect(model);
                case '2ch-exp-rect'
                    model = pred_runs_2ch_exp_rect(model);
                case '2ch-lin-quad'
                    model = pred_runs_2ch_lin_quad(model);
                case '2ch-exp-quad'
                    model = pred_runs_2ch_exp_quad(model);
                case '2ch-lin-cquad'
                    model = pred_runs_2ch_lin_cquad(model);
                case '2ch-exp-cquad'
                    model = pred_runs_2ch_exp_cquad(model);
                case '2ch-lin-sig'
                    model = pred_runs_2ch_lin_sig(model);
                case '2ch-exp-sig'
                    model = pred_runs_2ch_exp_sig(model);
                case '2ch-exp-crect'
                    model = pred_runs_2ch_exp_crect(model);
                case '3ch-exp-crect-crect'
                    model = pred_runs_3ch_exp_crect_crect(model);
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
                case '2ch-lin-lin'
                    model = pred_runs_2ch_lin_lin(model);
                case '2ch-pow-quad'
                    model = pred_runs_2ch_pow_quad(model);
                case '2ch-pow-rect'
                    model = pred_runs_2ch_pow_rect(model);
                case '2ch-div-quad'
                    model = pred_runs_2ch_div_quad(model);
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
                case '3ch-lin-quad-exp-opt'
                    model = pred_runs_3ch_lin_quad_exp_opt(model);
                case '3ch-lin-rect-exp-opt'
                    model = pred_runs_3ch_lin_rect_exp_opt(model);
                case '3ch-exp-quad-exp-opt'
                    model = pred_runs_3ch_exp_quad_exp_opt(model);
                case '2ch-lin-crect'
                    model = pred_runs_2ch_lin_crect(model);
                case '2ch-exp-dquad'
                    model = pred_runs_2ch_exp_dquad(model);
                case '2ch-lin-dquad'
                    model = pred_runs_2ch_lin_dquad(model);
                case '3ch-exp-cquad-rect'
                    model = pred_runs_3ch_exp_cquad_rect(model);
            end
        end
        
        % generate fMRI predictors for each session and trial type
        function model = pred_trials(model)
            switch model.type
                case '1ch-lin'
                    model = pred_trials_1ch_lin(model);
                case '1ch-exp'
                    model = pred_trials_1ch_exp(model);
                case '1ch-rect'
                    model = pred_trials_1ch_rect(model);
                case '1ch-quad'
                    model = pred_trials_1ch_quad(model);
                case '1ch-cquad'
                    model = pred_trials_1ch_cquad(model);
                case '1ch-sig'
                    model = pred_trials_1ch_sig(model);
                case '2ch-lin-rect'
                    model = pred_trials_2ch_lin_rect(model);
                case '2ch-exp-rect'
                    model = pred_trials_2ch_exp_rect(model);
                case '2ch-lin-quad'
                    model = pred_trials_2ch_lin_quad(model);
                case '2ch-exp-quad'
                    model = pred_trials_2ch_exp_quad(model);
                case '2ch-lin-cquad'
                    model = pred_trials_2ch_lin_cquad(model);
                case '2ch-exp-cquad'
                    model = pred_trials_2ch_exp_cquad(model);
                case '2ch-lin-sig'
                    model = pred_trials_2ch_lin_sig(model);
                case '2ch-exp-sig'
                    model = pred_trials_2ch_exp_sig(model);
                case '2ch-exp-crect'
                    model = pred_trials_2ch_exp_crect(model);
                case '3ch-exp-crect-crect'
                    model = pred_trials_3ch_exp_crect_crect(model);
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
                case '2ch-lin-lin'
                    model = pred_trials_2ch_lin_lin(model);
                case '2ch-pow-quad'
                    model = pred_trials_2ch_pow_quad(model);
                case '2ch-pow-rect'
                    model = pred_trials_2ch_pow_rect(model);
                case '2ch-div-quad'
                    model = pred_trials_2ch_div_quad(model);
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
                case '3ch-lin-quad-exp-opt'
                    model = pred_trials_3ch_lin_quad_exp_opt(model);
                case '3ch-lin-rect-exp-opt'
                    model = pred_trials_3ch_lin_rect_exp_opt(model);
                case '3ch-exp-quad-exp-opt'
                    model = pred_trials_3ch_exp_quad_exp_opt(model);
                case '2ch-exp-cquad-opt'
                    model = pred_trials_2ch_exp_cquad_opt(model);
                case '3ch-exp-quad-crect-opt'
                    model = pred_trials_3ch_exp_quad_crect_opt(model);
                case '2ch-exp-dquad'
                    model = pred_trials_2ch_exp_dquad(model);
                case '2ch-lin-dquad'
                    model = pred_trials_2ch_lin_dquad(model);
                case '3ch-exp-cquad-rect'
                    model = pred_trials_3ch_exp_cquad_rect(model);
            end
        end
        
    end
    
end
