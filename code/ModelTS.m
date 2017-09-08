
% ModelTS: Code for modeling fMRI responses to time-varying visual stimuli.
% 
% CONSTRUCTOR INPUTS
%   1) type: which model to use
%      Linear models:
%        'standard' -- standard general linear model
%        'htd'      -- hemodynamic temporal derivative (HTD; Henson 2002)
%        '2ch-lin'  -- linear version of 2 temporal-channel model
%      Nonlinear single-channel models:
%        'cts-pow'  -- CTS with power law (Zhou et al., 2017)
%        'cts-div'  -- CTS with divisive normalization (Zhou et al., 2017)
%        'dcts'     -- dynamic CTS (Zhou et al., 2017)
%        'balloon'  -- hemodynamic balloon model (Buxton et al., 1998)
%      Multi-channel models:
%        '2ch'      -- 2 temporal-channel model (Stigliani et al., 2017)
%        '2ch-rect' -- 2 temporal-channel model without offset responses
%        '2ch-div'  -- 2 temporal-channel model with CTS-div on sustained
%        '2ch-pow'  -- 2 temporal-channel model with CTS-pow on sustained
%        '2ch-dcts' -- 2 temporal-channel model with dCTS on sustained
%        '2ch-opt'  -- 2 temporal-channel model with dCTS on sustained
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
%   model = ModelTS('2ch', exps, ROI('V1', exps).sessions);
%   model = code_stim(model);
%   model = pred_runs(model);
%   model = pred_trials(model);
% 
% AS 2/2017

classdef ModelTS
    
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
    end

    properties (Constant, Hidden)
        % path to project directory
        project_dir = fileparts(fileparts(which(mfilename, 'class')));
        % descriptors for each model implemented
        types = {'standard' 'htd' '2ch-lin' ...
            'cts-pow' 'cts-div' 'dcts' 'balloon'...
            '2ch' '2ch-rect' '2ch-pow' '2ch-div' '2ch-dcts' '2ch-opt'};
        % experimental parameters
        tr = 1;         % fMRI TR (s)
        gap_dur = 1/60; % forced gap between stimuli (s)
        pre_dur = 4;    % pre-stimulus baseline period (s)
        post_dur = 12;  % post-stimulus baseline duration (s)
        fs = 1000;      % stimulus/neural sampling rate (Hz)
    end
    
    properties (Dependent, Hidden)
        num_runs  % numer of runs per experiment per session
        stimfiles % paths to stimulus history files for each run
    end
    
    methods
        
        % ModelTS class constructor
        function model = ModelTS(type, exps, sessions)
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
            [params, irfs] = init_params(type, nsess, model.fs);
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
                
        % code onset, offset, and category of each stimulus in experiments
        function model = code_stim(model)
            % get session and run parameters
            stimfiles = model.stimfiles; [nruns_max, nsess] = size(stimfiles);
            fs = model.fs; gd = model.gap_dur;
            % get image onsets/offsets and condition orders from stimfiles
            [on, off, c, ims, ton, toff, tc, rd, cl] = cellfun(@stimfileTS, stimfiles, 'uni', false);
            % store stimulus parameters in modelTS object
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
            if ~sum(strcmp(model.type, {'standard' 'cts-pow' 'cts-div' 'dcts'}))
                if custom_norm == 1
                    % code run_preds for all experiments in example subject
                    roi_list = dir(fullfile(model.sessions{1}, 'ROIs'));
                    roi = tcROI(roi_list(3).name, model.experiments);
                    imodel = modelTS(model.type, model.experiments, roi.sessions);
                    imodel.normT = 1;
                    imodel = code_stim(imodel);
                    imodel = pred_runs(imodel);
                    for ss = 1:length(imodel.sessions)
                        % get predictors for all runs in this session
                        iconcat = cell2mat(imodel.run_preds(:, ss));
                        npreds = size(iconcat, 2);
                        % find max of predictors
                        maxS = max(max(iconcat(:, 1:npreds / 2)));
                        maxT = max(max(iconcat(:, npreds / 2 + 1:npreds)));
                        % compute scalars to normalize max heights
                        normTs{ss} = maxS / maxT;
                    end
                    model.normT = mean([normTs{:}]);
                end
            end
        end
        
        % generate fMRI predictors for each session and run
        function model = pred_runs(model)
            switch model.type
                case 'standard'
                    model = pred_runs_standard(model);
                case 'htd'
                    model = pred_runs_htd(model);
                case '2ch-lin'
                    model = pred_runs_2ch_lin(model);
                case 'cts-pow'
                    model = pred_runs_cts_pow(model);
                case 'cts-div'
                    model = pred_runs_cts_div(model);
                case 'dcts'
                    model = pred_runs_dcts(model);
                case 'balloon'
                    model = pred_runs_balloon(model);
                case '2ch'
                    model = pred_runs_2ch(model);
                case '2ch-rect'
                    model = pred_runs_2ch_rect(model);
                case '2ch-pow'
                    model = pred_runs_2ch_pow(model);
                case '2ch-div'
                    model = pred_runs_2ch_div(model);
                case '2ch-dcts'
                    model = pred_runs_2ch_dcts(model);
                case '2ch-opt'
                    model = pred_runs_2ch_opt(model);
            end
        end
        
        % generate fMRI predictors for each session and trial type
        function model = pred_trials(model)
            switch model.type
                case 'standard'
                    model = pred_trials_standard(model);
                case 'htd'
                    model = pred_trials_htd(model);
                case 'balloon'
                    model = pred_trials_balloon(model);
                case '2ch-lin'
                    model = pred_trials_2ch_lin(model);
                case 'cts-pow'
                    model = pred_trials_cts_pow(model);
                case 'cts-div'
                    model = pred_trials_cts_div(model);
                case 'dcts'
                    model = pred_trials_dcts(model);
                case '2ch-pow'
                    model = pred_trials_2ch_pow(model);
                case '2ch'
                    model = pred_trials_2ch(model);
                case '2ch-rect'
                    model = pred_trials_2ch_rect(model);
                case '2ch-div'
                    model = pred_trials_2ch_div(model);
                case '2ch-dcts'
                    model = pred_trials_2ch_dcts(model);
                case '2ch-opt'
                    model = pred_trials_2ch_opt(model);
            end
        end
        
    end
    
end
