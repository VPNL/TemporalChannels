# TemporalChannels
Code for modeling fMRI responses to time-varying stimuli using a temporal channels framework.
* * *
*Technical notes:*
- The code in this repository is compatible with MATLAB R2015a and later versions. 
- An example dataset is available here: http://vpnl.stanford.edu/TemporalChannels/TemporalChannels_data.tar.gz
* * *
*Table of contents:*
1. [Instructions](#instructions)
    1. [Oranizing the data directory](#organizing-the-data-directory)
    2. [Generating stimulus timing parameter files](#generating-stimulus-timing-parameter-files)
    3. [Fitting and validating a model in a region of interest](#fitting-and-validating-a-model-in-a-region-of-interest)
        1. [Inputs](#inputs)
        2. [Outputs](#outputs)
2. [Example code](#example-code)
    1. [Using the model_roi function](#using-the-model_roi-function)
    2. [Using ModelTS class methods](#using-modelts-class-methods)
    3. [Using ROI class methods](#using-roi-class-methods)
* * *
## Instructions

### Organizing the data directory

To work with the example dataset, download the [archive](http://vpnl.stanford.edu/TemporalChannels/TemporalChannels_data.tar.gz) and extract all session directories in `~/TemporalChannels/data/` of your local branch of the repository.

To work with your own dataset, create a separate directory for each experimental session in  `~/TemporalChannels/data/`. Each session directory should contain the following subdirectories: 

1. *ROIs* – stores `*.mat` files sorted by experiment and labeled by run number (e.g., `~/TemporalChannels/data/*/ROIs/V1/Exp1/Run1.mat`) that contain the raw fMRI time series of each voxel in a predefined region of interest. Each run of data is stored as a matrix `tSeries` with each row indexing a TR (sorted in ascending order) and each column indexing a voxel (sorted in arbitrary order). 

2. *Stimuli* – stores the stimulus timing parameters for each run as `*.txt` files labeled by experiment and run number (e.g., `~/TemporalChannels/data/*/Stimuli/Exp1_Run1.txt`). 

### Generating stimulus timing parameter files

An example stimulus timing parameter file is available here: https://github.com/VPNL/TemporalChannels/tree/master/examples/Exp1_Run1.txt

#### Header information

Line 1 contains the name of the experiment followed by a list of all trial types (e.g., `Exp1: TrialA, TrialB, TrialC`). 
  - “Trials” are structured periods of stimulus presentation within a run, which are separated by prolonged baseline periods (ideally ≥12 s) that mark transitions between different conditions that repeat throughout the experiment and across sessions. 
  - The overall duration and number of stimuli in a trial can vary across trial types, but the timing of stimulus presentations and sequence of stimulus categories must be identical across trials of the same type. 
  - Note that trial segmentation does not affect the fit of the model but will affect noise ceiling estimation and visualization. 

Line 2 contains the the total duration of the run in seconds (e.g., `Run duration (s): 300`). 

#### Stimulus information

Lines 5 and below contain the following information about each stimulus delimited by spaces: 

1. *Trial* — trial number that increments from 1 to the total number of trials in the run. 

2. *Condition* — label indicating the type of trial to which the stimulus belongs (same as in header line). 

3. *Onset* — stimulus onset time in seconds relative to the beginning of the run. 

4. *Duration* — stimulus presentation duration in seconds. 

5. *Filename* — descriptor composed of a category label and stimulus-specific identifier delimited by a dash (e.g., `face-1.jpg`). Note that a separate array of predictors is coded for each unique category label. 

### Fitting and validating a model in a region of interest

The `model_roi.m` wrapper function is used to fit and validate various temporal models using the procedure described below: 

1. An object of the class ROI is generated that loads, preprocesses, and stores the fMRI time series of voxels in a region of interest for a group of experimental sessions. 

2. A object of the class ModelTS is generated that creates predictors for each run of data being used to fit the model for each session in the ROI object. 

3. Response amplitudes (β weights) for each predictor in the model are estimated using a general linear model (GLM).

4. Fitted β weights are used to predict responses to each trial type in the data being used to fit the model. 

5. Model validation is performed using these fitted β weights to predict responses in independent data. 

#### Inputs

Fitting a model using `model_roi.m` requires passing at least three input arguments: 

1. *name* — name of a region of interest (e.g., `'V1'`) in the session ROI directories (i.e., `~/TemporalChannels/data/*/ROIs/`). 

2. *type* — label indicating which model to use for predicting responses. 
    1. `‘standard’` — standard linear systems approach for analyzing fMRI data
    2. `‘htd’` — hemodynamic temporal derivative model (HTD) proposed by Henson et al. (2002)
    3. `‘cts’` — compressive temporal summation model (CTS) proposed by Zhou et al. (2017)
    4. `‘2ch’` — 2 temporal-channel model proposed by Stigliani et al. (2017)

3. *fit_exps* — which experiment/s to use for fitting the model (e.g., `{'Exp1' 'Exp2'}`) with experiment names matching the stems of filenames in the session Stimuli directories (i.e., `~/TemporalChannels/data/*/Stimuli/`). 

4. *val_exps* - optional arguement specifying which experiment/s to use for validating the model (e.g., `'Exp 3'`).

#### Outputs

After fitting the model in each session individually, the function plots the average measured vs. predicted responses for each trial type and returns two outputs: 

1. *roi* — object of the class `ROI` that stores fMRI data and model predictions for the region of interest.
    1. `roi(1).runs` — stores run time series for *fit_exps* (`roi(2).runs` for *val_exps* if applicable)
    2. `roi(1).trials` —  stores trial time series for *fit_exps* (`roi(2).trials` for *val_exps* if applicable)
    3. `roi(1).model` — stores model fit and quantification of performance
    4. See properties in `ROI.m` class file for more details
2. *model* — object of the class `ModelTS` that stores run and trial predictors for each session.
    1. `model(1).run_preds` — stores run predictors for *fit_exps* (`model(2).run_preds` for *val_exps* if applicable)
    2. `model(1).trial_preds` — stores trial predictors for *fit_exps* (`model(2).trial_preds` for *val_exps* if applicable)
    3. `model(1).irfs` — stores impulse response functions
    4. `model(1).params` — stores model parameters
    5. See properties in `ModelTS.m` class file for more details

By default, outputs are saved in `~/TemporalChannels/results/`. 

## Example code

### Using the model_roi function

Example of fitting the 2 temporal-channel model using V1 data from Exp1 & Exp2:

    [roi, model] = model_roi(‘V1', '2ch', {‘Exp1’ 'Exp2'});

Example of fitting the standard model using V1 data from Exp1 & Exp2 and validating on V1 data from Exp3:

    [roi, model] = model_roi(‘V1', 'standard', {‘Exp1’ 'Exp2'}, 'Exp3');

### Using ModelTS class methods

The `ModelTS.m` class file defines a class of objects that store predictors for each run and trial type using the methods listed in the example below.

Example of creating predictors using a ModelTS object (`model`): 

    model_type = '2ch';                              % specify the type of model to use
    fit_exps = {'Exp1' 'Exp2'};                      % list of experiments for fitting
    sessions = roi.sessions;                         % list of sessions to model (see below)
    model = ModelTS(model_type, fit_exps, sessions); % setup ModelTS object
    model = code_stim(model);                        % code the timing of stimuli 
    model = pred_runs(model);                        % generate run predictors
    model = pred_trials(model);                      % generate trial predictors

### Using ROI class methods

The `ROI.m` class file defines a class of objects that store and operate on fMRI time series for a given region of interest using the methods listed in the example below.

Example of fitting a model to a ROI object (`roi`): 

    roi_name = 'V1';                                 % name of region to model
    fit_exps = {'Exp1' 'Exp2'};                      % list of experiments for fitting
    roi = ROI(roi_name, fit_exps);                   % setup ROI object
    roi = tc_runs(roi);                              % preprocess run time series
    roi = tc_trials(roi, model);                     % compile responses to each trial
    roi = tc_fit(roi, model);                        % fit model for each session
    roi = tc_pred(roi, model);                       % predict responses for each trial type

Plot mean response to each trial type across all sessions:

    fig = plot_exps(roi);

Plot response vs. model prediction for each run in individual sesssions:

    fig = plot_model(roi);
