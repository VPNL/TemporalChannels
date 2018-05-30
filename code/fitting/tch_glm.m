function model = tch_glm(Y, X)
% Apply a general linear model (GLM) to data in Y using design matrix in X
%
% INPUTS
%   Y: TR by voxel matrix of fMRI time series data
%   X: model design matrix
%
% OUTPUT FIELDS
%   betas: fitted linear model weights
%   residual: error of model prediction for each time point
%   stdevs: estimated standard deviation for each beta weight
%   dof: degrees of freedom of the fitting
%
% Adapted from vistaoft (http://github.com/vistalab/vistasoft/)
% AS 2/2017

% check inputs
if nargin < 2
    error('not enough input arguments');
end
if size(X,1) ~= size(Y, 1)
    error('rows in data (Y) and design matrix (X) do not match');
else
    X = double(X);
    Y = double(Y);
end

% initialize model struct
model = struct('design_mat', [], 'dof', [], 'betas', [], 'residual', [], 'var_covar', []);

% compute degrees of freedom and store design matrix
model.design_mat = X;
model.dof = size(Y, 1) - rank(X);

% get number of TRs, voxels, and predictors
[num_trs, num_vox] = size(Y); num_preds = size(X, 2);

% estimate beta weights for each predictor in the design matrix
betas = X \ Y; model.betas = reshape(betas, [1 num_preds num_vox]);

% compute residual error at each time point
model.residual = Y - X * betas;

% compute residual variance and variance-covariance matrix
model.resid_var = sum(model.residual .^ 2) / model.dof;
vin = inv(X' * eye(num_trs) * X); model.var_covar = X' * X;

% compute standard deviations and standard errors of betas
stdevs = sqrt(diag(vin) .* diag(model.var_covar) * model.resid_var);
sems = stdevs ./ sqrt(repmat(round(sum(X))', [1 num_vox]));
model.stdevs = reshape(stdevs, [1 num_preds num_vox]);
model.sems = reshape(sems, [1 num_preds num_vox]);

end
