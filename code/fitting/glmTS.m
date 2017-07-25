function model = glmTS(Y,X)
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
%   DoF: degrees of freedom of the fitting
%
% adapted from vistaoft (http://github.com/vistalab/vistasoft/)
% AS 2/2017

% check inputs
if nargin < 2
    error('not enough input arguments');
end
if size(X,1) ~= size(Y,1)
    error('rows in data (Y) and design matrix (X) do not match');
else
    X = double(X);
    Y = double(Y);
end

% initialize model struct
model.designMat = [];
model.betas = [];
model.residual = [];
model.stdevs = [];
model.DoF = [];

% get number of TRs, voxels, and predictors
[nTRs,nVoxels] = size(Y);
nPredictors = size(X,2);
   
% compute degrees of freedom and store design matrix
model.DoF = size(Y,1)-rank(X);
model.designMat = X;

% estimate beta weights for each predictor in the design matrix
model.betas = X\Y;

% compute residual error at each time point
model.residual = Y-X*model.betas;

% estimate the standard deviation of each beta
indepNoise = inv(X'*eye(nTRs)*X);
errVar = sqrt(sum(model.residual.^2)/model.DoF);
model.stdevs = sqrt((diag(indepNoise).*diag(X'*X))*errVar.^2);

% reshape betas and stadard deviations
model.betas = reshape(model.betas,[1 nPredictors nVoxels]);
model.stdevs = reshape(model.stdevs,[1 nPredictors nVoxels]);

end
