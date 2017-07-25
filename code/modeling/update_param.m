function model = update_param(model, param, step_size)
% Updates a specified parameter in ModelTS object
% 
% INPUTS
%   1) model: inital ModelTS object
%   2) param: name of parameter to update
%   3) step_size: size and direction of update step
% 
% OUTPUT
%   model: updated ModelTS object
% 
% AS 2/2017

if nargin ~= 3
    error('Unexpected input arguements.');
end

% update the specified parameter
param_init = model.params.(param);
switch param
    case 'e'
        param_min = 1e-4;
    case 'tau1'
        param_min = 1;
    case 'tau2'
        param_min = 1;
    case 'sigma'
        param_min = 1e-4;
end
model.params.(param) = cellfun(@(X) max([param_min X + step_size]), param_init, 'uni', false);

% update IRFs accordingly
switch model.type
    case 'cts'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X / sum(X), 1, 1000 / model.fs)', nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
    case 'cts-norm'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X/sum(X), 1, 1000/model.fs)', nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
    case 'dcts'
        lpf = cellfun(@(X) exp(-(0:999) / X), model.params.tau2, 'uni', false);
        lpf = cellfun(@(X) X / sum(X), lpf, 'uni', false);
        model.irfs.lpf = lpf;
    case '2ch-dcts'
        lpf = cellfun(@(X) exp(-(0:999) / X), model.params.tau2, 'uni', false);
        lpf = cellfun(@(X) X / sum(X), lpf, 'uni', false);
        model.irfs.lpf = lpf;
    case '3ch'
        lpf = cellfun(@(X) exp(-(0:999) / X), model.params.tau2, 'uni', false);
        lpf = cellfun(@(X) X / sum(X), lpf, 'uni', false);
        model.irfs.lpf = lpf;
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X/sum(X), 1, 1000/model.fs)', nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
end

end
