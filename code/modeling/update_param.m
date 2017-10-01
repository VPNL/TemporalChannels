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
    case 'epsilon'
        param_min = 1e-4;
    case 'tau1'
        param_min = 1;
    case 'tau2'
        param_min = 1;
    case 'sigma'
        param_min = 1e-4;
    case 'tau_s'
        param_min = 1;
    case 'tau_t'
        param_min = 1;
    case 'tau_d'
        param_min = 1;
end
model.params.(param) = cellfun(@(X) max([param_min X + step_size]), param_init, 'uni', false);

% update IRFs accordingly
switch model.type
    case 'cts-pow'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X / sum(X), 1, 1000 / model.fs)', nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
    case 'cts-div'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X/sum(X), 1, 1000/model.fs)', nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
    case 'dcts'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X/sum(X), 1, 1000/model.fs)', nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
        lpf = cellfun(@(X) exp(-(0:999) / X), model.params.tau2, 'uni', false);
        lpf = cellfun(@(X) X / sum(X), lpf, 'uni', false);
        model.irfs.lpf = lpf;
    case '2ch-dcts-quad'
        lpf = cellfun(@(X) exp(-(0:999) / X), model.params.tau2, 'uni', false);
        lpf = cellfun(@(X) X / sum(X), lpf, 'uni', false);
        model.irfs.lpf = lpf;
    case '2ch-opt'
        model.irfs.nrfS = cellfun(@(X) tch_irfs('S', X, model.fs), ...
            model.params.tau_s, 'uni', false);
        model.irfs.nrfT = cellfun(@(X) tch_irfs('T', X, model.fs), ...
            model.params.tau_t, 'uni', false);
    case '3ch-lin-quad'
        nrfD = cellfun(@(X) tch_irfs('D', X, model.fs), model.params.tau_d, 'uni', false);
        model.irfs.nrfD = nrfD;
    case '3ch-lin-rect'
        nrfD = cellfun(@(X) tch_irfs('D', X, model.fs), model.params.tau_d, 'uni', false);
        model.irfs.nrfD = nrfD;
    case '3ch-pow-quad'
        nrfD = cellfun(@(X) tch_irfs('D', X, model.fs), model.params.tau_d, 'uni', false);
        model.irfs.nrfD = nrfD;
    case '3ch-pow-rect'
        nrfD = cellfun(@(X) tch_irfs('D', X, model.fs), model.params.tau_d, 'uni', false);
        model.irfs.nrfD = nrfD;
    case '3ch-opt'
        model.irfs.nrfS = cellfun(@(X) tch_irfs('S', X, model.fs), ...
            model.params.tau_s, 'uni', false);
        model.irfs.nrfT = cellfun(@(X) tch_irfs('T', X, model.fs), ...
            model.params.tau_t, 'uni', false);
        model.irfs.nrfD = cellfun(@(X) tch_irfs('D', X, model.fs), ...
            model.params.tau_d, 'uni', false);
end

end
