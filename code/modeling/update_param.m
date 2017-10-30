function model = update_param(model, param, param_step)
% Updates a specified parameter in tchModel object
% 
% INPUTS
%   1) model: inital tchModel object
%   2) param: name of parameter to update
%   3) param_step: step size for coordinate update
% 
% OUTPUT
%   model: updated tchModel object
% 
% AS 2/2017

if nargin ~= 3; error('Unexpected input arguements.'); end

% get minimum values for appropriate parameter
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
    case 'kappa'
        param_min = 1;
    case 'tau_ae'
        param_min = 100;
    case 'tau_de'
        param_min = 10;
    otherwise
        error('Input param is not recognized as valid parameter name.');
end

% update parameter respecting minimum values
model.params.(param) = cellfun(@(X) max([param_min X + param_step]), ...
    param_init, 'uni', false);

% update IRFs struct depending on model
switch model.type
    case '1ch-pow'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), ...
            model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X / sum(X), 1, 1000 / model.fs)', ...
            nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
    case '1ch-div'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), ...
            model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X/sum(X), 1, 1000/model.fs)', ...
            nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
    case '1ch-dcts'
        nrfS = cellfun(@(X) (0:999) .* exp(-(0:999) / X), ...
            model.params.tau1, 'uni', false);
        nrfS = cellfun(@(X) resample(X/sum(X), 1, 1000/model.fs)', ...
            nrfS, 'uni', false);
        model.irfs.nrfS = nrfS;
        lpf = cellfun(@(X) exp(-(0:999) / X), ...
            model.params.tau2, 'uni', false);
        lpf = cellfun(@(X) X / sum(X), lpf, 'uni', false);
        model.irfs.lpf = lpf;
    case '1ch-exp'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model);
    case '1ch-cexp'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model, 'cexp');
    case '2ch-dcts-quad'
        lpf = cellfun(@(X) exp(-(0:999) / X), ...
            model.params.tau2, 'uni', false);
        lpf = cellfun(@(X) X / sum(X), lpf, 'uni', false);
        model.irfs.lpf = lpf;
    case '2ch-exp-quad'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model);
    case '2ch-exp-rect'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model);
    case '2ch-cexp-quad'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model, 'cexp');
    case '2ch-cexp-rect'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model, 'cexp');
    case '3ch-lin-quad-exp'
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-lin-rect-exp'
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-pow-quad-exp'
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-pow-rect-exp'
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-exp-quad-exp'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model);
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-exp-rect-exp'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model);
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-cexp-quad-exp'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model);
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-cexp-rect-exp'
        adapt_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_ae, 'uni', false);
        model.irfs.adapt_exp = adapt_exps; model = code_adapt_decay(model);
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '2ch-lin-quad-opt'
        model.irfs.nrfS = cellfun(@(X, Y) tch_irfs('S', X, Y, model.fs), ...
            model.params.tau_s, model.params.kappa, 'uni', false);
        model.irfs.nrfT = cellfun(@(X, Y) tch_irfs('T', X, Y, model.fs), ...
            model.params.tau_s, model.params.kappa, 'uni', false);
    case '3ch-lin-quad-exp-opt'
        model.irfs.nrfS = cellfun(@(X, Y) tch_irfs('S', X, Y, model.fs), ...
            model.params.tau_s, model.params.kappa, 'uni', false);
        model.irfs.nrfT = cellfun(@(X, Y) tch_irfs('T', X, Y, model.fs), ...
            model.params.tau_s, model.params.kappa, 'uni', false);
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
    case '3ch-lin-rect-exp-opt'
        model.irfs.nrfS = cellfun(@(X, Y) tch_irfs('S', X, Y, model.fs), ...
            model.params.tau_s, model.params.kappa, 'uni', false);
        model.irfs.nrfT = cellfun(@(X, Y) tch_irfs('T', X, Y, model.fs), ...
            model.params.tau_s, model.params.kappa, 'uni', false);
        delay_exps = cellfun(@(X) exp(-(1:60000) / X), ...
            model.params.tau_de, 'uni', false);
        model.irfs.delay_exp = delay_exps; model = code_delay_decay(model);
end

end
