function optimize_flag = check_model_type(model_type)
% Function that checks if the model type argument passed as input to 
% tchModel object is a valid model implemented in the code. 
% INPUT
%   model_type: model name passed to tchModel (string)
% 
% OUTPUT
%   Error if input does not match a valid_model string, otherwise:
%   optimize_flag: 1 if model supports nonlinear optimization and 0 if not
% 
% AS 10/2017

valid_models = {'1ch-lin' '1ch-exp' '1ch-rect' '1ch-quad' '1ch-cquad' '1ch-sig'...
    '2ch-lin-rect' '2ch-exp-rect' '2ch-lin-quad' '2ch-exp-quad' ...
    '2ch-lin-cquad' '2ch-exp-cquad' '2ch-lin-sig' '2ch-exp-sig' ...
    '1ch-balloon' '1ch-pow' '1ch-div' '1ch-dcts'  '2ch-lin-lin' '2ch-lin-htd' ...
    '2ch-pow-quad' '2ch-pow-rect' '2ch-div-quad' '3ch-lin-quad-exp' ...
    '3ch-lin-rect-exp' '3ch-pow-quad-exp' '3ch-pow-rect-exp' '3ch-exp-quad-exp' ...
    '3ch-exp-rect-exp' '2ch-lin-crect' '2ch-exp-crect' '2ch-exp-dquad' ...
    '2ch-lin-dquad' '3ch-exp-cquad-rect' '3ch-exp-crect-crect'};
if sum(strcmp(model_type, valid_models)) == 0
    error('Invalid model');
end
no_opt_strs = {'balloon'};
if contains(model_type, no_opt_strs)
    optimize_flag = 0;
else
    optimize_flag = 1;
end

end
