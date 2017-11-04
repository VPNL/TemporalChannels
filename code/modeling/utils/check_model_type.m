function optimize_flag = check_model_type(model_type)

valid_types = {'1ch-lin' '1ch-balloon' '1ch-pow' '1ch-div' '1ch-dcts' '1ch-exp' ...
    '2ch-lin-lin' '2ch-lin-htd' '2ch-lin-quad' '2ch-lin-rect' ...
    '2ch-pow-quad' '2ch-pow-rect' '2ch-div-quad' '2ch-exp-quad' '2ch-exp-rect' ...
    '3ch-lin-quad-exp' '3ch-lin-rect-exp' '3ch-pow-quad-exp' '3ch-pow-rect-exp' ...
    '3ch-exp-quad-exp' '3ch-exp-rect-exp' ...
    '2ch-lin-quad-opt' '3ch-lin-quad-exp-opt' '3ch-lin-rect-exp-opt'};
if sum(strcmp(model_type, valid_types)) == 0; error('Invalid model'); end

optimize_flag = 0;
if contains(model_type, '-opt'); optimize_flag = 1; end
if contains(model_type, '-exp'); optimize_flag = 1; end
if contains(model_type, '-pow'); optimize_flag = 1; end
if contains(model_type, '-quad'); optimize_flag = 1; end
if contains(model_type, '-div'); optimize_flag = 1; end
if contains(model_type, '-dcts'); optimize_flag = 1; end

end
