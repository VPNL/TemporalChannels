function stats = tch_contrast(model, active, control)

stats = struct('contrast', [], 'ces', [], 'sem', [], 't', [], 'p', []);
cvecs = cellfun(@(X, Y) zeros(1, length(X) + length(Y)), ...
    model.betas, model.rbetas, 'uni', false);
cvecs = cellfun(@(X) code_stim_vec(X', active, 1, 1 / length(active))', ...
    cvecs, 'uni', false);
cvecs = cellfun(@(X) code_stim_vec(X', control, 1, -1 / length(control))', ...
    cvecs, 'uni', false); stats.contrast = cvec;
stats.ces = cellfun(@(C, X, Y) C * [X Y]', cvecs, model.betas, model.rbetas);
stats.sem = cellfun(@(R, C, N) sqrt(R * (C * N * C')), ...
    model.resid_var, cvecs, cellfun(@inv, model.var_covar, 'uni', false));
stats.t = stats.ces ./ stats.sem;
stats.p = cellfun(@(D, T) betainc(D ./ (D + T .^ 2), D / 2, 0.5), ...
    model.dof, num2cell(stats.t));

end
