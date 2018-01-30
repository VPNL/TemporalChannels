function fig = plot_roi_runs(roi, save_flag)
% plot measurement vs. prediction for runs in each session

% check inputs
if length(roi) > 1; roi = roi(1); end
if nargin < 2; save_flag = 0; end

% setup figure
fig_name = [roi.nickname ' run timecourses'];
fig = tch_fig(fig_name, [.1 .1 .8 length(roi.sessions) * .1]);

% plot run time series and predictors for each session
for ss = 1:length(roi.sessions)
    subplot(length(roi.sessions), 1, ss); hold on; tch_set_axes;
    plot(roi.model.run_tcs{ss}, 'k');
    plot(roi.model.run_preds{ss}, 'r');
    if ss == 1; ylabel('% signal'); end
    session_id = strrep(roi.session_ids{ss}, '_', '-'); R2 = roi.model.varexp{ss};
    leg = {[session_id ': ' num2str(R2 * 100, 2) '%'] 'pred'};
    legend(leg, 'Location', 'NorthWestOutside'); legend boxoff; axis tight;
    ylims = get(gca, 'YLim'); ylim([ylims(1) ceil(ylims(2))]);
    set(gca, 'XColor', 'w', 'YTick', [0 ceil(ylims(2))]);
end

% save to figures directory if applicable
if save_flag
    fpath = fullfile(roi.project_dir, 'figures');
    fname = [roi.nickname '_runs_' roi.model.type ...
        '_fit' [roi.model.fit_exps{:}] ...
        '_val' [roi.experiments{:}] ...
        '_' date '.fig'];
    saveas(fig, fullfile(fpath, fname), 'fig');
end

end
