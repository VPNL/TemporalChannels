function idxs = idx_trials(roi)
% Index trials of the same type across different experiments.

nexps = length(roi.experiments); cond_list = roi.model.cond_list;
all_conds = unique([cond_list{:}], 'stable'); nconds = length(all_conds);
if nexps > 1
    idxs = zeros(nconds, nexps);
    for cc = 1:nconds
        cond = all_conds{cc};
        cond_idx = cellfun(@(X) find(strcmp(cond, X)), cond_list, 'uni', false);
        for ee = 1:nexps
            if ~isempty(cond_idx{ee})
                idxs(cc, ee) = cond_idx{ee};
            end
        end
    end
else
    idxs = [1:length(cond_list{1})]';
end

end
