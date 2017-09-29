function groups = group_trials(roi)

nexps = size(roi.experiments, 2); groups = cell(1, nexps);
for ee = 1:nexps
    cstems = cellfun(@(X) X(1:4), roi.model.cond_list{ee}, 'uni', false);
    ustems = unique(cstems, 'stable');
    groups{ee} = cell(1, length(ustems));
    for cc = 1:length(ustems)
        groups{ee}{cc} = find(strcmp(ustems{cc}, cstems));
    end
end

end
