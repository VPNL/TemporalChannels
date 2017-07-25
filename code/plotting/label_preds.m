function labels = label_preds(model)
% Generates list of predictor labels for a given model. 

if sum(strcmp(model.type, {'standard' 'cts' 'dcts'}))
    labels = model.cat_list;
else
    xlabsS = cellfun(@(X) [X '-S'], model.cat_list, 'uni', false);
    xlabsT = cellfun(@(X) [X '-T'], model.cat_list, 'uni', false);
    if strcmp(model.type, '3ch')
        xlabsD = cellfun(@(X) [X '-D'], model.cat_list, 'uni', false);
        labels = [xlabsS xlabsT xlabsD];
    else
        labels = [xlabsS xlabsT];
    end
end

end
