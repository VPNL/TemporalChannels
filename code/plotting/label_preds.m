function labels = label_preds(model)
% Generates list of predictor labels for a given model. 

if sum(strcmp(model.type, {'glm' 'balloon' 'cts-pow' 'cts-div' 'dcts'}))
    labels = model.cat_list;
else
    xlabsS = cellfun(@(X) [X '-S'], model.cat_list, 'uni', false);
    xlabsT = cellfun(@(X) [X '-T'], model.cat_list, 'uni', false);
    labels = [xlabsS xlabsT];
end

if strcmp(model.type(1:3), '3ch')
    xlabsP = cellfun(@(X) [X '-P'], model.cat_list, 'uni', false);
    labels = [labels xlabsP];
end

end
