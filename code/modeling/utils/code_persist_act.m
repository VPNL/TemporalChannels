function persist_act = code_persist_act(stim)
% Helper function for coding persistent activity step functions.
% 
% INPUT
%   stim: matrix of stimulus step functions (TR x category)
% 
% OUTPUT
%   persist_act: matrix of persistent activity step functions (TR x category)
% 
% AS 10/2017

stim_diffs = diff(stim); stim_idxs = find(diff(sum(stim, 2)) == 1);
persist_act = zeros(size(stim));
for cc = 1:size(stim, 2)
    donsets = find(stim_diffs(:, cc) == -1);
    for oo = 1:length(donsets)
        next_stim_idx = find(stim_idxs > donsets(oo), 1);
        if ~isempty(next_stim_idx)
            doffset = stim_idxs(find(stim_idxs > donsets(oo), 1));
        else
            doffset = size(stim, 1);
        end
        persist_act(donsets(oo):doffset, cc) = 1;
    end
end

end
