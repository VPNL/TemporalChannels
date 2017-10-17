function delay_act = code_delay_act(stim)
% Helper function for coding delay activity.

delay_act = double(stim == 0);
for pp = 1:size(stim, 2)
    if sum(stim(:, pp)) > 0
        stim_on = find(stim(:, pp) == 1, 1);
    else
        stim_on = size(stim, 1);
    end
    delay_act(1:stim_on, pp) = 0;
end
stim_frames = sum(stim, 2) > 0;
delay_act(stim_frames, :) = 0;

end