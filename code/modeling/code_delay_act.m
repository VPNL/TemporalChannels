function delay_act = code_delay_act(stim)
% Helper function for coding delay activity.

delay_act = double(stim == 0);
for pp = 1:size(delay_act, 2)
    stim_on = find(stim(:, pp) == 1, 1);
    delay_act(1:stim_on, pp) = 0;
end
stim_frames = sum(stim, 2) > 0;
delay_act(stim_frames, :) = 0;

end