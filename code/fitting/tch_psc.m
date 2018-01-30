function tc = tch_psc(tc, detrend_option)
% Preprocess voxel time series and convert to precent signal change.
%
% INPUTS
%   1) tc: TR by voxel data matrix of fMRI time series
%   2) detrend_option: stages of detrending to perform on each voxel
%      1 = normaliztion without trend removal
%      2 = normalization with linear trend removal
%      3 = normalization with linear + quadratic trend removal
%      4 = normalization with low frequency baseline drift removal
%
% OUTPUT
%   tc: preprocessed TR by voxel data matrix
%
% Adapted from vistaoft (http://github.com/vistalab/vistasoft/)
% AS 2/2017

if nargin == 1; detrend_option = 3; end

if isempty(tc)
    tc = [];
else
    [num_frames, num_vox] = size(tc);
    % divide by mean of each voxel
    if detrend_option >= 1
        dc = nanmean(tc);
        dc(dc == 0 | isnan(dc)) = Inf;
        if sum(dc ~= Inf) > 0
            tc = tc ./ (ones(num_frames, 1) * dc);
        else
            tc = zeros(size(tc));
        end
    end
    % remove linear trends
    if detrend_option > 1
        model = [(1:num_frames); ones(1, num_frames)]';
        model = bsxfun(@rdivide, model, max(model));
        w = model \ tc;
        b = model * w;
        tc = tc - b;
    end
    % remove quadratic trends
    if detrend_option > 2
        fl = 1:num_frames;
        model = [fl .* fl; fl; ones(1, num_frames)]';
        model = bsxfun(@rdivide, model, max(model));
        w = model \ tc;
        b = model * w;
        tc = tc - b;
    end
    % remove low frequency baseline drifts
    if detrend_option > 3
        frame_win = 20; k = ones(frame_win, 1) / frame_win;
        pad_frames = num_frames + 2 * frame_win;  niter = 2;
        % initialize baseline array for single period
        baseline = zeros(pad_frames, num_vox);
        fwm  = mean(tc(1:frame_win, :));
        for ff = 1:frame_win
            baseline(ff, :) = fwm;
        end
        baseline(frame_win + 1:frame_win + num_frames, :) = tc;
        lwm = mean(tc(num_frames - frame_win + 1:num_frames, :));
        for ff = frame_win + num_frames + 1:pad_frames
            baseline(ff, :) = lwm;
        end
        % find all window indices and smooth with boxcar to remove drifts
        idxs  = niter * (frame_win - 1);
        start = floor(idxs / 2) + 1; stop = pad_frames + floor(idxs/2);
        for ii = 1:niter
            baseline = conv2(baseline, k);
        end
        baseline = baseline(start:stop, :);
        baseline = baseline(frame_win + 1:frame_win + num_frames, :);
        tc = tc - baseline;
    end
    % subtract mean and convert to percent signal change
    tc = (tc - ones(num_frames, 1) * mean(tc)) * 100;
    % remove spikes
    tc = medfilt1(tc, 5);
end

end
