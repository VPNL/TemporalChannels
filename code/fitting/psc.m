function tc = psc(tc, detrend_option)
% Preprocess voxel time series and convert to precent signal change
%
% INPUTS
% tc: TR by voxel matrix of fMRI time series data
% detrend_option:
%   1 = normaliztion without trend removal
%   2 = normalization with linear trend removal
%   3 = normalization with linear + quadratic trend removal
%
% OUTPUT
% tc: preprocessed TR by voxel data matrix
%
% adapted from vistaoft (http://github.com/vistalab/vistasoft/)
% AS 2/2017

if nargin == 1
    detrend_option = 3;
end

if isempty(tc)
    tc = [];
else
    nFrames = size(tc,1);

    % divide by mean of each voxel
    if detrend_option >= 1
        dc = nanmean(tc);
        dc(dc == 0 | isnan(dc)) = Inf;
        if sum(dc ~= Inf) > 0
            tc = tc ./ (ones(nFrames, 1) * dc);
        else
            tc = zeros(size(tc));
        end
    end
    
    % remove linear trend
    if detrend_option >= 2
        model = [(1:nFrames); ones(1,nFrames)]';
        model = bsxfun(@rdivide ,model, max(model));
        w = model \ tc;
        b = model * w;
        tc = tc - b;
    end
    
    % remove quadratic trend
    if detrend_option >= 3
        model = [(1:nFrames) .* (1:nFrames); (1:nFrames); ones(1, nFrames)]';
        model = bsxfun(@rdivide,model,max(model));
        w = model \ tc;
        b = model * w;
        tc = tc - b;
    end
    
    % subtract mean and convert to percent signal change
    tc = tc - ones(nFrames, 1) * mean(tc);
    tc = tc * 100;
    
    % remove spikes
    tc = medfilt1(tc, 5);
    
end

end
