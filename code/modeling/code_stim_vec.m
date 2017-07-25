function stim_out = code_stim_vec(stim_in, inds, col, val)
% Helper function for coding stimulus step functions.
% 
% INPUTS
%   1) stim_in: input stimulus time vector
%   2) inds: indicdes to code as val
%   3) col: column indices to code
%   4) val: value to code inds in col
% 
% OUTPUT
%   stim_out: stimulus time vector with indices in col set to val
% 
% AS 2/2017

stim_out = stim_in;
stim_out(inds, col) = val;

end
