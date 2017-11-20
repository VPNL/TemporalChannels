function data = tch_load(filename, field)
% Helper function that loads .mat files containing time series data with 
% error handling when data is not found.
% 
% INPUTS
% filename: path to *.mat file containing data (can be [])
% field: field in data structure with time series data
% 
% OUTPUT
% data: data loaded from *.mat file (or [] if not found)
% 
% AS 4/2017

if ~isempty(filename)
    try
        data = double(getfield(load(filename, field), field));
    catch
        data = [];
    end
else
    data = [];
end

end
