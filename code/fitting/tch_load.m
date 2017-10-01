function data = tch_load(filename,field)
% Helper function that loads .mat files with error handling.
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
        oo = load(filename, field);
        data = getfield(oo, field);
    catch
        data = [];
    end
else
    data = [];
end

end
