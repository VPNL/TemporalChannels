function output = loadTS(filename,field)
% Helper function that handles errors when loading files
% 
% INPUTS
% filename: path to *.mat file containing data (can be [])
% field: field in data structure with time series data
% 
% OUTPUT
% output: data loaded from *.mat file (or [] if not found)
% 
% AS 4/2017

if ~isempty(filename)
    try
        oo = load(filename,field);
        output = getfield(oo,field);
    catch
        output = [];
    end
else
    output = [];
end

end
