function all_sessions = find_sessions(project_dir)
% Finds paths to all sessions in project data directory and returns an
% error if none are found.
% 
% INPUT: path to project directory
% OUTPUT: paths to all session directories in project data directory
% 
% AS 5/2017

data_dirs = dir(fullfile(project_dir, 'data'));
all_sessions = {data_dirs([data_dirs.isdir]).name};
all_sessions(ismember(all_sessions, {'.' '..'})) = [];

if isempty(all_sessions)
    error('No sessions found in data directory.');
end
    
end

