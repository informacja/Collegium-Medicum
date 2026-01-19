
function [Files, Data] = load_all_mat_structs(dataDir)
%LOAD_ALL_MAT_STRUCTS Load all .mat files from a folder into memory (structures).
%
% USAGE:
%   [Files, Data] = load_all_mat_structs();                % uses pwd
%   [Files, Data] = load_all_mat_structs('C:\my\folder');  % specific folder
%
% OUTPUTS:
%   Files : table with variables:
%           - Name     : file name with extension
%           - Folder   : folder path
%           - FullPath : full absolute path to file
%   Data  : cell array, same length as Files:
%           - If the .mat contains exactly one structure variable, Data{k}
%             is that structure.
%           - Otherwise, Data{k} is a struct containing ALL variables
%             from the file (so nothing is lost).
%
% NOTES:
%   - Does NOT assign variables in the base workspace.
%   - Safe for mixed content: different files can hold different variable names.
%   - For very large v7.3 files, consider switching to `matfile` (later step).

    % -----------------------
    % 1) Resolve target path
    % -----------------------
    if nargin < 1 || isempty(dataDir)
        dataDir = pwd;  % default: current working directory
    end

    if ~isfolder(dataDir)
        error('The provided folder does not exist: %s', dataDir);
    end

    % -----------------------
    % 2) Enumerate .mat files
    % -----------------------
    fileList = dir(fullfile(dataDir, '*.mat'));

    if isempty(fileList)
        warning('No .mat files found in folder: %s', dataDir);
        Files = table(string.empty, string.empty, string.empty, ...
                      'VariableNames', {'Name','Folder','FullPath'});
        Data  = {};
        return;
    end

    % Build a table with file paths for traceability
    Names     = string({fileList.name}).';
    Folders   = string({fileList.folder}).';
    FullPaths = fullfile(Folders, Names);
    Files     = table(Names, Folders, FullPaths, ...
                      'VariableNames', {'Name','Folder','FullPath'});

    % -----------------------
    % 3) Load loop (safe)
    % -----------------------
    n = height(Files);
    Data = cell(n,1);

    % (Optional) progress text
    fprintf('Loading %d .mat file(s) from: %s\n', n, dataDir);

    for k = 1:n
        fp = Files.FullPath{k};
        try
            % Inspect variables without loading into workspace:
            info = whos('-file', fp);

            % Load everything into a temporary struct:
            tmp = load(fp);

            % If the file contains exactly ONE structure variable,
            % return that structure directly for convenience:
            if numel(info) == 1 && strcmp(info(1).class, 'struct')
                singleVarName = info(1).name;
                Data{k} = tmp.(singleVarName);
            else
                % Otherwise, preserve ALL variables under one struct:
                % (This avoids losing anything and keeps the interface uniform)
                Data{k} = tmp;
            end

            fprintf('  [%3d/%3d] OK  - %s\n', k, n, Files.Name{k});

        catch ME
            % If a file fails to load, store the error message
            % but do not abort the whole batch.
            Data{k} = struct('load_error', ME.message);
            fprintf('  [%3d/%3d] ERR - %s\n      -> %s\n', ...
                k, n, Files.Name{k}, ME.message);
        end
    end

    fprintf('Done.\n');
end
