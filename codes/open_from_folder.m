function [lower_lampshade, upper_spike, stats] = open_from_folder(file_folder_containing_data, x, hasHeader, cases, timesteps)
% OPEN_FROM_FOLDER  Load target/impactor CSVs from a given "upright_porosity_X" folder.
%   [lower_lampshade, upper_spike, stats] = open_from_folder(folder, x, hasHeader, cases, timesteps)
%
% Inputs
%   file_folder_containing_data : char/string, e.g. 'input/upright_porosity_0'
%   x           : (optional) integer X in 'upright_porosity_X'. If omitted, parsed from folder name.
%   hasHeader   : (optional) default true
%   cases       : (optional) default setdiff(1:8, 4)
%   timesteps   : (optional) default 0:99
%
% Outputs
%   lower_lampshade : matrix for kinds={'target'}
%   upper_spike     : matrix for kinds={'impactor'}
%   stats           : struct with fields .opened_target, .opened_impactor, .inputDir
%
% Notes
% - Assumes you have a function `load_porosity_arrays(x, hasHeader, cases, timesteps, kinds)`
%   available on the MATLAB path (we add ./codes to the path below).
% - This function lives in your project root alongside the 'codes' folder.

    % ---------------- defaults ----------------
    if nargin < 2 || isempty(x)
        % try to parse X from folder name '...upright_porosity_X'
        tok = regexp(char(file_folder_containing_data), 'upright_porosity_(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            x = str2double(tok{1});
        else
            x = 0;  % fallback
        end
    end
    if nargin < 3 || isempty(hasHeader),  hasHeader = true; end
    if nargin < 4 || isempty(cases),      cases = setdiff(1:8, 4); end
    if nargin < 5 || isempty(timesteps),  timesteps = 0:99; end

    % --------------- add ./codes ----------------
    % --- locate and add the codes/ folder exactly once ---
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    [maybeParent, last] = fileparts(thisDir);
    
    if strcmpi(last, 'codes')
        % This function is already inside codes/
        codesDir = thisDir;
    else
        % This function is in project root (or elsewhere) → expect a sibling codes/
        codesDir = fullfile(thisDir, 'codes');
        if ~isfolder(codesDir)
            % Try parent/codes (in case this function is in a subfolder of project root)
            codesDir = fullfile(maybeParent, 'codes');
        end
    end
    
    if isfolder(codesDir)
        if ~contains(path, [codesDir pathsep])   % avoid duplicate warning
            addpath(codesDir);
        end
    else
        warning('Could not find a "codes" folder near: %s', thisDir);
    end


    % --------------- sanity check folder --------
    % (We don't pass the folder directly to load_porosity_arrays because your
    %  implementation builds the path from projectRoot and x; we just verify.)
    if ~isfolder(file_folder_containing_data)
        warning('Folder not found: %s (proceeding anyway using x=%d)', file_folder_containing_data, x);
    end

    % --------------- load TARGET (lower_lampshade) ---------------
    [lower_lampshade, opened_target, inputDir_target] = load_porosity_arrays( ...
        x, hasHeader, cases, timesteps, {'target'});

    fprintf('Opened %d files from %s\n', opened_target, inputDir_target);
    if ~isempty(lower_lampshade)
        fprintf('Lower Lampshade Size: [%d rows, %d cols]\n', size(lower_lampshade,1), size(lower_lampshade,2));
    else
        fprintf('No data loaded for TARGET.\n');
    end

    % --------------- load IMPACTOR (upper_spike) ---------------
    [upper_spike, opened_impactor, inputDir_impactor] = load_porosity_arrays( ...
        x, hasHeader, cases, timesteps, {'impactor'});

    fprintf('Opened %d files from %s\n', opened_impactor, inputDir_impactor);
    if ~isempty(upper_spike)
        fprintf('Upper Spike Size: [%d rows, %d cols]\n', size(upper_spike,1), size(upper_spike,2));
    else
        fprintf('No data loaded for IMPACTOR.\n');
    end

    % --------------- pack stats ----------------
    stats = struct( ...
        'opened_target',   opened_target, ...
        'opened_impactor', opened_impactor, ...
        'inputDir_target', inputDir_target, ...
        'inputDir_impactor', inputDir_impactor, ...
        'x', x, ...
        'hasHeader', hasHeader, ...
        'cases', cases, ...
        'timesteps', timesteps ...
    );
end
