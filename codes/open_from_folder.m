function [lower_lampshade, upper_spike, stats] = open_from_folder(file_folder_containing_data, x, hasHeader, cases, timesteps, orientation)
% OPEN_FROM_FOLDER  Load target/impactor CSVs from a "<orientation>_porosity_X" folder.
%   [lower_lampshade, upper_spike, stats] = open_from_folder(folder, x, hasHeader, cases, timesteps, orientation)
%
% Inputs
%   file_folder_containing_data : char/string, e.g. 'input/upright_porosity_0' or 'input/flipped_porosity_5'
%   x           : (optional) integer X in '<orientation>_porosity_X'. If omitted, parsed from folder name.
%   hasHeader   : (optional) default true
%   cases       : (optional) default setdiff(1:8, 4)
%   timesteps   : (optional) default 0:99
%   orientation : (optional) char/string prefix before '_porosity_X', default 'upright'.
%
% Outputs
%   lower_lampshade : matrix for kinds={'target'}
%   upper_spike     : matrix for kinds={'impactor'}
%   stats           : struct with fields .opened_target, .opened_impactor, .inputDir_*, .orientation, .x

    % ---------------- defaults & parsing ----------------
    if nargin < 3 || isempty(hasHeader),  hasHeader = true; end
    if nargin < 4 || isempty(cases),      cases = setdiff(1:8, 4); end
    if nargin < 5 || isempty(timesteps),  timesteps = 0:99; end

    % Try to parse "<orientation>_porosity_<X>" from the provided folder path
    tok = regexp(char(file_folder_containing_data), '([A-Za-z]+)_porosity_(\d+)', 'tokens', 'once');

    if nargin < 2 || isempty(x)
        if ~isempty(tok), x = str2double(tok{2}); else, x = 0; end
    end

    if nargin < 6 || isempty(orientation)
        if ~isempty(tok), orientation = tok{1}; else, orientation = 'upright'; end
    end

    % --------------- add ./codes to path (once) --------
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    [maybeParent, last] = fileparts(thisDir);

    if strcmpi(last, 'codes')
        codesDir = thisDir;
    else
        codesDir = fullfile(thisDir, 'codes');
        if ~isfolder(codesDir)
            codesDir = fullfile(maybeParent, 'codes');
        end
    end
    if isfolder(codesDir) && ~contains(path, [codesDir pathsep])
        addpath(codesDir);
        addpath(genpath(codesDir));
        rehash;
    end

    % --------------- sanity check folder ----------------
    if ~isfolder(file_folder_containing_data)
        warning('Folder not found: %s (proceeding anyway with orientation="%s", x=%d)', ...
                 file_folder_containing_data, orientation, x);
    end

    % --------------- load TARGET (lower_lampshade) ------
    [lower_lampshade, opened_target, inputDir_target] = load_porosity_arrays( ...
        x, orientation, hasHeader, cases, timesteps, {'target'});

    fprintf('Opened %d files from %s\n', opened_target, inputDir_target);
    if ~isempty(lower_lampshade)
        fprintf('Lower Lampshade Size: [%d rows, %d cols]\n', size(lower_lampshade,1), size(lower_lampshade,2));
    else
        fprintf('No data loaded for TARGET.\n');
    end

    % --------------- load IMPACTOR (upper_spike) --------
    [upper_spike, opened_impactor, inputDir_impactor] = load_porosity_arrays( ...
        x, orientation, hasHeader, cases, timesteps, {'impactor'});

    fprintf('Opened %d files from %s\n', opened_impactor, inputDir_impactor);
    if ~isempty(upper_spike)
        fprintf('Upper Spike Size: [%d rows, %d cols]\n', size(upper_spike,1), size(upper_spike,2));
    else
        fprintf('No data loaded for IMPACTOR.\n');
    end

    % --------------- pack stats -------------------------
    stats = struct( ...
        'opened_target',       opened_target, ...
        'opened_impactor',     opened_impactor, ...
        'inputDir_target',     inputDir_target, ...
        'inputDir_impactor',   inputDir_impactor, ...
        'orientation',         orientation, ...
        'x',                   x, ...
        'hasHeader',           hasHeader, ...
        'cases',               cases, ...
        'timesteps',           timesteps ...
    );
end