function [all_data, openedCount, inputDir] = load_porosity_arrays(x, hasHeader, cases, timesteps, kinds)
%GATHER_POROSITY_ALL Load and vertically concatenate CSVs into one array.
%   [ALL_DATA, OPENEDCOUNT, INPUTDIR] = GATHER_POROSITY_ALL(X, HASHEADER, CASES, TIMESTEPS, KINDS)
%   - X: numeric (default 0) -> uses input/upright_porosity_X/
%   - HASHEADER: logical (default true) -> skip first line when reading
%   - CASES: vector (default 1..8, skipping 4)
%   - TIMESTEPS: vector (default 0..99)
%   - KINDS: cellstr (default {'impactor','target'})
%
%   ALL_DATA is a numeric matrix (rows stacked). If column counts differ
%   across files, shorter ones are padded with NaNs before concatenation.
%
%   This function assumes it lives in the `codes/` folder one level below
%   your project root (i.e., iSALE_to_FMMC/codes/thisfile.m).

    % -------------------- defaults --------------------
    if nargin < 1 || isempty(x), x = 0; end
    if nargin < 2 || isempty(hasHeader), hasHeader = true; end
    if nargin < 3 || isempty(cases), cases = setdiff(1:8, 4); end
    if nargin < 4 || isempty(timesteps), timesteps = 0:99; end
    if nargin < 5 || isempty(kinds), kinds = {'impactor','target'}; end

    % ----------------- resolve paths ------------------
    thisFile = mfilename('fullpath');
    codesDir = fileparts(thisFile);
    projectRoot = fileparts(codesDir);           % assumes this file is in ./codes
    inputDir = fullfile(projectRoot, 'input', sprintf('upright_porosity_%d', x));

    if ~isfolder(inputDir)
        error('Input folder not found: %s', inputDir);
    end

    % ------------------- read loop --------------------
    all_data = [];
    openedCount = 0;

    for y = cases
        for k = 1:numel(kinds)
            kind = kinds{k};
            for t = timesteps
                fname = sprintf('upright_porosity_%d_case_%d_%s_%d.csv', x, y, kind, t);
                fpath = fullfile(inputDir, fname);

                if ~isfile(fpath), continue; end

                try
                    if hasHeader
                        A = readmatrix(fpath, 'NumHeaderLines', 1);
                    else
                        A = readmatrix(fpath);
                    end
                catch
                    % unreadable -> skip
                    continue;
                end

                if isempty(A), continue; end

                % Ensure 2D row shape for vectors
                if isvector(A)
                    A = A(:).';  % row
                end

                % Attach to running array with NaN padding if needed
                all_data = append_to_array(all_data, A);
                openedCount = openedCount + 1;
            end
        end
    end


    % divide by the number of cases, as this is an ensemble average
    N = numel(cases);              % e.g., setdiff(1:8,4) -> N = 7
    all_data(:,6) = all_data(:,6) ./ N;


end

% ---------- local helper (pads columns, then appends rows) ----------
function B = append_to_array(B, A)
    if isempty(A), return; end
    if isempty(B)
        B = A;
        return;
    end
    c1 = size(B, 2);
    c2 = size(A, 2);
    if c1 < c2
        B(:, end+1:c2) = NaN;
    elseif c2 < c1
        A(:, end+1:c1) = NaN;
    end
    B = [B; A];
end
