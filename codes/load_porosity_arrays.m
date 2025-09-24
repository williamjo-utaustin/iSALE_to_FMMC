function [all_data, openedCount, inputDir] = load_porosity_arrays(x, orientation, hasHeader, cases, timesteps, kinds)
%LOAD_POROSITY_ARRAYS Load and vertically concatenate CSVs into one array.
%   [ALL_DATA, OPENEDCOUNT, INPUTDIR] = LOAD_POROSITY_ARRAYS(X, ORIENTATION, HASHEADER, CASES, TIMESTEPS, KINDS)
%   - X: numeric (default 0) -> uses input/<orientation>_porosity_X/
%   - ORIENTATION: char/string (default 'upright')
%   - HASHEADER: logical (default true) -> skip first line when reading
%   - CASES: vector (default 1..8, skipping 4)
%   - TIMESTEPS: vector (default 0..99)
%   - KINDS: cellstr (default {'impactor','target'})
%
%   ALL_DATA is a numeric matrix (rows stacked). If column counts differ,
%   shorter ones are padded with NaNs before concatenation.

    % -------------------- defaults --------------------
    if nargin < 1 || isempty(x),           x = 0; end
    if nargin < 2 || isempty(orientation), orientation = 'upright'; end
    if nargin < 3 || isempty(hasHeader),   hasHeader = true; end
    if nargin < 4 || isempty(cases),       cases = setdiff(1:8, 4); end
    if nargin < 5 || isempty(timesteps),   timesteps = 0:99; end
    if nargin < 6 || isempty(kinds),       kinds = {'impactor','target'}; end

    % ----------------- resolve paths ------------------
    thisFile    = mfilename('fullpath');
    codesDir    = fileparts(thisFile);
    projectRoot = fileparts(codesDir);  % assumes this file is in ./codes
    inputDir    = fullfile(projectRoot, 'input', sprintf('%s_porosity_%d', orientation, x));

    if ~isfolder(inputDir)
        error('Input folder not found: %s', inputDir);
    end

    % ------------------- read loop --------------------
    all_data    = [];
    openedCount = 0;

    for y = cases
        for k = 1:numel(kinds)
            kind = kinds{k};
            for t = timesteps
                fname = sprintf('%s_porosity_%d_case_%d_%s_%d.csv', orientation, x, y, kind, t);
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
                if isvector(A), A = A(:).'; end  % row shape

                all_data    = append_to_array(all_data, A);
                openedCount = openedCount + 1;
            end
        end
    end

    % divide by the number of cases, as this is an ensemble average
    if ~isempty(all_data) && size(all_data,2) >= 6
        N = numel(cases);
        all_data(:,6) = all_data(:,6) ./ N;
    end
end

% ---------- local helper (pads columns, then appends rows) ----------
function B = append_to_array(B, A)
    if isempty(A), return; end
    if isempty(B), B = A; return; end
    c1 = size(B, 2); c2 = size(A, 2);
    if c1 < c2
        B(:, end+1:c2) = NaN;
    elseif c2 < c1
        A(:, end+1:c1) = NaN;
    end
    B = [B; A];
end