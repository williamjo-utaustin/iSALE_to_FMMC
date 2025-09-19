function S = collapse_struct_array(A, fields)
%COLLAPSE_STRUCT_ARRAY Vertically concat fields across a struct array.
%   S = COLLAPSE_STRUCT_ARRAY(A)    % uses all fields from A(1)
%   S = COLLAPSE_STRUCT_ARRAY(A, {'pos','vel','acc'})  % only these fields
%
% Assumes each chosen field is numeric with consistent number of columns.
% Concatenates along the first dimension (rows).

    if isempty(A)
        error('Input struct array A is empty.');
    end
    if nargin < 2 || isempty(fields)
        fields = fieldnames(A(1));
    end

    S = struct();
    for i = 1:numel(fields)
        f = fields{i};
        % Collect pieces from each element
        parts = arrayfun(@(x) x.(f), A, 'UniformOutput', false);

        % Basic validation: numeric + same # of columns
        assert(all(cellfun(@isnumeric, parts)), ...
               'Field "%s" must be numeric in all elements.', f);
        ncols = cellfun(@(m) size(m,2), parts);
        assert(numel(unique(ncols)) == 1, ...
               'Field "%s" must have consistent number of columns.', f);

        % Efficient preallocation + fill (faster than vertcat on many parts)
        rows = cellfun(@(m) size(m,1), parts);
        total_rows = sum(rows);
        S.(f) = zeros(total_rows, ncols(1), class(parts{1}));

        idx = 1;
        for k = 1:numel(parts)
            r = rows(k);
            S.(f)(idx:idx+r-1, :) = parts{k};
            idx = idx + r;
        end
    end
end
