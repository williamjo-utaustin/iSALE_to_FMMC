function offsets = sliceOffsets(n, start_offset)
%SLICE_OFFSETS Return the starting angle (radians) for each slice.
%   offsets = SLICE_OFFSETS(n) returns n offsets from 0 to <2π.
%   offsets = SLICE_OFFSETS(n, start_offset) starts at start_offset (radians).

    if nargin < 2, start_offset = 0; end
    delta = 2*pi / n;
    offsets = start_offset + (0:n-1) * delta;
    offsets = mod(offsets, 2*pi);   % wrap to [0, 2π)
end
