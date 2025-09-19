function plume = full_plume_from_slice(particles_slice, clonesPer, pct_radial, pct_vertical, slice_width, baseSeed)
% FULL_PLUME_FROM_SLICE  Replicate a wedge (slice) into a full 2π plume.
%   plume = full_plume_from_slice(particles_slice, clonesPer, pct_radial, pct_vertical, slice_width, baseSeed)
%
% Inputs
%   particles_slice.pos : (N x 3) positions
%   particles_slice.vel : (N x 3) velocities [vx vy vz]
%   clonesPer           : clones per particle per slice (default 10)
%   pct_radial          : % 1σ multiplicative jitter for horizontal speed (default 5)
%   pct_vertical        : % 1σ multiplicative jitter for vy (default 5)
%   slice_width         : nominal wedge width [rad] (e.g., pi/16)
%   baseSeed            : base RNG seed ([] = don't set per-slice seeds)
%
% Output
%   plume.pos, plume.vel        : ((N*clonesPer*S) x 3)
%   plume.slice_idx             : slice index [0..S-1] for each clone
%   plume.parent_idx            : original particle index (before cloning)
%
% Notes
%   - Uses uniform azimuth within each wedge.
%   - Ensures exact 2π coverage by using Δφ = 2π / S, where S = round(2π/slice_width).
%   - If baseSeed is provided, each slice uses rng(baseSeed + s).

    if nargin < 2 || isempty(clonesPer),    clonesPer = 10; end
    if nargin < 3 || isempty(pct_radial),   pct_radial = 5;  end
    if nargin < 4 || isempty(pct_vertical), pct_vertical = 5; end
    if nargin < 5 || isempty(slice_width),  slice_width = pi/16; end
    if nargin < 6, baseSeed = []; end

    % # of wedges and exact per-wedge width to tile 0..2π
    S = max(1, round((2*pi) / slice_width));
    dphi = 2*pi / S;   % exact tile width
    N = size(particles_slice.pos, 1);

    POS = cell(S,1);
    VEL = cell(S,1);
    SLICE = cell(S,1);
    PARENT = cell(S,1);

    % Use uniform φ in [offset, offset + dphi) for each wedge
    phiUniform = true;
    angJitterDeg = 0;       % irrelevant when phiUniform=true
    clampPositive = true;
    keep_in_slice = true;

    for s = 0:S-1
        rngSeed = [];
        if ~isempty(baseSeed)
            rngSeed = double(baseSeed) + s;   % deterministic per slice
        end

        % Delegate cloning/jitter + uniform φ selection within this wedge
        pj = split_and_jitter_particles( ...
                particles_slice, ...
                clonesPer, ...
                pct_radial, ...
                pct_vertical, ...
                rngSeed, ...
                clampPositive, ...
                angJitterDeg, ...   % ignored when phiUniform=true
                phiUniform, ...
                dphi, ...           % slice_width for this wedge
                s*dphi, ...         % slice_offset for this wedge
                keep_in_slice);

        POS{s+1}    = pj.pos;
        VEL{s+1}    = pj.vel;
        SLICE{s+1}  = repmat(s, size(pj.pos,1), 1);
        PARENT{s+1} = pj.parent_idx;  % original index (before cloning)
    end

    plume.pos        = vertcat(POS{:});
    plume.vel        = vertcat(VEL{:});
    plume.slice_idx  = vertcat(SLICE{:});
    plume.parent_idx = vertcat(PARENT{:});
end