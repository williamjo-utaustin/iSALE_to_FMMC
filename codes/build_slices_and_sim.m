function out = build_slices_and_sim( ...
    particles_base, slice_width, slice_offset0, S, ...
    clonesPer, pct_radial, pct_vertical, angJitterDeg, phiUniform, keep_in_slice, rngSeed0, ...
    dt, steps, g, times_list)
% BUILD_SLICES_AND_SIM
% Replicate a wedge around 2π into S slices, simulate each, and concatenate.
%
% Inputs
%   particles_base : struct with .pos (N x 3), .vel (N x 3)
%   slice_width    : wedge width [rad] (e.g., pi/8)
%   slice_offset0  : base offset [rad] for the first slice (usually 0)
%   S              : number of slices (if empty, inferred as round(2*pi/slice_width))
%   clonesPer, pct_radial, pct_vertical, angJitterDeg, phiUniform, keep_in_slice, rngSeed0
%                   -> passed to split_and_jitter_particles (rng per slice = rngSeed0 + s)
%   dt, steps, g, times_list  -> passed to move_particles (times_list can be [] for dense)
%
% Output (struct)
%   .traj            : (Trec x Ntot x 3) concatenated positions [x y z]
%   .landed_mask     : (Trec x Ntot) logical
%   .frameTimesSec   : (1 x Trec) recorded times (s)
%   .particles_per_slice : 1xS array (# particles in each slice)
%   .slice_offsets   : 1xS array of offsets actually used
%   .meta            : per-slice meta (1xS cell), each with .frame_times_sec

    if nargin < 4 || isempty(S)
        S = max(1, round((2*pi)/slice_width));  % integer tiling
    end

    % storage
    parts = cell(1,S);
    trajs = cell(1,S);
    lands = cell(1,S);
    metas = cell(1,S);
    Ns    = zeros(1,S);
    slice_offsets = zeros(1,S);

    % build + sim each slice
    for s = 1:S
        offset_s = slice_offset0 + (s-1)*slice_width;
        slice_offsets(s) = offset_s;

        % unique seed per slice (reproducible; use a big stride)
        rng_s = [];
        if ~isempty(rngSeed0), rng_s = double(rngSeed0) + 104729*(s-1); end  % 104729 is prime

        % make a slice from the same base particles, but constrained/offset into this wedge
        p_s = split_and_jitter_particles( ...
            particles_base, clonesPer, pct_radial, pct_vertical, rng_s, true, ...
            angJitterDeg, phiUniform, slice_width, offset_s, keep_in_slice);

        Ns(s) = size(p_s.pos,1);
        parts{s} = p_s;

        % simulate this slice (sparse or dense)
        [p_s2, traj_s, land_s, meta_s] = move_particles(p_s, dt, steps, g, true, times_list); %#ok<ASGLU>
        trajs{s} = traj_s;
        lands{s} = land_s;
        metas{s} = meta_s;
    end

    % sanity: all slices should have the same recorded frame times
    frameTimesSec = metas{1}.frame_times_sec(:).';
    for s = 2:S
        if numel(metas{s}.frame_times_sec) ~= numel(frameTimesSec) || ...
           any(abs(metas{s}.frame_times_sec(:).' - frameTimesSec) > 1e-9)
            error('Recorded times differ across slices. Ensure the same times_list is used for all slices.');
        end
    end

    % concat along particle dimension
    Trec = size(trajs{1},1);
    Ntot = sum(Ns);
    traj_all   = zeros(Trec, Ntot, 3, 'like', trajs{1});
    land_all   = false(Trec, Ntot);

    col0 = 0;
    for s = 1:S
        cols = (col0+1):(col0+Ns(s));
        traj_all(:, cols, :) = trajs{s};
        land_all(:, cols)    = lands{s};
        col0 = col0 + Ns(s);
    end

    % pack output
    out = struct();
    out.traj              = traj_all;
    out.landed_mask       = land_all;
    out.frameTimesSec     = frameTimesSec;
    out.particles_per_slice = Ns;
    out.slice_offsets     = slice_offsets;
    out.meta              = metas;
end