function [particles_full, particles_combined_full] = generate_jitters_and_collapse( ...
    particles_10, particles, clonesPer, pct_radial, pct_vertical, ...
    rngSeed, angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice, verbose)
%GENERATE_JITTERS_AND_COLLAPSE
% Runs split_and_jitter_particles for every slice offset, adjusts mass per slice,
% and collapses the 1×N struct array into a single 1×1 struct (fields stacked).
%
% Inputs:
%   particles_10    : 1×1 struct (template from a prior call)
%   particles       : input particles struct with .pos, .vel, .mass
%   clonesPer       : clones per particle
%   pct_radial      : % 1σ multiplicative jitter for horizontal speed
%   pct_vertical    : % 1σ multiplicative jitter for vy
%   rngSeed         : base RNG seed
%   angJitterDeg    : azimuth jitter (deg) if phiUniform=false
%   phiUniform      : true => redraw φ uniformly in slice
%   slice_width     : slice width (radians)
%   slice_offset    : slice starting angle (radians)
%   keep_in_slice   : wrap φ back into [offset, offset+width)
%   verbose         : (optional) true to print sanity info (default: false)
%
% Outputs:
%   particles_full          : 1×N struct array (per-slice results)
%   particles_combined_full : 1×1 struct with fields vertically concatenated
%
% Requires helper functions on path:
%   numSlices, sliceOffsets, split_and_jitter_particles, collapse_struct_array

    if nargin < 12 || isempty(keep_in_slice), keep_in_slice = true; end
    if nargin < 13 || isempty(verbose),       verbose = false;      end

    % --- How many jitters to create ---
    totalSlices = numSlices(slice_width);                 % count
    sliceOffsetsAll = sliceOffsets(totalSlices, slice_offset);  % 1×N offsets (rad)

    % --- Preallocate from your existing 1×1 struct template ---
    particles_full = repmat(particles_10, 1, totalSlices);
    particles_full(1) = particles_10;

    % --- Generate per-slice jittered particles ---
    baseSeed = rngSeed;
    for idx = 2:totalSlices
        currentSeed = baseSeed + (idx - 1);               % unique seed per slice
        slice_offset = sliceOffsetsAll(idx);

        particles_full(idx) = split_and_jitter_particles( ...
            particles, clonesPer, pct_radial, pct_vertical, currentSeed, true, ...
            angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice);
    end

    % --- Divide mass by the number of slices (per-struct element) ---
    if isfield(particles_full, 'mass')
        for k = 1:totalSlices
            particles_full(k).mass = particles_full(k).mass / totalSlices;
        end
    else
        warning('Field "mass" not found in particles_full; skipping mass normalization.');
    end

    % --- Collapse 1×N struct array -> 1×1 struct by vertical concatenation ---
    particles_combined_full = collapse_struct_array(particles_full);   % all fields

    % --- Optional sanity checks ---
    if verbose
        try
            fprintf('Original N = %d, After first cloning (template) = %d\n', ...
                size(particles.pos,1), size(particles_10.pos,1));
            disp(particles_combined_full.vel(1:min(5,size(particles_combined_full.vel,1)),:));
        catch
            % keep silent if fields not present
        end
    end
end
