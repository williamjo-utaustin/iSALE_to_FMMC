clc;
clear all; 

% Locate repo root and add src to path (recursive)
thisFile = mfilename('codes');
repoRoot = fileparts(thisFile);                  % folder containing this script
srcDir   = fullfile(repoRoot, 'codes');            % adjust if your folder is named differently

if ~isfolder(srcDir)
    error('Could not find src folder at: %s', srcDir);
end

addpath(srcDir);
addpath(genpath(srcDir));    % include any subfolders under src
rehash;                      % refresh function/cache

% If the function is inside a package folder like +utils/open_from_folder.m,
% call it as utils.open_from_folder(...)

folder = fullfile(pwd, 'input', 'upright_porosity_0');
[lower_lampshade, upper_spike, stats] = open_from_folder(folder);

% Or override defaults:
% [lower_lampshade, upper_spike, stats] = open_from_folder(folder, 0, true, setdiff(1:8,4), 0:99);

%% 

% keep only rows where the 4th column > 200
lower_lampshade_nsp2 = lower_lampshade(lower_lampshade(:,4) >= 100 & lower_lampshade(:,5) < 2500, :);

% optional: show how many rows survived
fprintf('Filtered down to %d rows (from %d)\n', size(lower_lampshade,1), size(lower_lampshade,1));

%% 

% ---- choose your slice here (no global) ----
slice_width  = pi/8;   % wedge size
slice_offset = 0;       % start angle (rotate wedge if you want)

% Assuming you already have filtered_data
% Choose which columns hold radial & vertical speeds (edit these!)
col_radial   = 3;   % e.g., column with radial (horizontal) speed
col_vertical = 4;   % e.g., column with vertical speed

particles = make_particles_from_data(lower_lampshade_nsp2, col_radial, col_vertical, ...
                                     42, slice_width, slice_offset);
% Quick sanity check
N = size(lower_lampshade_nsp2,1);
fprintf('Created %d particles. Example velocity rows:\n', N);
disp(particles.vel(1:min(N,5), :))   % show first few [vx vy vz]

% for each particle created, create 10 sub particles where velocities are
% jittered from the actual launch velocity by 1%, 5% idk we shall see

% Clone & jitter (keep clones inside the same slice)
clonesPer     = 5;
pct_radial    = 5; 
pct_vertical  = 5;     
rngSeed       = 123;
angJitterDeg  = 3;      % small angle jitter within slice
phiUniform    = false;  % set true for uniform φ within slice (max spread)
keep_in_slice = true;

% Jitter particles to maintain spread
particles_10 = split_and_jitter_particles( ...
    particles, clonesPer, pct_radial, pct_vertical, rngSeed, true, ...
    angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice);

% --- How many jitters to create ---
totalSlices = numSlices(slice_width);                 % count
sliceOffsetsAll = sliceOffsets(totalSlices, slice_offset);  % 1×N offsets (rad)

% --- Preallocate from your existing 1×1 struct template ---
particles_full = repmat(particles_10, 1, totalSlices);
particles_full(1) = particles_10;

% --- Generate per-slice jittered particles ---
baseSeed = rngSeed;                                   % keep your original seed
for idx = 2:totalSlices
    currentSeed = baseSeed + (idx - 1);               % unique seed per slice
    slice_offset = sliceOffsetsAll(idx);              % DO NOT rename this var
                                                     % (function expects it)
    particles_full(idx) = split_and_jitter_particles( ...
        particles, clonesPer, pct_radial, pct_vertical, currentSeed, true, ...
        angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice);
end

% --- Collapse 1×N struct array -> 1×1 struct by vertical concatenation ---
particles_combined_full = collapse_struct_array(particles_full);   % all fields
% or: particles_combined_full = collapse_struct_array(particles_full, {'pos','vel','phi'});


% Sanity checks
fprintf('Original N = %d, After cloning = %d\n', size(particles.pos,1), size(particles_10.pos,1));
disp(particles_10.vel(1:5,:));   % peek first few jittered velocities

%% 
% --- simulate from t=0 to t=250 s, output every 1 s ---
dt = 0.1;                % seconds
steps = 2500;           % so t = 0:1:250
g = 1.62;              % m/s^2
recordTraj = true;

% Record every 1 s from 0..250
%times_list = 240:1:250;

% OR: every 1 s 0..250 plus every 0.5 s 230..250
times_list = unique([0:50:250, 240:1:250]);


[particles,   traj,   landed_mask,   meta ] = move_particles(particles,   dt, steps, g, true, times_list);
[particles_10,traj_10,landed_mask_10,meta2] = move_particles(particles_10, dt, steps, g, true, times_list);
[particles_combined_full, traj_combined_full, landed_mask_combined_full, meta3] = move_particles(particles_combined_full, dt, steps, g, true, times_list);

%% Plot the slice:

opts_10 = struct('show_base', false, 'show_overlay', true, ...
              'col_overlay', [0.7 0.7 0.7], 'sz_overlay', 2.5, ...
              'alpha_overlay', 0.5, 'viewAZEL', [0 0], 'zcap', [], ... % [] => auto vertical range
              'frameTimesSec',  meta.frame_times_sec, ...
              'frameTimesSec2', meta2.frame_times_sec);

% e.g., 50..250 every 50 s
plot_plumes_times(traj, traj_10, dt, struct('start',50,'stop',250,'step',50), ...
                  landed_mask, landed_mask_10, opts_10);

% or a custom list:
plot_plumes_times(traj, traj_10, dt, [240 242 244 246 248 250], ...
                   landed_mask, landed_mask_10, opts_10);

%% Plot the full plume

opts_10 = struct('show_base', false, 'show_overlay', true, ...
              'col_overlay', [0.7 0.7 0.7], 'sz_overlay', 2.5, ...
              'alpha_overlay', 0.5, 'viewAZEL', [0 0], 'zcap', [], ... % [] => auto vertical range
              'frameTimesSec',  meta.frame_times_sec, ...
              'frameTimesSec2', meta3.frame_times_sec);

plot_plumes_times(traj, traj_combined_full, dt, [240 242 244 246 248 250], ...
                   landed_mask, landed_mask_combined_full, opts);

plot_plumes_times(traj, traj_combined_full, dt, struct('start',0,'stop',250,'step',50), ...
                  landed_mask, landed_mask_combined_full, opts);



