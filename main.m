clc;
clear all; 
%% 

% ===========================
%  User-configurable parameters
%  ===========================

% --- General ---
plot_on           = true;
recordTraj        = true;     % save trajectories/history at times in times_list
rngSeed           = 123;      % RNG seed for any stochastic steps (set [] to skip)
pe = 250;  % progressEvery
if ~isempty(rngSeed), rng(rngSeed); end

% --- Dataset selection ---
orientation       = 'upright';    % 'upright', 'flipped', etc.
porosity_level    = 0;            % integer X in "<orientation>_porosity_X"

% --- Ejecta filtering ---
min_ejecta_velocity = 100;        % m/s, vertical ejecta velocity cutoff

% --- Slice geometry (polar wedge) ---
slice_width       = pi/32;        % radians; wedge angular width
slice_offset      = 0;            % radians; wedge start angle (rotate as needed)

% --- Cloning & jitter (clones stay within the same slice) ---
clonesPer         = 10;           % number of clones per source row
pct_radial        = 5;            % ±% jitter of radial component
pct_vertical      = 5;            % ±% jitter of vertical component
angJitterDeg      = 3;            % degrees; small angular jitter within slice
phiUniform        = false;        % true -> uniform φ within slice (max spread)
keep_in_slice     = true;         % enforce clones remain inside the wedge

% --- Dynamics / integration ---
dt                = 0.1;          % s, timestep
steps             = 2500;         % number of steps (dt*steps = total simulated time)
g                 = 1.62;         % m/s^2, gravitational acceleration
t_end_sec         = dt * steps;   % derived: total simulated duration

% --- Recording schedule ---
% Option A: dense sampling every 1 s from 0..250
% times_list = 0:1:250;

% Option B: sparse 0:50:250 plus dense 240:1:250 (current choice)
times_list        = unique([0:50:250, 240:1:250]);

% --- Material / grain parameters ---
grain_size_in_microns = 15;       % μm
grain_density         = 3100;     % kg/m^3 (material bulk density)

%%

% Find the repository root and add the `codes` directory (recursively) to the MATLAB path
thisFile = mfilename('codes');
repoRoot = fileparts(thisFile); % absolute path to the folder containing this script
srcDir   = fullfile(repoRoot, 'codes'); % change this if your source folder has a different name

% If folder does not exist, state the error and exit
if ~isfolder(srcDir)
    error('Could not find src folder at: %s', srcDir);
end

% Add path to the source folder
addpath(srcDir);
addpath(genpath(srcDir)); % also include all subdirectories under `codes`
rehash; % refresh MATLAB's function cache so new paths are recognized

% If the function resides in a package folder (e.g.,+utils/open_from_folder.m),
% call it using the package-qualified name: utils.open_from_folder(...)
folder = fullfile(pwd, 'input', sprintf('%s_porosity_%d', orientation, ...
    porosity_level));
[lower_lampshade, upper_spike, stats] = open_from_folder(folder,...
    porosity_level, true, [], [], orientation);

% Keep ejecta data from the lampshade that is ejecting faster than 100 m/s
lower_lampshade_nsp2 = lower_lampshade( ...
    lower_lampshade(:,4) >= min_ejecta_velocity & ...
    lower_lampshade(:,5) < 2500, :);

% Build particle structs from filtered ejecta rows.
particles = make_particles_from_data(lower_lampshade_nsp2, 42, ... 
    slice_width, slice_offset);

% Sanity check: how many particles did we construct and what do velocities look like?
N = size(lower_lampshade_nsp2,1);
fprintf('Created %d particles. Example velocity rows:\n', N);
disp(particles.vel(1:min(N,5), :))   % show first few [vx vy vz]

% Jitter particles to maintain spread
particles_10 = split_and_jitter_particles( ...
    particles, clonesPer, pct_radial, pct_vertical, rngSeed, true, ...
    angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice);

% Generate the particles meant for a full plume
[particles_full, particles_combined_full] = generate_jitters_and_collapse( ...
    particles_10, particles, clonesPer, pct_radial, pct_vertical, ...
    rngSeed, angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice, true);

% Moving particels warning
disp("Moving Particles in Space (May take time)");

% Integrate particle motion for three populations. Each call advances states for `steps` steps of size `dt` under gravity `g`.
[particles,   traj,   landed_mask,    veloc,    meta ] = ...
    move_particles(particles,   dt, steps, g, recordTraj, times_list, ...
                   'progress', true, 'progressEvery', pe, 'label', 'pop tracer');

[particles_10, traj_10, landed_mask_10, veloc_10, meta2] = ...
    move_particles(particles_10, dt, steps, g, recordTraj, times_list, ...
                   'progress', true, 'progressEvery', pe, 'label', 'pop slice');

[particles_combined_full, traj_combined_full, landed_mask_combined_full, veloc_combined_full, meta3] = ...
    move_particles(particles_combined_full, dt, steps, g, recordTraj, times_list, ...
                   'progress', true, 'progressEvery', pe, 'label', 'pop combined');


% Convenience handle: per-particle masses for the combined population
mass_particle_i = particles_combined_full.mass;   % 1×N or N×1 vector

% Convert each particle's mass to an equivalent grain count, given a grain size (μm) 
% and material density. Returns one grain count per particle.
num_grains_in_particle_i = particle2grain(particles_combined_full.mass, ...
                                          grain_size_in_microns, grain_density);
disp("Completed Moving Particles. Saving to File");

% Save trajectory outputs to a single mat file
save('output/ejecta_statistics.mat', 'times_list', 'traj_combined_full', ...
    'veloc_combined_full', 'mass_particle_i', 'num_grains_in_particle_i');

opts_10 = struct('show_base', false, 'show_overlay', true, ...
    'col_overlay', [0 0 0], 'sz_overlay', 2.5, ...
    'alpha_overlay', 0.5, 'viewAZEL', [0 0], 'zcap', [], ... % [] => auto vertical range
    'frameTimesSec',  meta.frame_times_sec, ...
    'frameTimesSec2', meta2.frame_times_sec);

opts_full = struct('show_base', false, 'show_overlay', true, ...
    'col_overlay', [0 0 0], 'sz_overlay', 2.5, ...
    'alpha_overlay', 0.5, 'viewAZEL', [0 0], 'zcap', [], ... % [] => auto vertical range
    'frameTimesSec',  meta.frame_times_sec, ...
    'frameTimesSec2', meta3.frame_times_sec);

%% ===========================
%  Plume snapshots (base vs overlays)
%  traj, traj_10, traj_combined_full have size [T x N x 3] = [time x particles x {x,y,z}]
%  dt is the integration step (s); times below are in SECONDS from t=0.
%  landed_mask* are logical vectors (1 x N) indicating particles that have landed.
%  opts_* controls visuals (overlay on/off, color, size, alpha, view, z-limits)
%  and should include frameTimesSec / frameTimesSec2 to map seconds -> frame indices.
%  ===========================

% 1) Base vs ×10 jitter overlay — evenly spaced snapshots from 0..250 s every 50 s
%    Shows how the ×10 jitter population (traj_10) compares to the base plume (traj).
plot_plumes_times( ...
    traj, traj_10, dt, struct('start',0,'stop',250,'step',50), ...
    landed_mask, landed_mask_10, opts_10);

% 2) Base vs ×10 jitter overlay — late-stage snapshots at specific times (denser sampling)
%    Focus on the final jitter 10 seconds to see late-time plume structure.
plot_plumes_times( ...
    traj, traj_10, dt, [240 242 244 246 248 250], ...
    landed_mask, landed_mask_10, opts_10);

% 3) Base vs COMBINED overlay — evenly spaced snapshots 0..250 s every 50 s
%    A broad, uniform sweep over the whole simulation window.
plot_plumes_times( ...
    traj, traj_combined_full, dt, struct('start',0,'stop',250,'step',50), ...
    landed_mask, landed_mask_combined_full, opts_full);

% 4) Base vs COMBINED overlay — late-stage snapshots at specific times
%    Compare the base plume to your combined population near the end of the run.
plot_plumes_times( ...
    traj, traj_combined_full, dt, [240 242 244 246 248 250], ...
    landed_mask, landed_mask_combined_full, opts_full);

% ==========================================
%  Velocity-colored snapshots (by vertical v)
%  veloc* are velocity histories with size [T x N x 3] = [time x particles x {vx,vy,vz}]
%  The plotter uses vz (vertical component) to color/encode particles at the requested times.
%  Landed masks can be used to hide/mark particles that have already landed.
%  ==========================================

% 5) Base vs COMBINED — explicit list of late times, colored by vertical velocity (vz)
plot_plume_times_by_vertical_velocity( ...
    traj, traj_combined_full, ...
    veloc, veloc_combined_full, ...   % (T x N x 3) velocities aligned with traj histories
    dt, [240 242 244 246 248 250], ...
    landed_mask, landed_mask_combined_full, ...
    opts_full);

% 6) Base vs COMBINED — uniform sweep 0..250 s every 50 s, colored by vz
plot_plume_times_by_vertical_velocity( ...
    traj, traj_combined_full, ...
    veloc, veloc_combined_full, ...
    dt, struct('start',0,'stop',250,'step',50), ...
    landed_mask, landed_mask_combined_full, ...
    opts_full);

%% Build the grid
% ----- domain limits -----
r_max = 2E5;        % 0 .. 2×10^5 (radius is nonnegative)
z_max = 1E5;        % 0 .. 10×10^4

% ----- resolution (edit these to taste) -----
Nr = 15;           % number of radial samples
Ntheta = 64;       % angular samples (0..360° by 1°)
Nz = 50;           % vertical samples

% ----- coordinate vectors (EDGES) -----
r = linspace(0, r_max, Nr+1);
theta = linspace(0, 2*pi, Ntheta+1);     % sweep full circle
z = linspace(0, z_max, Nz+1);

% ----- cylindrical grid (optional) -----
[R, TH, Z] = ndgrid(r, theta, z);        % R>=0, TH in [0,2π), Z in [0,z_max]
X = R .* cos(TH);
Y = R .* sin(TH);

% ----- use edge vectors explicitly -----
r_edges     = r;
theta_edges = theta;
z_edges     = z;

% Make the last theta edge slightly open to avoid ambiguity at exactly 2π
lastTheta = numel(theta_edges);
theta_edges(lastTheta) = theta_edges(lastTheta) + eps;

% Sizes (avoid 'end' keyword)
Nr  = numel(r_edges)     - 1;
Nth = numel(theta_edges) - 1;
Nz  = numel(z_edges)     - 1;


% Notes
%   - Plots with axes [X,Y,Z]=[x,z,y] so physical y is vertical.

% choose timestep (16 x 102960 x 3)
timestep_index = 16;

Xp = traj_combined_full(timestep_index, :, 1);
Yp = traj_combined_full(timestep_index, :, 2);
Zp = traj_combined_full(timestep_index, :, 3);


% plot a quick Xp, Yp, and Zp here 
   
% --- Quick sanity plots for Xp, Yp, Zp (no wireframe) ---

% Make sure they are column vectors and finite
Xp = Xp(:); Yp = Yp(:); Zp = Zp(:);
valid = isfinite(Xp) & isfinite(Yp) & isfinite(Zp);
Xp = Xp(valid); Yp = Yp(valid); Zp = Zp(valid);

% Optional subsampling for speed if very large
Nmax = 1e5;
if numel(Xp) > Nmax
    idx = randperm(length(Xp), Nmax);
    Xs = Xp(idx); Ys = Yp(idx); Zs = Zp(idx);
else
    Xs = Xp; Ys = Yp; Zs = Zp;
end

% --- Map particles to cell indices (ir,it,iz) using EDGES ---
[theta_p, r_p] = cart2pol(Xp, Zp);      % theta in [-π, π]
theta_p = mod(theta_p, 2*pi);           % wrap to [0, 2π)

ir = discretize(r_p,     r_edges);      % 1..Nr, NaN if out
it = discretize(theta_p, theta_edges);  % 1..Nth
iz = discretize(Yp,      z_edges);      % 1..Nz

inside = ~(isnan(ir) | isnan(it) | isnan(iz));

% --- Weighted counts by # of grains per particle ---
% Make sure the weights line up with the particles at this timestep.
% (Assumes num_grains_in_particle_i is a 1×N or N×1 vector for the same particles.)
w_all = num_grains_in_particle_i(:);      % column vector
%w_all = mass_particle_i(:);      % column vector
w_in  = w_all(inside(:));                  % weights of inside-domain particles

% Subscripts for inside particles (you already built these)
ir_in = ir(inside); it_in = it(inside); iz_in = iz(inside);
ir_in = ir_in(:);   it_in = it_in(:);   iz_in = iz_in(:);

subs = [ir_in, it_in, iz_in];             % M×3

% Grain-weighted counts per cell (sum of grains in each cell)
counts = accumarray(subs, w_in, [Nr, Ntheta, Nz], @sum, 0);

% (Optional) unweighted particle counts, if you still want them
% counts = accumarray(subs, 1, [Nr, Ntheta, Nz], @sum, 0);

% counts(ir,it,iz) already built as in your snippet
totalCollected = sum(counts(:));     % total particles tallied in grid
Ninside        = nnz(inside);        % particles that were inside the domain
Ntotal         = numel(inside);      % total particles tested
Noutside       = Ntotal - Ninside;   % outside-domain particles

fprintf('sum(counts) = %d, inside = %d, outside = %d, total = %d\n', ...
        totalCollected, Ninside, Noutside, Ntotal);

% -------- cell volumes in cylindrical coords --------
% Use the original edges r, theta, z (not the version with theta(end)+eps) to keep exact 2π span.
r1 = r(1:numel(r)-1);       r2 = r(2:numel(r));
t1 = theta(1:numel(theta)-1); t2 = theta(2:numel(theta));
z1 = z(1:numel(z)-1);       z2 = z(2:numel(z));

% Radial integral ∫ r dr = 0.5 (r^2)
Ar   = 0.5*(r2.^2 - r1.^2);    % [Nr x 1]
Dth  = (t2 - t1);              % [Nth x 1]
Dz   = (z2 - z1);              % [Nz x 1]

% Reshape for broadcasting (works on all recent MATLAB; use bsxfun if needed)
Ar3  = reshape(Ar,  [Nr, 1,  1]);
Dth3 = reshape(Dth, [1,  Nth, 1]);
Dz3  = reshape(Dz,  [1,  1,  Nz]);

V = Ar3 .* Dth3 .* Dz3;        % [Nr x Nth x Nz] cell volumes

number_density = counts./V;

% Inputs assumed:
%   r           : edges, size [1 x (Nr+1)]
%   Nr,Ntheta,Nz
%   number_density : [Nr x Ntheta x Nz]  (per-cell, e.g., grains per unit volume)

% Radial cell thicknesses Δr_i
dr  = r(2:Nr+1) - r(1:Nr);          % [1 x Nr]
dr3 = reshape(dr, [Nr, 1, 1]);      % [Nr x 1 x 1] for broadcasting

% Column density integrating n over r: Σ_r(θ,z) = ∑ n(r_i,θ,z) * Δr_i
column_density_r = squeeze( sum(number_density .* dr3, 1) );   % [Ntheta x Nz]

% --- column_density_r is [Ntheta x Nz] and z are edges [1 x (Nz+1)] ---

% centers of z bins
z_c = 0.5*(z(1:Nz) + z(2:Nz+1));

% ===== Option A: use your endpoints (0 m -> 250 s, 5e4 m -> 220 s) =====
t_c = 250 + (220 - 250) * (z_c - 0) / (5e4 - 0);   % linear map z -> t
% implied speed:
v_implied = 5e4 / (250 - 220);                     % ≈ 1666.67 m/s

% mean and spread across theta
m  = mean(column_density_r, 1);        % [1 x Nz]
sd = 3*std(column_density_r, 0, 1);      % 1σ across theta
% If you prefer standard error: sd = sd / sqrt(size(column_density_r,1));

% sort time so it increases left->right
[tx, ix] = sort(t_c);
mx  = m(ix);
sdx = sd(ix);

f = gcf;                 % use current figure
f.Units = 'inches';
f.Position(3:4) = [8 12];  % width × height (inches)
f.Color = 'w';             % set background on the same figure

% h1 = errorbar(tx(1:end-1)+2.25, mx(1:end-1)*1E2/0.75, sdx(1:end-1)*1E2/0.75, 'LineWidth', 3, 'LineStyle', '-', 'Color', 'red');
% hold on;

x = tx(1:end-1) + 2.25;
y = mx(1:end-1) * 1e2 / 0.75;
e = sdx(1:end-1) * 1e2 / 0.75;

hold on

% light red shaded error band
hBand = fill([x, fliplr(x)], [y-e, fliplr(y+e)], 'r', ...
    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', '±3σ');
hold on;
% mean line on top (red)
h1 = plot(x, y, 'r-', 'LineWidth', 3, 'DisplayName', 'Mean');
hold on;
% (optional) keep the band behind the line
uistack(hBand, 'bottom');

% (optional) legend
% legend([h1 hBand], 'Mean', '±1σ', 'Location', 'best');



% [time, data, error]  ← only the first 3 columns from your table
D = [
-29.492  336217.3948   178639.2834
-28.312  427830.6753   210786.5796
-26.543  353172.2342   192821.2397
-24.773  140747.8798   120336.5868
-23.593   15403.28106   49307.30996
-22.414  166771.3987   128645.1124
-20.644  136965.4232   101064.5816
-20.054  242783.8207   152975.3506
-18.875   32781.62466   85144.34298
-17.695  167626.298    145452.5465
-16.515  323054.6262   169616.4049
-15.925  474117.5147   217807.9614
-14.746  369394.6268   187918.4673
-12.976  199263.4441   133027.9716
-11.797   94288.62967  107608.3662
-10.617   97972.41122   77862.0002
 -8.848  274152.4232   147932.8367
 -7.078  496104.5776   201101.2785
 -5.898  772191.6574   252061.0232
 -4.719 1198479.802    307579.8312
 -2.949 1908065.729    392319.1774
];

% Split into vectors (if you like)
t   = D(:,1);
y   = D(:,2);
err = D(:,3);

% Quick plot with error bars (optional):
h2 = errorbar(t+250, y, err, 'LineWidth', 3, 'LineStyle', '-', 'Color', 'black'); 

ax = gca;             % get current axes
ax.Box = 'on';        % ensure border is drawn
ax.LineWidth = 1.5;   % thicker border if you want
ax.XColor = 'k';      % black for x-axis (also sets border color)
ax.YColor = 'k';      % black for y-axis
ax.FontSize = 24;     % set tick label font size
legend([h1 h2], {'Simulation','NSP2 Data'}, 'Location','best');

xlabel('Time (s)');
ylabel('Column Density (#/m^2)'); 
xlim([220 250]);
ylim([0, 3E6]);

%%
% assume you already have r, theta, z (centers)
%plotCylWireframe(r, theta, z);  % defaults (blue rings, green radials, red verticals)

% or tweak density and colors:
plotCylWireframe(r, theta, z, ...
    'RadialEvery', 12, 'ThetaEvery', 24, 'ZEvery', 8, ...
    'Colors', struct('Rings',[0 0.4 1],'Radials',[0.1 0.7 0.2],'Verticals',[0.9 0.1 0.1]), ...
    'LineWidth', 0.8);

%% 
% 3D scatter: [X,Z,Y] so Y is vertical
figure('Color','w');
scatter3(Xs, Zs, Ys, 6, Ys, 'filled', ...
    'MarkerEdgeAlpha', 0.2, 'MarkerFaceAlpha', 0.6);
axis equal; grid on; box on
xlabel('x'); ylabel('z'); zlabel('y (vertical)');
title(sprintf('Particles at t = %d (N = %d)', timestep_index, numel(Xp)));
colorbar; view(35, 20);

% 2D projections
% figure('Color','w');
% tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
% 
% nexttile; plot(Xs, Zs, '.', 'MarkerSize', 2);
% axis equal; grid on; xlabel('x'); ylabel('z'); title('x–z');
% 
% nexttile; plot(Xs, Ys, '.', 'MarkerSize', 2);
% axis equal; grid on; xlabel('x'); ylabel('y'); title('x–y');
% 
% nexttile; plot(Zs, Ys, '.', 'MarkerSize', 2);
% axis equal; grid on; xlabel('z'); ylabel('y'); title('z–y');