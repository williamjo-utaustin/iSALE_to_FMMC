function particles_out = split_and_jitter_particles( ...
    particles_in, clonesPer, pct_radial, pct_vertical, rngSeed, clampPositive, ...
    angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice)

% SPLIT_AND_JITTER_PARTICLES
% Clone and jitter velocities (radial magnitude + azimuth + vertical),
% optionally constraining azimuths to a given wedge (slice).
% Also clones mass and divides it equally among clones (mass conservation).
%
% Inputs
%   particles_in.pos : (N x 3)
%   particles_in.vel : (N x 3) [vx vy vz]
%   particles_in.mass: (N x C) scalar mass per particle (C>=1; commonly C=1)
%   clonesPer        : clones per particle (default 10)
%   pct_radial       : % 1σ multiplicative jitter for horizontal speed (default 5)
%   pct_vertical     : % 1σ multiplicative jitter for vy (default 5)
%   rngSeed          : optional integer seed (default [])
%   clampPositive    : clamp multiplicative scales to >= eps (default true)
%   angJitterDeg     : 1σ azimuth jitter in degrees (default 2) if phiUniform=false
%   phiUniform       : true => redraw φ ~ U(slice), ignoring original φ (default false)
%   slice_width      : wedge width [rad] (default 2*pi for no constraint)
%   slice_offset     : wedge start angle [rad] (default 0)
%   keep_in_slice    : true => wrap φ back into [offset, offset+width) (default true)
%
% Output
%   particles_out.pos, .vel, .mass, .parent_idx

    if nargin < 2 || isempty(clonesPer),     clonesPer = 10; end
    if nargin < 3 || isempty(pct_radial),    pct_radial = 5; end
    if nargin < 4 || isempty(pct_vertical),  pct_vertical = 5; end
    if nargin < 5, rngSeed = []; end
    if nargin < 6 || isempty(clampPositive), clampPositive = true; end
    if nargin < 7 || isempty(angJitterDeg),  angJitterDeg = 2; end
    if nargin < 8 || isempty(phiUniform),    phiUniform = false; end
    if nargin < 9 || isempty(slice_width),   slice_width = 2*pi; end
    if nargin < 10 || isempty(slice_offset), slice_offset = 0; end
    if nargin < 11 || isempty(keep_in_slice), keep_in_slice = true; end

    if ~isempty(rngSeed), rng(rngSeed); end

    % --- Required inputs ---
    pos = particles_in.pos;   % N x 3
    vel = particles_in.vel;   % N x 3
    if ~isfield(particles_in,'mass')
        error('particles_in.mass is required.');
    end
    mass = particles_in.mass; % N x C (C>=1)

    [N, ~] = size(pos);

    % --- Repeat rows for clones ---
    pos_rep  = repelem(pos,  clonesPer, 1);
    vel_rep  = repelem(vel,  clonesPer, 1);
    mass_rep = repelem(mass, clonesPer, 1) / clonesPer;   % split mass equally

    % --- Work with velocities in horizontal polar coords ---
    vx = vel_rep(:,1);
    vy = vel_rep(:,2);
    vz = vel_rep(:,3);

    r   = hypot(vx, vz);
    phi = atan2(vz, vx);

    % --- Jitter scales (multiplicative) ---
    sr = pct_radial   / 100;
    sz = pct_vertical / 100;
    fr = 1 + sr * randn(N*clonesPer, 1);   % radial magnitude scale
    fz = 1 + sz * randn(N*clonesPer, 1);   % vertical scale

    if clampPositive
        tiny = eps;
        fr = max(fr, tiny);
        fz = max(fz, tiny);
    end

    % --- Azimuth update ---
    if phiUniform
        phi2 = slice_offset + slice_width * rand(N*clonesPer,1);
    else
        sigma_phi = deg2rad(angJitterDeg);
        phi2 = phi + sigma_phi * randn(N*clonesPer,1);
        if keep_in_slice && slice_width < 2*pi
            % wrap into the slice interval
            phi2 = slice_offset + mod(phi2 - slice_offset, slice_width);
        end
    end

    % --- Apply magnitude/vertical jitters ---
    r2  = r  .* fr;
    vy2 = vy .* fz;

    % --- Back to Cartesian ---
    vx2 = r2 .* cos(phi2);
    vz2 = r2 .* sin(phi2);

    vel_out      = vel_rep;
    vel_out(:,1) = vx2;
    vel_out(:,2) = vy2;
    vel_out(:,3) = vz2;

    % --- Assemble output ---
    particles_out.pos        = pos_rep;
    particles_out.vel        = vel_out;
    particles_out.mass       = mass_rep;                     % cloned + scaled
    particles_out.parent_idx = repelem((1:N).', clonesPer, 1);

    % (Optional) quick mass-conservation sanity check:
    % if abs(sum(mass,'all') - sum(particles_out.mass,'all')) > 1e-9
    %     warning('Total mass changed by %.3g', ...
    %         sum(particles_out.mass,'all') - sum(mass,'all'));
    % end
end





% function particles_out = split_and_jitter_particles( ...
%     particles_in, clonesPer, pct_radial, pct_vertical, rngSeed, clampPositive, ...
%     angJitterDeg, phiUniform, slice_width, slice_offset, keep_in_slice)
% 
% % SPLIT_AND_JITTER_PARTICLES
% % Clone and jitter velocities (radial magnitude + azimuth + vertical),
% % optionally constraining azimuths to a given wedge (slice).
% %
% % Inputs
% %   particles_in.pos : (N x 3)
% %   particles_in.vel : (N x 3) [vx vy vz]
% %   clonesPer        : clones per particle (default 10)
% %   pct_radial       : % 1σ multiplicative jitter for horizontal speed (default 5)
% %   pct_vertical     : % 1σ multiplicative jitter for vy (default 5)
% %   rngSeed          : optional integer seed (default [])
% %   clampPositive    : clamp multiplicative scales to >= eps (default true)
% %   angJitterDeg     : 1σ azimuth jitter in degrees (default 2) if phiUniform=false
% %   phiUniform       : true => redraw φ ~ U(slice), ignoring original φ (default false)
% %   slice_width      : wedge width [rad] (default 2*pi for no constraint)
% %   slice_offset     : wedge start angle [rad] (default 0)
% %   keep_in_slice    : true => wrap φ back into [offset, offset+width) (default true)
% %
% % Output
% %   particles_out.pos, .vel, .parent_idx
% 
%     if nargin < 2 || isempty(clonesPer),     clonesPer = 10; end
%     if nargin < 3 || isempty(pct_radial),    pct_radial = 5; end
%     if nargin < 4 || isempty(pct_vertical),  pct_vertical = 5; end
%     if nargin < 5, rngSeed = []; end
%     if nargin < 6 || isempty(clampPositive), clampPositive = true; end
%     if nargin < 7 || isempty(angJitterDeg),  angJitterDeg = 2; end
%     if nargin < 8 || isempty(phiUniform),    phiUniform = false; end
%     if nargin < 9 || isempty(slice_width),   slice_width = 2*pi; end
%     if nargin < 10 || isempty(slice_offset), slice_offset = 0; end
%     if nargin < 11 || isempty(keep_in_slice), keep_in_slice = true; end
% 
%     if ~isempty(rngSeed), rng(rngSeed); end
% 
%     pos = particles_in.pos;   % N x 3
%     vel = particles_in.vel;   % N x 3
%     [N, ~] = size(pos);
% 
%     % Repeat rows for clones
%     pos_rep = repelem(pos, clonesPer, 1);
%     vel_rep = repelem(vel, clonesPer, 1);
% 
%     vx = vel_rep(:,1);
%     vy = vel_rep(:,2);
%     vz = vel_rep(:,3);
% 
%     % Horizontal polar
%     r   = hypot(vx, vz);
%     phi = atan2(vz, vx);
% 
%     % Jitter scales (multiplicative)
%     sr = pct_radial   / 100;
%     sz = pct_vertical / 100;
%     fr = 1 + sr * randn(N*clonesPer, 1);   % radial magnitude
%     fz = 1 + sz * randn(N*clonesPer, 1);   % vertical
% 
%     if clampPositive
%         tiny = eps;
%         fr = max(fr, tiny);
%         fz = max(fz, tiny);
%     end
% 
%     % Azimuth update
%     if phiUniform
%         phi2 = slice_offset + slice_width * rand(N*clonesPer,1);
%     else
%         sigma_phi = deg2rad(angJitterDeg);
%         phi2 = phi + sigma_phi * randn(N*clonesPer,1);
%         if keep_in_slice && slice_width < 2*pi
%             % wrap into the slice interval
%             phi2 = slice_offset + mod(phi2 - slice_offset, slice_width);
%         end
%     end
% 
%     % Apply magnitude/vertical jitters
%     r2  = r  .* fr;
%     vy2 = vy .* fz;
% 
%     % Back to Cartesian
%     vx2 = r2 .* cos(phi2);
%     vz2 = r2 .* sin(phi2);
% 
%     vel_out           = vel_rep;
%     vel_out(:,1)      = vx2;
%     vel_out(:,2)      = vy2;
%     vel_out(:,3)      = vz2;
% 
%     particles_out.pos         = pos_rep;
%     particles_out.vel         = vel_out;
%     particles_out.parent_idx  = repelem((1:N).', clonesPer, 1);
% end
