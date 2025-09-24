function particles = make_particles_from_data(D, rngSeed, slice_width, slice_offset)
% MAKE_PARTICLES_FROM_DATA
%   particles = make_particles_from_data(D, col_radial, col_vertical, rngSeed, slice_width, slice_offset)
%   D            : filtered_data (N x K)
%   col_radial   : column index for radial (horizontal) speed [m/s]
%   col_vertical : column index for vertical speed [m/s]
%   rngSeed      : (optional) integer for reproducible random azimuths
%   slice_width  : (optional) wedge width [rad], default = 2*pi (full circle)
%   slice_offset : (optional) wedge start angle [rad], default = 0
%
% Output
%   particles.pos : (N x 3), all zeros (launch from origin)
%   particles.vel : (N x 3), [vx vy vz]
%   particles.phi : (N x 1), azimuths within [slice_offset, slice_offset + slice_width)


    arguments
        D double
        rngSeed (1,1) {mustBeInteger} = []
        slice_width (1,1) double {mustBePositive} = 2*pi
        slice_offset (1,1) double = 0

    end

    if ~isempty(rngSeed), rng(rngSeed); end

    col_radial   = 3;   % e.g., column with radial (horizontal) speed
    col_vertical = 4;   % e.g., column with vertical speed
    col_mass = 6;       % e.g., column with mass

    N  = size(D,1);
    vr = D(:, col_radial);   % horizontal speed (can be signed)
    vy = D(:, col_vertical); % vertical speed

    % Random azimuths uniformly within the slice
    phi = slice_offset + slice_width * rand(N,1);

    % Carry sign of vr into horizontal direction
    sgn = sign(vr);  sgn(sgn==0) = 1;
    mag = abs(vr);

    % 3D velocity: horizontal in X–Z plane, vertical in Y
    vx = sgn .* mag .* cos(phi);
    vz = sgn .* mag .* sin(phi);

    pos = zeros(N,3);
    vel = [vx, vy, vz];

    mass = D(:,col_mass);

    particles = struct('pos', pos, 'vel', vel, 'phi', phi, 'mass', mass);
end
