function [traj_full, landed_full, frameTimesSec_out] = tile_traj_around_y( ...
    traj_slice, landed_slice, slice_width, frameTimesSec_in, slice_offset, anchor, jitter_deg)
% TILE_TRAJ_AROUND_Y  Tile a trajectory wedge around vertical (y) to fill 2π with alignment.
%
% Inputs
%   traj_slice        : (T x N x 3) positions [x y z] for one wedge over time
%   landed_slice      : (T x N) logical (or numeric) mask; [] allowed (treated as all false)
%   slice_width       : wedge width in radians, 0 < slice_width <= 2*pi
%   frameTimesSec_in  : (optional) 1xM recorded times (s) aligned with T rows; [] ok
%   slice_offset      : (optional) base rotation (radians), default 0
%   anchor            : (optional) 'edge' (default) or 'center'
%   jitter_deg        : (optional) small per-wedge random jitter in degrees (default 0)
%
% Outputs
%   traj_full         : (T x (N*S) x 3) positions tiled into S wedges
%   landed_full       : (T x (N*S)) logical
%   frameTimesSec_out : passthrough of frameTimesSec_in (or [] if not provided)

    if ndims(traj_slice) ~= 3 || size(traj_slice,3) ~= 3
        error('traj_slice must be (T x N x 3).');
    end
    if nargin < 2 || isempty(landed_slice)
        landed_slice = false(size(traj_slice,1), size(traj_slice,2));
    end
    if ~islogical(landed_slice), landed_slice = landed_slice ~= 0; end

    if nargin < 3 || isempty(slice_width) || ~isscalar(slice_width) || ~isfinite(slice_width) ...
            || ~(slice_width > 0 && slice_width <= 2*pi)
        error('slice_width must be a finite scalar with 0 < slice_width <= 2*pi.');
    end
    if nargin < 4, frameTimesSec_in = []; end
    if nargin < 5 || isempty(slice_offset), slice_offset = 0; end
    if nargin < 6 || isempty(anchor), anchor = 'edge'; end
    if nargin < 7 || isempty(jitter_deg), jitter_deg = 0; end

    [T, N, ~] = size(traj_slice);

    % Compute an integer number of wedges and use the *effective* width to avoid drift
    S = max(1, round((2*pi) / slice_width));
    dphi = 2*pi / S;  % effective slice width that exactly tiles the circle

    % Pre-allocate outputs
    traj_full   = zeros(T, N*S, 3, 'like', traj_slice);
    landed_full = false(T, N*S);

    % Split source once
    X = traj_slice(:,:,1);   % (T x N)
    Y = traj_slice(:,:,2);
    Z = traj_slice(:,:,3);

    % Optional per-wedge jitter (in radians)
    jitter_rad = deg2rad(jitter_deg);

    for s = 0:S-1
        % Base wedge angle (aligned to edges)
        phi = slice_offset + s * dphi;

        % If the input slice is centered at 0 (i.e., spans [-dphi/2, +dphi/2]),
        % shift so that the source "center" maps to each wedge center.
        if strcmpi(anchor, 'center')
            phi = phi + dphi/2;
        end

        % Optional tiny dither (kept very small if used)
        if jitter_rad > 0
            phi = phi + (2*rand-1) * jitter_rad;  % uniform in [-jitter, +jitter]
        end

        c  = cos(phi);
        si = sin(phi);

        % Rotate about y: [x', z'] = [x c + z si,  -x si + z c]
        Xr =  X*c + Z*si;
        Zr = -X*si + Z*c;

        cols = (s*N+1):((s+1)*N);
        traj_full(:, cols, 1) = Xr;
        traj_full(:, cols, 2) = Y;
        traj_full(:, cols, 3) = Zr;

        landed_full(:, cols) = landed_slice;
    end

    if nargout >= 3
        frameTimesSec_out = frameTimesSec_in;  % passthrough
    end
end












% function [traj_full, landed_full] = tile_traj_around_y(traj_slice, landed_slice, slice_width)
% % TILE_TRAJ_AROUND_Y  Rotate a simulated wedge/slice around +y to cover 2π.
% %   Inputs:
% %     traj_slice   : (T x N x 3), positions [x y z] over time for one slice
% %     landed_slice : (T x N), logical landed tracker for that slice
% %     slice_width  : nominal wedge width [rad], e.g., pi/16
% %   Outputs:
% %     traj_full    : (T x (N*S) x 3) tiled positions
% %     landed_full  : (T x (N*S))     tiled landed mask
% 
%     if nargin < 3 || isempty(slice_width), slice_width = pi/16; end
%     [T, N, ~] = size(traj_slice);
% 
%     % exact tiling so wedges butt perfectly
%     S    = max(1, round((2*pi) / slice_width));
%     dphi = 2*pi / S;
% 
%     traj_full   = zeros(T, N*S, 3, 'like', traj_slice);
%     if ~isempty(landed_slice)
%         landed_full = false(T, N*S);
%     else
%         landed_full = [];
%     end
% 
%     % pull once
%     X = squeeze(traj_slice(:,:,1));   % (T x N)
%     Y = squeeze(traj_slice(:,:,2));   % (T x N)
%     Z = squeeze(traj_slice(:,:,3));   % (T x N)
% 
%     for s = 0:S-1
%         c = cos(s*dphi);  sn = sin(s*dphi);
% 
%         % rotate around +y: [x'; z'] = [ c  sn; -sn  c ] * [x; z]
%         Xr =  X*c + Z*sn;     % (T x N)
%         Zr = -X*sn + Z*c;     % (T x N)
% 
%         idx = (s*N + 1) : ((s+1)*N);
%         traj_full(:, idx, 1) = Xr;
%         traj_full(:, idx, 2) = Y;     % y unchanged
%         traj_full(:, idx, 3) = Zr;
% 
%         if ~isempty(landed_slice)
%             landed_full(:, idx) = landed_slice;
%         end
%     end
% end
