function h = plot_plume_times_by_vertical_velocity(traj, traj_overlay, vel, vel_overlay, ...
    dt, timesSec, landed_mask, landed_mask_overlay, opts)
% PLOT_PLUME_TIMES_BY_VERTICAL_VELOCITY
% Plot particle clouds at specified times (<=6 panels) colored by vertical velocity (vy).
%
% h = plot_plume_times_by_vertical_velocity(traj, traj_overlay, vel, vel_overlay, ...
%       dt, timesSec, landed_mask, landed_mask_overlay, opts)
%
% Inputs
%   traj               : (T x N x 3) base positions [x y z]        (can be [])
%   traj_overlay       : (T2 x N2 x 3) overlay positions           (can be [])
%   vel                : (T x N x 3) base velocities [vx vy vz]    (can be [])
%   vel_overlay        : (T2 x N2 x 3) overlay velocities          (can be [])
%   dt                 : scalar time step (s)
%   timesSec           : numeric vector (e.g. [240 242 ...]) OR struct with fields .start .stop .step
%   landed_mask        : (T x N) logical for base (can be []; inferred by y==0)
%   landed_mask_overlay: (T2 x N2) logical for overlay (can be []; inferred by y==0)
%   opts (optional)    : struct
%       .show_base       (default: ~isempty(traj))
%       .show_overlay    (default: ~isempty(traj_overlay))
%       .show_plane      (default: true)
%       .sz_base         (default: 5)
%       .sz_overlay      (default: 2.5)
%       .alpha_overlay   (default: 0.5)
%       .col_landed      (default: [1 0 0])
%       .frameTimesSec   (default: []) % recorded times for traj
%       .frameTimesSec2  (default: []) % recorded times for traj_overlay
%       .vlim            (default: auto via robust 2–98 percentiles)
%       .cmap            (default: red->purple)
%       .viewAZEL        (default: [0 0])
%       .zcap            (default: []) % [] => auto vertical upper bound
%
% Output
%   h : struct of handles and meta (axes, tiles, indices, chosen times)
%
% Notes
%   - Plot axes are [X,Y,Z] = [x, z, y], so physical y is vertical (Z-axis).
%   - Color uses vertical velocity vy = vel(:,:,2).
%   - At most 6 times are shown (2x3). Extra requested times are dropped.

    if nargin < 9, opts = struct; end
    if isempty(dt) || ~isscalar(dt) || dt <= 0
        error('dt must be a positive scalar.');
    end

    % ---- defaults
    def.show_base      = ~isempty(traj);
    def.show_overlay   = ~isempty(traj_overlay);
    def.show_plane     = true;
    def.sz_base        = 5;
    def.sz_overlay     = 2.5;
    def.alpha_overlay  = 0.5;
    def.col_landed     = [1 0 0];
    def.frameTimesSec  = [];
    def.frameTimesSec2 = [];
    def.vlim           = [];                 % auto from data
    def.cmap           = redPurpleMap(256);  % red -> purple
    def.viewAZEL       = [0 0];
    def.zcap           = [];                 % [] => auto
    % merge defaults
    fn = fieldnames(def);
    for i = 1:numel(fn)
        f = fn{i}; if ~isfield(opts,f) || isempty(opts.(f)), opts.(f)=def.(f); end
    end

    % ---- derive available time (prefer explicitly recorded times)
    T  = size(traj,         1);
    T2 = size(traj_overlay, 1);
    ft  = opts.frameTimesSec(:).';
    ft2 = opts.frameTimesSec2(:).';

    if ~isempty(ft) || ~isempty(ft2)
        t_cap = max([0, ft, ft2]);
    else
        t_max  = (max(T,  1)-1)*dt;
        t_max2 = (max(T2, 1)-1)*dt;
        t_cap  = max([t_max, t_max2, 0]);
    end

    % ---- resolve list of requested seconds
    if isstruct(timesSec)
        s = timesSec;
        if ~all(isfield(s, {'start','stop','step'}))
            error('timesSec struct must have fields: start, stop, step');
        end
        ts = s.start : s.step : s.stop;
    else
        ts = timesSec(:).';
    end
    % keep valid times only, cap to 6
    ts = ts(ts >= 0 & ts <= t_cap);
    if isempty(ts)
        error('No valid times to plot within available trajectory (t_cap=%.3f s).', t_cap);
    end
    if numel(ts) > 6
        warning('plot_plume_times_by_vertical_velocity:TooManyTimes', ...
                'Requested %d times; displaying first 6 only.', numel(ts));
        ts = ts(1:6);
    end
    K = numel(ts);

    % ---- seconds -> frame indices (use recorded times if provided)
    idx = []; idx2 = [];
    if ~isempty(traj)
        if ~isempty(ft)
            idx = arrayfun(@(t) find(abs(ft - t) == min(abs(ft - t)), 1, 'first'), ts);
        else
            idx = 1 + floor(ts./dt + 1e-9);
        end
        idx = min(size(traj,1), max(1, idx));
    end
    if ~isempty(traj_overlay)
        if ~isempty(ft2)
            idx2 = arrayfun(@(t) find(abs(ft2 - t) == min(abs(ft2 - t)), 1, 'first'), ts);
        else
            idx2 = 1 + floor(ts./dt + 1e-9);
        end
        idx2 = min(size(traj_overlay,1), max(1, idx2));
    end

    % ---- gather positions & vertical velocities (force column vectors)
    poses   = cell(K,1);
    poses2  = cell(K,1);
    vvert   = cell(K,1);   % base vy
    vvert2  = cell(K,1);   % overlay vy
    for k = 1:K
        if ~isempty(traj)
            P = squeeze(traj(idx(k),:,:));     % [N x 3]
            poses{k} = P;
            if ~isempty(vel)
                % Always as column vector:
                vvert{k} = reshape(vel(idx(k),:,2), [], 1);
            end
        end
        if ~isempty(traj_overlay)
            Q = squeeze(traj_overlay(idx2(k),:,:));
            poses2{k} = Q;
            if ~isempty(vel_overlay)
                vvert2{k} = reshape(vel_overlay(idx2(k),:,2), [], 1);
            end
        end
    end

    % ---- bounds in plot coords [x z y] across whichever sets exist
    mins = [ inf  inf  inf];  maxs = [-inf -inf -inf];
    for k = 1:K
        stacks = [];
        if ~isempty(poses{k}),  P = poses{k};  stacks = [stacks; [P(:,1) P(:,3) P(:,2)]]; end
        if ~isempty(poses2{k}), Q = poses2{k}; stacks = [stacks; [Q(:,1) Q(:,3) Q(:,2)]]; end
        if ~isempty(stacks)
            mins = min(mins, min(stacks,[],1));
            maxs = max(maxs, max(stacks,[],1));
        end
    end
    pad = 0.05*(maxs - mins + eps);
    xlim_all = [mins(1)-pad(1), maxs(1)+pad(1)];
    ylim_all = [mins(2)-pad(2), maxs(2)+pad(2)];
    zmin_auto = max(0, mins(3)-pad(3));
    zmax_auto = maxs(3)+pad(3);
    if ~isempty(opts.zcap) && isfinite(opts.zcap)
        zmax_auto = min(zmax_auto, opts.zcap);
    end
    if zmax_auto <= zmin_auto, zmax_auto = zmin_auto + 1; end
    zlim_all = [zmin_auto, zmax_auto];

    % ---- global color limits (vy) for consistent scaling
    if ~isempty(opts.vlim)
        vmin = opts.vlim(1); vmax = opts.vlim(2);
    else
        vcat = [];
        for k = 1:K
            if ~isempty(vvert{k}),  vcat = [vcat; vvert{k}(:)];  end %#ok<AGROW>
            if ~isempty(vvert2{k}), vcat = [vcat; vvert2{k}(:)]; end %#ok<AGROW>
        end
        vcat = vcat(isfinite(vcat));   % drop NaN/Inf
        if isempty(vcat)
            vmin = -1; vmax = 1;       % fallback
        else
            lo = prctile(vcat, 2);
            hi = prctile(vcat, 98);
            if lo == hi, lo = lo - 1; hi = hi + 1; end
            vmin = lo; vmax = hi;
        end
    end

    % ---- figure (2x3) and tiles
    h.fig = figure('Name','Plumes by vertical velocity @ specified times', ...
                   'Units','normalized','Position',[0.04 0.06 0.92 0.84]);
    h.tl  = tiledlayout(2, 3, 'Padding','compact','TileSpacing','compact');

    h.base  = gobjects(K,1); h.over = gobjects(K,1);
    h.plane = gobjects(K,1); h.axes = gobjects(K,1);

    for k = 1:K
        nexttile;
        ax = gca; h.axes(k) = ax; hold(ax,'on');

        % Base
        if opts.show_base && ~isempty(poses{k})
            P = poses{k}; X = P(:,1); Yp = P(:,3); Z = P(:,2);
            if ~isempty(landed_mask)
                landed = squeeze(landed_mask(idx(k),:)).';
            else
                landed = (Z == 0);
            end
            air = ~landed;

            if any(air)
                if ~isempty(vvert{k})
                    C = vvert{k}(air);
                    h.base(k) = scatter3(X(air), Yp(air), Z(air), opts.sz_base, C, ...
                                         'filled','Parent',ax,'MarkerEdgeColor','none');
                    colormap(ax, opts.cmap); caxis(ax, [vmin vmax]);
                else
                    h.base(k) = scatter3(X(air), Yp(air), Z(air), opts.sz_base, ...
                                         'filled','Parent',ax);
                end
            end
            if any(landed)
                scatter3(X(landed), Yp(landed), Z(landed), opts.sz_base, ...
                         'filled', 'Parent', ax, 'MarkerFaceColor', opts.col_landed);
            end
        end

        % Overlay
        if opts.show_overlay && ~isempty(poses2{k})
            Q = poses2{k}; X2 = Q(:,1); Yp2 = Q(:,3); Z2 = Q(:,2);
            if ~isempty(landed_mask_overlay)
                landed2 = squeeze(landed_mask_overlay(idx2(k),:)).';
            else
                landed2 = (Z2 == 0);
            end
            air2 = ~landed2;

            if any(air2)
                if ~isempty(vvert2{k})
                    C2 = vvert2{k}(air2);
                    h.over(k) = scatter3(X2(air2), Yp2(air2), Z2(air2), opts.sz_overlay, C2, ...
                                         'filled','Parent',ax,'MarkerEdgeColor','none', ...
                                         'MarkerFaceAlpha', opts.alpha_overlay);
                    colormap(ax, opts.cmap); caxis(ax, [vmin vmax]);
                else
                    h.over(k) = scatter3(X2(air2), Yp2(air2), Z2(air2), opts.sz_overlay, ...
                                         'filled','Parent',ax,'MarkerFaceAlpha', opts.alpha_overlay);
                end
            end
            if any(landed2)
                scatter3(X2(landed2), Yp2(landed2), Z2(landed2), opts.sz_overlay, ...
                         'filled','Parent',ax,'MarkerFaceColor', opts.col_landed, ...
                         'MarkerFaceAlpha', 1.0);
            end
        end

        % Ground plane y=0 -> Z=0
        if opts.show_plane
            [XX,YY] = meshgrid(linspace(xlim_all(1),xlim_all(2),2), ...
                               linspace(ylim_all(1),ylim_all(2),2));
            s = surf(XX,YY,zeros(size(XX)), 'Parent', ax);
            s.FaceAlpha = 0.06; s.EdgeColor = 'none';
            h.plane(k) = s;
        end

        view(ax, opts.viewAZEL(1), opts.viewAZEL(2)); grid(ax,'on'); axis(ax,'equal');
        xlim(ax, xlim_all); ylim(ax, ylim_all); zlim(ax, zlim_all);
        xlabel(ax,'x'); ylabel(ax,'z'); zlabel(ax,'y (vertical)');
        title(ax, sprintf('t = %g s', ts(k)), 'FontWeight','bold');

        % colorbar per tile (comment if you want one global colorbar elsewhere)
        cb = colorbar(ax); cb.Label.String = 'vertical velocity v_y (m/s)';
        hold(ax,'off');
    end

    % meta out
    h.ts  = ts;    h.idx = idx;    h.idx2 = idx2;
    h.vlim = [vmin vmax];
end

% --- helper colormap: red (low) -> purple (high)
function cmap = redPurpleMap(n)
    if nargin<1, n=256; end
    r = linspace(1.0, 0.6, n).';
    g = linspace(0.0, 0.0, n).';
    b = linspace(0.0, 0.6, n).';
    cmap = [r g b];
end
