function h = plot_plumes_times(traj, traj_10, dt, timesSec, landed_mask, landed_mask_10, opts)
% PLOT_PLUMES_TIMES  Plot particle clouds at specified times in a 2x3 grid (<=6 panels).
%
% h = plot_plumes_times(traj, traj_10, dt, timesSec, landed_mask, landed_mask_10, opts)
%
% Supports sparse recording via:
%   opts.frameTimesSec  : 1xM vector of recorded times (s) for traj
%   opts.frameTimesSec2 : 1xM2 vector of recorded times (s) for traj_10
%
% Notes
%   - Plots with axes [X,Y,Z]=[x,z,y] so physical y is vertical.
%   - At most 6 times are plotted (2x3). Extra requested times are dropped with a warning.
%   - If landed masks are empty, "landed" is inferred by (y==0) at each sampled frame.

    if nargin < 7, opts = struct; end
    if isempty(dt) || ~isscalar(dt) || dt <= 0
        error('plot_plumes_times: dt must be a positive scalar.');
    end

    % ---- defaults
    def.show_base     = ~isempty(traj);
    def.show_overlay  = ~isempty(traj_10);
    def.show_plane    = true;
    def.sz_base       = 5;
    def.sz_overlay    = 2.5;
    def.col_overlay   = [0.7 0.7 0.7];
    def.alpha_overlay = 0.5;
    def.zcap          = [];          % [] => auto vertical range
    def.viewAZEL      = [0 0];
    def.frameTimesSec  = [];
    def.frameTimesSec2 = [];
    % merge defaults
    fn = fieldnames(def);
    for i = 1:numel(fn)
        f = fn{i}; if ~isfield(opts,f) || isempty(opts.(f)), opts.(f)=def.(f); end
    end

    % ---- derive available time (prefer explicit recorded times)
    T  = size(traj,   1);  T2 = size(traj_10,1);
    ft  = opts.frameTimesSec(:).';   % may be empty
    ft2 = opts.frameTimesSec2(:).';

    if ~isempty(ft) || ~isempty(ft2)
        t_cap = max([0, ft, ft2]);
    else
        t_max  = (max(T,1)-1)*dt;
        t_max2 = (max(T2,1)-1)*dt;
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
        warning('plot_plumes_times:TooManyTimes', ...
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
    if ~isempty(traj_10)
        if ~isempty(ft2)
            idx2 = arrayfun(@(t) find(abs(ft2 - t) == min(abs(ft2 - t)), 1, 'first'), ts);
        else
            idx2 = 1 + floor(ts./dt + 1e-9);
        end
        idx2 = min(size(traj_10,1), max(1, idx2));
    end

    % ---- gather positions
    poses   = cell(K,1);
    poses10 = cell(K,1);
    for k = 1:K
        if ~isempty(traj),    poses{k}   = squeeze(traj(   idx(k), :,:)); end
        if ~isempty(traj_10), poses10{k} = squeeze(traj_10(idx2(k),:,:)); end
    end

    % ---- bounds in plot coords [x z y] across whichever sets exist
    mins = [ inf  inf  inf];  maxs = [-inf -inf -inf];
    for k = 1:K
        stacks = [];
        if ~isempty(poses{k}),   P = poses{k};   stacks = [stacks; [P(:,1) P(:,3) P(:,2)]]; end
        if ~isempty(poses10{k}), Q = poses10{k}; stacks = [stacks; [Q(:,1) Q(:,3) Q(:,2)]]; end
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

    % ---- figure (2x3) and tiles
    h.fig = figure('Name','Plumes @ specified times', ...
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
                h.base(k) = scatter3(X(air), Yp(air), Z(air), opts.sz_base, 'filled', 'Parent', ax);
            end
            if any(landed)
                scatter3(X(landed), Yp(landed), Z(landed), opts.sz_base, ...
                         'filled', 'Parent', ax, 'MarkerFaceColor',[1 0 0]);
            end
        end

        % Overlay
        if opts.show_overlay && ~isempty(poses10{k})
            Q = poses10{k}; X2 = Q(:,1); Yp2 = Q(:,3); Z2 = Q(:,2);
            if ~isempty(landed_mask_10)
                landed10 = squeeze(landed_mask_10(idx2(k),:)).';
            else
                landed10 = (Z2 == 0);
            end
            air10 = ~landed10;
            if any(air10)
                h.over(k) = scatter3(X2(air10), Yp2(air10), Z2(air10), opts.sz_overlay, ...
                                     'filled','Parent',ax, ...
                                     'MarkerFaceColor',opts.col_overlay, ...
                                     'MarkerFaceAlpha',opts.alpha_overlay);
            end
            if any(landed10)
                scatter3(X2(landed10), Yp2(landed10), Z2(landed10), opts.sz_overlay, ...
                         'filled','Parent',ax,'MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',1.0);
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
        hold(ax,'off');
    end

    % meta out
    h.ts  = ts;    h.idx = idx;    h.idx2 = idx2;
end