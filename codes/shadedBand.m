function h = shadedBand(x, y, e, varargin)
% Plot a shaded band y ± e. Works with linear or log Y axes.
% Cleans NaN/Inf, sorts x, and clamps the lower band for log-scale axes.

    ax = gca;
    isLogY = strcmpi(ax.YScale,'log');

    % force columns
    x = x(:); y = y(:); e = e(:);

    % valid (finite) points
    good = isfinite(x) & isfinite(y) & isfinite(e);

    % if log-y, ensure band stays > 0
    yU = y + e;
    yL = y - e;
    if isLogY
        % small positive floor based on current positive data
        posRef = min(yU(good & yU>0));
        if isempty(posRef), posRef = 1e-12; end
        tiny = max(1e-12, 1e-9*posRef);
        yL = max(yL, tiny);
        yU = max(yU, tiny);  % in case y is tiny too
    end

    % re-evaluate good after adjustments
    good = good & isfinite(yL) & isfinite(yU) & yU>=yL;

    % bail if nothing to draw
    if ~any(good)
        warning('shadedBand:NoData','No finite/positive points to fill.');
        h = gobjects(0);
        return;
    end

    % sort by x to avoid self-crossing polygons
    [xs, ix] = sort(x(good));
    yUs = yU(good); yUs = yUs(ix);
    yLs = yL(good); yLs = yLs(ix);

    % build polygon
    X = [xs; flipud(xs)];
    Y = [yUs; flipud(yLs)];

    % draw (use patch – same as fill but returns a Patch object)
    h = patch('XData',X,'YData',Y, ...
              'FaceColor',[1 0 0], 'FaceAlpha',0.12, ...
              'EdgeColor','none', varargin{:});
end
