function plotCylWireframe(r, theta, z, varargin)
% plotCylWireframe(r,theta,z, ...)
% Draw a colored cylindrical wireframe given coordinate *centers*:
%   r     : [Nr x 1] radial centers (>=0)
%   theta : [Ntheta x 1] angular centers (0..2π)
%   z     : [Nz x 1] vertical centers
%
% Optional name-value args:
%   'RadialEvery' : stride along r (default: round(numel(r)/12))
%   'ThetaEvery'  : stride along theta (default: round(numel(theta)/24))
%   'ZEvery'      : stride along z (default: round(numel(z)/8))
%   'Colors'      : struct with fields Rings, Radials, Verticals (RGB or color char)
%                   default: Rings='b', Radials='g', Verticals='r'
%   'LineWidth'   : scalar line width (default: 0.7)
%   'MakeFigure'  : true/false create new figure (default: true)
%
% Example:
%   plotCylWireframe(r,theta,z,'RadialEvery',10,'ThetaEvery',15,'ZEvery',5);

% ---- defaults ----
radialEvery = max(1, round(numel(r)/12));
thetaEvery  = max(1, round(numel(theta)/24));
zEvery      = max(1, round(numel(z)/8));
colors.Rings     = 'b';
colors.Radials   = 'g';
colors.Verticals = 'r';
lw = 0.7;
mkFig = true;

% ---- parse varargin (lightweight) ----
for k = 1:2:numel(varargin)
    name = varargin{k};
    val  = varargin{k+1};
    switch lower(name)
        case 'radialevery', radialEvery = max(1, round(val));
        case 'thetaevery',  thetaEvery  = max(1, round(val));
        case 'zevery',      zEvery      = max(1, round(val));
        case 'colors',      colors = val;
        case 'linewidth',   lw = val;
        case 'makefigure',  mkFig = logical(val);
        otherwise, error('Unknown parameter: %s', name);
    end
end

if mkFig, figure('Color','w'); end
hold on

% ---- rings at constant z ----
for iz = 1:zEvery:numel(z)
    for ir = 1:radialEvery:numel(r)
        x = r(ir)*cos(theta);
        y = r(ir)*sin(theta);
        plot3(x, y, z(iz)*ones(size(theta)), '-', 'Color', colors.Rings, 'LineWidth', lw);
    end
end

% ---- radial lines at constant z ----
for iz = 1:zEvery:numel(z)
    for it = 1:thetaEvery:numel(theta)
        x = r.*cos(theta(it));
        y = r.*sin(theta(it));
        plot3(x, y, z(iz)*ones(size(r)), '-', 'Color', colors.Radials, 'LineWidth', lw);
    end
end

% ---- vertical lines ----
for ir = 1:radialEvery:numel(r)
    for it = 1:thetaEvery:numel(theta)
        x = r(ir)*cos(theta(it))*ones(size(z));
        y = r(ir)*sin(theta(it))*ones(size(z));
        plot3(x, y, z, '-', 'Color', colors.Verticals, 'LineWidth', lw);
    end
end

axis equal; box on; grid on
xlabel('x'); ylabel('y'); zlabel('z');
title('Cylindrical grid (colored wireframe)');
view(35, 20);
hold off
end
