function [zgnew, xgnew, ygnew] = regrid(xg, yg, zg, nx, ny)
%
% REGRID
%
%   Regrid a set of 2-D gridded data onto an (nx x ny) grid using bicubic
%   Hermite interpolation.
%
% USAGE: regrid.m
%
% INPUTS:
%
%   xg......coordinates of grid points in first dimension  (n x 1)
%
%   yg......coordinates of grid points in second dimension (m x 1)
%
%   zg......values of the underlying function at the grid points (m x n)
%
%   nx......number of grid points on new grid in first dimension
%
%   ny......number of grid points on new grid in second dimension
%
% OUTPUTS
%
%   z0......interpolated value at the query point
%
%   xgnew...coordinates of points on new grid in first dimension
%
%   ygnew...coordinates of points on new grid in second dimension
%
% METHOD: 
%
% AUTHOR: Patrick J. Vail
%
% DATE: 10/23/2018
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 10/23/2018
%
%...............................................................................

% Calculate grid points for the new (nx x ny) grid

xgnew = linspace(xg(3), xg(end-2), nx)';
ygnew = linspace(yg(3), yg(end-2), ny)';

% Interpolate gridded data on the new grid

zgnew = zeros(length(xgnew), length(ygnew));

for ii = 1:length(ygnew)
    for jj = 1:length(xgnew)
        zgnew(ii,jj) = bicubicHermite(xg, yg, zg, xgnew(jj), ygnew(ii));
    end
end

end
