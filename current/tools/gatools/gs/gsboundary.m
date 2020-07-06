% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   b = gsboundary(c,psizr,ra,za,nax,dnx);
% 
%   PURPOSE: Find boundary and related quantities needed to 
%            calculate flux from plasma
% 
%   INPUTS:  c, configuration data made by gsconfig.m
%            psizr, flux on the grid c.rg, c.zg
%            ra, za, [m] pick axis closest to ra, za (default nans)
%            nax, flag to reject axis > 2 grid units from ra, za (default false)
%            dnx, use nearx method when contour is within dnx of an
%              x-point [max=1, default=0.5, floating index units]
%              Nearx keeps boundary speed finite near x-points using
%              an interpolation that gives nonzero field at x-points
%            dpb, (default 0.005) moves boundary inward
%              Creates smooth transitions between bdef points
%   OUTPUTS: b, information about plasma boundary
%