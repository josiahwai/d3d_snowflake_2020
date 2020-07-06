% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   y = gs_interp2(rg, zg, yg, r, z)
%            Derivatives also available with:
%            [y, yr, yz, yrr, yrz, yzz, yrrr, yrrz, yrzz, yzzz] = ...
%            gs_interp2(rg, zg, yg, r, z)
%            To return indices and weights:
%            [i, w, wr, wz, wrr, wrz, wzz, wrrr, wrrz, wrzz, wzzz] = ...
%            gs_interp2(rg, zg, yg, r, z, 'WEIGHTS')
%            where y = sum(w'.*yg(i)'), yr = sum(wr'.*yg(i)'), etc.
% 
%   PURPOSE: Interpolate values of yg on grid rg, zg to get values at r, z
%            using cubic Hermite splines (which are used by gs codes)
%            http://en.wikipedia.org/wiki/Bicubic_interpolation
% 
%   INPUTS: rg, zg,  grid point coordinates (default 1:nr, 1:nz)
%           yg,      values at the nz x nr points on the grid
%           r, z,    interpolation points (default rg(:)', zg(:))
% 
%   OUTPUTS:  y,       interpolated values at points r, z
%             yr, etc, derivative of y w.r.t. to r, etc at points r, z
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%