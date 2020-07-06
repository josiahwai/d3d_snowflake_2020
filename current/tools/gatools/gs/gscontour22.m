% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE: [r,z,n,closed,around,P] = gscontour22(y,r0,z0,np,ra,za,vin)
% 
%   PURPOSE: Find closed contour in matrix y passing through r0, z0
% 
%   INPUTS:  y, matrix to find contours in
%            r0, z0, starting point (units are floating index)
%            np, length of output vectors
%            ra, za, contour must enclose ra,za if included and ~nan
%            vin, initial contour direction*vin will not be negative
% 
%   OUTPUTS: r,z, vectors of length np, beginning with r0, z0,
%            then n-1 points along contour back to r0, z0, then nans
%            closed, true if the contour is closed
%            around, true if the contour wraps around ra, za
%            P(m) = floor(2*z-1)+(2*nz-1)*floor(2*r-2) for r,z(m-1:m)
%            gs_interp2(y,r,z) = gs_interp2(y,r0,z0)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%