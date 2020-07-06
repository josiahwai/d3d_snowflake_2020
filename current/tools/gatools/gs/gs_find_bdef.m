%   USAGE:   gs_find_bdef
% 
%   PURPOSE: Find boundary-defining point (x or touch point)
% 
%   INPUTS: psizr, the flux at nz vertical positions * nr radial positions
%           ilimgg, flags for grid, -1 = outside vessel, >-1 = inside
%           rgg, zgg, dr, dz, nr, nz (grid variables)
%           For cubic interpolation on the grid:
%             mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%             neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
% 
%   OUTPUTS: rbdef, zbdef, position that defines the boundary (touch or x point)
%            rx1, zx1, psix1, position and flux of most important x-point below axis
%            rx2, zx2, psix2, position and flux of most important x-point below axis
%            rtl, ztl, psitl, position and flux of most important limiter-point
%            wb, iib, weights and indices such that psibry = wb*psizr(iib)
%            drbefdpsi, dzbdefdpsi, weights such that drbdef = drbdefdpsi*dpsizr(iib)
%            Weights and indices for x1, x2, tl are called:
%              wx1, iix1, drx1dpsi, dzx1dpsi, wx2, iix2, drx2dpsi, dzx2dpsi, wtl, iitl
%            ix1, flag for lower x-point defines the boundary
%            ix2, flag for upper x-point defines the boundary
%            itl, flag for plasma touches limiter
%            zbot, ztop, limits for z position of plasma
%            zbotg, ztopg, if zgg(k)==zbotg+dz/2 then zbot is in range zgg(k)+[0 dz]
% 	
%   METHOD: 
%