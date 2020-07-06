%   USAGE:   gs_trace_boundaryx
% 
%   PURPOSE: Trace the plasma boundary, this x-version also finds the x-point
% 
%   INPUTS: psizr, the flux at nz vertical positions * nr radial positions
%           psibry, the flux at the boundary (see gs_find_bdef)
%           ilimgg, flags for grid, -1 = outside vessel, >-1 = inside
%           nbbbs_max, maximum number of boundary points
%           zbot, ztop, min and max z position of plasma
%           zbotg = (floor((zbot-zg(1))/dz)+0.5)*dz+zg(1);
%           ztopg = (floor((ztop-zg(1))/dz)+0.5)*dz+zg(1);
%           rbdef, zbdef, position that defines the boundary (touch or x point)
%           ix1, flag which is true if lower x-point defines the boundary
%               if ix1 is true then these derivates of the flux are needed:
%               psix1_rr, psix1_rz, psix1_zz (see gs_find_bdef)
%           ix2, flag which is true if upper x-point defines the boundary
%               if ix2 is true then these derivates of the flux are needed:
%               psix2_rr, psix2_rz, psix2_zz (see gs_find_bdef)
%           rmaxis, zmaxis, position of axis (used to order points w.r.t theta)
%           rgg, zgg, dr, dz, nr, nz (grid variables)
%           R13 = (1+1i*sqrt(3))^2/4; % for solving cubic equations
%           For cubic interpolation on the grid:
%             mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%             neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
% 
%   OUTPUTS: rbbbs, zbbbs, nbbbs, R, Z of boundary and number of points
%            rhobbbs, thbbbs, distance to magnetic axis and poloidal angle
% 	
%   METHOD:  The grid is searched for points that are on either side of the boundary
%            The exact locations are then found by solving 1-dim cubic equations
%            The boundary defining point is added and also extra points if needed
%            to make all poloidal angles between adjacent points < dthbbbs_threshold 
%