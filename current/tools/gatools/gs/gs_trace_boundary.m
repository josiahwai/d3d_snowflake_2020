%   USAGE:   gs_trace_boundary
% 
%   PURPOSE: Trace plasma boundary
% 
%   INPUTS: rbdef, zbdef, position that defines the boundary (touch or x point)
%           psizr, the flux at nz vertical positions * nr radial positions
%           psibry, the flux at the boundary (see gs_find_bdef)
%           rmaxis, zmaxis, position of axis (used to order points w.r.t theta)
% 
%   OUTPUTS: rbbbs, zbbbs, nbbbs, R, Z of boundary and number of points
%            rhobbbs, thbbbs, distance to magnetic axis and poloidal angle
% 	
%   METHOD:  The grid is searched for points that are on either side of the boundary
%            The exact locations are then found by solving 1-dim cubic equations
%            The boundary defining point is added and also extra points if needed
%            to make all poloidal angles between adjacent points < dthbbbs_threshold 
%