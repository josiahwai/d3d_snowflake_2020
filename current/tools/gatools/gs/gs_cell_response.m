%   USAGE:   gs_cell_response
% 
%   PURPOSE: Calculate cell coverage and how it responds to changes in rbbbs, zbbbs
% 
%   INPUTS: data for points where the boundary intersects cell edges:
%             redge, zedge, nedge, iedge, fedge, xedge (see gs_trace_edge)
%           data for boundary points:
%             rbbbs, zbbbs, (see gs_trace_boundary)
%           Other precomputed variables:	
%             rgg, zgg, dr, dz, nr, nz (grid variables)
% 
%   OUTPUTS: Acell, plasma-covered area of each cell on the grid
% 	    RAcell, surface integral of R over the plasma in cell
%            ARcell, surface integral of 1/R over plasma in cell
%             dAcelldrbbbs, dAcelldzbbbs, response of  Acell to rbbbs, zbbbs
% 	    dRAcelldrbbbs,dRAcelldzbbbs, response of RAcell to rbbbs, zbbbs
% 	    dARcelldrbbbs,dARcelldzbbbs, response of ARcell to rbbbs, zbbbs
% 
%   METHOD:  Areas of polygons are computed by integrating R*dZ from corner to corner
%            going ccw around the cell. When integrating along boundary: R(x)*Z'(x)dx
%            Does all that gs_cell_coverage does *and* also response calculations
%