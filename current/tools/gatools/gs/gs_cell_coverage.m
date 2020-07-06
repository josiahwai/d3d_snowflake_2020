%   USAGE:   gs_cell_coverage
% 
%   PURPOSE: Calculate how much area of each grid cell is covered by plasma
%            A "cell" is a rectangle with a coordinate rg,zg at its center
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
% 
%   METHOD:  Areas of polygons are computed by integrating R*dZ from corner to corner
%            going ccw around the cell. When integrating along boundary: R(x)*Z'(x)dx
%