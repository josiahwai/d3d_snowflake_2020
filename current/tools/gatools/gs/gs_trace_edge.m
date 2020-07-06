%   USAGE:   gs_trace_edge
% 
%   PURPOSE: Find coordinates where the boundary intersects edges of grid cells
%            (used in calculation of how much of cells are covered by plasma)
% 
%   INPUTS: rbbbs, zbbbs, nbbbs: R, Z of boundary and number of points
%           rmaxis, zmaxis, position of axis (used to order points w.r.t theta)
%           Other precomputed variables:	
%             rg, zg, dr, dz (grid variables)
% 
%   OUTPUTS: redge, zedge, nedge: R, Z, and number of points that intersect grid cells
%            rhoedge, thedge, distance to magnetic axis and poloidal angle
%            fedge, direction when entering cells going ccw -1=down, 1=up, -nz=in, nz=out
%            xedge, values of x to use in interpolation formula to obtain redge, zedge
%            iedge, indices to the first of two bbbs points (j in the formula)
%            gedge, indices to cells entered at redge, zedge, going along bbbs
% 	
%   METHOD:  The linear interpolation formula:
%            R = rbbbs(j) + (rbbbs(j+1)-rbbbs(j))*x
%            Z = zbbbs(j) + (zbbbs(j+1)-zbbbs(j))*x
%            is solved for x-values that make either R or Z
%            intersect a cell edge. All solutions with x between 0 and 1 are
%            returned in redge, zedge. The points are sorted according to thedge.
%