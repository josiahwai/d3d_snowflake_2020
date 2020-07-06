% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   [jx,ix,m,vjx,vix] = gs_trace_x_step(jr,iz,j1,j2,i1,i2,y44,vj,vi);
% 
%   PURPOSE: Find where the contour bz = 0 exits the rectangle j1, j2, i1, i2
% 
%   INPUTS: jr, iz, entry point for the contour bz = 0 (both values >=0 & <=1)
%           j1, j2, min and max r of the rectangle (both values >=0 & <=1)
%           i1, i2, min and max z of the rectangle (both values >=0 & <=1)
%           y44, values of flux for jr,iz = -1, 0, 1, 2
% 
%   OUTPUTS:  jx, ix, r, z of exit point
%             m, where exit is, 0=x, -1=bottom, +1=top, -4=in, +4=out
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%