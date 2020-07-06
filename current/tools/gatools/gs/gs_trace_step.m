% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   [jx,ix,m,hit] = gs_trace_step(...
%   jr,iz,j1,j2,i1,i2,y44,yb,xhl,xvi,xhu,xvo, rzhit);
% 
%   PURPOSE: Find where contour y = yb exits rectangle j1, j2, i1, i2
% 
%   INPUTS: jr, iz, entry point for the contour y = yb (both values >=0 & <=1)
%           j1, j2, min and max r of the rectangle (both values >=0 & <=1)
%           i1, i2, min and max z of the rectangle (both values >=0 & <=1)
%           y44, values of y for jr,iz = -1, 0, 1, 2
%           yb, value for which contour is sought
%           xhl, (optional) 3 values of j where y = yb for i = j1 (or nan)
%           xvi, (optional) 3 values of i where y = yb for j = j1 (or nan)
%           xhu, (optional) 3 values of j where y = yb for i = j2 (or nan)
%           xvo, (optional) 3 values of i where y = yb for j = j2 (or nan)
%           rzhit, (optional) [r1 r2 z1 z2] of region for hit test
% 
%   OUTPUTS:  jx, ix, r, z of exit point
%             m, where exit is, 0=nowhere, -1=bottom, +1=top, -4=in, +4=out
%             hit, boolean true if contour goes through rzhit region
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%