% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   [jx,ix,m,vjx,vix] = gs_trace_step_v(jr,iz,j1,j2,i1,i2,y44,yb,vj,vi);
% 
%   PURPOSE: Find where the contour y = yb exits the rectangle j1, j2, i1, i2
% 
%   INPUTS: jr, iz, entry point for the contour y = yb (both values >=0 & <=1)
%           j1, j2, min and max r of the rectangle (both values >=0 & <=1)
%           i1, i2, min and max z of the rectangle (both values >=0 & <=1)
%           y44, values of y for jr,iz = -1, 0, 1, 2
%           yb, value for which contour is sought
% 
%   OUTPUTS:  jx, ix, r, z of exit point
%             m, where exit is, 0=nowhere, -1=bottom, +1=top, -4=in, +4=out
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%