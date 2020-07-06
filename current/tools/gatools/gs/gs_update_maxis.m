%   USAGE:   gs_update_maxis
% 
%   PURPOSE: Zoom in on magnetic axis
% 
%   INPUTS: rmaxis, zmaxis, a best estimate of the position such as found
%             by a linear prediction of change since previous analysis
%           psizr, the flux at nz vertical positions * nr radial positions
%           rgg, zgg, dr, dz, nr, nz (grid variables)
%           For bicubic interpolation on the grid:
%             mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%             neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
%           The function isinpoly must have been called with limiter information:
%           isinpoly([],[],Rlim,Zlim)
% 
%   OUTPUTS: maxis, a boolean true if axis still inside limiter after update
%            rmaxis, zmaxis, position of magnetic axis
%            wa, iia, weights and indices such that psimag = wa*psizr(iia)
%            drmaxisdpsi, weights such that drmaxis = drmaxisdpsi*dpsizr(iia)
%            dzmaxisdpsi, weights such that dzmaxis = dzmaxisdpsi*dpsizr(iia)
% 	
%   METHOD: interpolation with bicubic Hermite splines, Newton-Rhapson to zoom
%