%   USAGE:   gs_trace_contours
% 
%   PURPOSE: Trace contours of constant flux in the plasma
%            This version iterates to find the contours
% 
%   INPUTS: psibarzr, psibarzr = (psizr-psimag)/(psibry-psimag);
%           psibarc, normalized fluxes at contours (default psibar)
%           npola, number of poloidal angles (default 2*nr-1)
%           rbdef, zbdef, position that defines the boundary
%           anglesoption, for poloidal angles of contour points:
%             1. start at bdef point, equally spaced angles (default)
%             2. Like 1 for boundary but along gradient toward axis
%           rmaxis, zmaxis, position of axis
%           rg, zg, dr, dz, nr, nz (grid variables)
%           For cubic interpolation on the grid:
%             mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%             neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
% 
%   OUTPUTS: rcont, zcont, ncont, = R, Z, number of contours
% 	
%   METHOD: Newton-Rhapson method used on all contour-points in parallel
%