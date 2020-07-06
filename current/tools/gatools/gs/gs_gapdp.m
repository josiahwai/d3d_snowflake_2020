%   USAGE:   gs_gapdp
% 
%   PURPOSE: Find gaps between boundary and diagnostic points rdp, zdp
%            Gap distances measured along vector from axis to dp points
% 
%   INPUTS: psibarzr, psibarzr = (psizr-psimag)/(psibry-psimag);
%           rdp, zdp, diagnostic points
%           rmaxis, zmaxis, position of axis
%           rg, zg, dr, dz, nr, nz (grid variables)
%           For cubic interpolation on the grid:
%             mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%             neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
% 
%   OUTPUTS: gapdp, distance from dp to boundary in direction toward axis
%            rdpb, zdpb, = R, Z, of boundary points
%            idpb, indices to 16 grid points around each of rdpb, zdpb
%            wdpb, weights for 16 grid points around each of rdpb, zdpb
%            dgapdpdpsi, how gapdp respondes to (dpsizr(idpb) - dpsibry)
% 	
%   METHOD: Newton-Rhapson method used on all contour-points in parallel
%