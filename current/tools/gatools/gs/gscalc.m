%   USAGE: [psizr_pla, br, bz] = gscalc(jphi,rg,zg,mpp,nrb,nzb)
% 
%   PURPOSE: Calculate contribution to flux and fields from plasma
% 
%   INPUTS: jphi, current density in plasma [MA/m2]
%             rg, vector of R coordinates for the grid [m]
%             zg, vector of Z coordinates for the grid [m]
%            mpp, (optional) supply (nr*nz,nr) mutuals [H] if available.
%        nrb,nzb, (optional) number of boundary calculations per edge.
%                 Values from 2 to [nr,nz]. Default: sqrt([nr,nz]).
% 
%   OUTPUTS: psizr_pla, plasma flux [Wb]
%                   br, plasma radial field [T]
%                   bz, plasma vertical field [T]
% 
%   METHOD: The Grad-Shafranov equation is solved on finite-difference form
%           with boundary calculated by mutual inductances.
%           Tne FD equations are solved by sine transforms in Z
%           and tridiagonal matrix solutions in R.
%