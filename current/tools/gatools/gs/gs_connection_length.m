%   USAGE:   gs_connection_length
% 
%   PURPOSE: Calculate connection lengths from all grid points inside limiter
% 
%   INPUTS: psizr, flux at grid points rg, zg
% 
%   OUTPUTS: lconnfzr, connection lengths from grid points to limiter along phi
%            lconnbzr, connection lengths going opposite to phi direction
% 
%   METHOD:  countourc finds contours for fluxes psizr inside limiter
%