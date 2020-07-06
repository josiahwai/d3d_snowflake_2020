% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE: p = gsprofiles(e)
%          p = gsprofiles(e,b,c)
% 
%   PURPOSE: Calculate profiles for Grad-Shafranov equilibrium
% 
%   INPUTS: e, structure with fields:
%              psizr = flux on grid [Wb]
%              rg, zg = grid coordinates [m]
%              rbbbs, zbbbs, nbbbs = boundary coordinates
%              rmaxis, zmaxis = axis coordinate
%              fpol, ffprim, pres, pprime
%              xlim, ylim = limiter coordinates (R,Z values)
%           b, output from gsboundary (for response calculation)
%           c, output from gsconfig (for connection lengths calculation)
% 
%   OUTPUTS: p, structure with profiles versus normalized flux
%               p.info contains detailed information
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%