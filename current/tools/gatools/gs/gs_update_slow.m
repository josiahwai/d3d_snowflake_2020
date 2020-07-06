% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_update_slow
% 
%   PURPOSE: Update equilibrium: psizr, ic, iv, sp, sf and output y
%            
%   INPUTS:  Must have run gs_update_fast
% 	
%   OUTPUTS:  Updated:
%             psizr, flux on grid [Wb]
%             ic, coil currents [A]
%             iv, vessel currents [A]
%             sp, spline parameters for pressure versus normalized flux
%             sf, spline parameters for fpol^2/2-(rzero*bzero)^2/2
%             y, diagnostic outputs
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%