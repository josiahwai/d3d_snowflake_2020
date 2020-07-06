% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_update_eq
% 
%   PURPOSE: Update equilibrium: psizr, ic, iv, sp, sf and output y
%            Update response when dpsibar_limit_for_linear is reached
%            
%   INPUTS:  The gs workspace, and:
%            xc or xs, structure with inputs 
%                 ic: coil currents such that psizr_app = mpc*ic
%                 iv: vessel currents such that psizr_app = mpv*iv
%                 ip: total plasma current
%                 li: normalized inductance
%              betap: poloidal beta
%             plotit: Special flag for plotting equilibrium (default 0)
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