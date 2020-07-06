% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_update_fast
% 
%   PURPOSE: Update equilibrium: psizr, ic, iv, sp, sf and output y
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
%             g, largest fraction of dx(nx) that can be done
%             equilibrium_update_is_complete, if false run gs_update_slow
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%