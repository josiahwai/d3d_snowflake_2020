%   USAGE:   gs_eq_analysis
% 
%   PURPOSE: Find magnetic axis and boundary, calculate the plasma currents on the grid
%            Find the total flux from plasma and from external conductors
% 
%   INPUTS: several variables generated by gs_initialize are needed
% 
%   OUTPUTS: Some of the quantities are:
%            rmaxis, zmaxis, psimag, position of axis and its flux
%            rbdef, zbdef, psibry, boundary-defining point (x or touch) and its flux
%            rx1, zx1, psix1, position and flux of most important x-point below axis
%            rx2, zx2, psix2, position and flux of most important x-point below axis
%            rtl, ztl, psitl, position and flux of most important limiter-point
%            rbbbs, zbbbs, nbbbs, nbbbs boundary positions
%            pcurrt, current within grid cells (A)
%            Acell, RAcell, ARcell, amount of plasma-coverage in grid cells at boundary
%            psizr_pla, psizr_app, flux from plasma and from conductors
%            psizr_err = psizr-psizr_pla-psizr_app
%            cpasma = sum(pcurrt(:)), total plasma current
%            psipla, plasma flux
%            Vtot, total plasma volume
%            Wth, total thermal energy
%            betap, poloidal beta
%            li, normalized inductance
%