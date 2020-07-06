%   USAGE:   gsdesign_plot_progress
% 
%   PURPOSE: Plot progress of Grad-Shafranov (gs) calculations
% 
%   INPUTS: pprime, ffprim
%           cpasma, total plasma current
%           betap, poloidal beta
%           li, normalized inductance
%           psimag, the flux at the magnetix axis
%           psibry, the flux at the boundary
%           nulls, a structure with nulls that may affect boundary
%           psizr, the flux at nz vertical positions * nr radial positions
%           rbdef, zbdef, point that defines the boundary (x or touch)
%           ic, iv, currents in coils and vessel
% 
%   OUTPUTS: plots of conductor currents, pprime & ffprim, scalars,
%              contours, weights*errors, flux error
% 	
%   FEATURES: A resize function (ResizeFcn) scales fonts with window size 
%             Zooming in on errors reveals their individual names 
%