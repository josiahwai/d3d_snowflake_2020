%   USAGE:   gs_plot_progress
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
%   OUTPUTS: plots of pprime, ffprim, scalar values,
%              contours, conductor currents, flux error
% 	
%   METHOD:  A resize function (ResizeFcn) scales fonts with window size 
%