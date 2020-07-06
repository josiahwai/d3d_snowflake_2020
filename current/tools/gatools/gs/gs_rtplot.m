% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_rtplot
% 
%   PURPOSE: Real-time plot of plasma during simulation
% 
%   INPUTS:  plot_settings (optional)
%            Rlim, Zlim, vvdata, fcdata, ecdata (optional)
%            time-dependent data produced by gs_rtplot_prepare:
%              rbbbs, zbbbs, psizr, cpasma, li, betap
% 
%   OUTPUTS: Fast plots in figure window which is current on first call
%            Other figure windows can be used while gs_rtplot is running
% 
%   METHOD:  Handles to all objects remembered to allow speedy updates 
%