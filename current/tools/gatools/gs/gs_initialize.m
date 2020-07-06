% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_initialize
% 
%   PURPOSE: Initialize gseq, load an equilibrium
% 
%   INPUTS: init, structure with initial equilibrium quantities.
%                 These are used if available:
%                 rg, zg, psizr, psimag, rmaxis, zmaxis, psibry, 
%                 pprime, ffprim, rzero, bzero
%                 The grid in init does not need to match the grid in config
% 
%   OUTPUTS: Initialized equilibrium
% 
%   RESTRICTIONS: The initial equilibrium in gs will only approximate init
%                 if any of the following occurs:
%                 1. the grids in init and config don't match
%                 2. spline knots for pres, fpol in init don't match gs spec
%                 3. edge currents exist (gs considers partial cell coverage)
%                 4. init isn't converged
%                 Normally none of these will produce a serious discrepancy
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%