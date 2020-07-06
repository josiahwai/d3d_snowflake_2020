% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE: init = gsinit(c,equ)
% 
%   PURPOSE: Create initial data for gseq, gsdesign, gsevolve
% 
%   INPUTS: c, the output from gsconfig.m
%           equ, structure with initial equilibrium quantities:
%                 rg, zg required if equ and c grids differ
%                 psizr, psimag, psibry required if a plasma exists
%                 x is made with first available items in this list:
%                 ic (TokSys coil currents) or cc (EFIT coil currents)
%                 iv (TokSys vessel currents) or vc (EFIT vessel currents)
%                 fpol or ffprim,rzero,bzero
%                 pres or pprime
% 
%   OUTPUTS: init, initial equilibrium containing fields x and psizr
%                  where x is a state vector describing ic, iv, fpol, pres
% 
%   RESTRICTIONS: The initial equilibrium likely won't be perfectly 
%                 converged unless made by gsdesign
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%