% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_evolve_conductors
% 
%   PURPOSE: Calculate the time derivative of conductor currents
% 
%   INPUTS:  new_response_was_calculated, flag to update Amat, Bmat
%            constraints, flag for constraints on plasma profiles
%            lae.ys, flux at conductors from plasma for last analyzed equilibrium
%            dysdx, response of flux at conductors from plasma
%            dx, change of state vector since last analyzed equilibrium
%            f, fraction of dx to apply (normally 1)
%            mss, mutual inductances between conductors when no plasma
%            ic, iv, coil and vessel current vectors
%            u1, voltages applied to conductors
% 
%   OUTPUTS: x1, state of conductor system
%            x1dot, time derivative of x1
%            v1, voltage induced at conductors due to x1dot
%            y1, total flux at conductors including flux error corrections
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%