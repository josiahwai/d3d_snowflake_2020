% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   gs_find_dx
% 
%   PURPOSE: Find change of x (dx) since last analyzed equilibrium (lae)
% 
%   INPUTS:  One of either xc or xs
%            xc, array or structure, may contain:
%                ic, iv, sp, sf, cpasma, li, betap
%                (content depending on constraints)
%            xs, state vector (content controlled by evolve_option)
% 
%   OUTPUTS: dx, the change of the state vector, x (=[ic;iv;sp;sf;er])
%            dcpasma, dli, dbetap, dWth, changes of state plasma quantities
% 	
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%