%   USAGE:   gs_cp_analysis_response
% 
%   PURPOSE: Analyze a plasma described by circle model in same manner as
%            gs_eq_analysis does for grid model
%            Also calculate response immediately for each quantity
% 
%   INPUTS: r0, z0, a0 = geometric center and minor radius
%           All r, z that are derived from r0,z0,a0 including rbdef, zbdef
%           psih = total flux at rh, z0
%           sp, sf = parameters for pres and fpol
%           ic, iv = coil and vessel currents
% 
%   OUTPUTS: Many equilibrium quantities and how they respond to p,
%            where p = [psih r0 z0 a0]
%