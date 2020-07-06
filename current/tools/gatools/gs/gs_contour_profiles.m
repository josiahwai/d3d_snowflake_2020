%   USAGE:   gs_contour_profiles
% 
%   PURPOSE: Calculate profile functions using contours rcont, zcont
% 
%   INPUTS: output from gs_trace_contours
% 
%   OUTPUTS: Vc, volume within contours
%            Ac, area within contours
%            Lc, area-integral of 1/R within contours
%            All contour-derived profiles of size [ncont,1] have suffix c
% 	
%   METHOD: 
%