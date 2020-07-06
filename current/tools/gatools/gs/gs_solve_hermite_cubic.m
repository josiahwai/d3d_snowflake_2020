% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE:   [xb, rootin01, x01min, x01max] = gs_solve_hermite_cubic(y, yb)
% 
%   PURPOSE: Find where the function y = yb, with y given at x = -1, 0, 1, 2,
%            and interpolated between points using cubic Hermite splines
% 
%   INPUTS: y,  values at x = -1, 0, 1, 2
%           yb, value for which the roots xb are sought
% 
%   OUTPUTS:  xb, either 1 or 3 roots are real, first 1 always real
%             rootin01, true if at least one root is between 0 and 1
%             x01min, smallest root in the interval 0 to 1, if existent
%             x01max, biggest root in the interval 0 to 1, if existent
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%