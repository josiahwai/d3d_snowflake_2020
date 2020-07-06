% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%   USAGE: [c,h] = hdcontourf(X,Y,Z,N)
%          [c,h] = hdcontourf(X,Y,Z,V)
% 
%   PURPOSE: High-definition contouring. Works like contourf 
%            but the matrix is made 8x8 times larger with
%            gs_interp2 before contouring
% 
%   INPUTS:  X,Y = optional vectors with x and y coordinates
%            Z, a matrix with values to contour
%            N = number of contours, or V = values to contour
%            nd or ndx,ndy = extra resolution (default 8)
%            property and value, e.g. 'LineStyle','none'
% 
%   OUTPUTS: c = contours
%            h = plot handle
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%