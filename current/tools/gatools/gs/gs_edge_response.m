%   USAGE:   gs_edge_response
% 
%   PURPOSE: Calculate how xedge responds to changes in rbbbs, zbbbs
% 
%   INPUTS: nedge, fedge, (from gs_trace_edge)
%           DX = diff(rbbbs)', DY = diff(zbbbs)'
%           nz, grid variable
% 
%   OUTPUTS: dxedgedrbbbs0, dxedgedzbbbs0, derivatives w.r.t. point where xedge = 0
%            dxedgedrbbbs1, dxedgedzbbbs1, derivatives w.r.t. point where xedge = 1
%