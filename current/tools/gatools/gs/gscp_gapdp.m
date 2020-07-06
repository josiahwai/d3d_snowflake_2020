%   USAGE:   gscp_gapdp
% 
%   PURPOSE: Return gapdp for circular plasma model
%            By definition gapdp points toward the magnetic axis but
%            here the circle center is used instead since it is easier
% 
%   INPUTS: r0, z0, a0
% 
%   OUTPUTS: gapdp, distance from dp to boundary in direction toward axis
%            rdpb, zdpb, = R, Z, of boundary points
%            idpb, indices to 16 grid points around each of rdpb, zdpb
%            wdpb, weights for 16 grid points around each of rdpb, zdpb
%            dgapdpdpsi, how gapdp respondes to (dpsizr(idpb) - dpsibry)
% 	
%