function [L, r, z] = calcConnectionLength(r0, z0, psizr, rg, zg, ...
    bzero, rzero, limdata, iSign)
%
% CALCCONNECTIONLENGTH
%
%   Calculate the connection length for a fieldline in the scrape-off layer.
%
% USAGE: calcConnectionLength.m
%
% INPUTS:
%
%   r0........initial radial coordinate for point in SOL   [m]
%
%   z0........initial vertical coordinate for point in SOL [m] 
%
%   psizr.....matrix with dimensions (nz x nr) containing the magnetic flux
%             at nz vertical by nr radial grid points
%
%   rg........array containing the nr radial grid points
%
%   zg........array containing the nz vertical grid points
%
%   bzero.....applied toroidal field strength at radius rzero [T]
%
%   rzero.....radius at which the toroidal field strength bzero is given [m]
%
%   limdata...limiter (r,z) vertices as defined in tok_data_struct
%
%    
% OUTPUTS: 
%
%   L.........fieldline connection length [m]
%                       
%
%
%...............................................................................

zlim = limdata(1,:);
rlim = limdata(2,:);

% Define the fieldline length integration variable

s0 = 0;
L = 500;    % maximum allowable connection length of 100 m

ds = 0.01;  % integration stepsize of 1 cm
ns = L/ds;

for ii = 1:ns

    [~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, r0, z0);
    
    Br = -1/(2*pi*r0)*dpsidz;
    Bz =  1/(2*pi*r0)*dpsidr;
    Bphi = rzero*bzero/r0;
    
    B = norm([Br Bphi Bz]);
    
    r0 = r0 + iSign*ds/B*Br;
    z0 = z0 + iSign*ds/B*Bz;
    
    in = inpolygon(r0, z0, rlim, zlim);
    
    if ~in
        L = s0 + (ii-1)*ds;
        r = r0;
        z = z0;
        return
    end
    
end







