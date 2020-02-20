function [L, r, z] = calcConnectionLength_ToCurve(r0, z0, psizr, rg, zg, ...
    bzero, rzero, ipsign, rCurve, zCurve)
%
% CALCCONNECTIONLENGTH_TOCURVE
%
%   Calculate the connection length for a fieldline in the scrape-off layer.
%
%   The connection length is computed from the point (r0,z0) to the poloidal 
%   curve defined by the points (rCurve,zCurve).
%
% USAGE: calcConnectionLength_ToCurve.m
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
%   ipsign....direction of plasma current (+1/-1)
%
%   rCurve....array containing the radial coordinates of points which define
%             the poloidal curve curve [m]
%
%   zCurve....array containing the vertical coordinates of points which define
%             the poloidal curve curve [m]
%
% OUTPUTS: 
%
%   L.........fieldline connection length [m]
%
%   r.........radial coordinate of point at which connection length ends   [m]
%
%   z.........vertical coordinate of point at which connection length ends [m]
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 06/04/2018
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/04/2018
%
%...............................................................................

% Define the fieldline length integration variable

s0 = 0;
sE = 100;    % maximum allowable connection length of 100 m

ds = 1e-4;   % integration stepsize of 0.1 mm
ns = sE/ds;

phi0 = 0; 

for ii = 1:ns
    
    [~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, r0, z0);
    
    Br = ipsign*(-1/(2*pi*r0))*dpsidz;
    Bz = ipsign*( 1/(2*pi*r0))*dpsidr;
    
    Bphi = (rzero*bzero)/r0;
    
    kn1_r   = Br;
    kn1_phi = Bphi/r0;
    kn1_z   = Bz;
    
    kn2_r   = Br      + (ds/2)*kn1_r;
    kn2_phi = Bphi/r0 + (ds/2)*kn1_phi;
    kn2_z   = Bz      + (ds/2)*kn1_z;
    
    kn3_r   = Br      + (ds/2)*kn2_r;
    kn3_phi = Bphi/r0 + (ds/2)*kn2_phi;
    kn3_z   = Bz      + (ds/2)*kn2_z;
    
    kn4_r   = Br      + ds*kn3_r;
    kn4_phi = Bphi/r0 + ds*kn3_phi;
    kn4_z   = Bz      + ds*kn3_z;
        
    dr   = (ds/6)*(kn1_r   + 2*kn2_r   + 2*kn3_r   + kn4_r);
    dphi = (ds/6)*(kn1_phi + 2*kn2_phi + 2*kn3_phi + kn4_phi);
    dz   = (ds/6)*(kn1_z   + 2*kn2_z   + 2*kn3_z   + kn4_z);
    
    r0   = r0   + dr;
    phi0 = phi0 + dphi;
    z0   = z0   + dz;
    
    if (r0 >= min(rCurve)) && (r0 <= max(rCurve))
        
        zq = interp1(rCurve, zCurve, r0); 
        
        if abs(zq - z0) < 1e-4
            L = s0 + ii*ds;
            r = r0;
            z = z0;
            return
        end
        
    end
    
end

L = 0;
r = 0;
z = 0;

end
