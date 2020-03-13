function [L, r1, z1] = calcConnectionLength(r0, z0, psizr, rg, zg, ...
    bzero, rzero, limdata, iSign, startin)
%
% CALCCONNECTIONLENGTH
%
%   Calculate the connection length for a fieldline in the scrape-off layer.
%
% USAGE: calcConnectionLength.m
%
% INPUTS:
%
%   r0........initial radial coordinate vector for points in SOL   [m]
%
%   z0........initial vertical coordinate vector for points in SOL [m] 
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
%   iSign... direction to integrate magnetic field +/- 1
%   
%   startin... (default=true(size(r0)) vector indicating,
%              Did the point start inside the limiter? This will also 
%              change the integration direction for that point
%
% OUTPUTS: 
%
%   L.........fieldline connection length [m]
%                       
%
%
%...............................................................................

if ~exist('startin','var')
    startin = true(size(r0));
end

    
zlim = limdata(1,:);
rlim = limdata(2,:);
n = length(r0);

L = 200*ones(n,1);  % maximum allowable connection length of 200 m
r1 = zeros(n,1);
z1 = zeros(n,1);
ds = 0.01;  % integration stepsize of 1 cm
ns = L(1)/ds;


% Define the fieldline length integration variable
for k = 1:n
    
    r = r0(k);
    z = z0(k);
    s0 = 0; 

    for ii = 1:ns

        [~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, r, z);

        Br = -1/(2*pi*r)*dpsidz;
        Bz =  1/(2*pi*r)*dpsidr;
        Bphi = rzero*bzero/r;

        B = norm([Br Bphi Bz]);

        if startin(k)
            dir = iSign;
        else
            dir = -iSign;
        end
        
        r = r + dir*ds/B*Br;
        z = z + dir*ds/B*Bz;        
        
        in = inpolygon(r, z, rlim, zlim);               
        
        % crossed the limiter?
        if in ~= startin(k)               
            L(k) = s0 + (ii-1)*ds;
            
            % r,z are close to the limiter now
            % project so that they are literally on the limiter           
            [rz] = distance2curve([rlim' zlim'], [r z]);
            r1(k) = rz(1);
            z1(k) = rz(2);

            break
        end

    end

end






