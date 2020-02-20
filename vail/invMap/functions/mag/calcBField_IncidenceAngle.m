function thetaB = calcBField_IncidenceAngle(r, z, psizr, rg, zg, ...
    bzero, rzero, limdata)
%
% CALCBFIELD_INCIDENCEANGLE
%
%   Compute the angle-of-incidence of the total magnetic field (Bp + BT) at a 
%   point along the limiter with coordinates (r,z).
%
% USAGE: calcBField_IncidenceAngle.m
%
% INPUTS:
%
%   r.........radial coordinate of point on the limiter   [m]
%
%   z.........vertical coordinate of point on the limiter [m]
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
% OUTPUTS: 
%
%   thetaB....field line angle-of-incidence [rad]
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 06/04/2018
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/04/2018
%
%...............................................................................

zlim = limdata(1,:);
rlim = limdata(2,:);

rgeo = mean(rlim);
zgeo = mean(zlim);

% Determine on which limiter segment the point lies

dr = r - rlim;
dz = z - zlim;

[~, idx] = min(sqrt(dr.*dr + dz.*dz));

if z > 0
    
    if z > zlim(idx)
        idxs = [idx idx+1];
    elseif r > rlim(idx)
        idxs = [idx idx+1];
    else
        idxs = [idx-1 idx];
    end
    
else
    
    if zlim(idx) > z
        idxs = [idx-1 idx];
    elseif r > rlim(idx)
        idxs = [idx-1 idx];
    else 
        idxs = [idx idx+1];
    end
    
end

% Compute (unit) normal vector to limiter segment (pointing into device)

if rlim(idxs(1)) == rlim(idxs(2))    % vertical segment
    
    if rlim(idxs(1)) > rgeo
        nhatr = -1;
        nhatz =  0;
    else
        nhatr =  1;
        nhatz =  0;
    end
          
elseif zlim(idxs(1)) == zlim(idxs(2)) % horizontal segment
    
    if zlim(idxs(1)) > zgeo
        nhatr =  0;
        nhatz = -1;
    else
        nhatr =  0;
        nhatz =  1;
    end
    
else
    
    m    = (zlim(idxs(2)) - zlim(idxs(1)))/(rlim(idxs(2)) - rlim(idxs(1)));
    
    if zlim(idxs(1)) > zgeo
        
        if (-1/m) > 0
            nhatr = -1/(1*sqrt(1 + (1/m)*(1/m)));
            nhatz =  1/(m*sqrt(1 + (1/m)*(1/m)));
        else
            nhatr =  1/(1*sqrt(1 + (1/m)*(1/m)));
            nhatz = -1/(m*sqrt(1 + (1/m)*(1/m)));
        end
        
    else
        
        if (-1/m) > 0
            nhatr =  1/(1*sqrt(1 + (1/m)*(1/m)));
            nhatz = -1/(m*sqrt(1 + (1/m)*(1/m)));
        else
            nhatr = -1/(1*sqrt(1 + (1/m)*(1/m)));
            nhatz =  1/(m*sqrt(1 + (1/m)*(1/m)));
        end
        
    end
    
end

% Compute poloidal magnetic field components Br and Bz at the point (r,z)

[~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr, r, z);

Br = -1/(2*pi*r)*(dpsidz);
Bz =  1/(2*pi*r)*(dpsidr);

% Compute toroidal magnetic field BT at the point (r,z)

BT = (rzero*bzero)/r;

% Compute (unit) magnetic field vector

B = sqrt(Br*Br + Bz*Bz + BT*BT);

bhatr = Br/B;
bhatz = Bz/B;
bhatT = BT/B;

% Compute the angle-of-incidence

thetaB = abs(asin(bhatr*nhatr + bhatz*nhatz + bhatT*0));

end
