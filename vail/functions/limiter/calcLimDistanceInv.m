function [r, z] = calcLimDistanceInv(s, limdata)
%
% CALCLIMDISTANCEINV
%
%   Compute the (r,z) coordinates of a point which is a distance of s-units 
%   along the limiter (measured in the clockwise direction relative to the 
%   centerstack midplane).
%
%   This function performs the inverse operation of calcLimDistance.
%
% USAGE: calcLimDistanceInv.m
%
% INPUTS:
%
%   s.........distance along the limiter [m]
%
%   limdata...limiter (r,z) vertices as defined in tok_data_struct
%
% OUTPUTS: 
%
%   r.........radial coordinate of point on the limiter   [m]
%
%   z.........vertical coordinate of point on the limiter [m]
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 06/11/2018
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/11/2018
%
%...............................................................................

zlim = limdata(1,:);
rlim = limdata(2,:);

stest = 0;

ii = 1;

while stest < s
    
    dR = rlim(ii+1) - rlim(ii);
    dZ = zlim(ii+1) - zlim(ii);
    
    stest = stest + sqrt(dR*dR + dZ*dZ);
    
    ii = ii + 1;
    
end

ii = ii - 1;
stest = stest - sqrt(dR*dR + dZ*dZ);

if rlim(ii) == rlim(ii+1)      % vertical segment
    
    r = rlim(ii);
    
    if zlim(ii+1) > zlim(ii)
        z = zlim(ii) + (s - stest);
    else
        z = zlim(ii) - (s - stest);
    end
        
elseif zlim(ii) == zlim(ii+1)  % horizontal segment
    
    z = zlim(ii);
    
    if rlim(ii+1) > rlim(ii)
        r = rlim(ii) + (s - stest);
    else
        r = rlim(ii) - (s - stest);
    end
     
else
    
    lengthRatio = sqrt(dR*dR + dZ*dZ)/(s - stest);
    
    dr = abs(dR/lengthRatio);
    dz = abs(dZ/lengthRatio);
    
    if rlim(ii+1) > rlim(ii)
        r = rlim(ii) + dr;
    else 
        r = rlim(ii) - dr;
    end
    
    if zlim(ii+1) > zlim(ii)
        z = zlim(ii) + dz;
    else
        z = zlim(ii) - dz;
    end
    
end

end
