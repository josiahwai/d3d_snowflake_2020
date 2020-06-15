function [r, z] = calcLimDistanceInv(s, limdata)
%
% CALCLIMDISTANCEINV
%
%   Compute the (r(k),z(k)) coordinates of a point which is a distance of s(k)-units
%   along the limiter (measured in the clockwise direction relative to the
%   centerstack midplane).
%
%   This function performs the inverse operation of calcLimDistance.
%
% USAGE: calcLimDistanceInv.m
%
% INPUTS:
%
%   s(k).........distance along the limiter [m]
%
%   limdata...limiter (r(k),z(k)) vertices as defined in tok_data_struct
%
% OUTPUTS:
%
%   r(k).........radial coordinate of point on the limiter   [m]
%
%   z(k).........vertical coordinate of point on the limiter [m]
%
% AUTHOR: Patrick J. Vail
%
% DATE: 06/11/2018
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/11/2018
%   Josiah Wai: loop over multiple s(k) values
%...............................................................................


zlim = limdata(1,:);
rlim = limdata(2,:);

r = [];
z =[];
for k = 1:length(s)
  try
    stest = 0;
    
    ii = 1;
    
    while stest < s(k)
      
      dR = rlim(ii+1) - rlim(ii);
      dZ = zlim(ii+1) - zlim(ii);
      
      stest = stest + sqrt(dR*dR + dZ*dZ);
      
      ii = ii + 1;
      
    end
    
    ii = ii - 1;
    stest = stest - sqrt(dR*dR + dZ*dZ);
    
    if rlim(ii) == rlim(ii+1)      % vertical segment
      
      r(k) = rlim(ii);
      
      if zlim(ii+1) > zlim(ii)
        z(k) = zlim(ii) + (s(k) - stest);
      else
        z(k) = zlim(ii) - (s(k) - stest);
      end
      
    elseif zlim(ii) == zlim(ii+1)  % horizontal segment
      
      z(k) = zlim(ii);
      
      if rlim(ii+1) > rlim(ii)
        r(k) = rlim(ii) + (s(k) - stest);
      else
        r(k) = rlim(ii) - (s(k) - stest);
      end
      
    else
      
      lengthRatio = sqrt(dR*dR + dZ*dZ)/(s(k) - stest);
      
      dr = abs(dR/lengthRatio);
      dz = abs(dZ/lengthRatio);
      
      if rlim(ii+1) > rlim(ii)
        r(k) = rlim(ii) + dr;
      else
        r(k) = rlim(ii) - dr;
      end
      
      if zlim(ii+1) > zlim(ii)
        z(k) = zlim(ii) + dz;
      else
        z(k) = zlim(ii) - dz;
      end
      
    end
  catch
    r(k) = nan;
    z(k) = nan;
  end
end

