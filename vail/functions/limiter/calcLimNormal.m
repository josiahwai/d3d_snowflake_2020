function [nhatr, nhatz] = calcLimNormal(r, z, limdata)
%
% CALCLIMNORMAL
%
%   Calculate the unit normal vector at a point along the limiter with
%   coordinates (r,z).
%
%   The unit normal pointing into the device is computed.
%
% USAGE: calcLimNormal.m
%
% INPUTS:
%
%   r.........radial coordinate of point on the limiter   [m]
%
%   z.........vertical coordinate of point on the limiter [m]
%
%   limdata...limiter (r,z) vertices as defined in tok_data_struct
%
% OUTPUTS: 
%
%   nhatr.....radial component of unit normal
%
%   nhatz.....vertical component of unit normal
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 10/16/2018
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 10/16/2018
%
%...............................................................................

zlim = limdata(1,:);
rlim = limdata(2,:);

rgeo = mean(rlim);
zgeo = mean(zlim);

% Compute the distance along limiter for each limiter vertex 

slim = zeros(length(zlim),1);

for ii = 1:length(slim)
    
    slim(ii) = calcLimDistance(rlim(ii), zlim(ii), limdata);
    
end

% Compute the distance along limiter for point with coordinates (r,z)

s = calcLimDistance(r, z, limdata);

% Determine on which limiter segment the point lies

[~, idx] = min(abs(s - slim));
 
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

end
