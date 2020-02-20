function s = calcLimDistance(r, z, limdata)
%
% CALCLIMDISTANCE
%
%   Compute the distance along the limiter (measured in the clockwise direction 
%   relative to the centerstack midplane) of a point on the limiter defined by 
%   the coordinates (r,z).
%
% USAGE: calcLimDistance.m
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
%   s.........distance along the limiter [m]
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 05/29/2018
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 05/29/2018
%
%...............................................................................

zlim = limdata(1,:);
rlim = limdata(2,:);

nlim = length(zlim);

% Ensure that the limiter boundary is a closed curve

if zlim(1) ~= zlim(end)
    
    rlim = [rlim rlim(1)];
    zlim = [zlim zlim(1)];
    
end

% Return either half or full distance if point lies on geometric midplane

if z == 0
    
    drs = diff(rlim);
    dzs = diff(zlim);
    
    s = sum(sqrt(drs.*drs + dzs.*dzs));
    
    if r == rlim(1) % full distance
        return
    else            % half distance
        s = s/2;
        return
    end
    
end

for ii = 1:(nlim-1)
    
    % Compute slope of the limiter segment
    
    rS = rlim(ii);
    zS = zlim(ii);
    rE = rlim(ii+1);
    zE = zlim(ii+1);
    
    % Determine if the point lies on current limiter segment
    
    crossproduct = (z - zS)*(rE - rS) - (r - rS)*(zE - zS);
    
    if abs(crossproduct) < eps
        
        dotproduct = (r - rS)*(rE - rS) + (z - zS)*(zE - zS);
        
        if dotproduct < 0
            onflag = 0;
        else
            slength = (rE - rS)*(rE - rS) + (zE - zS)*(zE - zS);
            if dotproduct > slength
                onflag = 0;
            else
                onflag = 1;
            end
        end
        
    else
        onflag = 0;
    end
    
    if onflag == 1
        
        % Determine total distance up to (but not including) the active segment
        
        drs = diff(rlim(1:ii));
        dzs = diff(zlim(1:ii));
        
        s = sum(sqrt(drs.*drs + dzs.*dzs));
        
        % Determine the distance along the active segment
        
        if rlim(ii) == rlim(ii+1)      % vertical segment
            
            s = s + abs(z - zlim(ii));
            
        elseif zlim(ii) == zlim(ii+1)  % horizontal segment
            
            s = s + abs(r - rlim(ii));
            
        else
            
            dr = r - rlim(ii);
            dz = z - zlim(ii);
            
            s = s + sqrt(dr*dr + dz*dz);
            
        end
        
        return
        
    end
   
end

end
