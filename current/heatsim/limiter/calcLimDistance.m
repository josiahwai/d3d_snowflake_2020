function s = calcLimDistance(r,z,limdata)
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
% AUTHOR: Josiah Wai
%
%...............................................................................

zlim = [limdata(1,:) limdata(1,1)];
rlim = [limdata(2,:) limdata(2,1)];
slim = cumsum([0 sqrt(diff(rlim).^2 + diff(zlim).^2)]);

rlim0 = rlim(1:end-1);
rlim1 = rlim(2:end);
zlim0 = zlim(1:end-1);
zlim1 = zlim(2:end);

slack = sqrt(eps);

s=[];
% figure()
% hold on
% plot(rlim,zlim)
for k = 1:length(r)
    
    % check if the point is literally on a limiter defining pt
    [mindist, dum] = min((r(k)-rlim1).^2 + (z(k)-zlim1).^2);
    
    if mindist < 10*slack
        s(k) = slim(dum+1);
    else
        % index of segment the point lies on
        i = find(...
            r(k) > min(rlim0, rlim1) - slack & ...
            r(k) < max(rlim0, rlim1) + slack & ...
            z(k) > min(zlim0, zlim1) - slack & ...
            z(k) < max(zlim0, zlim1) + slack);
                
        s(k) = slim(i) + norm([r(k)-rlim(i), z(k)-zlim(i)]);
    end

%     scatter(rlim(i:i+1), zlim(i:i+1),'k','filled')
%     scatter(r(k),z(k),'g','filled')
        
end

















