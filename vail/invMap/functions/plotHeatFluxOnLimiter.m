%...............................................................................
%
% USAGE: plotHeatFluxOnLimiter.m
%
% AUTHOR: Patrick J. Vail
%
% DATE: 10/16/2018
%
% PURPOSE:
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 10/16/2018
%
%...............................................................................

scaleHF = 0.1;

dataHF = scaleHF*dataHF/maxQ;

% Calculate the unit normal vector at each heat flux point on the limiter

nhats = zeros(length(rHF),2);

for ii = 1:length(rHF)
    
    [nhatr, nhatz] = calcLimNormal(rHF(ii), zHF(ii), limdata);
    
    nhats(ii,1) = nhatr;
    nhats(ii,2) = nhatz;
    
end

% Calculate the (r,z) coordinates of the heat flux profile

rzplotHF = zeros(length(rHF),2);

for ii = 1:length(rHF)
    
    nhatr = nhats(ii,1);
    nhatz = nhats(ii,2);
    
    theta = atan(nhatz/nhatr);
    
    r = rHF(ii) + dataHF(ii)*cos(theta);
    z = zHF(ii) + dataHF(ii)*sin(theta);
    
    rzplotHF(ii,1) = r;
    rzplotHF(ii,2) = z;
    
end

plot(rzplotHF(:,1), rzplotHF(:,2), '-r')
    